package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Remco's "blue-region" CCD regularisation with <em>full support</em> (every tree
 * gets nonzero probability), parameterised by the reserve depth {@code k}.
 *
 * <p>
 * The split set is first widened by {@link RegCCD}'s CCD split expansion and
 * additive-{@code alpha} smoothing, so a split among observed clades is always
 * representable ("red", priced by its smoothed CCP). The only genuinely unseen
 * structures are then <em>novel clades</em>. To score a tree we colour each
 * internal node red (split in the expanded CCD) or blue (not — i.e. it involves a
 * novel clade), cluster blue nodes into maximal regions, and price:
 * <ul>
 *   <li>each red node at clade {@code C} by {@code (1 - mu) * p(split|C)};</li>
 *   <li>each blue region with top clade {@code C} and {@code j} novel clades by
 *       {@code eps(C)^j / pathcount} — one {@code log eps} per novel clade.
 *       {@code pathcount} is the number of all-novel resolutions of {@code C} into
 *       the region's boundary subclades (a cheap subset-DP over the boundary).</li>
 * </ul>
 * A region with a boundary of {@code m} observed subclades has {@code m - 2} novel
 * clades. Regions of <em>any</em> depth are scored, so any tree gets nonzero
 * probability.
 *
 * <p>
 * The single hyperparameter is {@code mu}, the per-clade escape probability.
 * {@code eps(C)} is the root of {@code R(C) = sum_{j>=1} N_j(C) eps^j = mu}, where
 * {@code N_j(C)} counts the region boundaries of {@code C} with {@code j} novel
 * clades. Computing the exact reserve is #P-hard (a set-partition count over a
 * non-laminar clade family), so we <em>decouple coverage from normalisation</em>:
 * {@code eps} is solved from only the cheap low orders {@code j = 1..k} (boundaries
 * {@code 3..k+2}; {@code N_1} is {@code O(m^2)}, {@code N_2} is {@code O(m^3)} in
 * the number {@code m} of observed subclades), while regions of any depth are still
 * scored with that {@code eps}. The omitted deep orders contribute {@code O(mu^2/m)}
 * to the reserve — a rounding error — so this barely affects probabilities while
 * keeping the model tractable and full-support. Reserve depth {@code k = 2} (the
 * default) is recommended; an op-budget ({@code -Dkreg.enumOps}, default {@code 1e9})
 * bounds even the low-order enumeration on pathological clades, and a clade with
 * only deep regions climbs to its first nonzero order so {@code eps} stays positive.
 *
 * @author Claude (CCD-Sophie)
 */
public class KRegCCD extends RegCCD {

    /** Per-clade enumeration-op budget: a region-top clade's k is capped at the
     * largest boundary whose N_j enumeration stays within this many steps. Based
     * on actual work, so the disjointness pruning keeps most large clades uncapped;
     * only genuine combinatorial blow-ups are capped. Tunable, default 1e9. */
    private static final long OPS_BUDGET =
            Long.parseLong(System.getProperty("kreg.enumOps", "1000000000"));

    /** Thrown internally when a single clade's N_j enumeration exceeds OPS_BUDGET. */
    private static final class BudgetExceeded extends RuntimeException {
        BudgetExceeded() {
            super(null, null, false, false);
        }
    }

    private static final BudgetExceeded BUDGET_EXCEEDED = new BudgetExceeded();

    /** Enumeration steps spent on the current clade's N_j computation. */
    private long enumOps;

    /** Penalty exponent reduction: a boundary-{@code m} region (with {@code m-2}
     * novel clades) is priced {@code eps^(m - EXP_REDUCTION)} = {@code eps^(#novel
     * clades)}, i.e. one factor of {@code eps} per novel clade. */
    private static final int EXP_REDUCTION = 2;

    /** Per-clade escape probability (reserved mass for unseen resolutions). */
    private final double mu;

    private final double log1mMu;

    /** Reserve depth: include novel-clade orders 1..k (boundaries 3..k+2) when
     * solving eps. Scoring is NOT truncated at this depth — every region is scored
     * (full support), so coverage is independent of k. */
    private final int reserveBoundary; // = k + 2

    /**
     * How to correct for the omitted reserve tail (orders &gt; reserve depth) when
     * discounting red splits by {@code 1 - mu - tail(C)}:
     * <ul>
     *   <li>{@link #NONE}: no correction ({@code tail = 0}); slightly
     *       super-normalised (inflates held-out scores by the tail).</li>
     *   <li>{@link #BOUND}: geometric upper bound on the tail; provably
     *       sub-normalised (never inflates).</li>
     *   <li>{@link #SAMPLED}: Knuth estimate of the actual tail; near-exactly
     *       normalised (up to Monte-Carlo noise).</li>
     * </ul>
     */
    public enum TailMode {NONE, BOUND, SAMPLED}

    private final TailMode tailMode;

    /** mu above which the un-corrected inflation becomes material (&gt; ~0.01 nat). */
    private static final double MU_WARN_THRESHOLD = 0.05;

    /** Knuth samples per order when estimating the tail (SAMPLED mode). */
    private static final int TAIL_SAMPLES = 1000;

    /** Below this, the bounded tail is negligible, so SAMPLED skips the sampling. */
    private static final double TAIL_SAMPLE_THRESHOLD = 1e-5;

    private final java.util.Random rng = new java.util.Random(12345);

    /** Per-clade reserve info: whether C reserves mu, log eps(C) (NEGATIVE_INFINITY
     * if nothing computable), and the tail correction to subtract from (1 - mu). */
    private record CladeReg(boolean reservable, double logEps, double tail) {
    }

    private final Map<Clade, CladeReg> regCache = new HashMap<>();

    /**
     * Default variant: reserve depth k = 2 (eps from the one- and two-novel-clade
     * orders), full-support scoring.
     *
     * @param trees  training trees (post burn-in)
     * @param burnin fraction discarded as burn-in (0.0 if already removed)
     * @param mu     per-clade escape probability, in (0, 1)
     * @param alpha  additive-smoothing pseudocount for the split-expanded backbone
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha) {
        this(trees, burnin, mu, alpha, 2);
    }

    /**
     * @param k reserve depth (see below); tail bound correction on.
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha, int k) {
        this(trees, burnin, mu, alpha, k, TailMode.BOUND);
    }

    /**
     * @param trees    training trees (post burn-in)
     * @param burnin   fraction discarded as burn-in (0.0 if already removed)
     * @param mu       per-clade escape probability, in (0, 1)
     * @param alpha    additive-smoothing pseudocount for the split-expanded backbone
     * @param k        reserve depth: include novel-clade orders {@code 1..k} when
     *                 solving {@code eps} (>= 1). Scoring is NOT truncated at
     *                 {@code k}, so coverage is full regardless; {@code k} only
     *                 affects accuracy.
     * @param tailMode how to correct for the omitted reserve tail (see {@link TailMode})
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha, int k,
                   TailMode tailMode) {
        super(trees, burnin, false); // plain CCD1; split expansion + smoothing applied below
        if (mu <= 0 || mu >= 1) {
            throw new IllegalArgumentException("mu must be in (0, 1), got " + mu);
        }
        if (k < 1) {
            throw new IllegalArgumentException("reserve depth k must be >= 1, got " + k);
        }
        if (Double.isNaN(alpha)) {
            throw new IllegalArgumentException("KRegCCD requires split expansion (finite alpha)");
        }
        if (mu > MU_WARN_THRESHOLD) {
            System.err.println("WARNING: KRegCCD mu = " + mu + " > " + MU_WARN_THRESHOLD
                    + "; the reserve-tail inflation grows ~mu^3 and becomes material here "
                    + "(~" + String.format("%.2g", Math.pow(mu / 0.05, 3) * 0.01)
                    + " nat/tree at this mu). Use a smaller mu, or TailMode.BOUND/SAMPLED.");
        }
        this.mu = mu;
        this.log1mMu = Math.log(1.0 - mu);
        this.reserveBoundary = k + 2; // boundary size = novel clades + 2
        this.tailMode = tailMode;
        expandRegCCD();
        regulariseRegCCD(alpha);
    }

    @Override
    public String toString() {
        return "KRegCCD [mu = " + mu + ", reserve depth k = " + (reserveBoundary - 2)
                + ", tail=" + tailMode + ", full support, split-expanded]";
    }

    /* ----------------------------------------------------------------------
     * Tree scoring
     * ------------------------------------------------------------------- */

    @Override
    public double getLogProbabilityOfTree(Tree tree) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(tree.getRoot(), bits);
        double[] logp = new double[]{0.0};
        scoreFresh(tree.getRoot(), bits, logp);
        return logp[0];
    }

    /* ----------------------------------------------------------------------
     * Calibration: total model mass on trees with at least one novel clade
     * ------------------------------------------------------------------- */

    /**
     * The model's total probability mass on trees containing at least one novel
     * clade, {@code P(novel) = 1 - M(root)}, where
     * {@code M(c) = disc(c) * sum_{observed (L,R)} ccp(L,R) * M(L) * M(R)} is the
     * mass on subtrees of {@code c} that stay entirely within the observed clade
     * set, and {@code disc(c) = 1 - mu} for clades that reserve escape mass (the
     * tail correction is negligible for this aggregate). Compare against the
     * empirical fraction of held-out trees that {@link #containsNovelClade}.
     *
     * @return model probability that a tree contains a novel clade
     */
    public double getNovelCladeProbability() {
        return 1.0 - massNoNovelClade(getRootClade(), new HashMap<>());
    }

    private double massNoNovelClade(Clade c, Map<Clade, Double> memo) {
        if (c.isLeaf()) {
            return 1.0;
        }
        Double cached = memo.get(c);
        if (cached != null) {
            return cached;
        }
        double sum = 0.0;
        for (CladePartition p : c.getPartitions()) {
            Clade[] ch = p.getChildClades();
            sum += p.getCCP() * massNoNovelClade(ch[0], memo) * massNoNovelClade(ch[1], memo);
        }
        double m = (cheapReservable(c) ? (1.0 - mu) : 1.0) * sum;
        memo.put(c, m);
        return m;
    }

    /** Cheap reservability test (early-exit over the reserve-depth orders), used
     * for the aggregate {@link #getNovelCladeProbability} where the exact tail is
     * not needed. */
    private boolean cheapReservable(Clade c) {
        if (c.size() < 3) {
            return false;
        }
        List<Clade> subs = observedSubclades(c);
        enumOps = 0;
        try {
            for (int j = 3; j <= reserveBoundary && j <= c.size(); j++) {
                if (countNj(c, subs, j, true) > 0) {
                    return true;
                }
            }
        } catch (BudgetExceeded e) {
            // boundary 3 (O(m^2)) always completes
        }
        return false;
    }

    /** Whether the given tree contains a clade that is not in the (expanded) CCD. */
    public boolean containsNovelClade(Tree tree) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(tree.getRoot(), bits);
        for (Map.Entry<Node, BitSet> e : bits.entrySet()) {
            if (!e.getKey().isLeaf() && getClade(e.getValue()) == null) {
                return true;
            }
        }
        return false;
    }

    private BitSet computeBits(Node v, Map<Node, BitSet> bits) {
        BitSet b = BitSet.newBitSet(leafArraySize);
        if (v.isLeaf()) {
            b.set(v.getNr());
        } else {
            b.or(computeBits(v.getChildren().get(0), bits));
            b.or(computeBits(v.getChildren().get(1), bits));
        }
        bits.put(v, b);
        return b;
    }

    private void scoreFresh(Node v, Map<Node, BitSet> bits, double[] logp) {
        if (v.isLeaf() || logp[0] == Double.NEGATIVE_INFINITY) {
            return;
        }
        Node c1 = v.getChildren().get(0);
        Node c2 = v.getChildren().get(1);
        BitSet vb = bits.get(v);

        if (isRed(vb, bits.get(c1), bits.get(c2))) {
            Clade c = getClade(vb);
            CladeReg reg = computeReg(c);
            if (reg.reservable()) {
                // tail = 0 for NONE -> plain (1 - mu); otherwise (1 - mu - tail)
                logp[0] += Math.log(1.0 - mu - reg.tail());
            }
            logp[0] += rawLogCCP(c, bits.get(c1), bits.get(c2));
            scoreFresh(c1, bits, logp);
            scoreFresh(c2, bits, logp);
        } else {
            // v is the top of a maximal blue region. No truncation: every region
            // is scored, so any tree gets nonzero probability (full support).
            List<Node> boundary = new ArrayList<>();
            collectRegion(v, bits, boundary);
            int m = boundary.size();
            Clade c = getClade(vb); // region top is always an observed clade
            BitSet[] parts = new BitSet[m];
            for (int i = 0; i < m; i++) {
                parts[i] = bits.get(boundary.get(i));
            }
            int pathcount = countAllNovelResolutions(vb, parts);
            double logEps = computeReg(c).logEps();
            if (logEps == Double.NEGATIVE_INFINITY) {
                // C reserves nothing computable: give this region fallback mass
                // mu/pathcount (treat its own boundary as the sole reservation)
                logp[0] += Math.log(mu) - Math.log(pathcount);
            } else {
                // m - 2 = number of novel clades in this region: one log eps each
                logp[0] += (m - EXP_REDUCTION) * logEps - Math.log(pathcount);
            }
            for (Node b : boundary) {
                scoreFresh(b, bits, logp);
            }
        }
    }

    private void collectRegion(Node v, Map<Node, BitSet> bits, List<Node> boundary) {
        for (Node child : v.getChildren()) {
            if (child.isLeaf()) {
                boundary.add(child);
            } else if (isRed(bits.get(child), bits.get(child.getChildren().get(0)),
                    bits.get(child.getChildren().get(1)))) {
                boundary.add(child);
            } else {
                collectRegion(child, bits, boundary);
            }
        }
    }

    /* ----------------------------------------------------------------------
     * Observed-split / observed-clade queries
     * ------------------------------------------------------------------- */

    private boolean isSplitObserved(BitSet parentBits, BitSet c1Bits, BitSet c2Bits) {
        Clade parent = getClade(parentBits);
        if (parent == null) {
            return false;
        }
        Clade a = getClade(c1Bits);
        if (a == null) {
            return false;
        }
        Clade b = getClade(c2Bits);
        if (b == null) {
            return false;
        }
        return parent.getCladePartition(a, b) != null;
    }

    private boolean isRed(BitSet parentBits, BitSet c1Bits, BitSet c2Bits) {
        return isSplitObserved(parentBits, c1Bits, c2Bits);
    }

    private double rawLogCCP(Clade parent, BitSet c1Bits, BitSet c2Bits) {
        Clade a = getClade(c1Bits);
        Clade b = getClade(c2Bits);
        CladePartition p = parent.getCladePartition(a, b);
        return p.getLogCCP();
    }

    /* ----------------------------------------------------------------------
     * pathcount: number of all-novel resolutions of C into the given parts,
     * via a subset DP over the parts (handles any boundary size).
     * ------------------------------------------------------------------- */

    private int countAllNovelResolutions(BitSet cBits, BitSet[] parts) {
        int k = parts.length;
        if (k == 1) {
            return 1;
        }
        int full = (1 << k) - 1;
        BitSet[] unionOf = new BitSet[1 << k];
        unionOf[0] = BitSet.newBitSet(leafArraySize);
        for (int mask = 1; mask <= full; mask++) {
            int low = Integer.numberOfTrailingZeros(mask);
            BitSet u = (BitSet) unionOf[mask & (mask - 1)].clone();
            u.or(parts[low]);
            unionOf[mask] = u;
        }
        int[] f = new int[1 << k];
        for (int mask = 1; mask <= full; mask++) {
            if (Integer.bitCount(mask) == 1) {
                f[mask] = 1;
                continue;
            }
            int low = mask & (-mask); // lowest set bit, fixed in S1 to avoid double counting
            int rest = mask ^ low;
            int count = 0;
            // iterate non-empty proper submasks S1 of mask that contain `low`
            for (int sub = rest; ; sub = (sub - 1) & rest) {
                int s1 = sub | low;
                int s2 = mask ^ s1;
                if (s2 != 0 && !isSplitObserved(unionOf[mask], unionOf[s1], unionOf[s2])) {
                    count += f[s1] * f[s2];
                }
                if (sub == 0) {
                    break;
                }
            }
            f[mask] = count;
        }
        return f[full];
    }

    /* ----------------------------------------------------------------------
     * Per-clade reservation: N_j and eps
     * ------------------------------------------------------------------- */

    /* Per-clade reserve: solve eps from R(C) = sum_{j>=1} N_{j+2} eps^j = mu using
     * orders up to the reserve depth (boundaries 3..reserveBoundary; climbing to the
     * first nonzero order if those are all zero, so eps stays positive), and bound
     * the omitted tail (orders > k) geometrically from the two computed orders:
     *   rho = (N_2/N_1) eps,  tail <= N_2 eps^2 * rho/(1-rho)
     * (rigorous if the N_j are log-concave; clamped to [0, mu]). Computed once per
     * clade and cached. */
    private CladeReg computeReg(Clade c) {
        CladeReg cached = regCache.get(c);
        if (cached != null) {
            return cached;
        }
        CladeReg reg;
        if (c.size() < 3) {
            reg = new CladeReg(false, Double.NEGATIVE_INFINITY, 0.0);
            regCache.put(c, reg);
            return reg;
        }
        List<Clade> subs = observedSubclades(c);
        int[] n = new int[c.size() + 1];
        int last = 2;
        boolean anyNonzero = false;
        enumOps = 0;
        for (int j = 3; j <= c.size(); j++) {
            if (j > reserveBoundary && anyNonzero) {
                break; // reached reserve depth and have a usable (positive) reserve
            }
            try {
                n[j] = countNj(c, subs, j, false); // accumulates enumOps across j
            } catch (BudgetExceeded e) {
                break;
            }
            last = j;
            if (n[j] > 0) {
                anyNonzero = true;
            }
        }
        double logEps = anyNonzero ? solveLogEps(n, last, mu) : Double.NEGATIVE_INFINITY;

        double tail = 0.0;
        int n1 = n[3];
        int n2 = (n.length > 4) ? n[4] : 0;
        if (tailMode != TailMode.NONE && logEps != Double.NEGATIVE_INFINITY && n1 > 0 && n2 > 0) {
            double eps = Math.exp(logEps);
            // geometric upper bound from the two computed orders (rigorous if log-concave)
            double rho = ((double) n2 / n1) * eps;
            double tailUB = (rho > 0 && rho < 1) ? n2 * eps * eps * rho / (1 - rho) : mu;
            tailUB = Math.min(tailUB, mu);

            if (tailMode == TailMode.BOUND || tailUB < TAIL_SAMPLE_THRESHOLD) {
                // BOUND mode, or a negligible bound: use the (cheap) bound directly
                tail = tailUB;
            } else {
                // SAMPLED mode with a non-negligible bound: Knuth-estimate the actual
                // tail (orders > reserve depth), summing the two leading orders
                tail = 0.0;
                for (int m = reserveBoundary + 1; m <= reserveBoundary + 2 && m <= c.size(); m++) {
                    double nm = estimateNj(c, subs, m, TAIL_SAMPLES);
                    tail += nm * Math.pow(eps, m - EXP_REDUCTION);
                }
                tail = Math.min(tail, mu);
            }
        }

        reg = new CladeReg(anyNonzero, logEps, tail);
        regCache.put(c, reg);
        return reg;
    }

    /* Knuth (1975) backtrack-tree estimator of N_m (number of valid boundary-m
     * partitions): average over random root-to-leaf descents of the product of
     * branching factors, times the leaf-valid indicator. Each descent is O(m * |subs|);
     * unbiased. */
    private double estimateNj(Clade c, List<Clade> subs, int m, int samples) {
        BitSet cBits = c.getCladeInBits();
        double sum = 0.0;
        for (int s = 0; s < samples; s++) {
            sum += knuthDescent(cBits, subs, m);
        }
        return sum / samples;
    }

    private double knuthDescent(BitSet cBits, List<Clade> subs, int m) {
        BitSet used = BitSet.newBitSet(leafArraySize);
        double product = 1.0;
        int startIdx = 0;
        BitSet prev = null;
        BitSet[] chosen = new BitSet[m];
        for (int level = 0; level < m - 1; level++) {
            // children = candidates at index >= startIdx disjoint from `used`
            int d = 0;
            for (int i = startIdx; i < subs.size(); i++) {
                if (!subs.get(i).getCladeInBits().intersects(used)) {
                    d++;
                }
            }
            if (d == 0) {
                return 0.0; // dead end
            }
            product *= d;
            int target = rng.nextInt(d);
            int pickIdx = -1, seen = 0;
            for (int i = startIdx; i < subs.size(); i++) {
                if (!subs.get(i).getCladeInBits().intersects(used)) {
                    if (seen == target) {
                        pickIdx = i;
                        break;
                    }
                    seen++;
                }
            }
            BitSet pb = subs.get(pickIdx).getCladeInBits();
            chosen[level] = pb;
            used = (BitSet) used.clone();
            used.or(pb);
            prev = pb;
            startIdx = pickIdx + 1;
        }
        BitSet last = (BitSet) cBits.clone();
        last.andNot(used);
        if (last.cardinality() == 0 || getClade(last) == null) {
            return 0.0;
        }
        if (compareBitSets(prev, last) >= 0) {
            return 0.0; // canonical order violated
        }
        chosen[m - 1] = last;
        return countAllNovelResolutions(cBits, chosen) > 0 ? product : 0.0;
    }

    /* Observed clades strictly contained in C, sorted in canonical (lexicographic) order. */
    private List<Clade> observedSubclades(Clade c) {
        List<Clade> subs = new ArrayList<>();
        int cSize = c.size();
        for (Clade x : getClades()) {
            int s = x.size();
            if (s >= 1 && s < cSize && c.contains(x.getCladeInBitsTaxaOnly())) {
                subs.add(x);
            }
        }
        subs.sort((a, b) -> compareBitSets(a.getCladeInBits(), b.getCladeInBits()));
        return subs;
    }

    /* Count m-part boundaries of C into observed subclades admitting an all-novel
     * resolution. earlyExit returns 1 as soon as one is found. */
    private int countNj(Clade c, List<Clade> subs, int m, boolean earlyExit) {
        BitSet cBits = c.getCladeInBits();
        BitSet used = BitSet.newBitSet(leafArraySize);
        List<BitSet> chosen = new ArrayList<>(m);
        return enumerate(cBits, subs, m, 0, used, chosen, earlyExit);
    }

    /* Recursive canonical enumeration: pick (m-1) parts at strictly increasing
     * indices (= increasing canonical order), then derive the last part as the
     * complement and require it to be the canonical-largest. */
    private int enumerate(BitSet cBits, List<Clade> subs, int m, int startIdx,
                          BitSet used, List<BitSet> chosen, boolean earlyExit) {
        if (++enumOps > OPS_BUDGET) {
            throw BUDGET_EXCEEDED;
        }
        if (chosen.size() == m - 1) {
            BitSet last = (BitSet) cBits.clone();
            last.andNot(used);
            if (last.cardinality() == 0 || getClade(last) == null) {
                return 0;
            }
            BitSet prev = chosen.get(chosen.size() - 1);
            // canonical: last must be strictly greater than the previous (largest) chosen
            if (compareBitSets(prev, last) >= 0) {
                return 0;
            }
            BitSet[] parts = new BitSet[m];
            for (int i = 0; i < m - 1; i++) {
                parts[i] = chosen.get(i);
            }
            parts[m - 1] = last;
            return countAllNovelResolutions(cBits, parts) > 0 ? 1 : 0;
        }
        int count = 0;
        for (int i = startIdx; i < subs.size(); i++) {
            BitSet pb = subs.get(i).getCladeInBits();
            if (pb.intersects(used)) {
                continue;
            }
            chosen.add(pb);
            BitSet newUsed = (BitSet) used.clone();
            newUsed.or(pb);
            count += enumerate(cBits, subs, m, i + 1, newUsed, chosen, earlyExit);
            chosen.remove(chosen.size() - 1);
            if (earlyExit && count > 0) {
                return count;
            }
        }
        return count;
    }

    /* Monotone bisection solve of sum_{m>=3} n[m] eps^(m-2) = mu on eps >= 0
     * (eps^(m-2) = eps^(#novel clades)). */
    private static double solveLogEps(int[] n, int maxDeg, double mu) {
        boolean any = false;
        for (int j = 3; j <= maxDeg; j++) {
            if (n[j] > 0) {
                any = true;
                break;
            }
        }
        if (!any) {
            return Double.NEGATIVE_INFINITY; // not reservable; eps unused
        }
        double lo = 0.0;
        double hi = 1.0;
        while (evalPoly(n, maxDeg, hi) < mu) {
            hi *= 2.0;
        }
        for (int it = 0; it < 100; it++) {
            double mid = 0.5 * (lo + hi);
            if (evalPoly(n, maxDeg, mid) < mu) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return Math.log(0.5 * (lo + hi));
    }

    /* Consistent total order on clade bitsets: compare sorted set-bit indices
     * lexicographically (a proper prefix sorts before its extension). */
    private static int compareBitSets(BitSet a, BitSet b) {
        int ia = a.nextSetBit(0);
        int ib = b.nextSetBit(0);
        while (ia >= 0 && ib >= 0) {
            if (ia != ib) {
                return Integer.compare(ia, ib);
            }
            ia = a.nextSetBit(ia + 1);
            ib = b.nextSetBit(ib + 1);
        }
        // whichever still has set bits is "greater"; -1 (exhausted) sorts first
        return Integer.compare(ia, ib);
    }

    /* sum_{m>=3} n[m] * x^(m - 2). */
    private static double evalPoly(int[] n, int maxDeg, double x) {
        double s = 0.0;
        for (int m = 3; m <= maxDeg; m++) {
            if (n[m] != 0) {
                s += n[m] * Math.pow(x, m - EXP_REDUCTION);
            }
        }
        return s;
    }
}
