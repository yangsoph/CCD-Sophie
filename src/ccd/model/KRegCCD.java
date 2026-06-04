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

    /** log eps(C) per region-top clade; NEGATIVE_INFINITY if C reserves nothing
     * computable (then scoring falls back to per-region mass mu/pathcount). */
    private final Map<Clade, Double> logEpsCache = new HashMap<>();
    private final Map<Clade, Boolean> reservableCache = new HashMap<>();

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
     * @param trees  training trees (post burn-in)
     * @param burnin fraction discarded as burn-in (0.0 if already removed)
     * @param mu     per-clade escape probability, in (0, 1)
     * @param alpha  additive-smoothing pseudocount for the split-expanded backbone
     * @param k      reserve depth: include novel-clade orders {@code 1..k} when
     *               solving {@code eps} (>= 1; {@code k = 1} uses only the NNI
     *               order). Scoring is NOT truncated at {@code k}, so coverage is
     *               full regardless; {@code k} only affects {@code eps} accuracy.
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha, int k) {
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
        this.mu = mu;
        this.log1mMu = Math.log(1.0 - mu);
        this.reserveBoundary = k + 2; // boundary size = novel clades + 2
        expandRegCCD();
        regulariseRegCCD(alpha);
    }

    @Override
    public String toString() {
        return "KRegCCD [mu = " + mu + ", reserve depth k = " + (reserveBoundary - 2)
                + ", full support, split-expanded]";
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
            if (reservable(c)) {
                logp[0] += log1mMu;
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
            double logEps = cladeLogEps(c);
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

    /** Does clade C reserve mu? (admits a shallow novel-clade region). Used only
     * for the (1 - mu) discount on red splits; checks boundaries up to the reserve
     * depth with early exit, which is cheap. */
    private boolean reservable(Clade c) {
        Boolean cached = reservableCache.get(c);
        if (cached != null) {
            return cached;
        }
        boolean result = false;
        if (c.size() >= 3) {
            List<Clade> subs = observedSubclades(c);
            enumOps = 0;
            try {
                for (int j = 3; j <= reserveBoundary && j <= c.size(); j++) {
                    if (countNj(c, subs, j, true) > 0) {
                        result = true;
                        break;
                    }
                }
            } catch (BudgetExceeded e) {
                // keep result as found so far (boundary 3 is O(m^2), always completes)
            }
        }
        reservableCache.put(c, result);
        return result;
    }

    /* log eps(C), solved from the reserve R(C) = sum_{j>=1} N_{j+2} eps^j = mu using
     * orders up to the reserve depth (boundaries 3..reserveBoundary). If those are
     * all zero (a clade with only deeper novel regions), climb to the first nonzero
     * order so eps stays positive. Returns NEGATIVE_INFINITY only if nothing is
     * found within the op budget (then scoring uses a per-region fallback). */
    private double cladeLogEps(Clade c) {
        Double cached = logEpsCache.get(c);
        if (cached != null) {
            return cached;
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
        logEpsCache.put(c, logEps);
        return logEps;
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
