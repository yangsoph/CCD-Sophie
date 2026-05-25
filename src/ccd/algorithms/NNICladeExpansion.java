package ccd.algorithms;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.bitsets.BitSet;
import ccd.tools.CCDToolUtil;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Clade-level NNI expansion of a CCD: the analogue of the CCD0 split expansion
 * ({@link ccd.model.CCD0#expand()}), but operating on the <em>clades</em>
 * themselves rather than on clade partitions.
 *
 * <p>
 * The CCD0 expansion adds clade partitions {@code P -> {L, R}} that were never
 * observed but whose parent {@code P} and children {@code L}, {@code R} all are.
 * This class instead enumerates novel <em>clades</em> that a single
 * nearest-neighbour-interchange (NNI) move can produce from the splits already
 * stored in the CCD, without ever building a tree.
 * </p>
 *
 * <p>
 * <b>The move at the clade level.</b> An internal edge of a tree is a stored
 * split {@code P -> {L, R}}: the parent clade {@code P} branches into children
 * {@code L} and {@code R}. An NNI at that edge needs the grandchildren: pick a
 * child to descend into, say {@code L -> {A, B}}, with sibling {@code R}. The
 * two NNI moves recombine {@code {A, B, R}} so that the child of {@code P} is no
 * longer {@code L = A ∪ B} but either {@code A ∪ R} or {@code B ∪ R}. Those
 * unions are exactly the novel clades; each comes with the NNI-implied seed
 * split into {@code {kept, sibling}} (e.g. {@code (A ∪ R) -> {A, R}}), whose two
 * constituents are themselves existing clades. We read {@code A}, {@code B},
 * {@code R} straight out of the split table.
 * </p>
 *
 * <p>
 * <b>Pairing modes.</b> The CCD graph mixes splits drawn from many posterior
 * trees, so pairing a split {@code P -> {L, R}} with a grandchild split
 * {@code L -> {A, B}} can combine splits that never co-occurred in one tree.
 * The {@link PairingMode} controls which pairs are allowed:
 * <ul>
 * <li>{@link PairingMode#GRAPH_WIDE} pairs <em>any</em> stored split on
 * {@code P} with <em>any</em> stored split on its child {@code L}. This is the
 * cheapest, broadest reading; it generates clades that no single sampled tree is
 * one NNI move away from.</li>
 * <li>{@link PairingMode#CO_OCCURRING} only recombines grandchild splits that
 * actually appeared together under {@code P} in some single base tree, i.e. the
 * true NNI neighbourhood of the posterior sample. It walks the stored base trees
 * and therefore requires the CCD to have been built with stored base trees
 * (currently non-sampled-ancestor CCDs only).</li>
 * </ul>
 * </p>
 *
 * <p>
 * The pass is read-only: it does not modify the CCD. Cost is proportional to the
 * size of the CCD itself ({@code O(|splits| × d̄)}, with {@code d̄} the average
 * number of splits per clade), not to the number of trees the CCD represents.
 * </p>
 *
 * @author Claude (CCD-Sophie)
 */
public class NNICladeExpansion {

    /** How splits may be paired to form NNI recombinations. */
    public enum PairingMode {
        /** Pair any stored split on P with any stored split on its child L. */
        GRAPH_WIDE,
        /** Only pair splits that co-occurred under P in some single base tree. */
        CO_OCCURRING
    }

    /**
     * One way a novel clade can be produced by an NNI move. The novel clade's
     * taxa are {@code kept ∪ sibling}; its NNI-implied seed split is
     * {@code {kept, sibling}}; it appears as a new child of {@code parent} in
     * place of the original child {@code kept ∪ dropped}.
     *
     * @param parent  the parent clade P at the rearranged edge
     * @param sibling the sibling clade R that joins the kept grandchild
     * @param kept    the grandchild (A or B) that joins the sibling
     * @param dropped the other grandchild, left behind under P
     */
    public record Provenance(Clade parent, Clade sibling, Clade kept, Clade dropped) {
    }

    /** A candidate clade not present in the CCD, with all the ways it arises. */
    public static final class NovelClade {
        /** Taxa of the novel clade (taxa-only bits). */
        public final BitSet taxa;
        /** Number of taxa in the novel clade. */
        public final int size;
        /** Every (parent, sibling, kept, dropped) that produces this clade. */
        public final List<Provenance> provenances = new ArrayList<>(2);

        private NovelClade(BitSet taxa) {
            this.taxa = taxa;
            this.size = taxa.cardinality();
        }
    }

    private final AbstractCCD ccd;
    private final PairingMode mode;

    private final Map<BitSet, NovelClade> novelClades = new HashMap<>();
    private long candidatesEmitted = 0;
    private long candidatesAlreadyPresent = 0;
    private boolean computed = false;

    /** Taxa-only bits of every clade in the CCD, for fast presence tests. */
    private Set<BitSet> existingTaxa;

    /**
     * @param ccd  the CCD whose clade set is expanded (not modified)
     * @param mode how splits may be paired into NNI recombinations
     */
    public NNICladeExpansion(AbstractCCD ccd, PairingMode mode) {
        this.ccd = ccd;
        this.mode = mode;
    }

    /** Runs the expansion if not already done; idempotent. */
    public NNICladeExpansion compute() {
        if (computed) {
            return this;
        }

        existingTaxa = new HashSet<>(ccd.getNumberOfClades() * 2);
        for (Clade clade : ccd.getClades()) {
            existingTaxa.add(clade.getCladeInBitsTaxaOnly());
        }

        switch (mode) {
            case GRAPH_WIDE -> expandGraphWide();
            case CO_OCCURRING -> expandCoOccurring();
        }

        computed = true;
        return this;
    }

    /* -- GRAPH-WIDE EXPANSION (pair any stored adjacent splits) -- */

    private void expandGraphWide() {
        for (Clade parent : ccd.getClades()) {
            // need at least 3 taxa to have a (cherry-or-larger child) + sibling
            if (parent.size() < 3) {
                continue;
            }
            for (CladePartition split : parent.getPartitions()) {
                Clade childA = split.getChildClades()[0];
                Clade childB = split.getChildClades()[1];
                // descend into each child, with the other child as the sibling
                recombineAllGrandchildSplits(parent, childA, childB);
                recombineAllGrandchildSplits(parent, childB, childA);
            }
        }
    }

    /* For parent P, descend-into child L and sibling R; recombine every L-split. */
    private void recombineAllGrandchildSplits(Clade parent, Clade descend, Clade sibling) {
        for (CladePartition grandSplit : descend.getPartitions()) {
            Clade a = grandSplit.getChildClades()[0];
            Clade b = grandSplit.getChildClades()[1];
            emit(parent, sibling, a, b); // novel = A ∪ R, drop B
            emit(parent, sibling, b, a); // novel = B ∪ R, drop A
        }
    }

    /* -- CO-OCCURRING EXPANSION (only splits seen together in one base tree) -- */

    private void expandCoOccurring() {
        List<Tree> baseTrees = ccd.getBaseTrees();
        if (baseTrees == null) {
            throw new IllegalStateException(
                    "CO_OCCURRING pairing needs all base trees stored; "
                            + "build the CCD with storeBaseTrees = true.");
        }
        for (Tree tree : baseTrees) {
            Map<Node, Clade> nodeToClade = new HashMap<>();
            mapNodesToClades(tree.getRoot(), nodeToClade);
            recombineTree(tree.getRoot(), nodeToClade);
        }
    }

    /* Map every node of a base tree to its CCD clade via taxa bits. */
    private BitSet mapNodesToClades(Node vertex, Map<Node, Clade> map) {
        BitSet bits = BitSet.newBitSet(ccd.getSizeOfLeavesArray());
        if (vertex.isLeaf()) {
            bits.set(vertex.getNr());
        } else {
            for (Node child : vertex.getChildren()) {
                bits.or(mapNodesToClades(child, map));
            }
        }
        Clade clade = ccd.getClade(bits);
        if (clade == null) {
            throw new UnsupportedOperationException(
                    "No clade found for a base-tree node; CO_OCCURRING pairing "
                            + "currently supports non-sampled-ancestor CCDs only.");
        }
        map.put(vertex, clade);
        return bits;
    }

    /* Walk internal edges of one base tree, emitting their NNI recombinations. */
    private void recombineTree(Node vertex, Map<Node, Clade> map) {
        if (vertex.isLeaf()) {
            return;
        }
        List<Node> children = vertex.getChildren();
        if (children.size() == 2) {
            Node left = children.get(0);
            Node right = children.get(1);
            Clade parent = map.get(vertex);
            // edge above an internal child = NNI site; the other child is the sibling
            if (!left.isLeaf()) {
                recombineAt(parent, left, right, map);
            }
            if (!right.isLeaf()) {
                recombineAt(parent, right, left, map);
            }
        }
        for (Node child : children) {
            recombineTree(child, map);
        }
    }

    /* The two NNI moves at the edge (parent, descend), with given sibling. */
    private void recombineAt(Clade parent, Node descend, Node siblingNode, Map<Node, Clade> map) {
        if (descend.getChildren().size() != 2) {
            return;
        }
        Clade sibling = map.get(siblingNode);
        Clade a = map.get(descend.getChildren().get(0));
        Clade b = map.get(descend.getChildren().get(1));
        emit(parent, sibling, a, b);
        emit(parent, sibling, b, a);
    }

    /* -- SHARED EMIT -- */

    /* Record the candidate clade kept ∪ sibling; novel iff not already in C. */
    private void emit(Clade parent, Clade sibling, Clade kept, Clade dropped) {
        BitSet taxa = (BitSet) kept.getCladeInBitsTaxaOnly().clone();
        taxa.or(sibling.getCladeInBitsTaxaOnly());

        candidatesEmitted++;
        if (existingTaxa.contains(taxa)) {
            candidatesAlreadyPresent++;
            return;
        }
        NovelClade novel = novelClades.computeIfAbsent(taxa, NovelClade::new);
        novel.provenances.add(new Provenance(parent, sibling, kept, dropped));
    }

    /* -- RESULTS -- */

    /** @return the distinct novel candidate clades (not present in the CCD) */
    public Collection<NovelClade> getNovelClades() {
        compute();
        return novelClades.values();
    }

    /** @return number of distinct novel clades found */
    public int getNumberOfNovelClades() {
        compute();
        return novelClades.size();
    }

    /** @return total candidate emissions, counting multiplicity */
    public long getNumberOfCandidatesEmitted() {
        compute();
        return candidatesEmitted;
    }

    /** @return emissions whose taxa already form a clade in the CCD */
    public long getNumberOfCandidatesAlreadyPresent() {
        compute();
        return candidatesAlreadyPresent;
    }

    /** @return novel clades found per size; index i holds the count of size i */
    public int[] novelSizeHistogram() {
        compute();
        int[] hist = new int[ccd.getSizeOfLeavesArray() + 1];
        for (NovelClade novel : novelClades.values()) {
            hist[novel.size]++;
        }
        return hist;
    }

    /* -- REPORTING -- */

    /** @return a human-readable summary of this single expansion */
    public String report() {
        compute();
        int numClades = ccd.getNumberOfClades();
        int numNovel = novelClades.size();
        double blowup = (numClades == 0) ? 0 : (numNovel / (double) numClades);
        double presentFrac = (candidatesEmitted == 0) ? 0
                : (candidatesAlreadyPresent / (double) candidatesEmitted);

        StringBuilder sb = new StringBuilder();
        sb.append("NNI clade expansion [").append(mode).append("]\n");
        sb.append(String.format("  taxa                 : %d%n", ccd.getNumberOfLeaves()));
        sb.append(String.format("  base trees           : %d%n", ccd.getNumberOfBaseTrees()));
        sb.append(String.format("  clades in C          : %d%n", numClades));
        sb.append(String.format("  candidates emitted   : %d%n", candidatesEmitted));
        sb.append(String.format("  already in C         : %d (%.1f%%)%n",
                candidatesAlreadyPresent, 100 * presentFrac));
        sb.append(String.format("  novel clades         : %d%n", numNovel));
        sb.append(String.format("  blow-up (novel/|C|)  : %.3f%n", blowup));
        sb.append("  novel by size        :");
        int[] hist = novelSizeHistogram();
        for (int s = 2; s < hist.length; s++) {
            if (hist[s] > 0) {
                sb.append(' ').append(s).append(':').append(hist[s]);
            }
        }
        sb.append('\n');
        return sb.toString();
    }

    /**
     * Runs both pairing modes on the given CCD and reports each, plus how the
     * graph-wide novel set relates to the strict co-occurring one.
     *
     * @param ccd CCD to analyse (must have base trees stored for the strict mode)
     * @return a combined report
     */
    public static String compareModes(AbstractCCD ccd) {
        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        StringBuilder sb = new StringBuilder();
        sb.append(broad.report());

        if (ccd.getBaseTrees() == null) {
            sb.append("NNI clade expansion [CO_OCCURRING]: skipped (base trees not stored)\n");
            return sb.toString();
        }

        NNICladeExpansion strict = new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();
        sb.append(strict.report());

        Set<BitSet> broadSet = broad.novelClades.keySet();
        Set<BitSet> strictSet = strict.novelClades.keySet();
        int onlyBroad = 0;
        for (BitSet b : broadSet) {
            if (!strictSet.contains(b)) {
                onlyBroad++;
            }
        }
        int onlyStrict = 0;
        for (BitSet b : strictSet) {
            if (!broadSet.contains(b)) {
                onlyStrict++;
            }
        }
        sb.append("comparison\n");
        sb.append(String.format("  shared novel clades  : %d%n", strictSet.size() - onlyStrict));
        sb.append(String.format("  graph-wide only      : %d%n", onlyBroad));
        sb.append(String.format("  co-occurring only    : %d%n", onlyStrict));
        return sb.toString();
    }

    /**
     * CLI: report NNI clade-expansion properties of a tree file.
     * <p>
     * Usage: {@code NNICladeExpansion <trees-file> [burnin-percentage]}
     *
     * @param args trees file path and optional burn-in percentage (default 10)
     * @throws Exception if the tree file cannot be read
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: NNICladeExpansion <trees-file> [burnin-percentage]");
            System.exit(1);
        }
        int burnin = (args.length >= 2) ? Integer.parseInt(args[1]) : 10;
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(args[0], burnin);
        CCD0 ccd = new CCD0(treeSet, true);
        System.out.print(compareModes(ccd));
    }
}
