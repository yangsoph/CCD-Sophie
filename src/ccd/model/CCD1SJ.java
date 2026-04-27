package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.List;

/**
 * CCD1 with extended-split (SJ) identity for sampled-ancestor trees.
 *
 * <p>Each internal clade carries SA flags for its direct leaf children only
 * (no propagation up the tree). Two trees that share the same base topology
 * but differ in which leaf is the sampled ancestor produce different
 * extended-clade identities at the cherry that contains the SA leaf.
 *
 * <p>Leaves are represented as plain (taxon-only) clades; their SA-ness is
 * carried by their parent's flags.
 *
 * <p>The placeholder {@code rootClade} holds one {@link CladePartition} per
 * observed extended root variant: each partition's children are
 * {@code (variant, emptyClade)}, where {@code emptyClade} is a shared sister
 * carrying no taxa. The marginal probability of a root variant is therefore
 * {@code rootClade.getCladePartition(variant, emptyClade).getCCP() =
 * variantCount / numBaseTrees}. This structure keeps the standard CCD graph
 * conventions (every real clade reachable from rootClade via partitions) so
 * algorithms such as clade-probability BFS and per-clade max-CCP caching work
 * on the variants without special casing.
 *
 * <p>Identity bitset layout (width 2n+1):
 * <ul>
 *   <li>bits [0, n): taxa membership of the clade
 *   <li>bits [n, 2n): for an internal clade, bit n+i is set iff leaf i is a
 *       direct leaf child of this clade and is non-SA. For a plain leaf, no
 *       bit in this range is set.
 *   <li>bit 2n: sentinel for the placeholder root clade (only the
 *       AbstractCCD-managed rootClade sets it; variants and other clades do
 *       not use it).
 * </ul>
 */
public class CCD1SJ extends CCD1 {

    /**
     * Shared empty sister used in rootClade's (variant, emptyClade) partitions.
     * Lazy-initialised: super constructor calls cladifyTree before subclass
     * field initialisers run, so this is populated in initializeRootClade.
     */
    protected Clade emptyClade;

    public CCD1SJ(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public CCD1SJ(TreeSet treeSet) {
        super(treeSet);
    }

    public CCD1SJ(TreeSet treeSet, int numTreesToUse) {
        super(treeSet, numTreesToUse);
    }

    public CCD1SJ(TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public CCD1SJ(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    /* -- ROOT INITIALIZATION -- */

    @Override
    protected void initializeRootClade(int numLeaves) {
        this.leafArraySize = numLeaves;

        // Placeholder root clade (sentinel bit 2n set, no taxa bits). Its
        // partitions are (variant, emptyClade) entries added during
        // cladification.
        BitSet rootBitSet = BitSet.newBitSet(2 * numLeaves + 1);
        rootBitSet.set(2 * numLeaves);
        this.rootClade = new Clade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);

        // Shared empty sister: distinct bitset key (no bits set) so it lives
        // in cladeMapping without colliding with any real clade.
        BitSet emptyKey = BitSet.newBitSet(2 * numLeaves + 1);
        this.emptyClade = addNewClade(emptyKey);
    }

    /* -- TREE INSERTION -- */

    @Override
    protected Clade cladifyVertex(Node vertex) {
        return cladifySJVertex(vertex, true);
    }

    private Clade cladifySJVertex(Node vertex, boolean isRoot) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex.getNr());
            Clade leafClade = cladeMapping.get(leafKey);
            if (leafClade == null) {
                leafClade = addNewClade(leafKey);
            }
            leafClade.increaseOccurrenceCount(vertex.getHeight());
            return leafClade;
        }

        Clade firstChildClade = cladifySJVertex(vertex.getChildren().get(0), false);
        Clade secondChildClade = cladifySJVertex(vertex.getChildren().get(1), false);

        BitSet selfKey = extendedSelfKey(vertex);
        Clade currentClade = cladeMapping.get(selfKey);
        if (currentClade == null) {
            currentClade = addNewClade(selfKey);
        }
        currentClade.increaseOccurrenceCount(vertex.getHeight());

        CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        if (partition == null) {
            partition = currentClade.createCladePartition(firstChildClade, secondChildClade);
        }
        partition.increaseOccurrenceCount(vertex.getHeight());

        if (isRoot) {
            // Attach this root variant under the placeholder rootClade via a
            // (variant, emptyClade) partition. rootClade's own occurrence
            // count ticks up once per tree so that
            // partition.getCCP() = variantCount / numBaseTrees.
            rootClade.increaseOccurrenceCount(vertex.getHeight());
            CladePartition rootPart = rootClade.getCladePartition(currentClade, emptyClade);
            if (rootPart == null) {
                rootPart = rootClade.createCladePartition(currentClade, emptyClade);
            }
            rootPart.increaseOccurrenceCount(vertex.getHeight());
        }

        return currentClade;
    }

    /**
     * @return whether the given clade has a sampled ancestor at its root
     */
    @Override
    public boolean isSampledAncestor(Clade clade) {
        BitSet bits = clade.getCladeInBits();
        // check SA region: [n, 2n)
        for (int i = leafArraySize; i < 2 * leafArraySize; i++) {
            if (bits.get(i)) {
                return true;
            }
        }
        return false;
    }

    /* -- PROBABILITY -- */

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        Clade variant = computeProbSJVertex(vertex, runningProbability, computeLog);
        if (variant == null) {
            return null;
        }
        // Multiply in the marginal probability of this root variant,
        // captured by the (variant, emptyClade) partition on rootClade.
        CladePartition rootPart = rootClade.getCladePartition(variant, emptyClade);
        if (rootPart == null) {
            setComputedNoProbability(runningProbability, computeLog);
            return null;
        }
        double ccp = rootPart.getCCP();
        if (computeLog) {
            runningProbability[0] += Math.log(ccp);
        } else {
            runningProbability[0] *= ccp;
        }
        return variant;
    }

    private Clade computeProbSJVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex.getNr());
            Clade leafClade = cladeMapping.get(leafKey);
            if (leafClade == null) {
                setComputedNoProbability(runningProbability, computeLog);
                return null;
            }
            // SJ leaves are plain (one per taxon, present in every tree);
            // marginal 1 — no contribution to the product.
            return leafClade;
        }

        Clade firstChildClade = computeProbSJVertex(vertex.getChildren().get(0), runningProbability, computeLog);
        Clade secondChildClade = computeProbSJVertex(vertex.getChildren().get(1), runningProbability, computeLog);

        if (firstChildClade == null || secondChildClade == null) {
            setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        BitSet selfKey = extendedSelfKey(vertex);
        Clade currentClade = cladeMapping.get(selfKey);
        if (currentClade == null) {
            setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        if (partition == null) {
            setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        double ccp = partition.getCCP();
        if (computeLog) {
            runningProbability[0] += Math.log(ccp);
        } else {
            runningProbability[0] *= ccp;
        }
        return currentClade;
    }

    /* -- KEY CONSTRUCTION -- */

    private BitSet leafKey(int taxonIndex) {
        BitSet key = BitSet.newBitSet(2 * leafArraySize + 1);
        key.set(taxonIndex);
        return key;
    }

    private BitSet extendedSelfKey(Node vertex) {
        BitSet key = BitSet.newBitSet(2 * leafArraySize + 1);
        collectTaxaBits(vertex, key);
        for (Node child : vertex.getChildren()) {
            if (child.isLeaf() && child.getLength() == 0) {
                key.set(leafArraySize + child.getNr());
            }
        }
        return key;
    }

    private void collectTaxaBits(Node vertex, BitSet key) {
        if (vertex.isLeaf()) {
            key.set(vertex.getNr());
        } else {
            for (Node child : vertex.getChildren()) {
                collectTaxaBits(child, key);
            }
        }
    }

    /* -- SAMPLING & MAP -- */

    @Override
    protected Tree getTreeBasedOnStrategy(SamplingStrategy samplingStrategy,
                                          HeightSettingStrategy heightStrategy) {
        tidyUpCacheIfDirty();

        CladePartition rootPart = pickRootPartition(samplingStrategy);
        Clade variant = variantChildOf(rootPart);

        int[] innerIdx = new int[]{this.getSizeOfLeavesArray()};
        Node root = buildSubtree(variant, samplingStrategy, innerIdx);
        return new Tree(root);
    }

    private Clade variantChildOf(CladePartition partition) {
        Clade[] kids = partition.getChildClades();
        return (kids[0] == emptyClade) ? kids[1] : kids[0];
    }

    private CladePartition pickRootPartition(SamplingStrategy strategy) {
        ArrayList<CladePartition> partitions = rootClade.getPartitions();
        if (partitions.isEmpty()) {
            throw new AssertionError("CCD1SJ: rootClade has no variant partitions.");
        }
        if (strategy == SamplingStrategy.MAP) {
            double bestScore = Double.NEGATIVE_INFINITY;
            CladePartition best = partitions.get(0);
            for (CladePartition p : partitions) {
                Clade variant = variantChildOf(p);
                double score = Math.log(p.getCCP()) + variant.getMaxSubtreeLogCCP();
                if (score > bestScore) {
                    bestScore = score;
                    best = p;
                }
            }
            return best;
        }
        // Sampling (and any non-MAP strategy): weighted draw by CCP.
        double total = 0;
        for (CladePartition p : partitions) total += p.getCCP();
        double r = random.nextDouble() * total;
        double cum = 0;
        for (CladePartition p : partitions) {
            cum += p.getCCP();
            if (r < cum) return p;
        }
        return partitions.get(partitions.size() - 1);
    }

    private Node buildSubtree(Clade clade, SamplingStrategy strategy, int[] innerIdx) {
        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = this.getSomeBaseTree().getTaxaNames()[leafNr];
            Node node = new Node(taxonName);
            node.setNr(leafNr);
            // Height set by parent, which knows whether this leaf is SA.
            return node;
        }

        CladePartition partition = pickInternalPartition(clade, strategy);
        if (partition == null) {
            throw new AssertionError("No partitions on clade: " + clade.getCladeInBits());
        }

        Node firstChild = buildSubtree(partition.getChildClades()[0], strategy, innerIdx);
        Node secondChild = buildSubtree(partition.getChildClades()[1], strategy, innerIdx);

        Node node = new Node();
        node.setNr(innerIdx[0]++);
        node.addChild(firstChild);
        node.addChild(secondChild);

        // "One" height semantics: parent = max(children)+1; SA leaves get
        // parent's height (branch length 0).
        double parentHeight = Math.max(safeHeight(firstChild), safeHeight(secondChild)) + 1.0;
        node.setHeight(parentHeight);

        for (Node child : node.getChildren()) {
            if (child.isLeaf()) {
                int leafNr = child.getNr();
                boolean isSA = clade.getCladeInBits().get(leafArraySize + leafNr);
                child.setHeight(isSA ? parentHeight : 0.0);
            }
        }
        return node;
    }

    private CladePartition pickInternalPartition(Clade clade, SamplingStrategy strategy) {
        ArrayList<CladePartition> partitions = clade.getPartitions();
        if (strategy == SamplingStrategy.MAP) {
            return clade.getMaxSubtreeCCPPartition();
        }
        double total = 0;
        for (CladePartition p : partitions) total += p.getCCP();
        double r = random.nextDouble() * total;
        double cum = 0;
        for (CladePartition p : partitions) {
            cum += p.getCCP();
            if (r < cum) return p;
        }
        return partitions.get(partitions.size() - 1);
    }

    private double safeHeight(Node n) {
        return n.isLeaf() ? 0.0 : n.getHeight();
    }

    /* -- MISC -- */

    @Override
    public String toString() {
        return "CCD1-SJ " + super.toString().replaceFirst("CCD1 ", "");
    }

    @Override
    protected String getSampledAncestorInfoString(Clade clade) {
        BitSet bits = clade.getCladeInBits();
        List<Integer> sampledAncestorTaxa = new ArrayList<>();
        for (int i = leafArraySize; i < 2 * leafArraySize; i++) {
            if (bits.get(i)) {
                sampledAncestorTaxa.add(i - leafArraySize);
            }
        }
        if (sampledAncestorTaxa.isEmpty()) {
            return "sampled ancestor = none";
        }
        return "sampled ancestor taxa = " + sampledAncestorTaxa;
    }

    @Override
    public AbstractCCD copy() {
        CCD1SJ copy = new CCD1SJ(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();
        AbstractCCD.buildCopy(this, copy);
        return copy;
    }
}
