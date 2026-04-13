package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

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
 * <p>Each extended root variant lives in {@code cladeMapping} like any other
 * extended clade. At the root, the conditional split probability is
 * {@code partition.count / numBaseTrees} (rather than the parent's own
 * count), matching the SJ semantics where each observed root variant
 * contributes 1/nTrees to the marginal at the root.
 *
 * <p>Identity bitset layout (width 2n+1):
 * <ul>
 *   <li>bits [0, n): taxa membership of the clade
 *   <li>bits [n, 2n): for an internal clade, bit n+i is set iff leaf i is a
 *       direct leaf child of this clade and is non-SA. For a plain leaf, no
 *       bit in this range is set.
 *   <li>bit 2n: sentinel for the placeholder root clade (only the
 *       AbstractCCD-managed rootClade sets it; SJ never uses it for splits).
 * </ul>
 */
public class CCD1SJ extends CCD1 {

    /** Extended root variants (clades whose taxa-bits cover all leaves).
     *  Lazy-initialised: super constructor calls cladifyTree before subclass
     *  field initialisers run, so this is created on first use. */
    protected Set<Clade> extendedRootVariants;

    private Set<Clade> rootVariants() {
        if (extendedRootVariants == null) extendedRootVariants = new LinkedHashSet<>();
        return extendedRootVariants;
    }

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

        // Placeholder root clade. SJ does not use it for probability or
        // sampling; extended root variants are stored in cladeMapping
        // alongside other extended clades. The placeholder exists only to
        // satisfy AbstractCCD invariants that expect a non-null rootClade.
        BitSet rootBitSet = BitSet.newBitSet(2 * numLeaves + 1);
        rootBitSet.set(2 * numLeaves);
        this.rootClade = new Clade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);
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
            if (isRoot) {
                rootVariants().add(currentClade);
            }
        }
        currentClade.increaseOccurrenceCount(vertex.getHeight());

        CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        if (partition == null) {
            partition = currentClade.createCladePartition(firstChildClade, secondChildClade);
        }
        partition.increaseOccurrenceCount(vertex.getHeight());

        return currentClade;
    }

    /* -- PROBABILITY -- */

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        return computeProbSJVertex(vertex, runningProbability, computeLog, true);
    }

    private Clade computeProbSJVertex(Node vertex, double[] runningProbability, boolean computeLog, boolean isRoot) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex.getNr());
            Clade leafClade = cladeMapping.get(leafKey);
            if (leafClade == null) {
                setComputedNoProbability(runningProbability, computeLog);
                return null;
            }
            // SJ leaves are plain (one per taxon, present in every tree); their
            // marginal probability is 1 and contributes nothing to the product.
            // Skip the parent-class BFS-based probability lookup, which is not
            // valid for SJ since computeCladeProbabilities walks down from the
            // unused placeholder rootClade and never visits leaves via splits.
            return leafClade;
        }

        Clade firstChildClade = computeProbSJVertex(vertex.getChildren().get(0), runningProbability, computeLog, false);
        Clade secondChildClade = computeProbSJVertex(vertex.getChildren().get(1), runningProbability, computeLog, false);

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

        double ccp;
        if (isRoot) {
            // Denominator at the root is the total tree count, summed across
            // all extended root variants. Computed locally rather than reading
            // numBaseTrees, which AbstractCCD does not populate when its
            // List<Tree> constructor is called with burnin == 0.
            ccp = partition.getNumberOfOccurrences() / (double) totalRootCount();
        } else {
            ccp = partition.getCCP();
        }
        if (computeLog) {
            runningProbability[0] += Math.log(ccp);
        } else {
            runningProbability[0] *= ccp;
        }
        return currentClade;
    }

    /* -- SAMPLING -- */

    @Override
    public Tree sampleTree(HeightSettingStrategy heightStrategy) {
        return sampleSJTree(heightStrategy);
    }

    private Tree sampleSJTree(HeightSettingStrategy heightStrategy) {
        // 1. Sample a root variant weighted by its occurrence count.
        Clade rootVariant = sampleRootVariant();

        // 2. Recursively sample partitions.
        int[] runningInnerIndex = new int[]{this.getSizeOfLeavesArray()};
        double[] heightCounter = new double[]{0.0};
        Node root = sampleVertex(rootVariant, runningInnerIndex, heightStrategy, heightCounter);
        return new Tree(root);
    }

    private int totalRootCount() {
        int total = 0;
        for (Clade rv : rootVariants()) {
            total += rv.getNumberOfOccurrences();
        }
        return total;
    }

    private Clade sampleRootVariant() {
        Set<Clade> variants = rootVariants();
        int total = totalRootCount();
        int draw = random.nextInt(Math.max(total, 1));
        int cum = 0;
        for (Clade rv : variants) {
            cum += rv.getNumberOfOccurrences();
            if (draw < cum) return rv;
        }
        return variants.iterator().next();
    }

    private Node sampleVertex(Clade clade, int[] runningInnerIndex,
                              HeightSettingStrategy heightStrategy, double[] heightCounter) {
        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = this.getSomeBaseTree().getTaxaNames()[leafNr];
            Node vertex = new Node(taxonName);
            vertex.setNr(leafNr);
            // height set by parent (parent knows whether this leaf is SA)
            return vertex;
        }

        CladePartition partition = sampleCladePartition(clade);
        if (partition == null) {
            throw new AssertionError("No partitions on clade: " + clade.getCladeInBits());
        }

        Node firstChild = sampleVertex(partition.getChildClades()[0], runningInnerIndex, heightStrategy, heightCounter);
        Node secondChild = sampleVertex(partition.getChildClades()[1], runningInnerIndex, heightStrategy, heightCounter);

        Node vertex = new Node();
        vertex.setNr(runningInnerIndex[0]++);
        vertex.addChild(firstChild);
        vertex.addChild(secondChild);

        // Set heights so SA leaves get branch length 0.
        // Default: parent gets height max(child)+1 (heightStrategy=One semantics);
        // a direct leaf child marked SA in this clade's flags is then lifted to
        // share the parent's height.
        double parentHeight = Math.max(safeHeight(firstChild), safeHeight(secondChild)) + 1.0;
        vertex.setHeight(parentHeight);

        for (Node child : vertex.getChildren()) {
            if (child.isLeaf()) {
                int leafNr = child.getNr();
                boolean nonSA = clade.getCladeInBits().get(leafArraySize + leafNr);
                if (nonSA) {
                    child.setHeight(0.0);
                } else {
                    // SA leaf: contemporaneous with parent
                    child.setHeight(parentHeight);
                }
            }
        }
        return vertex;
    }

    private double safeHeight(Node n) {
        return n.isLeaf() ? 0.0 : n.getHeight();
    }

    private CladePartition sampleCladePartition(Clade clade) {
        ArrayList<CladePartition> partitions = clade.getPartitions();
        int total = 0;
        for (CladePartition p : partitions) {
            total += p.getNumberOfOccurrences();
        }
        int draw = random.nextInt(Math.max(total, 1));
        int cum = 0;
        for (CladePartition p : partitions) {
            cum += p.getNumberOfOccurrences();
            if (draw < cum) return p;
        }
        return partitions.get(partitions.size() - 1);
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
            if (child.isLeaf() && child.getLength() != 0) {
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

    /* -- MISC -- */

    @Override
    public String toString() {
        return "CCD1-SJ " + super.toString().replaceFirst("CCD1 ", "");
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
