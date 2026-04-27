package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.bitsets.BitSet;

import java.util.List;
import java.util.Map;

public class CCD1CP extends CCD1 {

    public CCD1CP(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public CCD1CP(TreeAnnotator.TreeSet treeSet) {
        super(treeSet);
    }

    public CCD1CP(TreeAnnotator.TreeSet treeSet, int numTreesToUse) {
        super(treeSet, numTreesToUse);
    }

    public CCD1CP(TreeAnnotator.TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public CCD1CP(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    /* -- ROOT INITIALIZATION -- */
    @Override
    protected void initializeRootClade(int numLeaves) {
        this.leafArraySize = numLeaves;

        BitSet rootBitSet = BitSet.newBitSet(leafArraySize + 1);
        rootBitSet.set(0, numLeaves + 1); // root cannot be a sampled ancestor, so index numLeaves is set to 1

        this.rootClade = new Clade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);
    }

    /* -- TREE INSERTION -- */
    @Override
    protected Clade cladifyVertex(Node vertex) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex);
            Clade leafClade = cladeMapping.get(leafKey);
            if (leafClade == null) {
                leafClade = addNewClade(leafKey);
            }
            leafClade.increaseOccurrenceCount(vertex.getHeight());
            return leafClade;
        }

        Clade firstChildClade = cladifyVertex(vertex.getChildren().get(0));
        Clade secondChildClade = cladifyVertex(vertex.getChildren().get(1));

        BitSet selfKey = extendedSelfKey(firstChildClade, secondChildClade);
        Clade currentClade = cladeMapping.get(selfKey);
        if (currentClade == null) {
            currentClade = addNewClade(selfKey);
        }
        currentClade.increaseOccurrenceCount(vertex.getHeight());

        CladePartition currentPartition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        if (currentPartition == null) {
            currentPartition = currentClade.createCladePartition(firstChildClade, secondChildClade);
        }
        currentPartition.increaseOccurrenceCount(vertex.getHeight());

        return currentClade;
    }

    private BitSet leafKey(Node vertex) {
        BitSet key = BitSet.newBitSet(leafArraySize + 1);
        key.set(vertex.getNr());
        if (vertex.getLength() != 0) {
            key.set(leafArraySize);
        }
        return key;
    }

    private BitSet extendedSelfKey(Clade firstChildClade, Clade secondChildClade) {
        BitSet key = BitSet.newBitSet(leafArraySize + 1);
        key.or(firstChildClade.getCladeInBits());
        key.or(secondChildClade.getCladeInBits());
        return key;
    }

    /**
     * @return whether this clade is marked as a sampled ancestor
     */
    @Override
    public boolean isSampledAncestor(Clade clade) {
        // the last bit, i.e. at index = number of taxa is SA flag
        // Note: we use 0 to represent SA, and 1 to represent non-SA
        return !clade.getCladeInBits().get(leafArraySize);
    }

    /* -- SAMPLING & MAP -- */

    @Override
    protected double computeParentHeight(CladePartition partition, Node firstChild, Node secondChild) {
        if (isSampledAncestorPartition(partition)) {
            return Math.max(firstChild.getHeight(), secondChild.getHeight());
        } else {
            return Math.max(firstChild.getHeight(), secondChild.getHeight()) + 1.0;
        }
    }

    /**
     * @return whether one of the children is a sampled ancestor
     */
    protected boolean isSampledAncestorPartition(CladePartition partition) {
        return isSampledAncestor(partition.getChildClades()[0]) ||
                isSampledAncestor(partition.getChildClades()[1]);
    }

    /* -- HELPER METHODS -- */

    @Override
    protected Clade reduceCladeCount(Node vertex) {
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize + 1);

        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex);
            Clade leafClade = cladeMapping.get(leafKey);
            leafClade.decreaseOccurrenceCount(vertex.getHeight());
            return leafClade;
        }

        Clade firstChildClade = reduceCladeCount(vertex.getChildren().get(0));
        Clade secondChildClade = reduceCladeCount(vertex.getChildren().get(1));

        BitSet selfKey = extendedSelfKey(firstChildClade, secondChildClade);

        // retrieve clade and reduce count
        Clade currentClade = cladeMapping.get(selfKey);
        currentClade.decreaseOccurrenceCount(vertex.getHeight());

        // reduce counts for its clade partitions
        CladePartition currentPartition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        currentPartition.decreaseOccurrenceCount(vertex.getHeight());

        removeCladePartitionIfNecessary(currentClade, currentPartition);

        return currentClade;
    }

    // public int getNumberOfSampledAncestorTaxa() {
    //     int count = 0;
    //     for (Clade clade : this.getClades()) {
    //         if (clade.isSampledAncestor()) {
    //             count++;
    //         }
    //     }
    //     return count;
    // }

    @Override
    protected BitSet computeCladeToNodeMapping(Node vertex, Map<Clade, Node> map) {
        BitSet bits;
        if (vertex.isLeaf()) {
            bits = leafKey(vertex);
        } else {
            bits = computeCladeToNodeMapping(vertex.getLeft(), map);
            BitSet otherBits = computeCladeToNodeMapping(vertex.getRight(), map);
            bits.or(otherBits);
        }

        Clade clade = getClade(bits);
        if (clade == null) {
            System.err.println("No clade found in CCD for this vertex (" + bits + ").");
        }
        map.put(clade, vertex);

        return bits;
    }

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        if (vertex.isLeaf()) {
            // set bitSet
            BitSet leafKey = leafKey(vertex);
            Clade leafClade = cladeMapping.get(leafKey);

            // probability of the leaf
            if (computeLog) {
                runningProbability[0] += leafClade.getLogProbability();
            } else {
                runningProbability[0] *= leafClade.getProbability();
            }
            return leafClade;
        }

        Clade firstChildClade = computeProbabilityOfVertex(vertex.getChildren().get(0), runningProbability, computeLog);
        Clade secondChildClade = computeProbabilityOfVertex(vertex.getChildren().get(1), runningProbability, computeLog);

        if (computeLog && runningProbability[0] > 0) {
            return null;
        }

        if ((firstChildClade == null) || (secondChildClade == null)) {
            setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        BitSet selfKey = extendedSelfKey(firstChildClade, secondChildClade);
        Clade currentClade = cladeMapping.get(selfKey);
        if (currentClade != null) {
            CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
            if (partition != null) {
                if (computeLog) {
                    runningProbability[0] += partition.getLogCCP();
                } else {
                    runningProbability[0] *= partition.getCCP();
                }
            } else {
                setComputedNoProbability(runningProbability, computeLog);
            }
        } else {
            setComputedNoProbability(runningProbability, computeLog);
        }
        return currentClade;
    }

    /* -- MISC -- */
    @Override
    public String toString() {
        return "CCD1-CP " + super.toString().replaceFirst("CCD1 ", "");
    }

    @Override
    protected String getSampledAncestorInfoString(Clade clade) {
        BitSet bits = clade.getCladeInBits();
        boolean isNonSA = bits.get(leafArraySize);
        return isNonSA ? "sampled ancestor = false" : "sampled ancestor = true";
    }

    @Override
    public AbstractCCD copy() {
        CCD1CP copy = new CCD1CP(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();
        AbstractCCD.buildCopy(this, copy);
        return copy;
    }
}
