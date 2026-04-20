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
        // Clade tempClade = new Clade(this);
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize + 1); // plus 1 bit for sampled ancestor flag
        Clade firstChildClade = null;
        Clade secondChildClade = null;

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            // tempClade.addTaxon(index);
            cladeInBits.set(index);
            if (vertex.getLength() != 0) {
                // tempClade.markAsNonSampledAncestor();
                cladeInBits.set(leafArraySize);
                // } else {
                // tempClade.markAsSampledAncestor();
            }
        } else {
            firstChildClade = cladifyVertex(vertex.getChildren().get(0));
            secondChildClade = cladifyVertex(vertex.getChildren().get(1));

            // tempClade.combineClades(firstChildClade, secondChildClade);
            cladeInBits.or(firstChildClade.getCladeInBits());
            cladeInBits.or(secondChildClade.getCladeInBits());
        }
        // BitSet cladeInBits = tempClade.getCladeInBits();
        Clade currentClade = cladeMapping.get(cladeInBits);
        if (currentClade == null) {
            currentClade = addNewClade(cladeInBits);
        }
        currentClade.increaseOccurrenceCount(vertex.getHeight());

        if (!vertex.isLeaf()) {
            CladePartition currentPartition = currentClade.getCladePartition(firstChildClade, secondChildClade);
            if (currentPartition == null) {
                currentPartition = currentClade.createCladePartition(firstChildClade, secondChildClade);
            }
            currentPartition.increaseOccurrenceCount(vertex.getHeight());
        }
        return currentClade;
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
        // 1. build BitSet to retrieve clade and call recursion
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize + 1);
        Clade firstChildClade = null;
        Clade secondChildClade = null;

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);
            if (vertex.getLength() != 0) { // if vertex is not a sampled ancestor, set the last bit to 1
                cladeInBits.set(leafArraySize);
            }
        } else {
            firstChildClade = reduceCladeCount(vertex.getChildren().get(0));
            secondChildClade = reduceCladeCount(vertex.getChildren().get(1));

            cladeInBits.or(firstChildClade.getCladeInBits());
            cladeInBits.or(secondChildClade.getCladeInBits());
        }

        // 2. retrieve clade and reduce count
        Clade currentClade = this.cladeMapping.get(cladeInBits);
        currentClade.decreaseOccurrenceCount(vertex.getHeight());

        // 3. reduce counts for its clade partitions
        if (!vertex.isLeaf()) {
            CladePartition currentPartition = currentClade.getCladePartition(firstChildClade, secondChildClade);
            currentPartition.decreaseOccurrenceCount(vertex.getHeight());

            removeCladePartitionIfNecessary(currentClade, currentPartition);
        }

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
            bits = BitSet.newBitSet(getSizeOfLeavesArray() + 1);
            bits.set(vertex.getNr());
            if (vertex.getLength() != 0) { // if vertex is not a sampled ancestor, set the last bit to 1
                bits.set(leafArraySize);
            }
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
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize + 1);

        if (vertex.isLeaf()) {
            // set bitSet
            int index = vertex.getNr();
            cladeInBits.set(index);
            if (vertex.getLength() != 0) { // if vertex is not a sampled ancestor, set the last bit to 1
                cladeInBits.set(leafArraySize);
            }

            Clade currentClade = cladeMapping.get(cladeInBits);

            // probability of the leaf
            if (computeLog) {
                runningProbability[0] += currentClade.getLogProbability();
            } else {
                runningProbability[0] *= currentClade.getProbability();
            }
            return cladeMapping.get(cladeInBits);

        } else {
            Clade firstChildClade = computeProbabilityOfVertex(vertex.getChildren().get(0), runningProbability, computeLog);
            Clade secondChildClade = computeProbabilityOfVertex(vertex.getChildren().get(1), runningProbability, computeLog);

            if (computeLog && runningProbability[0] > 0) {
                return null;
            }

            if ((firstChildClade == null) || (secondChildClade == null)) {
                setComputedNoProbability(runningProbability, computeLog);
                return null;
            }

            cladeInBits.or(firstChildClade.getCladeInBits());
            cladeInBits.or(secondChildClade.getCladeInBits());

            Clade currentClade = cladeMapping.get(cladeInBits);
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
    }

    /* -- MISC -- */
    @Override
    public String toString() {
        return "CCD1-CP " + super.toString().replaceFirst("CCD1 ", "");
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
