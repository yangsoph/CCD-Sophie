package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.List;

public class CCD0CP extends CCD0 {

    public CCD0CP(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public CCD0CP(TreeAnnotator.TreeSet treeSet) {
        super(treeSet);
    }

    public CCD0CP(TreeAnnotator.TreeSet treeSet, int numTreesToUse) {
        super(treeSet, numTreesToUse);
    }

    public CCD0CP(TreeAnnotator.TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public CCD0CP(int numLeaves, boolean storeBaseTrees) {
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

    @Override
    public double setPartitionProbabilities(Clade clade, boolean useCladeParameters) {
        if (clade.getSumCladeCredibilities() > 0) {
            return clade.getSumCladeCredibilities();
        }

        double cladeValue = useCladeParameters ? clade.getCladeParameter() : clade.getCladeCredibility();

        if (clade.isLeaf()) {
            clade.setSumCladeCredibilities(cladeValue);
            return cladeValue;
        } else {
            double sumSubtreeProbabilities = 0.0;
            double[] sumPartitionSubtreeProbabilities = new double[clade.getPartitions().size()];

            // compute sum of probabilities over all partitions ...
            int i = 0;
            for (CladePartition partition : clade.getPartitions()) {
                sumPartitionSubtreeProbabilities[i] =
                        setPartitionProbabilities(partition.getChildClades()[0], useCladeParameters)
                                * setPartitionProbabilities(partition.getChildClades()[1], useCladeParameters);
                sumSubtreeProbabilities += sumPartitionSubtreeProbabilities[i];
                i++;
            }

            // ... and then normalize
            if (sumSubtreeProbabilities == 0) {
                // probability of this subtree is so small, that we have an underflow problem;
                // we can try to use log transformed probabilities to set CCPs
                // but the probability of the clade will still become zero

                // try log-sum-exp trick
                double logMax = Double.NEGATIVE_INFINITY;
                double[] logProbs = new double[clade.getPartitions().size()];
                i = 0;
                for (CladePartition partition : clade.getPartitions()) {
                    double left = partition.getChildClades()[0].getSumCladeCredibilities();
                    double right = partition.getChildClades()[1].getSumCladeCredibilities();
                    logProbs[i] = Math.log(left) + Math.log(right);
                    logMax = Math.max(logProbs[i], logMax);
                    i++;
                }

                // log(Σ exp(xᵢ)) = max + log(Σ exp(xᵢ - max))
                double intermSum = 0;
                for (int j = 0; j < logProbs.length; j++) {
                    if (Double.isFinite(logProbs[j])) {
                        intermSum += Math.exp(logProbs[j] - logMax);
                    }
                }
                double logSum = logMax + Math.log(intermSum);

                // normalizing with normalized log pi = log pi - log sum
                i = 0;
                double pSum = 0;
                for (CladePartition partition : clade.getPartitions()) {
                    if (!Double.isFinite(logProbs[i])) {
                        // we are still encountering underflows
                        throw new UnderflowException("An underflow has occurred.");
                    } else {
                        double logProbability = logProbs[i] - logSum;
                        double probability = Math.exp(logProbability);
                        pSum += probability;

                        if (Double.isNaN(probability)) {
                            out.println("NaN probability = " + probability);
                            out.println("logProbability = " + logProbability);
                            out.println("logProbs[i] = " + logProbs[i]);
                            out.println("logSum = " + logSum);
                        }

                        partition.setCCP(probability);
                    }
                    i++;
                }
            } else {
                i = 0;
                for (CladePartition partition : clade.getPartitions()) {
                    double probability = sumPartitionSubtreeProbabilities[i] / sumSubtreeProbabilities;
                    if (Double.isNaN(probability)) {
                        out.println("clade = " + clade);
                        out.println("partition = " + partition);
                        out.println("sumPartitionSubtreeProbabilities = " + sumPartitionSubtreeProbabilities[i]);
                        out.println("sumSubtreeProbabilities = " + sumSubtreeProbabilities);
                    }
                    partition.setCCP(probability);
                    i++;
                }
            }

            // combined with probability of clade, we get sum of all subtree probabilities
            double sumCladeCredibilities = sumSubtreeProbabilities * cladeValue;
            clade.setSumCladeCredibilities(sumCladeCredibilities);
            return sumCladeCredibilities;
        }
    }

    // @Override
    // protected List<List<Clade>> processCladeBuckets(List<Clade> clades) {
    // }
}
