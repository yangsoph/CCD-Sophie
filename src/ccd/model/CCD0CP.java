package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.bitsets.BitSet;

import java.util.*;

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
        CPSupport.initializeRoot(this, numLeaves);
    }

    @Override
    protected Clade cladifyVertex(Node vertex) {
        return CPSupport.cladifyVertex(this, vertex);
    }

    @Override
    public boolean isSampledAncestor(Clade clade) {
        return CPSupport.isSampledAncestor(clade, leafArraySize);
    }

    @Override
    protected String getSampledAncestorInfoString(Clade clade) {
        return CPSupport.saInfoString(clade, leafArraySize);
    }

    @Override
    protected double computeParentHeight(CladePartition partition, Node firstChild, Node secondChild) {
        return CPSupport.computeParentHeight(partition, firstChild, secondChild, leafArraySize);
    }

    @Override
    protected void expand() {
        List<Clade> realClades = new ArrayList<>();
        Map<BitSet, List<Clade>> taxaMaskToClades = new HashMap<>();
        for (Clade c : cladeMapping.values()) {
            realClades.add(c);
            BitSet taxaMask = c.getCladeInBitsTaxaOnly();
            taxaMaskToClades.computeIfAbsent(taxaMask, k -> new ArrayList<>()).add(c);
        }

        int n = realClades.size();
        for (int i = 0; i < n; i++) {
            Clade A = realClades.get(i);
            BitSet taxaA = A.getCladeInBitsTaxaOnly();
            for (int j = i + 1; j < n; j++) {
                Clade B = realClades.get(j);
                BitSet taxaB = B.getCladeInBitsTaxaOnly();

                // Taxa-disjointness check via base mask.
                if (taxaA.intersects(taxaB)) continue;

                // Compute taxa union and look up all parent clades sharing
                // that taxa-mask (any number of SA variants).
                BitSet unionMask = (BitSet) taxaA.clone();
                unionMask.or(taxaB);

                List<Clade> parents = taxaMaskToClades.get(unionMask);
                if (parents == null) continue;

                for (Clade parent : parents) {
                    if (parent == A || parent == B) continue;
                    if (parent.size() <= 2) continue;
                    if (parent.getCladePartition(A, B) != null) continue;
                    parent.createCladePartition(A, B);
                }
            }
        }
    }

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        return CPSupport.computeProbabilityOfVertex(vertex, runningProbability, computeLog, this, leafArraySize);
    }

    @Override
    public double setPartitionProbabilities(Clade clade, boolean useCladeParameters) {

        if (clade.getSumCladeCredibilities() > 0) {
            return clade.getSumCladeCredibilities();
        }

        double cladeValue = useCladeParameters ? clade.getCladeParameter() : clade.getCladeCredibility();

        // Leaf: still trivial, bcs all leaves, SA or not, all have probability 1 locally
        if (clade.isLeaf()) {
            clade.setSumCladeCredibilities(1.0);
            return 1.0;
        }

        // NO special cherry case, bcs cherry is not trivial anymore

        double totalOccurrences = clade.getNumberOfOccurrences();
        double sumSubtreeProbabilities = 0.0;

        for (CladePartition partition : clade.getPartitions()) {

            // --- CP normalization ---
            double occ = partition.getNumberOfOccurrences();
            double ccp = (totalOccurrences == 0) ? 0.0 : occ / totalOccurrences;
            partition.setCCP(ccp);

            // --- recurse ---
            double left = setPartitionProbabilities(partition.getChildClades()[0], useCladeParameters);
            double right = setPartitionProbabilities(partition.getChildClades()[1], useCladeParameters);

            sumSubtreeProbabilities += ccp * left * right;
        }

        double sumCladeCredibilities = sumSubtreeProbabilities * cladeValue;
        clade.setSumCladeCredibilities(sumCladeCredibilities);

        return sumCladeCredibilities;
    }

    @Override
    public String toString() {
        return "CCD0-CP " + super.toString().replaceFirst("CCD0 ", "");
    }

    @Override
    public AbstractCCD copy() {
        CCD0CP copy = new CCD0CP(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();
        AbstractCCD.buildCopy(this, copy);
        return copy;
    }

}
