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

}
