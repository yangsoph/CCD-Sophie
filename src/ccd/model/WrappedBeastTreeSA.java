package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.bitsets.BitSet;

public class WrappedBeastTreeSA extends WrappedBeastTree {

    public WrappedBeastTreeSA(Tree wrappedTree) {
        super(wrappedTree);
    }

    @Override
    protected BitSet initCladeBitSet(Node vertex) {
        BitSet cladeAsBitSet;
        if (vertex.isLeaf()) {
            cladeAsBitSet = leafKey(vertex.getNr());
        } else {
            for (Node child : vertex.getChildren()) {
                initCladeBitSet(child);
            }
            cladeAsBitSet = extendedSelfKey(vertex);
        }
        cladeOfVertex[vertex.getNr()] = cladeAsBitSet;
        return cladeAsBitSet;
    }

    private BitSet leafKey(int taxonIndex) {
        BitSet key = BitSet.newBitSet(2 * wrappedTree.getLeafNodeCount());
        key.set(taxonIndex);
        return key;
    }

    private BitSet extendedSelfKey(Node vertex) {
        BitSet key = BitSet.newBitSet(2 * wrappedTree.getLeafNodeCount());
        collectTaxaBits(vertex, key);
        for (Node child : vertex.getChildren()) {
            if (child.isLeaf() && child.getLength() == 0) {
                key.set(wrappedTree.getLeafNodeCount() + child.getNr());
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
}
