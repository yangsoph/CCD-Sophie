package ccd.experiments;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static ccd.algorithms.LoadOrStoreTrees.storeTrees;

public class SmallTreeSetTest {
    public static void main(String[] args) throws IOException {

        String name = "4taxa_alexei.trees";

        String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/" + name;
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        CCD1CP ccd = new CCD1CP(treeSet, false);
        WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd.getMAPTree());
        List<Tree> trees = new ArrayList<>();
        trees.add(mapTreeCCD1.getWrappedTree());
        File treeFile = new File("/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/map_" + name);
        storeTrees(trees, treeFile);
        System.out.println("CCD MAP tree log prob = " + ccd.getMaxLogTreeProbability());
        System.out.println("entropy = " + ccd.getEntropy());
        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            System.out.println(clade.getCladeInBits());
            System.out.println(ccd.isSampledAncestor(clade));
            for (CladePartition partition : clade.getPartitions()) {
                System.out.println(partition);
            }
        }
    }
}
