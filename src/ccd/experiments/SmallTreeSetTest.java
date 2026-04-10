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

        String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/example_ccd0_test.trees";
        // String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/example_4taxa.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        CCD1 ccd1 = new CCD1(treeSet, false);
        CCD0 ccd0 = new CCD0(treeSet, false);
        WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd1.getMAPTree());
        WrappedBeastTree mapTreeCCD0 = new WrappedBeastTree(ccd0.getMAPTree());
        List<Tree> trees = new ArrayList<>();
        trees.add(mapTreeCCD1.getWrappedTree());
        File treeFile = new File("/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/example_ccd0_map.trees");
        // File treeFile = new File("/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/example_4taxa_map.trees");
        storeTrees(trees, treeFile);
        System.out.println("CCD1 MAP tree log prob = " + ccd1.getMaxLogTreeProbability());
        System.out.println("CCD0 MAP tree log prob = " + ccd0.getMaxLogTreeProbability());
        System.out.println("entropy = " + ccd1.getEntropy());
        System.out.println("num of clades = " + ccd1.getNumberOfClades());
        int[] cladeBin = new int[90];
        int sampledAncestorCount = 0;
        for (Clade clade : ccd1.getClades()) {
            System.out.println(clade);
            for (CladePartition partition : clade.getPartitions()) {
                System.out.println(partition);
            }
            cladeBin[clade.size()] += 1;
            if (clade.isSampledAncestor()) {
                sampledAncestorCount++;
            }
        }
        // for (int i = 0; i < cladeBin.length; i++) {
        //     System.out.println(i + " taxa: " + cladeBin[i]);
        // }
        System.out.println("sampledAncestorCount " + sampledAncestorCount);
    }
}
