package ccd.experiments;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.*;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.bitsets.BitSet;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static ccd.algorithms.LoadOrStoreTrees.storeTrees;

public class TreeProbabilityTest {
    public static void main(String[] args) throws IOException {

        String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_4taxa_rootSA.trees";
        // String dataPath = "/Volumes/DYNABOOK/spike_model_Douglas_2024/cephalopod.spikes.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        CCD1 ccd1 = new CCD1(treeSet, false);
        WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd1.getMAPTree());
        List<Tree> trees = new ArrayList<>();
        trees.add(mapTreeCCD1.getWrappedTree());
        File treeFile = new File("/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_4taxa_rootSA_map.trees");
        storeTrees(trees, treeFile);
        System.out.println("CCD1 MAP tree log prob = " + ccd1.getMaxLogTreeProbability());
        // System.out.println("CCD1 MAP tree prob = " + ccd1.getMaxTreeProbability());
        System.out.println("entropy = " + ccd1.getEntropy());
        // System.out.println("root prob = " + ccd1.getRootClade().getProbability());
        System.out.println("num of clades = " + ccd1.getNumberOfClades());
        // System.out.println("num of splits  = " + ccd1.getNumberOfCladePartitions());
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
