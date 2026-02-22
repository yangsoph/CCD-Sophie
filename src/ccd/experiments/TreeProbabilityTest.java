package ccd.experiments;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.CCD1;
import ccd.model.WrappedBeastTree;
import ccd.algorithms.LoadOrStoreTrees;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class TreeProbabilityTest {
    public static void main(String[] args) throws IOException {
        // String dataPathSA = "/Volumes/DYNABOOK/spike_model_Douglas2024/cephalopod.spikes_short.trees";
        // TreeAnnotator.MemoryFriendlyTreeSet treeSetSA = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPathSA, 0);
        // File inputTreeFileSA = new File(dataPathSA);
        // List<Tree> treeListSA = LoadOrStoreTrees.loadTrees(inputTreeFileSA, 0);
        // CCD1 ccd1SA = new CCD1(treeSetSA, false);
        // WrappedBeastTree mapTree = new WrappedBeastTree(ccd1SA.getMAPTree());
        // System.out.println(mapTree);
        // System.out.println("MAP tree prob = " + ccd1SA.getMaxTreeProbability());
        // System.out.println("entropy = " + ccd1SA.getEntropy());
        // System.out.println("num of clades = " + ccd1SA.getNumberOfClades());

        String dataPath = "/Users/zyan598/Desktop/local_test/thinned3_yule-n50-1.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        File inputTreeFile = new File(dataPath);
        List<Tree> treeList = LoadOrStoreTrees.loadTrees(inputTreeFile, 0);
        CCD1 ccd1 = new CCD1(treeSet, false);
        WrappedBeastTree mapTree = new WrappedBeastTree(ccd1.getMAPTree());
        // for (int j = 0; j < treeList.size(); j++) { // loop through input tree set
        //     Tree testTree = treeList.get(j);
        //     System.out.println(ccd1.getProbabilityOfTree(testTree));
        // }
    }
}
