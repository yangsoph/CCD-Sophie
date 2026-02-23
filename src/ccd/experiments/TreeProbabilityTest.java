package ccd.experiments;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.WrappedBeastTree;
import ccd.algorithms.LoadOrStoreTrees;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class TreeProbabilityTest {
    public static void main(String[] args) throws IOException {
        String dataPath = "/Volumes/DYNABOOK/spike_model_Douglas2024/cephalopod.spikes.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        CCD1 ccd1 = new CCD1(treeSet, false);
        WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd1.getMAPTree());
        System.out.println("CCD1 MAP tree log prob = " + ccd1.getMaxLogTreeProbability());
        System.out.println("entropy = " + ccd1.getEntropy());
        System.out.println("num of clades = " + ccd1.getNumberOfClades());
        System.out.println("num of splits  = " + ccd1.getNumberOfCladePartitions());
        System.out.println("num of topologies  = " + ccd1.getRootClade().getNumberOfTopologies());

        // CCD0 ccd0 = new CCD0(treeSet, false);
        // WrappedBeastTree mapTreeCCD0 = new WrappedBeastTree(ccd0.getMAPTree());
        // System.out.println("CCD0 MAP tree log prob = " + ccd0.getMaxLogTreeProbability());
        // System.out.println("entropy = " + ccd0.getEntropy());

        // String dataPath = "/Users/zyan598/Desktop/local_test/thinned3_yule-n50-1.trees";
        // TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        // File inputTreeFile = new File(dataPath);
        // List<Tree> treeList = LoadOrStoreTrees.loadTrees(inputTreeFile, 0);
        // CCD1 ccd1 = new CCD1(treeSet, false);
        // WrappedBeastTree mapTree = new WrappedBeastTree(ccd1.getMAPTree());
        // for (int j = 0; j < treeList.size(); j++) { // loop through input tree set
        //     Tree testTree = treeList.get(j);
        //     System.out.println(ccd1.getProbabilityOfTree(testTree));
        // }
    }
}
