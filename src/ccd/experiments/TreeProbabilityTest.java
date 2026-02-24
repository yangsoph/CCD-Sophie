package ccd.experiments;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.Clade;
import ccd.model.WrappedBeastTree;
import ccd.algorithms.LoadOrStoreTrees;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TreeProbabilityTest {
    public static void main(String[] args) throws IOException {
        String dataPath = "/Volumes/DYNABOOK/spike_model_Douglas2024/cephalopod.spikes.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        CCD1 ccd1 = new CCD1(treeSet, false);
        WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd1.getMAPTree());
        System.out.println("CCD1 MAP tree log prob = " + ccd1.getMaxLogTreeProbability());
        System.out.println("entropy = " + ccd1.getEntropy());
        System.out.println("root prob = " + ccd1.getRootClade().getProbability());
        System.out.println("num of clades = " + ccd1.getNumberOfClades());
        // System.out.println("num of splits  = " + ccd1.getNumberOfCladePartitions());
        int[] cladeBin = new int[81];
        int sampledAncestorCount = 0;
        for (Clade clade : ccd1.getClades()) {
            cladeBin[clade.size()] += 1;
            if (clade.isSampledAncestor()) {
                sampledAncestorCount++;
            }
        }
        int totalCount = 0;
        for (int i = 0; i < cladeBin.length; i++) {
            totalCount += cladeBin[i];
            System.out.println(i + " taxa: " + cladeBin[i]);
        }
        System.out.println(totalCount);
        System.out.println("sampledAncestorCount " + sampledAncestorCount);

    }
}
