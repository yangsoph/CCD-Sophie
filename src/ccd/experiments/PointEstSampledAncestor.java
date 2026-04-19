package ccd.experiments;

import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.TreeDistances;
import ccd.model.CCD1;
import ccd.model.CCD1SJ;
import ccd.model.WrappedBeastTree;
import jam.panels.AddRemovePanel;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class PointEstSampledAncestor {

    public static int reps = 100;

    public static void main(String[] args) throws IOException {

        String outputPathName = "/nesi/nobackup/uoa04397/sophie/fossilBD/experiment/result_SJ_numSA_pointEst.csv";
        File outputFile = new File(outputPathName);
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);

        // header
        String separator = ",";
        StringBuilder sb = new StringBuilder();
        sb.append("rep").append(separator);
        sb.append("numTaxa").append(separator);
        sb.append("numSA_in_CCD1map").append(separator);
        sb.append("numSA_in_truth").append(separator);
        sb.append("RF_dist").append(separator);
        sb.append("SA_dist");
        writer.println(sb.toString());
        writer.flush();

        for (int i = 0; i < reps; i++) {

            // String treeSetDir = "/Volumes/DYNABOOK/wcss_sa/mcmc_tree_set/FBD-" + i + ".trees";
            String treeSetDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/beast_output/FBD-" + i + ".trees";
            // String trueTreeDir = "/Volumes/DYNABOOK/wcss_sa/true_tree/FBD-" + i + "_true_fossilTree.trees";
            String trueTreeDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/true_tree/FBD-" + i + "_true_fossilTree.trees";

            TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeSetDir, 10);
            TreeAnnotator.MemoryFriendlyTreeSet trueTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(trueTreeDir, 0);
            trueTreeSet.reset();
            WrappedBeastTree trueTree = new WrappedBeastTree(trueTreeSet.next());

            int numTaxa = trueTree.getWrappedTree().getLeafNodeCount();

            CCD1SJ ccd1 = new CCD1SJ(treeSet, false);
            WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd1.getMAPTree());
            int numSACCD1MAP = mapTreeCCD1.getNumberOfSampledAncestors();
            int numSACCD1Truth = trueTree.getNumberOfSampledAncestors();
            int distRFCCD1 = TreeDistances.robinsonsFouldDistance(trueTree, mapTreeCCD1);
            int distSACCD1 = TreeDistances.sampledAncestorDistance(trueTree, mapTreeCCD1);

            sb = new StringBuilder();
            sb.append(i).append(separator);
            sb.append(numTaxa).append(separator);
            sb.append(numSACCD1MAP).append(separator);
            sb.append(numSACCD1Truth).append(separator);
            sb.append(distRFCCD1).append(separator);
            sb.append(distSACCD1);
            writer.println(sb.toString());
            writer.flush();
        }
        writer.close();
    }

}
