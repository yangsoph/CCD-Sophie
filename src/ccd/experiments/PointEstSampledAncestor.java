package ccd.experiments;

import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.TreeDistances;
import ccd.model.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class PointEstSampledAncestor {

    public static int startRep = 0;
    public static int endRep = 100;

    public static void main(String[] args) throws IOException {

        // String outputPathName = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/experiment/local_test.csv";
        String outputPathName = "/nesi/nobackup/uoa04397/sophie/fossilBD/experiment/result_numSA_pointEst_CP.csv";
        File outputFile = new File(outputPathName);
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);

        // header. RF_dist and sym_SA_dist preserved for backward compatibility
        // with existing analysis scripts. topoRF_dist and SA_RF_dist are the
        // new columns: pure topology RF (taxa-only, halved) and the proper
        // sampled-ancestor RF defined in the paper (full symmetric difference
        // of cuttable clade sets).
        String separator = ",";
        StringBuilder sb = new StringBuilder();
        sb.append("rep").append(separator);
        sb.append("numTaxa").append(separator);
        sb.append("numSA_in_CCD1map").append(separator);
        sb.append("numSA_in_truth").append(separator);
        sb.append("RF_dist").append(separator);
        sb.append("sym_SA_dist").append(separator);
        sb.append("topoRF_dist").append(separator);
        sb.append("SA_RF_dist");
        writer.println(sb.toString());
        writer.flush();

        for (int i = startRep; i < endRep; i++) {

            // String treeSetDir = "/Volumes/DYNABOOK/wcss_sa/mcmc_tree_set/FBD-" + i + ".trees";
            String treeSetDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/beast_output/FBD-" + i + ".trees";
            // String trueTreeDir = "/Volumes/DYNABOOK/wcss_sa/true_tree/FBD-" + i + "_true_fossilTree.trees";
            String trueTreeDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/true_tree/FBD-" + i + "_true_fossilTree.trees";

            TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeSetDir, 10);
            TreeAnnotator.MemoryFriendlyTreeSet trueTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(trueTreeDir, 0);
            trueTreeSet.reset();
            WrappedBeastTreeWithSampledAncestor trueTree = new WrappedBeastTreeWithSampledAncestor(trueTreeSet.next());
            int numTaxa = trueTree.getWrappedTree().getLeafNodeCount();

            CCD1CP ccd1 = new CCD1CP(treeSet, false);
            // CCD1SJ ccd1 = new CCD1SJ(treeSet, false);
            WrappedBeastTreeWithSampledAncestor mapTreeCCD1 = new WrappedBeastTreeWithSampledAncestor(ccd1.getMAPTree());

            int numSACCD1MAP = mapTreeCCD1.getNumberOfSampledAncestors();
            int numSACCD1Truth = trueTree.getNumberOfSampledAncestors();
            int distRFCCD1 = TreeDistances.robinsonsFouldDistance(trueTree, mapTreeCCD1);
            double symDistSACCD1 = TreeDistances.sampledAncestorSymmetricDistance(trueTree, mapTreeCCD1);
            int distTopoRFCCD1 = TreeDistances.topologyRobinsonFoulds(trueTree, mapTreeCCD1);
            int distSARFCCD1 = TreeDistances.sampledAncestorRobinsonFoulds(trueTree, mapTreeCCD1);

            sb = new StringBuilder();
            sb.append(i).append(separator);
            sb.append(numTaxa).append(separator);
            sb.append(numSACCD1MAP).append(separator);
            sb.append(numSACCD1Truth).append(separator);
            sb.append(distRFCCD1).append(separator);
            sb.append(symDistSACCD1).append(separator);
            sb.append(distTopoRFCCD1).append(separator);
            sb.append(distSARFCCD1);
            writer.println(sb.toString());
            writer.flush();
        }
        writer.close();
    }

}
