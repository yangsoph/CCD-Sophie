package ccd.experiments;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class LOTOSampledAncestor {
    public static int startRep = 77;
    public static int endRep = 100;

    public static void main(String[] args) throws IOException {

        for (int i = startRep; i < endRep; i++) {

            // output file
            String outputPathName = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/experiment/LOTO_output/result_LOTO_rep" + i + ".csv";
            // String outputPathName = "/nesi/nobackup/uoa04397/sophie/fossilBD/experiment/result_LOTO_rep" + i + ".csv";
            File outputFile = new File(outputPathName);
            FileWriter fileWriter = new FileWriter(outputFile);
            PrintWriter writer = new PrintWriter(fileWriter);
            String separator = ",";
            // header
            StringBuilder sb = new StringBuilder();
            sb.append("heldOutTreeIndex").append(separator);
            sb.append("probEI").append(separator);
            sb.append("logProbEI").append(separator);
            sb.append("probEL").append(separator);
            sb.append("logProbEL");
            writer.println(sb.toString());
            writer.flush();

            // input file
            String treeSetDir = "/Volumes/DYNABOOK/wcss_sa/mcmc_tree_set/FBD-" + i + ".trees";
            // String treeSetDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/beast_output/FBD-" + i + ".trees";
            TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeSetDir, 10);

            CCD1SJ ccdSJ = new CCD1SJ(treeSet);
            CCD1CP ccdCP = new CCD1CP(treeSet);

            File inputTreeFile = new File(treeSetDir);
            List<Tree> treeList = LoadOrStoreTrees.loadTrees(inputTreeFile, 0.1);

            for (int j = 0; j < treeList.size(); j++) {
                Tree heldOutTree = treeList.get(j);
                double probSJ = ccdSJ.getProbOfHeldOutTree(heldOutTree, 0);
                double logProbSJ = ccdSJ.getLogProbOfHeldOutTree(heldOutTree, 0);
                double probCP = ccdCP.getProbOfHeldOutTree(heldOutTree, 0);
                double logProbCP = ccdCP.getLogProbOfHeldOutTree(heldOutTree, 0);

                sb = new StringBuilder();
                sb.append(j).append(separator);
                sb.append(probSJ).append(separator);
                sb.append(logProbSJ).append(separator);
                sb.append(probCP).append(separator);
                sb.append(logProbCP);
                writer.println(sb.toString());
                writer.flush();
            }
        }
    }
}