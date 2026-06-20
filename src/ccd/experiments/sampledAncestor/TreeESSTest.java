package ccd.experiments.sampledAncestor;

import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;
import ccd.model.CCD1SJ;
import ccd.model.WrappedBeastTreeWithSampledAncestor;
import ccd.algorithms.TreeDistances;
// import ccp.binomialess.TraceStatistics;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class TreeESSTest {

    // public static void main(String[] args) throws IOException {
    //
    //     String outputPathName = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/experiment/ess_0.csv";
    //     // String outputPathName = "/nesi/nobackup/uoa04397/sophie/fossilBD/experiment/ess.csv";
    //     File outputFile = new File(outputPathName);
    //     FileWriter fileWriter = new FileWriter(outputFile);
    //     PrintWriter writer = new PrintWriter(fileWriter);
    //     // header
    //     String separator = ",";
    //     StringBuilder sb = new StringBuilder();
    //     sb.append("rep").append(separator);
    //     sb.append("ess");
    //     writer.println(sb.toString());
    //     writer.flush();
    //
    //     File inputTreeSet;
    //     for (int i = 0; i < 100; i++) {
    //         System.out.println(i);
    //         inputTreeSet = new File("/Volumes/DYNABOOK/wcss_sa/mcmc_tree_set/FBD-" + i + ".trees");
    //         // inputTreeSet = new File("/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/beast_output/FBD-" + i + ".trees");
    //         int ess = estimateESS(inputTreeSet);
    //         sb = new StringBuilder();
    //         sb.append(i).append(separator);
    //         sb.append(ess);
    //         writer.println(sb.toString());
    //         writer.flush();
    //     }
    //     writer.close();
    // }
    //
    // public static int estimateESS(File inputTreesFile) throws IOException {
    //
    //     // calculate CCD0 based on trees
    //     List<Tree> inputTrees = loadTrees(inputTreesFile);
    //     System.out.println("ccd");
    //     // CCD1CP ccd = new CCD1CP(inputTrees.get(0).getLeafNodeCount(), false);
    //     CCD1SJ ccd = new CCD1SJ(inputTrees.get(0).getLeafNodeCount(), false);
    //
    //     for (Tree tree : inputTrees) {
    //         ccd.addTree(tree);
    //     }
    //
    //     ccd.initialize();
    //     WrappedBeastTreeWithSampledAncestor ccdMAP = new WrappedBeastTreeWithSampledAncestor(ccd.getMAPTree());
    //
    //     // calculate expected RF distance to CCD0 MAP
    //     double[] expectedRfTrace = new double[inputTrees.size()];
    //
    //     for (int i = 0; i < inputTrees.size(); i++) {
    //         Tree tree = inputTrees.get(i);
    //         WrappedBeastTreeWithSampledAncestor wrappedTree = new WrappedBeastTreeWithSampledAncestor(tree);
    //         expectedRfTrace[i] = TreeDistances.robinsonsFouldDistance(wrappedTree, ccdMAP);
    //     }
    //
    //     // estimate ESS
    //     TraceStatistics traceStatistics = new TraceStatistics(expectedRfTrace, 1);
    //     int ess = (int) Math.floor(traceStatistics.getESS());
    //
    //     return ess;
    //
    // }
    //
    // public static List<Tree> loadTrees(File treeFile) throws IOException {
    //     NexusParser parser = new NexusParser();
    //     parser.parseFile(treeFile);
    //     return parser.trees;
    // }
}
