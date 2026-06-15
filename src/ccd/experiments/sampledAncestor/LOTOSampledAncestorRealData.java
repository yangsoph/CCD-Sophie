package ccd.experiments.sampledAncestor;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.CCD1CP;
import ccd.model.CCD1SJ;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class LOTOSampledAncestorRealData {

    public static void main(String[] args) throws IOException {

        String[] dataList = {
                "douglas_2024_spike_model",
                "near_2021_phylogeny_sunfishes",
                "savelyev_2020_bayesian_turkic_language"
        };

        for (String dataName : dataList) {

            String basePath = "/Volumes/DYNABOOK/phyloData_SA/" + dataName + "/";
            File folder = new File(basePath);
            File[] treeFiles = folder.listFiles((dir, name) -> name.endsWith(".trees") && !name.startsWith("._"));

            if (treeFiles != null) {
                for (File file : treeFiles) {
                    System.out.println("Processing: " + file.getName());
                    String shortName = file.getName().replaceAll("\\.trees$", "");
                    System.out.println("shortName: " + shortName);

                    // output file
                    String outputPathName = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/experiment/LOTO_realData_output/" + shortName + ".csv";
                    File outputFile = new File(outputPathName);
                    FileWriter fileWriter = new FileWriter(outputFile);
                    PrintWriter writer = new PrintWriter(fileWriter);

                    // header
                    String separator = ",";
                    StringBuilder sb = new StringBuilder();
                    sb.append("heldOutTreeIndex").append(separator);
                    sb.append("probEI").append(separator);
                    sb.append("logProbEI").append(separator);
                    sb.append("probEL").append(separator);
                    sb.append("logProbEL");
                    writer.println(sb.toString());
                    writer.flush();

                    // input file
                    String inputTreePath = file.getAbsolutePath();
                    TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(inputTreePath, 0);

                    CCD1SJ ccdSJ = new CCD1SJ(treeSet);
                    CCD1CP ccdCP = new CCD1CP(treeSet);

                    File inputTreeFile = new File(inputTreePath);
                    List<Tree> treeList = LoadOrStoreTrees.loadTrees(inputTreeFile, 0);

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
    }
}
