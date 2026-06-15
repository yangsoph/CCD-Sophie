package ccd.experiments.sampledAncestor;

import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.TreeDistances;
import ccd.model.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class CCDGraphSizeTest {

    public static int startRep = 0;
    public static int endRep = 100;

    public static void main(String[] args) throws IOException {

        String outputPathName = "/nesi/nobackup/uoa04397/sophie/fossilBD/experiment/result_CCD_graph_size.csv";
        File outputFile = new File(outputPathName);
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);

        // header
        String separator = ",";
        StringBuilder sb = new StringBuilder();
        sb.append("rep").append(separator);
        sb.append("numTaxa").append(separator);
        sb.append("CCD1CP_numCLades").append(separator);
        sb.append("CCD1CP_numPartitions").append(separator);
        sb.append("CCD1SJ_numCLades").append(separator);
        sb.append("CCD1SJ_numPartitions").append(separator);
        sb.append("CCD0CP_numCLades").append(separator);
        sb.append("CCD0CP_numPartitions").append(separator);
        sb.append("CCD0SJ_numCLades").append(separator);
        sb.append("CCD0SJ_numPartitions");
        writer.println(sb.toString());
        writer.flush();

        for (int i = startRep; i < endRep; i++) {

            String treeSetDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/beast_output/FBD-" + i + ".trees";
            String trueTreeDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/true_tree/FBD-" + i + "_true_fossilTree.trees";

            TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeSetDir, 10);
            TreeAnnotator.MemoryFriendlyTreeSet trueTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(trueTreeDir, 0);
            trueTreeSet.reset();
            WrappedBeastTreeWithSampledAncestor trueTree = new WrappedBeastTreeWithSampledAncestor(trueTreeSet.next());
            int numTaxa = trueTree.getWrappedTree().getLeafNodeCount();

            CCD1CP ccd1cp = new CCD1CP(treeSet, false);
            CCD1SJ ccd1sj = new CCD1SJ(treeSet, false);
            CCD0CP ccd0cp = new CCD0CP(treeSet, false);
            CCD0SJ ccd0sj = new CCD0SJ(treeSet, false);

            sb = new StringBuilder();
            sb.append(i).append(separator);
            sb.append(numTaxa).append(separator);
            sb.append(ccd1cp.getNumberOfClades()).append(separator);
            sb.append(ccd1cp.getNumberOfCladePartitions()).append(separator);
            sb.append(ccd1sj.getNumberOfClades()).append(separator);
            sb.append(ccd1sj.getNumberOfCladePartitions()).append(separator);
            sb.append(ccd0cp.getNumberOfClades()).append(separator);
            sb.append(ccd0cp.getNumberOfCladePartitions()).append(separator);
            sb.append(ccd0sj.getNumberOfClades()).append(separator);
            sb.append(ccd0sj.getNumberOfCladePartitions());
            writer.println(sb.toString());
            writer.flush();
        }
        writer.close();
    }

}
