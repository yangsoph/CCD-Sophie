package ccd.experiments.regularisation;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.CCD1;
import ccd.model.KRegCCD;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class KRegLOTOSimData {

    public static String dataName;
    public static int rep;
    public static int subsampleSize = 3000;
    public static String separator = ",";

    // KRegCCD parameters
    private static final double MU = KRegCCD.DEFAULT_MU;
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws IOException {

        if (args.length > 1) { // if bash file pass arguments
            dataName = args[0];
            rep = Integer.parseInt(args[1]);
        } else { // default local test
            dataName = "yule50";
            rep = 1;
        }

        // output file
        // String outputPathName = "/Users/zyan598/Documents/GitHub/CCD_regularisation/data/KReg/kreg_LOTO_rep" + rep + ".csv";
        String outputPathName = "/nesi/nobackup/uoa04397/sophie/smoothing/KReg/LOTO_output/kreg_LOTO_" + dataName + "_rep" + rep + ".csv";
        File outputFile = new File(outputPathName);
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);
        // header
        StringBuilder sb = new StringBuilder();
        sb.append("heldOutTreeIndex").append(separator);
        sb.append("ccd1Prob").append(separator);
        sb.append("ccd1LogProb").append(separator);
        sb.append("kregProb").append(separator);
        sb.append("kregLogProb");
        writer.println(sb.toString());
        writer.flush();

        // String inputTreePath = "/Volumes/DYNABOOK/wcss/Yule50/treeSets_subsampleSize1000/" +
        // "thinned" + subsampleSize + "_" + dataName.substring(0, 4) + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
        String inputTreePath = "/nesi/nobackup/uoa04397/sophie/wcss/" + dataName + "/treeSets_subsampleSize" + subsampleSize + "/thinned"
                + subsampleSize + "_" + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(inputTreePath, 0);

        File inputTreeFile = new File(inputTreePath);
        List<Tree> treeList = LoadOrStoreTrees.loadTrees(inputTreeFile, 0);

        CCD1 ccd1 = new CCD1(treeSet);
        KRegCCD kreg = new KRegCCD(treeSet, MU, ALPHA, K, KRegCCD.TailMode.NONE);

        for (int j = 0; j < treeList.size(); j++) {
            Tree heldOutTree = treeList.get(j);
            double ccd1Prob = ccd1.getProbOfHeldOutTree(heldOutTree, ALPHA);
            double ccd1LogProb = ccd1.getLogProbOfHeldOutTree(heldOutTree, ALPHA);
            double kregProb = kreg.getProbabilityOfTree(heldOutTree);
            double kregLogProb = kreg.getLogProbabilityOfTree(heldOutTree);

            sb = new StringBuilder();
            sb.append(j).append(separator);
            sb.append(ccd1Prob).append(separator);
            sb.append(ccd1LogProb).append(separator);
            sb.append(kregProb).append(separator);
            sb.append(kregLogProb);
            writer.println(sb.toString());
            writer.flush();
        }
    }
}
