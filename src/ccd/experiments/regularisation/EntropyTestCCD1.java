package ccd.experiments.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccp.model.CCD1;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class EntropyTestCCD1 {
    public static String dataName;
    public static int subsampleSize;
    public static int startRep;
    public static int endRep;

    public static void main(String[] args) throws IOException {

        if (args.length > 1) { // if bash file pass arguments
            dataName = args[0];
            subsampleSize = Integer.parseInt(args[1]);
            startRep = Integer.parseInt(args[2]);
            endRep = Integer.parseInt(args[3]);
        } else { // default local test
            dataName = "Yule50";
            subsampleSize = 3000;
            startRep = 18;
            endRep = 18;
        }

        // output file
        String outputPathName = "/nesi/nobackup/uoa04397/sophie/entropy/entropyCCD1_" + dataName + "_sub" + subsampleSize + "_" + startRep + ".csv";
        // String outputPathName = "/Users/zyan598/Desktop/local_test/entropyCCD1_" + dataName + "_sub" + subsampleSize + "_" + startRep + ".csv";
        File outputFile = new File(outputPathName);
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);

        // header
        String separator = ",";
        StringBuilder sb = new StringBuilder();
        sb.append("rep").append(separator);
        sb.append("entropyCCD1");
        writer.println(sb.toString());
        writer.flush();

        for (int rep = startRep; rep <= endRep; rep++) {

            String treesDir = "/nesi/nobackup/uoa04397/sophie/wcss_full/" + dataName + "/rep" + rep + "/run1/"
                    + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
            // String treesDir = "/Volumes/DYNABOOK/wcss_full/" + dataName + "/rep" + rep + "/run1/"
            //         + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
            File treesFile = new File(treesDir);
            List<Tree> treesList = LoadOrStoreTrees.loadTrees(treesFile, 0, subsampleSize);

            CCD1 ccd = new CCD1(treesList, 0);

            sb = new StringBuilder();
            sb.append(rep).append(separator);
            sb.append(ccd.getEntropy());
            writer.println(sb.toString());
            writer.flush();
        }
    }
}
