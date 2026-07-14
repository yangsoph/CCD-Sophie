package ccd.experiments.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.algorithms.credibleSets.ProbabilityBasedCredibleSetComputer;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.KRegCCD;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class KRegCredibleSetExperiment {

    public static String dataName;
    public static int subsampleSize;
    public static String ccdType; // ccd0, ccd1, kreg
    public static int startRep; // the starting index of reps, used in the case when nesi killed job half way
    public static int endRep; // the ending index of reps, used in the case when nesi killed job half way
    public static int precision = 50; // the number of thresholds used to determine credible levels in credible set.

    public static void main(String[] args) throws IOException {

        if (args.length > 1) { // if bash file pass arguments
            dataName = args[0];
            subsampleSize = Integer.parseInt(args[1]);
            ccdType = args[2];
            startRep = Integer.parseInt(args[3]);
            endRep = Integer.parseInt(args[4]);
        } else { // default local test
            dataName = "Yule50";
            subsampleSize = 300;
            ccdType = "ccd1";
            startRep = 1;
            endRep = 1;
        }

        // output file
        // String outputPathName = "/nesi/nobackup/uoa04397/sophie/smoothing/KReg/credSet_fixedAlpha_output/credibleSet_" + dataName + "_sub" + subsampleSize + "_" + ccdType + "_" + startRep + ".csv";
        String outputPathName = "/Users/zyan598/Desktop/local_test/credibleSetCheck.csv";
        File outputFile = new File(outputPathName);
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);

        // header
        String separator = ",";
        StringBuilder sb = new StringBuilder();
        sb.append("rep").append(separator);
        for (int j = 1; j <= precision; j++) {
            double percentage = j / (double) precision;
            sb.append("ccdCredibleSet" + percentage);
            if (j != precision) sb.append(separator);
        }
        writer.println(sb.toString());
        writer.flush();

        for (int rep = startRep; rep <= endRep; rep++) {

            // String trainingTreesDir = "/nesi/nobackup/uoa04397/sophie/wcss_full/" + dataName + "/rep" + rep + "/run1/"
            //         + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
            // String testingTreesDir = "/nesi/nobackup/uoa04397/sophie/wcss_full/" + dataName + "/rep" + rep + "/run2/"
            //         + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";

            String trainingTreesDir = "/Volumes/DYNABOOK/wcss_full/" + dataName + "/rep" + rep + "/run1/"
                    + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
            String testingTreesDir = "/Volumes/DYNABOOK/wcss_full/" + dataName + "/rep" + rep + "/run2/"
                    + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";

            File trainingTreesFile = new File(trainingTreesDir);
            File testingTreesFile = new File(testingTreesDir);

            // TreeAnnotator.MemoryFriendlyTreeSet mcmcTreesHalf1 = new TreeAnnotator().new MemoryFriendlyTreeSet(treeSetDirHalf1, 0);
            List<Tree> trainingTrees = LoadOrStoreTrees.loadTrees(trainingTreesFile, 0, subsampleSize);
            List<Tree> testingTrees = LoadOrStoreTrees.loadTrees(testingTreesFile, 0, subsampleSize);

            ProbabilityBasedCredibleSetComputer probBasedCredSetComputer = null;

            if (ccdType.equals("ccd0")) {
                CCD0 ccd0 = new CCD0(trainingTrees, 0);
                probBasedCredSetComputer = new ProbabilityBasedCredibleSetComputer(ccd0);
            } else if (ccdType.equals("ccd1")) {
                CCD1 ccd1 = new CCD1(trainingTrees, 0);
                probBasedCredSetComputer = new ProbabilityBasedCredibleSetComputer(ccd1);
            } else {
                KRegCCD kreg = KRegCCD.withOptimisedMu(trainingTrees);
                probBasedCredSetComputer = new ProbabilityBasedCredibleSetComputer(kreg);
            }
            int[] countsPerCredibleLevel = getCountsPerCredibleLevel(probBasedCredSetComputer, testingTrees);

            // write output
            sb = new StringBuilder();
            sb.append(rep).append(separator);

            for (int k = 0; k < countsPerCredibleLevel.length; k++) {
                double fraction = (double) countsPerCredibleLevel[k] / subsampleSize;
                sb.append(fraction);
                if (k != countsPerCredibleLevel.length - 1) sb.append(separator);
            }
            writer.println(sb.toString());
            writer.flush();
        }
        writer.flush();
        writer.close();
    }

    public static int[] getCountsPerCredibleLevel(ProbabilityBasedCredibleSetComputer probBasedCredSetComputer, List<Tree> treeList) {
        int[] countsPerCredibleLevel = new int[precision];
        for (Tree tree : treeList) {
            double credibleLevel = probBasedCredSetComputer.getCredibleLevel(tree);
            if (credibleLevel < 0) continue; // skip trees with 0 probability: getCredibleLevel returns -1 if prob == 0
            // int index = (int) Math.floor(credibleLevel * precision) - 1;
            int index = (int) Math.ceil(credibleLevel * precision) - 1;
            for (int i = index; i < precision; i++) {
                countsPerCredibleLevel[i]++;
            }
        }
        return countsPerCredibleLevel;
    }
}
