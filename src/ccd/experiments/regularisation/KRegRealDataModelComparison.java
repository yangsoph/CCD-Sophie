package ccd.experiments.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.algorithms.credibleSets.ProbabilityBasedCredibleSetComputer;
import ccd.model.KRegCCD;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class KRegRealDataModelComparison {
    public static String dataName;

    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws IOException {

        if (args.length > 0) { // if bash file pass arguments
            dataName = args[0];
        }

        String dir = "/nesi/nobackup/uoa04397/sophie/real_data/" + dataName;
        File directory = new File(dir);
        File[] treeFiles = directory.listFiles((d, name) -> name.endsWith(".trees") && !name.startsWith("._"));

        String[] files = new String[treeFiles.length];
        String[] labels = new String[treeFiles.length];

        for (int i = 0; i < treeFiles.length; i++) {
            files[i] = treeFiles[i].getName();
            labels[i] = files[i].replaceFirst("\\.trees$", "").replaceFirst("^" + dataName + "_", "");
        }

        String outputDir = "/nesi/nobackup/uoa04397/sophie/real_data/model_comparison/" + dataName + ".tex";

        comparePosteriorTreeSets(dir, files, labels, outputDir);
    }

    public static void comparePosteriorTreeSets(String directory, String[] filenames, String[] labels, String outputTexFile) throws IOException {

        int n = filenames.length;

        KRegCCD[] ccds = new KRegCCD[n];
        Tree[] mapTrees = new Tree[n];

        // Build one CCD per posterior
        for (int i = 0; i < n; i++) {
            File file = new File(directory, filenames[i]);
            List<Tree> trees = LoadOrStoreTrees.loadTrees(file, 0.1);
            System.out.println("Building kreg CCD for " + filenames[i]);
            ccds[i] = KRegCCD.withOptimisedMu(trees);
            System.out.println("Getting MAP tree");
            mapTrees[i] = ccds[i].getMAPTree();
        }
        String optiMuDir = "/nesi/nobackup/uoa04397/sophie/real_data/model_comparison/optimised_mu.csv";
        writeOptimisedMuCSV(optiMuDir, filenames, ccds);

        double[][] result = new double[n][n];

        // Pairwise comparisons
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    result[i][j] = Double.NaN;
                    continue;
                }
                ProbabilityBasedCredibleSetComputer cred = new ProbabilityBasedCredibleSetComputer(ccds[j]);
                result[i][j] = cred.getCredibleLevel(mapTrees[i]);
                System.out.printf("%s -> %s : %.4f%n", labels[i], labels[j], result[i][j]);
            }
        }

        writeLatexTable(labels, result, outputTexFile);
    }

    private static void writeLatexTable(String[] labels, double[][] values, String outputFile) throws IOException {

        try (PrintWriter out = new PrintWriter(new FileWriter(outputFile))) {
            int n = labels.length;
            out.println("\\begin{tabular}{|l|" + "c|".repeat(n) + "}");
            out.println("\\hline");
            out.print("model");

            for (String label : labels) out.print(" & " + label);

            out.println(" \\\\");
            out.println("\\hline");

            for (int i = 0; i < n; i++) {
                out.print(labels[i]);
                for (int j = 0; j < n; j++) {
                    if (i == j)
                        out.print(" & ");
                    else {
                        double x = values[i][j];
                        String s;
                        if (x >= 0.95) s = "\\\\re{" + String.format("%.2f", x) + "}";
                        else s = String.format("%.2f", x);
                        out.print(" & " + s);
                    }
                }
                out.println(" \\\\");
                out.println("\\hline");
            }
            out.println("\\end{tabular}");
        }
    }

    private static void writeOptimisedMuCSV(String outputFile, String[] filenames, KRegCCD[] ccds) throws IOException {
        try (PrintWriter out = new PrintWriter(new FileWriter(outputFile))) {
            out.println("file,optimised_mu");
            for (int i = 0; i < filenames.length; i++) {
                out.printf("%s,%s%n", filenames[i], ccds[i].getMu());
            }
        }
        System.out.println("Optimised mu written to " + outputFile);
    }
}