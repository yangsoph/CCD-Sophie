package ccd.experiments;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.KRegCCD;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Head-to-head of the two novel-clade weightings on a large sample, where SHARED (clustering-
 * penalised: each boundary's escape mass is shared over its pathcount resolutions) is expected to
 * beat FLAT (clustering-neutral) once enough clades are directly sampled. On the 8192-tree block of
 * the mu-vs-sample-size experiment, for each weighting we select mu by cross-validated held-out tree
 * log-probability and then report the entropy and novel-clade probability of the whole-block model
 * at that mu. Higher held-out log-probability is the better fit.
 *
 * Usage: KRegCCDFlatVsShared <trees-file> [burnin=0.1] [folds=5]
 */
public class KRegCCDFlatVsShared {

    private static final int OFFSET = 128 + 512 + 2048; // = 2688
    private static final int SIZE = 8192;
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;
    private static final double MU_LO = 1e-4;
    private static final double MU_HI = KRegCCD.MU_RELIABLE_MAX;
    private static final int MU_GRID = 49;

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: KRegCCDFlatVsShared <trees-file> [burnin=0.1] [folds=5]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;
        int folds = args.length > 2 ? Integer.parseInt(args[2]) : 5;

        List<Tree> trees = LoadOrStoreTrees.loadTrees(treeFile, burnin);
        List<Tree> block = new ArrayList<>(trees.subList(OFFSET, OFFSET + SIZE));
        int nTaxa = block.get(0).getLeafNodeCount();
        int b = Math.max(2, Math.min(folds, SIZE));
        double[] muGrid = logspace(MU_LO, MU_HI, MU_GRID);

        // fixed strided folds (shared by both weightings)
        List<List<Tree>> trainByFold = new ArrayList<>(b);
        List<List<Tree>> testByFold = new ArrayList<>(b);
        int testCount = 0;
        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < SIZE; i++) {
                (i % b == f ? test : train).add(block.get(i));
            }
            trainByFold.add(train);
            testByFold.add(test);
            testCount += test.size();
        }

        System.out.printf("%n=== FLAT vs SHARED on the %d-tree block (trees %d-%d): %s ===%n",
                SIZE, OFFSET + 1, OFFSET + SIZE, treeFile.getName());
        System.out.printf("taxa=%d, alpha=%.2f, k=%d, folds=%d, held-out trees=%d%n%n",
                nTaxa, ALPHA, K, b, testCount);
        System.out.printf("%-8s %-10s %-16s %-10s %-10s%n",
                "mode", "mu*", "held-out nat/tr", "H(mu*)", "P(novel)");

        double[] heldOut = new double[2];
        KRegCCD.NovelMode[] modes = {KRegCCD.NovelMode.FLAT, KRegCCD.NovelMode.SHARED};
        for (int mi = 0; mi < modes.length; mi++) {
            KRegCCD.NovelMode mode = modes[mi];
            KRegCCD[] models = new KRegCCD[b];
            for (int f = 0; f < b; f++) {
                models[f] = new KRegCCD(trainByFold.get(f), 0.0, KRegCCD.DEFAULT_MU, ALPHA, K,
                        KRegCCD.TailMode.NONE, mode);
            }
            double bestMu = muGrid[0];
            double bestTotal = Double.NEGATIVE_INFINITY;
            for (double mu : muGrid) {
                double total = 0.0;
                for (int f = 0; f < b; f++) {
                    for (Tree t : testByFold.get(f)) {
                        total += models[f].getLogProbabilityOfTree(t, mu);
                    }
                }
                if (total > bestTotal) {
                    bestTotal = total;
                    bestMu = mu;
                }
            }
            heldOut[mi] = bestTotal;

            KRegCCD whole = new KRegCCD(block, 0.0, bestMu, ALPHA, K, KRegCCD.TailMode.NONE, mode);
            double entropy = whole.getEntropyRecursive();
            double pNovel = whole.getNovelCladeProbability();
            System.out.printf("%-8s %-10.5f %-16.4f %-10.4f %-10.4f%n",
                    mode, bestMu, bestTotal / testCount, entropy, pNovel);
        }

        double diff = heldOut[1] - heldOut[0]; // SHARED - FLAT
        System.out.printf("%n=> %s fits better: held-out logP higher by %.2f total (%.4f nat/tree over %d trees)%n",
                diff > 0 ? "SHARED" : "FLAT", Math.abs(diff), Math.abs(diff) / testCount, testCount);
    }

    private static double[] logspace(double lo, double hi, int n) {
        double[] grid = new double[n];
        double logLo = Math.log(lo), logHi = Math.log(hi);
        for (int i = 0; i < n; i++) {
            grid[i] = Math.exp(logLo + (logHi - logLo) * i / (n - 1));
        }
        return grid;
    }
}
