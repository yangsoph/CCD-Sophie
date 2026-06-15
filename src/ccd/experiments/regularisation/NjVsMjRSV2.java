package ccd.experiments.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.KRegCCD;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Head-to-head held-out comparison of the two KRegCCD novel-clade weightings on a posterior
 * tree set: SHARED (N_j boundary counts + per-tree /pathcount; clustering-penalised) vs FLAT
 * (M_j tree counts, no /pathcount; clustering-neutral). Both modes are scored on the SAME folds
 * over a shared mu grid (alpha fixed), so the only thing that varies is the weighting. Higher
 * held-out total log-probability fits the posterior better.
 * <p>
 * Usage: NjVsMjRSV2 <trees-file> [burninFraction] [folds] [maxTrees]
 */
public class NjVsMjRSV2 {

    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int RESERVE_DEPTH = KRegCCD.DEFAULT_RESERVE_DEPTH;
    // mu grid capped at the reliability ceiling (MU_RELIABLE_MAX = 0.05).
    private static final double[] MU_GRID =
            {0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05};

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: NjVsMjRSV2 <trees-file> [burnin=0.1] [folds=5] [maxTrees=all]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;
        int folds = args.length > 2 ? Integer.parseInt(args[2]) : 5;
        int maxTrees = args.length > 3 ? Integer.parseInt(args[3]) : Integer.MAX_VALUE;

        List<Tree> trees = LoadOrStoreTrees.loadTrees(treeFile, burnin);
        if (trees.size() > maxTrees) {
            // even thinning to maxTrees (keep the chain spread, not a contiguous block)
            List<Tree> thinned = new ArrayList<>(maxTrees);
            double step = (double) trees.size() / maxTrees;
            for (int i = 0; i < maxTrees; i++) {
                thinned.add(trees.get((int) (i * step)));
            }
            trees = thinned;
        }
        int n = trees.size();
        System.out.printf("Loaded %d trees (burnin %.2f) from %s%n", n, burnin, treeFile.getName());
        System.out.printf("alpha=%.3f, reserve depth k=%d, folds=%d%n%n", ALPHA, RESERVE_DEPTH, folds);

        int b = Math.max(2, Math.min(folds, n));
        // total[mode][muIndex] accumulated over folds
        double[] sharedTotal = new double[MU_GRID.length];
        double[] flatTotal = new double[MU_GRID.length];
        int testCount = 0;

        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                (i % b == f ? test : train).add(trees.get(i)); // strided folds
            }
            testCount += test.size();

            long t0 = System.currentTimeMillis();
            KRegCCD shared = new KRegCCD(train, 0.0, KRegCCD.DEFAULT_MU, ALPHA,
                    RESERVE_DEPTH, KRegCCD.TailMode.NONE, KRegCCD.NovelMode.SHARED);
            KRegCCD flat = new KRegCCD(train, 0.0, KRegCCD.DEFAULT_MU, ALPHA,
                    RESERVE_DEPTH, KRegCCD.TailMode.NONE, KRegCCD.NovelMode.FLAT);

            for (Tree t : test) {
                for (int g = 0; g < MU_GRID.length; g++) {
                    sharedTotal[g] += shared.getLogProbabilityOfTree(t, MU_GRID[g]);
                    flatTotal[g] += flat.getLogProbabilityOfTree(t, MU_GRID[g]);
                }
            }
            System.out.printf("fold %d/%d done (train=%d, test=%d) in %.1fs%n",
                    f + 1, b, train.size(), test.size(),
                    (System.currentTimeMillis() - t0) / 1000.0);
        }

        System.out.printf("%n%-8s | %14s | %14s | %14s%n", "mu", "SHARED (N_j)", "FLAT (M_j)", "FLAT - SHARED");
        System.out.println("---------+----------------+----------------+----------------");
        double bestShared = Double.NEGATIVE_INFINITY, bestFlat = Double.NEGATIVE_INFINITY;
        double bestSharedMu = 0, bestFlatMu = 0;
        for (int g = 0; g < MU_GRID.length; g++) {
            System.out.printf("%-8.3f | %14.2f | %14.2f | %+14.2f%n",
                    MU_GRID[g], sharedTotal[g], flatTotal[g], flatTotal[g] - sharedTotal[g]);
            if (sharedTotal[g] > bestShared) {
                bestShared = sharedTotal[g];
                bestSharedMu = MU_GRID[g];
            }
            if (flatTotal[g] > bestFlat) {
                bestFlat = flatTotal[g];
                bestFlatMu = MU_GRID[g];
            }
        }

        System.out.printf("%nBest SHARED (N_j): total=%.2f at mu=%.3f  (%.4f nat/tree)%n",
                bestShared, bestSharedMu, bestShared / testCount);
        System.out.printf("Best FLAT   (M_j): total=%.2f at mu=%.3f  (%.4f nat/tree)%n",
                bestFlat, bestFlatMu, bestFlat / testCount);
        double diff = bestFlat - bestShared;
        System.out.printf("%n=> %s fits better by %.2f nat total (%.4f nat/tree over %d held-out trees)%n",
                diff > 0 ? "FLAT (M_j)" : "SHARED (N_j)", Math.abs(diff),
                Math.abs(diff) / testCount, testCount);
    }
}
