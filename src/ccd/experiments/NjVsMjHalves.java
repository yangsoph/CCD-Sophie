package ccd.experiments;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.KRegCCD;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Temporal-split held-out comparison of KRegCCD SHARED (N_j) vs FLAT (M_j): train on the FIRST
 * half of an MCMC tree file, test on the SECOND half. Unlike strided CV (which interleaves the
 * chain, so test trees are near-duplicates of train trees and carry mostly single, dispersed
 * novel clades), chain drift across the halves puts genuinely CLUSTERED novelty into the test set
 * -- exactly the regime where FLAT (locally Boltzmann, no clustering penalty) and SHARED (clustering
 * penalised via /pathcount) diverge.
 *
 * Reports a novel-clade-count histogram of the test set, and the FLAT-SHARED held-out log-prob
 * margin split by novel-clade bucket, so the difference can be read off where it actually lives.
 *
 * Usage: NjVsMjHalves <trees-file> [burnin=0.1] [trainSize=n/2] [testTail=rest]
 *   trainSize : number of trees from the START (post-burnin) used to train.
 *   testTail  : number of trees from the END used to test (a gap between them widens the
 *               train/test separation, surfacing more novel/clustered structure).
 */
public class NjVsMjHalves {

    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int RESERVE_DEPTH = KRegCCD.DEFAULT_RESERVE_DEPTH;
    private static final double[] MU_GRID =
            {0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05};

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: NjVsMjHalves <trees-file> [burnin=0.1] [maxTrees=all]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;

        List<Tree> trees = LoadOrStoreTrees.loadTrees(treeFile, burnin);
        int n = trees.size();
        int trainSize = args.length > 2 ? Integer.parseInt(args[2]) : n / 2;
        int testTail = args.length > 3 ? Integer.parseInt(args[3]) : (n - trainSize);
        trainSize = Math.min(trainSize, n);
        testTail = Math.min(testTail, n);
        List<Tree> train = trees.subList(0, trainSize);
        List<Tree> test = trees.subList(n - testTail, n);
        int gap = (n - testTail) - trainSize;
        System.out.printf("Loaded %d trees (burnin %.2f) from %s%n", n, burnin, treeFile.getName());
        System.out.printf("TEMPORAL split: train=first %d, test=last %d, gap=%d; alpha=%.3f, k=%d%n",
                train.size(), test.size(), gap, ALPHA, RESERVE_DEPTH);
        if (gap < 0) {
            System.out.printf("  (WARNING: train and test overlap by %d trees)%n", -gap);
        }
        System.out.println();

        long t0 = System.currentTimeMillis();
        KRegCCD shared = new KRegCCD(new ArrayList<>(train), 0.0, KRegCCD.DEFAULT_MU, ALPHA,
                RESERVE_DEPTH, KRegCCD.TailMode.NONE, KRegCCD.NovelMode.SHARED);
        KRegCCD flat = new KRegCCD(new ArrayList<>(train), 0.0, KRegCCD.DEFAULT_MU, ALPHA,
                RESERVE_DEPTH, KRegCCD.TailMode.NONE, KRegCCD.NovelMode.FLAT);
        System.out.printf("built both models in %.1fs%n", (System.currentTimeMillis() - t0) / 1000.0);

        // Novel-clade-count bucket per test tree (mode-independent: same backbone).
        // bucket 0 -> 0 novel clades, 1 -> 1, 2 -> >=2 (the clustered/multi regime).
        int[] bucketOf = new int[test.size()];
        int[] hist = new int[8]; // raw counts 0..6, 7 = 7+
        int multi = 0;
        for (int i = 0; i < test.size(); i++) {
            int nc = shared.novelCladeCount(test.get(i));
            bucketOf[i] = nc == 0 ? 0 : (nc == 1 ? 1 : 2);
            hist[Math.min(nc, 7)]++;
            if (nc >= 2) {
                multi++;
            }
        }
        System.out.println("\nnovel-clade-count histogram over test trees:");
        for (int k = 0; k < hist.length; k++) {
            if (hist[k] > 0) {
                System.out.printf("  %s%d novel : %5d trees (%.1f%%)%n",
                        k == 7 ? ">=" : "", k, hist[k], 100.0 * hist[k] / test.size());
            }
        }
        System.out.printf("  -> %d/%d test trees (%.1f%%) have >=2 novel clades (clustered regime)%n",
                multi, test.size(), 100.0 * multi / test.size());

        // Per-mu totals, overall and split by bucket.
        System.out.printf("%n%-8s | %14s | %14s | %14s%n", "mu", "SHARED (N_j)", "FLAT (M_j)", "FLAT - SHARED");
        System.out.println("---------+----------------+----------------+----------------");
        double bestShared = Double.NEGATIVE_INFINITY, bestFlat = Double.NEGATIVE_INFINITY;
        double bestSharedMu = 0, bestFlatMu = 0;
        double[][] bucketDiffAtBestFlat = new double[3][1]; // filled at FLAT's best mu
        for (double mu : MU_GRID) {
            double sTot = 0, fTot = 0;
            double[] sB = new double[3], fB = new double[3];
            for (int i = 0; i < test.size(); i++) {
                double s = shared.getLogProbabilityOfTree(test.get(i), mu);
                double f = flat.getLogProbabilityOfTree(test.get(i), mu);
                sTot += s; fTot += f;
                sB[bucketOf[i]] += s; fB[bucketOf[i]] += f;
            }
            System.out.printf("%-8.4f | %14.2f | %14.2f | %+14.2f%n", mu, sTot, fTot, fTot - sTot);
            if (sTot > bestShared) { bestShared = sTot; bestSharedMu = mu; }
            if (fTot > bestFlat) {
                bestFlat = fTot; bestFlatMu = mu;
                for (int b = 0; b < 3; b++) bucketDiffAtBestFlat[b][0] = fB[b] - sB[b];
            }
        }

        System.out.printf("%nBest SHARED (N_j): total=%.2f at mu=%.4f  (%.4f nat/tree)%n",
                bestShared, bestSharedMu, bestShared / test.size());
        System.out.printf("Best FLAT   (M_j): total=%.2f at mu=%.4f  (%.4f nat/tree)%n",
                bestFlat, bestFlatMu, bestFlat / test.size());
        double diff = bestFlat - bestShared;
        System.out.printf("=> %s fits better by %.2f nat total (%.4f nat/tree over %d test trees)%n",
                diff > 0 ? "FLAT (M_j)" : "SHARED (N_j)", Math.abs(diff),
                Math.abs(diff) / test.size(), test.size());

        System.out.printf("%n(FLAT - SHARED) at mu=%.4f, split by novel-clade bucket:%n", bestFlatMu);
        String[] labels = {"0 novel ", "1 novel ", ">=2 novel"};
        for (int b = 0; b < 3; b++) {
            System.out.printf("  %s : %+10.2f nat  over %d trees%n",
                    labels[b], bucketDiffAtBestFlat[b][0], hist1(hist, b));
        }
    }

    private static int hist1(int[] hist, int bucket) {
        if (bucket == 0) return hist[0];
        if (bucket == 1) return hist[1];
        int s = 0; for (int k = 2; k < hist.length; k++) s += hist[k]; return s;
    }
}
