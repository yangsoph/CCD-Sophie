package ccd.experiments;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.KRegCCD;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Dumps the (alpha, mu) parameter surface for KRegCCD: the cross-validated held-out tree
 * log-probability -- the very objective {@link ccd.algorithms.regularisation.KRegCCDParameterOptimiser}
 * maximises -- evaluated on a 2D log-spaced grid so the ridge can be plotted. The tree set is thinned
 * evenly across the post-burnin chain to a fixed size (default 1024) so the surface is comparable
 * across datasets and cheap enough for a dense grid.
 *
 * <p>For each alpha the per-fold backbones are built once and the whole mu grid is scored by reusing
 * them (mu only re-solves eps from the cached counts), so the cost is one backbone rebuild per
 * (alpha, fold), not per grid point -- the same trick the optimiser uses. Tail correction is NONE,
 * matching the optimiser's objective; mu is capped at {@link KRegCCD#MU_RELIABLE_MAX} where the
 * depth-k reserve is still reliable.
 *
 * <p>Writes a TSV (alpha, mu, heldOutTotal, heldOutPerTree) to the output file, one row per grid
 * point, for plotting; the argmax is printed to stdout.
 *
 * Usage: KRegCCDParameterSurface &lt;trees-file&gt; [burnin=0.1] [thinTo=1024] [folds=5]
 *        [alphaGrid=31] [muGrid=31] [out=surface.tsv]
 */
public class KRegCCDParameterSurface {

    private static final double ALPHA_LO = 1e-2;
    private static final double ALPHA_HI = 2.0;
    private static final double MU_LO = 1e-3;
    private static final double MU_HI = KRegCCD.MU_RELIABLE_MAX;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: KRegCCDParameterSurface <trees-file> [burnin=0.1] [thinTo=1024] "
                    + "[folds=5] [alphaGrid=31] [muGrid=31] [out=surface.tsv]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;
        int thinTo = args.length > 2 ? Integer.parseInt(args[2]) : 1024;
        int folds = args.length > 3 ? Integer.parseInt(args[3]) : 5;
        int nAlpha = args.length > 4 ? Integer.parseInt(args[4]) : 31;
        int nMu = args.length > 5 ? Integer.parseInt(args[5]) : 31;
        String out = args.length > 6 ? args[6] : "surface.tsv";

        // Parse only the thinned subsample; the trees in between are skipped in the reader.
        List<Tree> trees = LoadOrStoreTrees.loadTrees(treeFile, burnin, thinTo);
        int n = trees.size();
        int nTaxa = trees.get(0).getLeafNodeCount();
        int b = Math.max(2, Math.min(folds, n));

        // Fixed strided fold partition (only the backbones change with alpha, not the split).
        List<List<Tree>> trainByFold = new ArrayList<>(b);
        List<List<Tree>> testByFold = new ArrayList<>(b);
        int testCount = 0;
        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                (i % b == f ? test : train).add(trees.get(i));
            }
            trainByFold.add(train);
            testByFold.add(test);
            testCount += test.size();
        }

        double[] alphaGrid = logspace(ALPHA_LO, ALPHA_HI, nAlpha);
        double[] muGrid = logspace(MU_LO, MU_HI, nMu);

        System.out.printf("%n=== (alpha, mu) held-out surface: %s ===%n", treeFile.getName());
        System.out.printf("taxa=%d, thinned to %d (burnin %.2f), folds=%d, held-out trees=%d, k=%d%n",
                nTaxa, n, burnin, b, testCount, K);
        System.out.printf("alpha in [%.3g, %.3g] x %d, mu in [%.3g, %.3g] x %d -> %d grid points%n%n",
                ALPHA_LO, ALPHA_HI, nAlpha, MU_LO, MU_HI, nMu, nAlpha * nMu);

        PrintWriter w = new PrintWriter(out);
        w.println("alpha\tmu\theldout_total\theldout_per_tree");

        double bestAlpha = alphaGrid[0], bestMu = muGrid[0], bestTotal = Double.NEGATIVE_INFINITY;
        long t0 = System.currentTimeMillis();
        for (int ai = 0; ai < nAlpha; ai++) {
            double alpha = alphaGrid[ai];
            // Build the per-fold backbones once for this alpha; mu only re-solves eps when scoring.
            KRegCCD[] models = new KRegCCD[b];
            for (int f = 0; f < b; f++) {
                models[f] = new KRegCCD(trainByFold.get(f), 0.0, KRegCCD.DEFAULT_MU, alpha,
                        K, KRegCCD.TailMode.NONE);
            }
            double rowBestTotal = Double.NEGATIVE_INFINITY, rowBestMu = muGrid[0];
            for (double mu : muGrid) {
                double total = 0.0;
                for (int f = 0; f < b; f++) {
                    for (Tree t : testByFold.get(f)) {
                        total += models[f].getLogProbabilityOfTree(t, mu);
                    }
                }
                w.printf("%.6g\t%.6g\t%.6f\t%.6f%n", alpha, mu, total, total / testCount);
                if (total > rowBestTotal) { rowBestTotal = total; rowBestMu = mu; }
                if (total > bestTotal) { bestTotal = total; bestAlpha = alpha; bestMu = mu; }
            }
            System.out.printf("alpha=%.5f : best mu=%.5f, held-out=%.4f nat/tree   [%.0fs]%n",
                    alpha, rowBestMu, rowBestTotal / testCount, (System.currentTimeMillis() - t0) / 1000.0);
        }
        w.close();

        System.out.printf("%nargmax: alpha=%.5f, mu=%.5f, held-out=%.4f nat/tree (%.2f total)%n",
                bestAlpha, bestMu, bestTotal / testCount, bestTotal);
        if (bestMu >= MU_HI * (1 - 1e-9)) {
            System.err.println("WARNING: held-out logP still rising in mu at the reliability ceiling "
                    + MU_HI + "; treat as a capped optimum.");
        }
        System.out.printf("surface written to %s (%d rows)%n", out, nAlpha * nMu);
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
