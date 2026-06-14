package ccd.experiments;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.KRegCCD;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Tracks how the held-out-optimal escape probability {@code mu*} changes with posterior sample
 * size. After burn-in, takes <em>disjoint, consecutive</em> blocks of growing size (128, 512, 2048,
 * 8192 by default) and, for each, selects {@code mu*} by maximising cross-validated held-out tree
 * log-probability over a log-spaced grid.
 *
 * <p>{@code alpha} (backbone smoothing) is held fixed: in a joint (alpha, mu) optimisation on this
 * data {@code mu*} was insensitive to {@code alpha} (the same grid point won for every {@code alpha}
 * Brent tried), so fixing it isolates the sample-size effect on {@code mu} and avoids the ~10x cost
 * of an inner alpha search. Each fold's backbone is built once and reused across the whole mu grid
 * (mu only re-solves {@code eps} from the cached counts), so the sweep is cheap.
 *
 * <p>Expected trend: more trees -> a richer observed backbone -> fewer genuinely novel clades in the
 * held-out folds -> less escape mass wanted -> {@code mu*} decreases.
 *
 * Usage: KRegCCDMuVsSampleSize <trees-file> [burnin=0.1] [folds=5]
 */
public class KRegCCDMuVsSampleSize {

    private static final int[] BLOCK_SIZES = {128, 512, 2048, 8192};
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int RESERVE_DEPTH = KRegCCD.DEFAULT_RESERVE_DEPTH;
    // mu grid spans well below the optimiser's usual floor (so a shrinking mu* stays interior) up to
    // the reliability ceiling MU_RELIABLE_MAX (above it the depth-k reserve approximation breaks).
    private static final double MU_LO = 1e-4;
    private static final double MU_HI = KRegCCD.MU_RELIABLE_MAX;
    private static final int MU_GRID = 49;

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: KRegCCDMuVsSampleSize <trees-file> [burnin=0.1] [folds=5]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;
        int folds = args.length > 2 ? Integer.parseInt(args[2]) : 5;

        List<Tree> trees = LoadOrStoreTrees.loadTrees(treeFile, burnin);
        int nTaxa = trees.get(0).getLeafNodeCount();
        double[] muGrid = logspace(MU_LO, MU_HI, MU_GRID);

        System.out.printf("%n=== mu* vs sample size: %s ===%n", treeFile.getName());
        System.out.printf("taxa=%d, available after burnin=%d, alpha=%.2f (fixed), k=%d, folds=%d%n",
                nTaxa, trees.size(), ALPHA, RESERVE_DEPTH, folds);
        System.out.printf("mu grid: %d log-spaced points in [%.0e, %.3f]%n%n", MU_GRID, MU_LO, MU_HI);

        System.out.printf("%-7s %-13s %-7s %-10s %-14s %-9s %-10s%n",
                "block", "range", "clades", "mu*", "heldOut nat/tr", "P(novel)", "H(mu*)");

        int start = 0;
        for (int size : BLOCK_SIZES) {
            if (start + size > trees.size()) {
                System.out.printf("(not enough trees for a block of %d starting at %d; stopping)%n",
                        size, start);
                break;
            }
            List<Tree> block = new ArrayList<>(trees.subList(start, start + size)); // disjoint, consecutive

            // strided cross-validation folds over the block
            int b = Math.max(2, Math.min(folds, size));
            List<List<Tree>> trainByFold = new ArrayList<>(b);
            List<List<Tree>> testByFold = new ArrayList<>(b);
            for (int f = 0; f < b; f++) {
                List<Tree> train = new ArrayList<>();
                List<Tree> test = new ArrayList<>();
                for (int i = 0; i < size; i++) {
                    (i % b == f ? test : train).add(block.get(i));
                }
                trainByFold.add(train);
                testByFold.add(test);
            }

            // build each fold's backbone once (NONE tail; scoring overrides mu)
            KRegCCD[] models = new KRegCCD[b];
            for (int f = 0; f < b; f++) {
                models[f] = new KRegCCD(trainByFold.get(f), 0.0, KRegCCD.DEFAULT_MU, ALPHA,
                        RESERVE_DEPTH, KRegCCD.TailMode.NONE, KRegCCD.NovelMode.FLAT);
            }
            int testCount = 0;
            for (List<Tree> test : testByFold) {
                testCount += test.size();
            }

            // sweep mu, reusing the backbones; pick the argmax of total held-out log-probability
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

            // model on the whole block at (alpha, mu*) for P(novel) and entropy
            KRegCCD whole = new KRegCCD(block, 0.0, bestMu, ALPHA,
                    RESERVE_DEPTH, KRegCCD.TailMode.NONE, KRegCCD.NovelMode.FLAT);
            double pNovel = whole.getNovelCladeProbability();
            double entropy = whole.getEntropyRecursive();

            String boundary = bestMu <= MU_LO * 1.0001 ? " (=grid floor)"
                    : bestMu >= MU_HI * 0.9999 ? " (=ceiling)" : "";
            System.out.printf("%-7d %-13s %-7d %-10.5f %-14.4f %-9.4f %-10.4f%s%n",
                    size, (start + 1) + "-" + (start + size), whole.getNumberOfClades(),
                    bestMu, bestTotal / testCount, pNovel, entropy, boundary);

            start += size;
        }
    }

    private static double[] logspace(double lo, double hi, int n) {
        double[] grid = new double[n];
        double logLo = Math.log(lo);
        double logHi = Math.log(hi);
        for (int i = 0; i < n; i++) {
            grid[i] = Math.exp(logLo + (logHi - logLo) * i / (n - 1));
        }
        return grid;
    }
}
