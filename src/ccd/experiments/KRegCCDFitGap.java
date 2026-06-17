package ccd.experiments;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.KRegCCD;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Isolates how much of the held-out-vs-entropy gap is genuine misfit rather than a comparison
 * artifact. The held-out average surprisal is the cross-entropy {@code H(Q,P) = H(Q) + KL(Q||P)} of
 * the true posterior Q against the fitted model P; comparing it to the entropy of a model fitted on
 * the <em>whole</em> block mixes in a train-vs-whole size difference. Here we compare it instead to
 * the self-entropy of the <em>same</em> cross-validation training models, removing that confound:
 * <pre>
 *   same-model gap = (-held-out) - mean_f H(P_train_f)   (cross-entropy vs the same models' entropy)
 *   confounded gap = (-held-out) - H(P_whole)            (the originally reported gap)
 * </pre>
 * Reproduces the 8192-tree block of the mu-vs-sample-size experiment (skip the 128+512+2048 earlier
 * blocks, take the next 8192) so the numbers line up with the reported -0.38 nat gap.
 *
 * Usage: KRegCCDFitGap <trees-file> [burnin=0.1] [mu=0.00061] [folds=5]
 */
public class KRegCCDFitGap {

    private static final int OFFSET = 128 + 512 + 2048; // = 2688
    private static final int SIZE = 8192;
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: KRegCCDFitGap <trees-file> [burnin=0.1] [mu=0.00061] [folds=5]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;
        double mu = args.length > 2 ? Double.parseDouble(args[2]) : 0.00061; // 8192-block mu*
        int folds = args.length > 3 ? Integer.parseInt(args[3]) : 5;

        List<Tree> trees = LoadOrStoreTrees.loadTrees(treeFile, burnin);
        if (OFFSET + SIZE > trees.size()) {
            System.err.printf("need %d trees after burnin, have %d%n", OFFSET + SIZE, trees.size());
            System.exit(1);
        }
        List<Tree> block = new ArrayList<>(trees.subList(OFFSET, OFFSET + SIZE));
        int nTaxa = block.get(0).getLeafNodeCount();
        int b = Math.max(2, Math.min(folds, SIZE));

        System.out.printf("%n=== fit gap on the %d-tree block (trees %d-%d): %s ===%n",
                SIZE, OFFSET + 1, OFFSET + SIZE, treeFile.getName());
        System.out.printf("taxa=%d, mu=%.5f, alpha=%.2f, k=%d, folds=%d%n%n", nTaxa, mu, ALPHA, K, b);

        double crossTotal = 0.0;
        int testCount = 0;
        double selfSum = 0.0;
        System.out.printf("%-6s %-8s %-8s %-12s%n", "fold", "train", "test", "H(P_train)");
        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < SIZE; i++) {
                (i % b == f ? test : train).add(block.get(i));
            }
            KRegCCD model = new KRegCCD(train, 0.0, mu, ALPHA, K, KRegCCD.TailMode.NONE,
                    KRegCCD.NovelMode.FLAT);
            for (Tree t : test) {
                crossTotal += model.getLogProbabilityOfTree(t);
            }
            testCount += test.size();
            double selfEnt = model.getEntropyRecursive();
            selfSum += selfEnt;
            System.out.printf("%-6d %-8d %-8d %-12.4f%n", f, train.size(), test.size(), selfEnt);
        }

        double crossEntropyPerTree = -crossTotal / testCount; // = -(held-out nat/tree)
        double selfEntropyAvg = selfSum / b;

        KRegCCD whole = new KRegCCD(block, 0.0, mu, ALPHA, K, KRegCCD.TailMode.NONE,
                KRegCCD.NovelMode.FLAT);
        double wholeEnt = whole.getEntropyRecursive();

        System.out.printf("%n%-34s %12.4f%n", "held-out cross-entropy (-held-out):", crossEntropyPerTree);
        System.out.printf("%-34s %12.4f%n", "fold-model self-entropy (mean):", selfEntropyAvg);
        System.out.printf("%-34s %12.4f%n", "whole-block self-entropy:", wholeEnt);
        System.out.printf("%n%-34s %+12.4f%n", "same-model gap (cross - foldH):",
                crossEntropyPerTree - selfEntropyAvg);
        System.out.printf("%-34s %+12.4f%n", "confounded gap (cross - wholeH):",
                crossEntropyPerTree - wholeEnt);
    }
}
