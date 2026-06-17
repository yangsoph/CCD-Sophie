package ccd.experiments;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.KRegCCD;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Tests whether the model over-reserves novel-clade mass: compares the model's closed-form
 * probability that a tree contains a novel clade, {@code P(novel)}, against the empirical held-out
 * rate (fraction of held-out trees that actually contain a clade absent from the training backbone).
 * If {@code P(novel) > } held-out rate, the held-out-optimal mu reserves more total novel mass than
 * the data needs --- one component of the model being flatter (higher entropy) than the posterior.
 *
 * Reproduces the 8192-tree block of the mu-vs-sample-size experiment (offset 128+512+2048, size
 * 8192). Cheap: no entropy recursion, just the closed-form P(novel) and a novelty count.
 *
 * Usage: KRegCCDNovelCalibration <trees-file> [burnin=0.1] [mu=0.00061] [folds=5]
 */
public class KRegCCDNovelCalibration {

    private static final int OFFSET = 128 + 512 + 2048; // = 2688
    private static final int SIZE = 8192;
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: KRegCCDNovelCalibration <trees-file> [burnin=0.1] [mu=0.00061] [folds=5]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;
        double mu = args.length > 2 ? Double.parseDouble(args[2]) : 0.00061;
        int folds = args.length > 3 ? Integer.parseInt(args[3]) : 5;

        List<Tree> trees = LoadOrStoreTrees.loadTrees(treeFile, burnin);
        List<Tree> block = new ArrayList<>(trees.subList(OFFSET, OFFSET + SIZE));
        int nTaxa = block.get(0).getLeafNodeCount();
        int b = Math.max(2, Math.min(folds, SIZE));

        System.out.printf("%n=== novel-clade calibration on the %d-tree block (trees %d-%d): %s ===%n",
                SIZE, OFFSET + 1, OFFSET + SIZE, treeFile.getName());
        System.out.printf("taxa=%d, mu=%.5f, alpha=%.2f, k=%d, folds=%d%n%n", nTaxa, mu, ALPHA, K, b);

        System.out.printf("%-6s %-8s %-8s %-14s %-16s%n",
                "fold", "train", "test", "model P(novel)", "held-out novel");
        double modelSum = 0.0;
        int novelCount = 0;
        int testCount = 0;
        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < SIZE; i++) {
                (i % b == f ? test : train).add(block.get(i));
            }
            KRegCCD model = new KRegCCD(train, 0.0, mu, ALPHA, K, KRegCCD.TailMode.NONE,
                    KRegCCD.NovelMode.FLAT);
            double modelPnovel = model.getNovelCladeProbability();
            int foldNovel = 0;
            for (Tree t : test) {
                if (model.containsNovelClade(t)) {
                    foldNovel++;
                }
            }
            modelSum += modelPnovel;
            novelCount += foldNovel;
            testCount += test.size();
            System.out.printf("%-6d %-8d %-8d %-14.5f %-16.5f%n",
                    f, train.size(), test.size(), modelPnovel, foldNovel / (double) test.size());
        }

        double modelAvg = modelSum / b;
        double heldOutRate = novelCount / (double) testCount;
        double se = Math.sqrt(heldOutRate * (1 - heldOutRate) / testCount);
        System.out.printf("%n%-26s %.5f%n", "model P(novel) (mean):", modelAvg);
        System.out.printf("%-26s %.5f  (+/- %.5f)%n", "held-out novel rate:", heldOutRate, se);
        System.out.printf("%-26s %+.5f  (%.1fx)%n", "model - held-out:",
                modelAvg - heldOutRate, modelAvg / heldOutRate);
    }
}
