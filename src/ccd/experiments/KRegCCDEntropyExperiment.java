package ccd.experiments;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.algorithms.regularisation.KRegCCDParameterOptimiser;
import ccd.model.CCD1;
import ccd.model.KRegCCD;

import java.io.File;
import java.util.List;
import java.util.Random;

/**
 * Runs the two KRegCCD entropy estimators -- the generalised Lewis recursion
 * ({@link KRegCCD#getEntropyRecursive()}) and the Monte-Carlo plug-in
 * ({@link KRegCCD#getEntropyMonteCarlo(int)}) -- on a real posterior tree set, where the topology
 * space is far too large to enumerate. Reports both alongside a plain CCD1 (observed-support-only)
 * entropy, so the contribution of KRegCCD's novel-clade / escape support is visible.
 *
 * <p>The hyperparameters can be fixed (a numeric {@code mu}, {@code alpha = DEFAULT_ALPHA}) or
 * selected per dataset by maximising cross-validated held-out tree log-probability (pass
 * {@code mu = opt}). Because the entropy is a function of {@code mu} -- it grows from the
 * observed-only backbone entropy at {@code mu -> 0} as escape mass increases -- a fixed default is
 * arbitrary; {@code opt} reports the entropy of the predictively best-fit model instead.
 *
 * Usage: KRegCCDEntropyExperiment <trees-file> [burnin=0.1] [maxTrees=all] [mu=0.01|opt] [mcSamples=20000] [folds=5]
 */
public class KRegCCDEntropyExperiment {

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: KRegCCDEntropyExperiment <trees-file> [burnin=0.1] "
                    + "[maxTrees=all] [mu=0.01|opt] [mcSamples=20000] [folds=5]");
            System.exit(1);
        }
        File treeFile = new File(args[0]);
        double burnin = args.length > 1 ? Double.parseDouble(args[1]) : 0.1;
        int maxTrees = args.length > 2 ? Integer.parseInt(args[2]) : Integer.MAX_VALUE;
        String muArg = args.length > 3 ? args[3] : String.valueOf(KRegCCD.DEFAULT_MU);
        boolean optimise = muArg.equalsIgnoreCase("opt");
        double mu = optimise ? Double.NaN : Double.parseDouble(muArg);
        int mcSamples = args.length > 4 ? Integer.parseInt(args[4]) : 20_000;
        int folds = args.length > 5 ? Integer.parseInt(args[5]) : KRegCCDParameterOptimiser.DEFAULT_FOLDS;

        // When a cap is given, parse only the evenly-thinned subsample (skipping the rest in the
        // reader); otherwise load every post-burnin tree in a single pass.
        List<Tree> trees = (maxTrees == Integer.MAX_VALUE)
                ? LoadOrStoreTrees.loadTrees(treeFile, burnin)
                : LoadOrStoreTrees.loadTrees(treeFile, burnin, maxTrees);
        int n = trees.size();
        int nTaxa = trees.get(0).getLeafNodeCount();
        System.out.printf("%n=== %s ===%n", treeFile.getName());

        // Hyperparameters: either fixed, or selected per dataset by cross-validated held-out
        // tree log-probability. The entropy depends on mu (it grows with escape mass), so a fixed
        // default is arbitrary; the optimised pair reports the entropy of the best-fit model.
        double alpha = KRegCCD.DEFAULT_ALPHA;
        if (optimise) {
            System.out.printf("optimising (alpha, mu) by %d-fold held-out log-probability...%n", folds);
            long tOpt = System.currentTimeMillis();
            KRegCCDParameterOptimiser.Params p = KRegCCDParameterOptimiser.optimise(trees, folds);
            alpha = p.alpha();
            mu = p.mu();
            System.out.printf("selected alpha=%.5f, mu=%.5f (held-out logP=%.2f) [%.1fs]%n",
                    alpha, mu, p.heldOutLogProb(), secs(tOpt));
        }
        System.out.printf("taxa=%d, trees used=%d (burnin %.2f), mu=%.5f, alpha=%.3f, k=%d%n%n",
                nTaxa, n, burnin, mu, alpha, KRegCCD.DEFAULT_RESERVE_DEPTH);

        // observed-support-only baseline
        long t0 = System.currentTimeMillis();
        CCD1 ccd1 = new CCD1(trees, 0.0);
        double ccd1Entropy = ccd1.getEntropy();
        System.out.printf("CCD1 (observed-only) entropy : %12.4f nats   [%d clades, %.1fs]%n",
                ccd1Entropy, ccd1.getNumberOfClades(), secs(t0));

        // SELF-CONSISTENT target: tail = 0 (NONE), escape mass exactly mu, blue regions only up to
        // the reserve depth k. This is a proper NORMALISED distribution (red 1-mu + escape mu = 1)
        // -- the one SamplingFidelity.SELF_CONSISTENT samples and the recursion scores -- but it is
        // NOT the full-support model: it omits regions deeper than k (the reserve "tail", mass
        // ~mu^3 per reserving clade). Both estimators agree here because they share this target.
        t0 = System.currentTimeMillis();
        KRegCCD ccd = new KRegCCD(trees, 0.0, mu, alpha,
                KRegCCD.DEFAULT_RESERVE_DEPTH, KRegCCD.TailMode.NONE, KRegCCD.NovelMode.FLAT);
        System.out.printf("KRegCCD built                : %12s        [%d clades, %.1fs]%n",
                "", ccd.getNumberOfClades(), secs(t0));

        double pNovel = ccd.getNovelCladeProbability();
        System.out.printf("model P(tree has novel clade): %12.4f%n%n", pNovel);

        // 1) generalised Lewis recursion (deterministic, exact for the self-consistent distribution)
        t0 = System.currentTimeMillis();
        double recursion = ccd.getEntropyRecursive();
        double recSecs = secs(t0);

        // 2) Monte-Carlo plug-in (SELF_CONSISTENT sampler -> same target as the recursion)
        ccd.setRandom(new Random(20260614L));
        t0 = System.currentTimeMillis();
        double[] mc = ccd.getEntropyMonteCarlo(mcSamples);
        double mcSecs = secs(t0);

        System.out.println("KRegCCD self-consistent entropy (escape=mu, regions <= k, normalised):");
        System.out.printf("  recursion                  : %12.4f nats   [%.1fs]%n", recursion, recSecs);
        System.out.printf("  Monte-Carlo (%d samples)%s: %12.4f nats  +/- %.4f   [%.1fs]%n",
                mcSamples, mcSamples >= 100000 ? "" : "   ", mc[0], mc[1], mcSecs);
        System.out.printf("  MC - recursion             : %+12.4f nats  (%.2f standard errors)%n",
                mc[0] - recursion, mc[1] > 0 ? (mc[0] - recursion) / mc[1] : 0.0);

        // 3) FULL-SUPPORT cross-check: tail-included model (SAMPLED) + tail-sampling sampler
        // (FULL_SUPPORT draws the leading tail orders k+1, k+2). This estimates the entropy of the
        // actual full-support model -- every tree, all region depths -- so the gap to the
        // self-consistent number above is the reserve tail's contribution.
        t0 = System.currentTimeMillis();
        KRegCCD full = new KRegCCD(trees, 0.0, mu, alpha,
                KRegCCD.DEFAULT_RESERVE_DEPTH, KRegCCD.TailMode.SAMPLED, KRegCCD.NovelMode.FLAT);
        full.setSamplingFidelity(KRegCCD.SamplingFidelity.FULL_SUPPORT);
        full.setRandom(new Random(20260614L));
        double[] mcFull = full.getEntropyMonteCarlo(mcSamples);
        System.out.printf("%nKRegCCD full-support entropy (all region depths, tail-corrected):%n");
        System.out.printf("  Monte-Carlo (%d samples)%s: %12.4f nats  +/- %.4f   [%.1fs]%n",
                mcSamples, mcSamples >= 100000 ? "" : "   ", mcFull[0], mcFull[1], secs(t0));
        System.out.printf("  full-support - self-consistent: %+9.4f nats  (the reserve tail)%n",
                mcFull[0] - recursion);

        System.out.printf("%n  escape adds over CCD1 backbone: %+.4f nats (%.1f%%)%n",
                recursion - ccd1Entropy, 100.0 * (recursion - ccd1Entropy) / ccd1Entropy);
    }

    private static double secs(long t0) {
        return (System.currentTimeMillis() - t0) / 1000.0;
    }
}
