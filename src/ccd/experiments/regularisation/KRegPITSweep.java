package ccd.experiments.regularisation;

import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import ccd.model.ITreeDistribution;
import ccd.model.KRegCCD;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Random;

/**
 * Sweeps the {@link KRegPITExperiment} PIT calibration test over several models, sample sizes and
 * replicates, pooling the randomized PIT across replicates into one histogram per (model, size).
 *
 * <p>Pooling is statistically sound: each test tree's PIT is computed under <em>its own</em>
 * replicate's fitted model, and if every replicate is calibrated then each PIT is Uniform(0,1), so
 * their mixture is Uniform(0,1) too. Pooling simply gives a smoother, higher-power calibration
 * check than any single replicate (the worst single-replicate noise is at the smallest sample
 * size). The chi-square and KS tests are run on the pooled values.</p>
 *
 * <p>Expects the figshare Yule-style layout: {@code <dataRoot>/rep<r>/run1} (training) and
 * {@code .../run2} (testing), each containing one {@code <prefix>*.trees} posterior file. The trees
 * file is found by globbing, not by reconstructing {@code <prefix>-<r>.trees}, so a mislabelled rep
 * number does not silently drop the replicate (e.g. rep100's run trees are named {@code yule-n50-0.trees}).
 * The single combined tidy TSV it writes feeds straight into {@code doc/plot_pit.py}.</p>
 *
 * <p>Usage: {@code KRegPITSweep <dataRoot> [prefix=yule-n50] [startRep=1] [endRep=10]
 * [sizes=300,1000,3000] [models=kreg,ccd1,ccd0] [numBins=20] [numSamples=100000] [seed=42]
 * [out=pit_sweep_hist.tsv]}</p>
 */
public class KRegPITSweep {

    public static void main(String[] args) throws IOException {
        if (args.length < 1) {
            System.err.println("Usage: KRegPITSweep <dataRoot> [prefix=yule-n50] [startRep=1] [endRep=10] "
                    + "[sizes=300,1000,3000] [models=kreg,ccd1,ccd0] [numBins=20] [numSamples=100000] "
                    + "[seed=42] [out=pit_sweep_hist.tsv]");
            return;
        }

        String dataRoot = args[0];
        String prefix = args.length > 1 ? args[1] : "yule-n50";
        int startRep = args.length > 2 ? Integer.parseInt(args[2]) : 1;
        int endRep = args.length > 3 ? Integer.parseInt(args[3]) : 10;
        int[] sizes = parseInts(args.length > 4 ? args[4] : "300,1000,3000");
        String[] models = (args.length > 5 ? args[5] : "kreg,ccd1,ccd0").split(",");
        int numBins = args.length > 6 ? Integer.parseInt(args[6]) : KRegPITExperiment.DEFAULT_NUM_BINS;
        int numSamples = args.length > 7 ? Integer.parseInt(args[7]) : KRegPITExperiment.DEFAULT_NUM_SAMPLES;
        long seed = args.length > 8 ? Long.parseLong(args[8]) : KRegPITExperiment.DEFAULT_SEED;
        String out = args.length > 9 ? args[9] : "pit_sweep_hist.tsv";
        String paramsOut = (out.endsWith(".tsv") ? out.substring(0, out.length() - 4) : out) + ".params.tsv";

        Randomizer.setSeed(seed); // make the model sampling (sampleTreeLogProbability, CV) reproducible
        Random rng = new Random(seed);
        long t0 = System.currentTimeMillis();

        try (PrintWriter w = new PrintWriter(new File(out));
             PrintWriter pw = new PrintWriter(new File(paramsOut))) {
            KRegPITExperiment.writeTidyHeader(w);
            pw.println("rep\tsampleSize\talpha\tmu"); // one row per KRegCCD fit (the only model with fitted params)

            for (String model : models) {
                for (int size : sizes) {
                    List<Double> uPool = new ArrayList<>();
                    List<Double> vPool = new ArrayList<>();
                    int zeroedPool = 0;

                    for (int rep = startRep; rep <= endRep; rep++) {
                        File trainDir = new File(dataRoot + "/rep" + rep + "/run1");
                        File testDir = new File(dataRoot + "/rep" + rep + "/run2");
                        File train = findTreesFile(trainDir, prefix);
                        File test = findTreesFile(testDir, prefix);
                        if (train == null || test == null) {
                            System.err.printf(Locale.US, "skip rep %d: no %s*.trees in %s%n", rep, prefix,
                                    train == null ? trainDir : testDir);
                            continue;
                        }

                        List<Tree> trainTrees = KRegPITExperiment.loadTrees(train, size);
                        List<Tree> testTrees = KRegPITExperiment.loadTrees(test, size);
                        ITreeDistribution m = KRegPITExperiment.buildModel(model, trainTrees);
                        if (m instanceof KRegCCD kreg) { // only KRegCCD has fitted (alpha, mu)
                            pw.printf(Locale.US, "%d\t%d\t%.6f\t%.6f%n", rep, size, kreg.getAlpha(), kreg.getMu());
                            pw.flush();
                        }
                        double[] sortedLogP = KRegPITExperiment.sortedSampleLogProbs(m, numSamples);
                        int z = KRegPITExperiment.accumulatePIT(m, testTrees, sortedLogP, rng, uPool, vPool);
                        zeroedPool += z;

                        System.out.printf(Locale.US, "%s size=%d rep=%d: scored %d, zeroed %d  (pool=%d, %.0fs)%n",
                                model, size, rep, testTrees.size() - z, z, vPool.size(),
                                (System.currentTimeMillis() - t0) / 1000.0);
                    }

                    if (vPool.isEmpty()) {
                        System.err.printf(Locale.US, "no scored trees for %s size=%d%n", model, size);
                        continue;
                    }

                    long[] observed = KRegPITExperiment.histogram(vPool, numBins);
                    double[] stats = KRegPITExperiment.uniformityTests(vPool, observed);
                    KRegPITExperiment.writeTidyGroup(w, model, size, observed,
                            stats[0], stats[1], stats[2], stats[3], stats[4], vPool.size(), zeroedPool);
                    w.flush();

                    System.out.printf(Locale.US, "=> POOLED %s size=%d: m=%d L1=%.3f chiP=%.3g ksP=%.3g zeroed=%d%n",
                            model, size, vPool.size(), stats[4], stats[1], stats[3], zeroedPool);
                }
            }
        }
        System.out.printf(Locale.US, "wrote %s and %s  (%.0fs total)%n",
                out, paramsOut, (System.currentTimeMillis() - t0) / 1000.0);
    }

    /**
     * Finds the single {@code <prefix>*.trees} posterior file in a run directory. We glob rather than
     * reconstruct {@code <prefix>-<rep>.trees} because the data pipeline mislabels some reps (e.g.
     * rep100's run trees are named yule-n50-0.trees, not yule-n50-100.trees); each run dir holds
     * exactly one matching trees file. Returns null if the directory has none.
     */
    private static File findTreesFile(File runDir, String prefix) {
        File[] hits = runDir.listFiles((dir, name) -> name.startsWith(prefix) && name.endsWith(".trees"));
        if (hits == null || hits.length == 0) return null;
        if (hits.length > 1) {
            Arrays.sort(hits); // deterministic pick, but flag the ambiguity
            System.err.printf(Locale.US, "warning: %d %s*.trees files in %s; using %s%n",
                    hits.length, prefix, runDir, hits[0].getName());
        }
        return hits[0];
    }

    private static int[] parseInts(String csv) {
        String[] parts = csv.split(",");
        int[] r = new int[parts.length];
        for (int i = 0; i < parts.length; i++) r[i] = Integer.parseInt(parts[i].trim());
        return r;
    }
}
