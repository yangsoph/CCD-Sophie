package ccd.algorithms.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.regularisation.NNIHeldOutComparison.FoldAssignment;
import ccd.model.KRegCCD;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import java.util.ArrayList;
import java.util.List;

/**
 * Selects {@link KRegCCD}'s hyperparameters by maximising <em>cross-validated held-out</em>
 * tree log-probability — the counterpart of {@link RegCCDParameterOptimiser} (which tunes the
 * backbone smoothing {@code alpha} alone) for the full-support model's two parameters.
 *
 * <p>
 * The objective is the total log-probability of held-out trees under batched leave-one-tree-out
 * (B-fold cross-validation): each fold is held out in turn while the backbone is rebuilt on the
 * rest and scores the fold. Scoring training trees instead would drive {@code mu -> 0} (the
 * backbone fits them exactly, so any reserved escape mass is pure loss), so the trees scored must
 * be out of training.
 *
 * <p>
 * The two parameters are optimised jointly but with different search strategies that exploit their
 * structure:
 * <ul>
 *   <li>{@code alpha} (backbone smoothing) changes the clade/split graph, so each value needs a
 *       full rebuild per fold; it has a sharp interior optimum and is found by Brent's method
 *       (as in {@link RegCCDParameterOptimiser}).</li>
 *   <li>{@code mu} (escape probability) does not touch the counts — only {@code eps(C)} re-solves
 *       from the cached {@code N_j} — so for a fixed {@code alpha} the whole {@code mu} sweep is
 *       evaluated on one set of trained backbones via
 *       {@link KRegCCD#getLogProbabilityOfTree(Tree, double)}. Its held-out curve is shallow, so
 *       it is searched on a (log-spaced) grid and taken at the argmax rather than by a
 *       derivative-based method that would wander on the flat ridge.</li>
 * </ul>
 * For each candidate {@code alpha}, {@code mu} is profiled out (best over the grid); Brent then
 * maximises that profiled objective over {@code alpha}.
 *
 * @author Claude (CCD-Sophie)
 */
public class KRegCCDParameterOptimiser {

    /** The selected hyperparameters and the cross-validated held-out log-probability at them. */
    public record Params(double alpha, double mu, double heldOutLogProb) {
    }

    /** Default number of cross-validation folds. */
    public static final int DEFAULT_FOLDS = 5;

    /** Reserve depth used for the candidate models (matches {@link KRegCCD#DEFAULT_RESERVE_DEPTH}). */
    private static final int RESERVE_DEPTH = KRegCCD.DEFAULT_RESERVE_DEPTH;

    /* alpha search interval (strictly positive; mirrors RegCCDParameterOptimiser's range). */
    private static final double ALPHA_LO = 1e-3;
    private static final double ALPHA_HI = 2.0;
    private static final double ALPHA_START = KRegCCD.DEFAULT_ALPHA;

    /* mu grid: log-spaced over a range spanning negligible to heavy escape. */
    private static final double MU_LO = 1e-3;
    private static final double MU_HI = 0.5;
    private static final int MU_GRID = 41;

    private static final int MAX_ALPHA_EVALS = 60;

    public static Params optimise(List<Tree> trees) {
        return optimise(trees, DEFAULT_FOLDS);
    }

    public static Params optimise(List<Tree> trees, int folds) {
        return optimise(trees, folds, FoldAssignment.STRIDED);
    }

    public static Params optimise(List<Tree> trees, int folds, FoldAssignment assignment) {
        int n = trees.size();
        if (n < 2) {
            throw new IllegalArgumentException("need at least 2 trees to cross-validate, got " + n);
        }
        int b = Math.max(2, Math.min(folds, n));

        // Fixed fold partition (only the backbones change with alpha, not the split).
        List<List<Tree>> trainByFold = new ArrayList<>(b);
        List<List<Tree>> testByFold = new ArrayList<>(b);
        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                boolean isTest = switch (assignment) {
                    case STRIDED -> (i % b == f);
                    case CONTIGUOUS -> (i >= (int) ((long) f * n / b)) && (i < (int) ((long) (f + 1) * n / b));
                };
                (isTest ? test : train).add(trees.get(i));
            }
            trainByFold.add(train);
            testByFold.add(test);
        }

        double[] muGrid = logspace(MU_LO, MU_HI, MU_GRID);

        UnivariateFunction profiled = alpha -> {
            double[] best = bestMuLogProb(trainByFold, testByFold, alpha, muGrid);
            System.out.println(String.format(
                    "  testing alpha = %.5f -> best mu = %.5f, held-out logP = %.5f",
                    alpha, best[0], best[1]));
            return best[1];
        };

        BrentOptimizer optimiser = new BrentOptimizer(1e-4, 1e-4);
        UnivariatePointValuePair sol = optimiser.optimize(
                new MaxEval(MAX_ALPHA_EVALS),
                new UnivariateObjectiveFunction(profiled),
                GoalType.MAXIMIZE,
                new SearchInterval(ALPHA_LO, ALPHA_HI, ALPHA_START));

        double alphaStar = sol.getPoint();
        double[] best = bestMuLogProb(trainByFold, testByFold, alphaStar, muGrid);
        System.out.println(String.format(
                "KRegCCD optimised: alpha = %.5f, mu = %.5f, held-out logP = %.5f",
                alphaStar, best[0], best[1]));
        return new Params(alphaStar, best[0], best[1]);
    }

    /**
     * The cross-validated held-out total log-probability at a fixed {@code (alpha, mu)} — the
     * objective the optimiser maximises. Rebuilds the per-fold backbones, so prefer the internal
     * grid path for a full search; this is for evaluating or comparing individual points.
     */
    public static double crossValidatedLogProb(List<Tree> trees, int folds, FoldAssignment assignment,
                                               double alpha, double mu) {
        int n = trees.size();
        int b = Math.max(2, Math.min(folds, n));
        double total = 0.0;
        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                boolean isTest = switch (assignment) {
                    case STRIDED -> (i % b == f);
                    case CONTIGUOUS -> (i >= (int) ((long) f * n / b)) && (i < (int) ((long) (f + 1) * n / b));
                };
                (isTest ? test : train).add(trees.get(i));
            }
            KRegCCD model = new KRegCCD(train, 0.0, KRegCCD.DEFAULT_MU, alpha,
                    RESERVE_DEPTH, KRegCCD.TailMode.NONE);
            for (Tree t : test) {
                total += model.getLogProbabilityOfTree(t, mu);
            }
        }
        return total;
    }

    /**
     * For a fixed {@code alpha}: builds one NONE-tail backbone per fold and returns
     * {@code {argmax mu over the grid, that held-out total log-probability}}. All grid values reuse
     * the same backbones — {@code mu} only re-solves {@code eps} from the cached counts.
     */
    private static double[] bestMuLogProb(List<List<Tree>> trainByFold, List<List<Tree>> testByFold,
                                          double alpha, double[] muGrid) {
        int b = trainByFold.size();
        KRegCCD[] models = new KRegCCD[b];
        for (int f = 0; f < b; f++) {
            // mu here is a placeholder: scoring overrides it via getLogProbabilityOfTree(tree, mu).
            models[f] = new KRegCCD(trainByFold.get(f), 0.0, KRegCCD.DEFAULT_MU, alpha,
                    RESERVE_DEPTH, KRegCCD.TailMode.NONE);
        }
        double bestMu = muGrid[0];
        double bestLogProb = Double.NEGATIVE_INFINITY;
        for (double mu : muGrid) {
            double total = 0.0;
            for (int f = 0; f < b; f++) {
                for (Tree t : testByFold.get(f)) {
                    total += models[f].getLogProbabilityOfTree(t, mu);
                }
            }
            if (total > bestLogProb) {
                bestLogProb = total;
                bestMu = mu;
            }
        }
        return new double[]{bestMu, bestLogProb};
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
