package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.KRegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Checks the closed-form (symbolic) probabilities of the 15 four-taxon trees against the
 * actual KRegCCD scorer. Sample: T1 = (((A,B),C),D) and T2 = (((D,C),B),A), each once.
 *
 * <p>With smoothing alpha and escape mu (NovelMode.FLAT, TailMode.NONE) the whole 15-tree
 * distribution reduces to three building blocks:
 * <ul>
 *   <li>root CCPs: cObs = (1+a)/(2+3a) and cBal = a/(2+3a), with 2*cObs + cBal = 1;</li>
 *   <li>root escape root eps solving 2 eps + 6 eps^2 = mu, i.e. eps = (sqrt(1+6 mu) - 1)/6;</li>
 *   <li>interior 3-clade escape root eps3 = mu/2.</li>
 * </ul>
 * The fifteen terms sum to (1 - mu)*(2 cObs + cBal) + mu = 1 exactly.
 */
public class KRegFourTaxaSymbolicTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");
    private static final double ALPHA = 0.4;
    private static final double MU = 0.05;
    private static final double TOL = 1e-12;

    private static Tree p(String nwk) {
        return new TreeParser(TAXA, nwk, 1, false);
    }

    private static List<Tree> sample() {
        return new ArrayList<>(Arrays.asList(
                p("(((A:1,B:1):1,C:1):1,D:1):0;"),
                p("(((D:1,C:1):1,B:1):1,A:1):0;")));
    }

    @Test
    public void symbolicMatchesScorer() {
        double a = ALPHA, mu = MU;
        double cObs = (1 + a) / (2 + 3 * a);
        double cBal = a / (2 + 3 * a);
        double eps = (Math.sqrt(1 + 6 * mu) - 1) / 6.0;
        double eps3 = mu / 2.0;

        double pObs = (1 - mu) * (1 - mu) * cObs;        // T1, T2
        double pBal = (1 - mu) * cBal;                   // T3 (alpha expansion)
        double pTriple = eps;                            // 1 blue clade over an observed cherry
        double pCherry = (1 - mu) * cObs * eps3;         // novel cherry under an observed triple
        double pTwoBlue = eps * eps;                     // 2 novel clades

        // newick -> expected closed form
        String[][] cases = {
                {"(((A:1,B:1):1,C:1):1,D:1):0;", "obs"},
                {"(((D:1,C:1):1,B:1):1,A:1):0;", "obs"},
                {"((A:1,B:1):1,(C:1,D:1):1):0;", "bal"},
                {"(((A:1,B:1):1,D:1):1,C:1):0;", "triple"},
                {"(((C:1,D:1):1,A:1):1,B:1):0;", "triple"},
                {"(((A:1,C:1):1,B:1):1,D:1):0;", "cherry"},
                {"(((B:1,C:1):1,A:1):1,D:1):0;", "cherry"},
                {"(((B:1,C:1):1,D:1):1,A:1):0;", "cherry"},
                {"(((B:1,D:1):1,C:1):1,A:1):0;", "cherry"},
                {"(((A:1,D:1):1,B:1):1,C:1):0;", "twoblue"},
                {"(((B:1,D:1):1,A:1):1,C:1):0;", "twoblue"},
                {"(((A:1,C:1):1,D:1):1,B:1):0;", "twoblue"},
                {"(((A:1,D:1):1,C:1):1,B:1):0;", "twoblue"},
                {"((A:1,C:1):1,(B:1,D:1):1):0;", "twoblue"},
                {"((A:1,D:1):1,(B:1,C:1):1):0;", "twoblue"},
        };

        KRegCCD ccd = new KRegCCD(sample(), 0, mu, a, 2, KRegCCD.TailMode.NONE);

        double sum = 0;
        for (String[] c : cases) {
            double expected = switch (c[1]) {
                case "obs" -> pObs;
                case "bal" -> pBal;
                case "triple" -> pTriple;
                case "cherry" -> pCherry;
                case "twoblue" -> pTwoBlue;
                default -> throw new IllegalStateException(c[1]);
            };
            double actual = ccd.getProbabilityOfTree(p(c[0]));
            assertEquals(expected, actual, TOL, "mismatch for " + c[0] + " (" + c[1] + ")");
            sum += actual;
        }

        // Backbone CCPs partition unity, and the escaped mass equals the reserve mu.
        assertEquals(1.0, 2 * cObs + cBal, TOL);
        assertEquals(mu, 2 * eps + 6 * eps * eps, TOL);
        assertEquals(1.0, sum, 1e-10, "the 15 trees must sum to 1");
    }
}
