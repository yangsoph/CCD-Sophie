package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.CCD1;
import ccd.model.CCD1CP;
import ccd.model.Clade;
import ccd.model.WrappedBeastTree;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Reproduces the four-taxon worked example from doc/sa-models.tex.
 * <p>
 * Three sampled trees on taxa {A,B,C,D}, each with count 1:
 * Tree 1: ((A,B:0),(C,D))   — B is SA
 * Tree 2: ((A:0,B),(C,D))   — A is SA
 * Tree 3: (((A,B),C:0),D)   — C is SA
 * <p>
 * Expected probabilities (in eighty-firsts):
 * CP:  18, 18, 9
 * SJ:  27, 27, 27
 */

public class CCD1CPTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:0,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:0):1,D:2):0;"));
        return trees;
    }

    private static TreeAnnotator.MemoryFriendlyTreeSet sampleTreeSet() throws IOException {
        String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/4taxa_alexei.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        return treeSet;
    }

    private static double rounded(double x) {
        return Math.round(x * 81 * 1e6) / 1e6 / 81.0;
    }

    @Test
    public void treeSetTest() throws IOException {
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = sampleTreeSet();
        // CCD1 ccd = new CCD1CP(treeSet, false);
        // WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd.getMAPTree());
        // System.out.println("CCD MAP tree log prob = " + ccd.getMaxLogTreeProbability());
    }

    @Test
    public void testCCD1CPProbabilities() {
        List<Tree> trees = sampleTreeList();
        CCD1CP cp = new CCD1CP(trees, 0.0);

        double p1 = cp.getProbabilityOfTree(trees.get(0));
        double p2 = cp.getProbabilityOfTree(trees.get(1));
        double p3 = cp.getProbabilityOfTree(trees.get(2));

        assertEquals(18.0 / 81.0, rounded(p1), 1e-9, "CP P(Tree 1)");
        assertEquals(18.0 / 81.0, rounded(p2), 1e-9, "CP P(Tree 2)");
        assertEquals(9.0 / 81.0, rounded(p3), 1e-9, "CP P(Tree 3)");

    }
}
