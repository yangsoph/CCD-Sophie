package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.CCD0;
import ccd.model.CCD0CP;
import ccd.model.Clade;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class CCD0Test {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampledTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:2):1,D:3):0;"));
        trees.add(parse("(((D:1,C:1):1,B:2):1,A:3):0;"));
        return trees;
    }

    private static TreeAnnotator.MemoryFriendlyTreeSet sampledTreeSet() throws IOException {
        String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/ccd0_nonSA.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        return treeSet;
    }

    private static List<Tree> unsampledTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        return trees;
    }

    @Test
    public void testCCD0NonSA() throws IOException {
        List<Tree> treeList = sampledTreeList();
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = sampledTreeSet();
        CCD0 ccd = new CCD0(treeSet);
        double p1 = ccd.getProbabilityOfTree(treeList.get(0));
        System.out.println("p1 = " + p1);
    }

    @Test
    public void testCCD0CP() throws IOException {
        List<Tree> treeList = sampledTreeList();
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = sampledTreeSet();
        CCD0CP ccd = new CCD0CP(treeSet);
        double p1 = ccd.getProbabilityOfTree(treeList.get(0));
        System.out.println("p1 = " + p1);
    }

    /**
     * The {@link ccd.model.AbstractCCD#AbstractCCD(List, double)} constructor
     * has two branches; the {@code burnin == 0} branch must still set
     * {@code numBaseTrees} so {@link Clade#getCladeCredibility()} is finite.
     * Before the fix, {@code burnin == 0} left {@code numBaseTrees = 0}, so
     * {@code numOccurrences / 0 = Infinity} for every clade and the CCD0
     * constructor blew up with NaN inside {@code setPartitionProbabilities}.
     */
    @Test
    public void testZeroBurninSetsBaseTreeCount() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,C:1):1,(B:1,D:1):1):0;"));

        CCD0 ccd = new CCD0(trees, 0.0);

        assertEquals(trees.size(), ccd.getNumberOfBaseTrees());
        for (Clade c : ccd.getClades()) {
            double cred = c.getCladeCredibility();
            assertTrue(Double.isFinite(cred), "clade credibility must be finite, got " + cred + " for " + c);
        }
    }

}
