package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.CCD0;
import ccd.model.CCD0CP;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

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
}
