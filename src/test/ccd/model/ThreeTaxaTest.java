package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.CCD1SJ;
import ccd.model.Clade;
import ccd.model.CladePartition;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class ThreeTaxaTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,C:0):0;"));
        trees.add(parse("((A:1,B:1):1,C:1):0;"));
        trees.add(parse("((A:1,B:0):1,C:0):0;"));
        return trees;
    }

    @Test
    public void testCCD1SJProbabilities() {
        List<Tree> trees = sampleTreeList();
        CCD1SJ sj = new CCD1SJ(trees, 0.0);

        for (Clade clade : sj.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }

        System.out.println("p(tree1) = " + sj.getProbabilityOfTree(trees.get(0)));
        System.out.println("p(tree2) = " + sj.getProbabilityOfTree(trees.get(1)));
        System.out.println("p(tree3) = " + sj.getProbabilityOfTree(trees.get(2)));
    }
}
