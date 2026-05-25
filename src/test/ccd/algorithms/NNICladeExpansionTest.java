package test.ccd.algorithms;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.algorithms.NNICladeExpansion;
import ccd.algorithms.NNICladeExpansion.NovelClade;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.algorithms.NNICladeExpansion.Provenance;
import ccd.model.CCD0;
import ccd.model.Clade;
import ccd.model.bitsets.BitSet;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for {@link NNICladeExpansion} against hand-computed NNI neighbourhoods.
 */
public class NNICladeExpansionTest {

    private static final List<String> TAXA5 = Arrays.asList("A", "B", "C", "D", "E");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA5, newick, 1, false);
    }

    /** Taxa of a novel clade as a sorted letter string, e.g. "ABD". */
    private static String name(BitSet taxa, Tree tree) {
        StringBuilder sb = new StringBuilder();
        String[] names = tree.getTaxaNames();
        for (int i = taxa.nextSetBit(0); i >= 0; i = taxa.nextSetBit(i + 1)) {
            sb.append(names[i]);
        }
        return sb.toString();
    }

    private static Set<String> novelNames(NNICladeExpansion exp, Tree tree) {
        Set<String> names = new TreeSet<>();
        for (NovelClade novel : exp.getNovelClades()) {
            names.add(name(novel.taxa, tree));
        }
        return names;
    }

    /**
     * Single tree (((A,B),C),D). The NNI neighbourhood of its three internal
     * edges produces exactly the novel clades {ABD, CD, AC, BC}; with one tree
     * the graph-wide and co-occurring modes must coincide.
     */
    @Test
    public void testSingleTreeExactNeighbourhood() {
        List<Tree> trees = new ArrayList<>();
        Tree tree = parse("(((A:1,B:1):1,C:2):1,D:3):0;");
        trees.add(tree);
        CCD0 ccd = new CCD0(trees, 0.0);

        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        NNICladeExpansion strict = new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();

        Set<String> expected = new TreeSet<>(Arrays.asList("ABD", "CD", "AC", "BC"));
        assertEquals(expected, novelNames(broad, tree), "graph-wide novel clades");
        assertEquals(expected, novelNames(strict, tree), "co-occurring novel clades");
    }

    /**
     * Every novel clade's NNI seed split is {kept, sibling}, and both of those
     * are clades that already exist in C - which is what makes the Phase 2
     * extension feasible. Check this on the single-tree case.
     */
    @Test
    public void testSeedSplitConstituentsExist() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:2):1,D:3):0;"));
        CCD0 ccd = new CCD0(trees, 0.0);

        NNICladeExpansion exp = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        for (NovelClade novel : exp.getNovelClades()) {
            assertTrue(!novel.provenances.isEmpty(), "novel clade has a provenance");
            for (Provenance prov : novel.provenances) {
                assertTrue(ccd.getClades().contains(prov.kept()), "kept piece is in C");
                assertTrue(ccd.getClades().contains(prov.sibling()), "sibling is in C");
                // novel taxa == kept ∪ sibling
                BitSet union = (BitSet) prov.kept().getCladeInBitsTaxaOnly().clone();
                union.or(prov.sibling().getCladeInBitsTaxaOnly());
                assertEquals(novel.taxa, union, "novel taxa equal kept ∪ sibling");
            }
        }
    }

    /**
     * Two 5-taxon trees that share the clade {A,B,C} but split it differently
     * ({AB,C} under {ABCD} in tree 1; {AC,B} under {ABCE} in tree 2). Graph-wide
     * pairing recombines split contexts that never co-occurred, so it strictly
     * exceeds the co-occurring (true sample-NNI) neighbourhood.
     */
    @Test
    public void testGraphWideExceedsCoOccurring() {
        List<Tree> trees = new ArrayList<>();
        Tree t1 = parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;");
        Tree t2 = parse("((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;");
        trees.add(t1);
        trees.add(t2);
        CCD0 ccd = new CCD0(trees, 0.0);

        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        NNICladeExpansion strict = new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();

        Set<String> broadNames = novelNames(broad, t1);
        Set<String> strictNames = novelNames(strict, t1);

        Set<String> expectedBroad = new TreeSet<>(Arrays.asList(
                "ABD", "CD", "ACD", "BD", "ABE", "CE", "ACE", "BE", "DE", "BC"));
        Set<String> expectedStrict = new TreeSet<>(Arrays.asList(
                "DE", "ABD", "CD", "BC", "ACE", "BE"));
        assertEquals(expectedBroad, broadNames, "graph-wide novel clades");
        assertEquals(expectedStrict, strictNames, "co-occurring novel clades");

        // strict neighbourhood is contained in the graph-wide one
        assertTrue(broadNames.containsAll(strictNames), "strict ⊆ broad");

        // the cross-tree recombinations are exactly the graph-wide extras
        Set<String> onlyBroad = new TreeSet<>(broadNames);
        onlyBroad.removeAll(strictNames);
        assertEquals(new TreeSet<>(Arrays.asList("ACD", "BD", "ABE", "CE")), onlyBroad,
                "graph-wide-only (cross-tree) clades");
    }

    /**
     * Co-occurring mode needs all base trees; a CCD built without storing them
     * should fail loudly rather than silently producing wrong results.
     */
    @Test
    public void testCoOccurringRequiresBaseTrees() {
        // empty CCD with storeBaseTrees = false: only the first tree is kept,
        // so getBaseTrees() returns null and co-occurring mode cannot run
        CCD0 ccd = new CCD0(5, false);
        ccd.addTree(parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;"));
        ccd.addTree(parse("((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;"));
        ccd.initialize();

        boolean threw = false;
        try {
            new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();
        } catch (IllegalStateException expected) {
            threw = true;
        }
        // graph-wide still works without stored base trees
        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        assertTrue(broad.getNumberOfNovelClades() >= 0);
        assertTrue(threw, "co-occurring mode should require stored base trees");
    }
}
