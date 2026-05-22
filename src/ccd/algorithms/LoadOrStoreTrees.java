package ccd.algorithms;

import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class LoadOrStoreTrees {

    public static List<Tree> loadTrees(File treeFile, double burninPercentage) throws IOException {
        NexusParser parser = new NexusParser();
        parser.parseFile(treeFile);
        System.out.println("loading tree file: " + treeFile);
        List<Tree> fullTreeList = parser.trees;
        List<Tree> treesToUse;
        if (burninPercentage == 0) {
            treesToUse = fullTreeList;
        } else {
            int numDiscardedTrees = (int) (fullTreeList.size() * burninPercentage);
            int numUsedTrees = fullTreeList.size() - numDiscardedTrees;
            treesToUse = new ArrayList<Tree>(numUsedTrees);
            treesToUse.addAll(fullTreeList.subList(numDiscardedTrees, fullTreeList.size()));
        }
        return treesToUse;
    }

    // TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.getAbsolutePath(), 10);
    // AbstractCCD initialCCD1 = new CCD1(treeSet, false);  // Create initial CCD1 object

    public static void storeTrees(List<Tree> trees, File treeFile) throws FileNotFoundException {
        try (PrintStream handle = new PrintStream(treeFile)) {

            trees.get(0).init(handle); // Write header info
            handle.println(); // Ensure the semicolon before the first tree is on its own line

            for (int i = 0; i < trees.size(); i++) {
                Tree tree = trees.get(i);

                int id = (tree.getID() == null) ? i : Integer.parseInt(tree.getID().replace("_", ""));
                tree.log(id, handle);

                handle.println(); // blank line between trees
            }

            trees.get(trees.size() - 1).close(handle);
        }
    }
}
