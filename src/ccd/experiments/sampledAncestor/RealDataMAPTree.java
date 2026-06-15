package ccd.experiments.sampledAncestor;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static ccd.algorithms.LoadOrStoreTrees.storeTrees;

public class RealDataMAPTree {

    public static void main(String[] args) throws IOException {

        String[] dataList = {
                "douglas_2024_spike_model",
                "near_2021_phylogeny_sunfishes",
                "savelyev_2020_bayesian_turkic_language"
        };

        for (String dataName : dataList) {

            String basePath = "/Volumes/DYNABOOK/phyloData_SA/" + dataName + "/";
            File folder = new File(basePath);
            File[] treeFiles = folder.listFiles((dir, name) -> name.endsWith(".trees") && !name.startsWith("._"));

            if (treeFiles != null) {
                for (File file : treeFiles) {
                    System.out.println("Processing: " + file.getName());
                    String shortName = file.getName().replaceAll("\\.trees$", "");
                    System.out.println("shortName: " + shortName);

                    // input file
                    String inputTreePath = file.getAbsolutePath();
                    TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(inputTreePath, 0);

                    CCD0 ccd = new CCD0(treeSet);
                    WrappedBeastTree mapTree = new WrappedBeastTree(ccd.getMAPTree());

                    // output file
                    List<Tree> trees = new ArrayList<>();
                    trees.add(mapTree.getWrappedTree());
                    File treeFile = new File(basePath + shortName + "_MAP_CCD0.trees");
                    storeTrees(trees, treeFile);

                }
            }
        }
    }
}
