package ccd.experiments;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Function;

import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.TreeDistances;
import ccd.model.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class PointEstSampledAncestor {

    public static int startRep = 0;
    public static int endRep = 100;

    public static void main(String[] args) throws IOException {

        // model name -> constructor
        Map<String, Function<TreeAnnotator.MemoryFriendlyTreeSet, AbstractCCD>> ccdModels = new LinkedHashMap<>();

        // ccdModels.put("CCD1CP", ts -> new CCD1CP(ts, false));
        // ccdModels.put("CCD1SJ", ts -> new CCD1SJ(ts, false));
        // ccdModels.put("CCD0CP", ts -> new CCD0CP(ts, false));
        // ccdModels.put("CCD0SJ", ts -> new CCD0SJ(ts, false));
        ccdModels.put("CCD1", ts -> new CCD1(ts, false));
        ccdModels.put("CCD0", ts -> new CCD0(ts, false));

        for (Map.Entry<String, Function<TreeAnnotator.MemoryFriendlyTreeSet, AbstractCCD>> entry : ccdModels.entrySet()) {

            String modelName = entry.getKey();
            Function<TreeAnnotator.MemoryFriendlyTreeSet, AbstractCCD> constructor = entry.getValue();

            // output file per model
            String outputPathName = "/nesi/nobackup/uoa04397/sophie/fossilBD/experiment/result_pointEst_" + modelName + ".csv";
            File outputFile = new File(outputPathName);

            try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {

                String separator = ",";
                // header
                StringBuilder sb = new StringBuilder();
                sb.append("rep").append(separator);
                sb.append("numTaxa").append(separator);
                sb.append("numSA_in_CCDmap").append(separator);
                sb.append("numSA_in_truth").append(separator);
                sb.append("SA_RF_dist").append(separator);
                sb.append("extended_RF_dist").append(separator);
                sb.append("topoRF_dist").append(separator);
                sb.append("SA_dist");
                writer.println(sb.toString());

                for (int i = startRep; i < endRep; i++) {

                    String treeSetDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/beast_output/FBD-" + i + ".trees";
                    String trueTreeDir = "/nesi/nobackup/uoa04397/sophie/fossilBD/wcss_sa/true_tree/FBD-" + i + "_true_fossilTree.trees";
                    TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeSetDir, 10);
                    TreeAnnotator.MemoryFriendlyTreeSet trueTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(trueTreeDir, 0);
                    trueTreeSet.reset();
                    WrappedBeastTreeWithSampledAncestor trueTree = new WrappedBeastTreeWithSampledAncestor(trueTreeSet.next());

                    int numTaxa = trueTree.getWrappedTree().getLeafNodeCount();

                    // instantiate CCD model
                    AbstractCCD ccd = constructor.apply(treeSet);

                    WrappedBeastTreeWithSampledAncestor mapTreeCCD = new WrappedBeastTreeWithSampledAncestor(ccd.getMAPTree());

                    int numSACCDMAP = mapTreeCCD.getNumberOfSampledAncestors();
                    int numSACCDTruth = trueTree.getNumberOfSampledAncestors();

                    int distExtendedRFCCD = TreeDistances.robinsonsFouldDistance(trueTree, mapTreeCCD);
                    double symDistSACCD = TreeDistances.sampledAncestorSymmetricDistance(trueTree, mapTreeCCD);
                    int distTopoRFCCD = TreeDistances.topologyRobinsonFoulds(trueTree, mapTreeCCD);
                    int distSARFCCD = TreeDistances.sampledAncestorRobinsonFoulds(trueTree, mapTreeCCD);

                    sb = new StringBuilder();
                    sb.append(i).append(separator);
                    sb.append(numTaxa).append(separator);
                    sb.append(numSACCDMAP).append(separator);
                    sb.append(numSACCDTruth).append(separator);
                    sb.append(distSARFCCD).append(separator);
                    sb.append(distExtendedRFCCD).append(separator);
                    sb.append(distTopoRFCCD).append(separator);
                    sb.append(symDistSACCD);
                    writer.println(sb.toString());
                }
                System.out.println("Finished: " + modelName);
            }
        }
    }
}
