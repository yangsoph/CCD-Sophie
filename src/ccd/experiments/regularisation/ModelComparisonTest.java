package ccd.experiments.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.algorithms.credibleSets.ProbabilityBasedCredibleSetComputer;
import ccd.model.KRegCCD;
import ccd.model.WrappedBeastTree;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class ModelComparisonTest {

    private static final double MU = KRegCCD.DEFAULT_MU;
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws IOException {

        // String filePath = "/Users/zyan598/Downloads/realData/output.txt";
        // File file = new File(filePath);
        // File parentDir = file.getParentFile();
        // if (parentDir != null && !parentDir.exists()) {
        //     parentDir.mkdirs();
        // }
        //
        // String[] models = new String[]{"dra-ctmc-est-strict-yule", "dra-ctmc-est-ucln-yule", "dra-ctmc-fixed-ucln-yule",
        //         "dra-ctmc4g-est-strict-yule", "dra-ctmc4g-est-ucln-yule", "dra-ctmc4g-fixed-strict-yule", "dra-ctmc4g-fixed-ucln-yule",
        //         "dra-sdollo-est-ucln-yule", "drav_cov_est_strict_yule", "drav_cov_est_ucln_yule", "drav_cov_fixed_strict_yule",
        //         "drav_cov_fixed_ucln_yule"};
        //
        // for (int i = 0; i < models.length; i++) {
        //     for (int j = 0; j < models.length; j++) {
        //         if (i != j) {
        //             String pointEstDataPath = "/Volumes/DYNABOOK/sophie/kolipakam_et_al2018-Dravidian/" + models[i] + ".trees";
        //             String credSetDataPath = "/Volumes/DYNABOOK/sophie/kolipakam_et_al2018-Dravidian/" + models[j] + ".trees";
        //             double credLevel = checkMAPTreeInCredSet(pointEstDataPath, credSetDataPath);
        //             System.out.println("pointEst = " + models[i] + " , credSet = " + models[j]);
        //             try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath, true))) {
        //                 writer.write("pointEst = " + models[i] + " , credSet = " + models[j] + " , cred level = " + credLevel);
        //                 writer.newLine();
        //             } catch (IOException e) {
        //                 e.printStackTrace();
        //             }
        //         }
        //     }
        // }

        // String pointEstDataPath = "/nesi/nobackup/uoa04397/sophie/real_data/kolipakam_2018_Dravidian_language/kolipakam_2018_ctmc-est-ucln.trees";
        String pointEstDataPath = "/Volumes/DYNABOOK/phyloData/kolipakam_2018_Dravidian_language/kolipakam_2018_ctmc-est-ucln.trees";
        // String credSetDataPath = "/nesi/nobackup/uoa04397/sophie/real_data/stervander_2019_origin_flightless_bird/stervander_2019_2Nc3Mt_1tree.trees";
        String credSetDataPath = "/Volumes/DYNABOOK/phyloData/kolipakam_2018_Dravidian_language/kolipakam_2018_ctmc-est-strict.trees";
        checkMAPTreeInCredSet(pointEstDataPath, credSetDataPath);
    }

    public static double checkMAPTreeInCredSet(String pointEstDataPath, String credSetDataPath) throws IOException {

        File pointEstFile = new File(pointEstDataPath);
        List<Tree> pointEstTreeList = LoadOrStoreTrees.loadTrees(pointEstFile, 0);

        File credSetFile = new File(credSetDataPath);
        List<Tree> credSetreeList = LoadOrStoreTrees.loadTrees(credSetFile, 0);

        // MAP tree
        KRegCCD kreg = KRegCCD.withOptimisedMu(pointEstTreeList);
        // double optiMu = 0.017747683;
        // KRegCCD kreg = new KRegCCD(pointEstTreeList, 0, MU, ALPHA, K);
        WrappedBeastTree mapTree = new WrappedBeastTree(kreg.getMAPTree());

        // Cred Set
        // double optiMuCredSet = 0.002016054;
        KRegCCD kregCredSet = KRegCCD.withOptimisedMu(credSetreeList);
        // KRegCCD kregCredSet = new KRegCCD(credSetreeList, 0, MU, ALPHA, K);
        ProbabilityBasedCredibleSetComputer credSetComputer = new ProbabilityBasedCredibleSetComputer(kregCredSet);

        // check if the MAP tree is in the credible set
        double credibleLevel = credSetComputer.getCredibleLevel(mapTree.getWrappedTree());
        System.out.println("The MAP tree has a credible level of " + credibleLevel + " in the credible set");

        return credibleLevel;
    }
}
