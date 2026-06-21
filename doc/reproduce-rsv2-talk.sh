#!/usr/bin/env bash
# Reproduce every RSV2 (and the Yule50 PIT) result used in the Nagoya talk
# ("Tractable time tree distributions for Bayesian phylogenetic inference").
#
# Run from the CCD-Sophie repo root:
#     bash doc/reproduce-rsv2-talk.sh
#
# Produces, all from ONE matched RSV2 chain and matched sample sizes
# (train = 1000 evenly-spread trees from the first half of the chain,
#  test  = 1000 evenly-spread trees from the second half; alpha = 0.4):
#   - the three-model coverage / log-prob table          -> stdout + doc/rsv2-threemodel.out
#   - the mu-calibration sweep (40 log-spaced mu)         -> doc/rsv2-musweep.{tsv,pdf}
#   - the normalisation (inflation vs tail bound) table   -> doc/rsv2-inflation.tsv
#   - the Yule50 PIT calibration grid (CCD1/regCCD/KReg)  -> doc/pit-yule50-3model.{tsv,pdf}
# and copies the two figures into the talk. The slide tables (RSV2 coverage,
# Normalisation) are transcribed by hand from rsv2-threemodel.out / rsv2-inflation.tsv.
set -euo pipefail

# --- paths (override via environment) ---
RSV2=${RSV2:-$HOME/Git/bayesianPhylogeneticInstability/data/zenodo/bayesian_phylogenetic_instability_posterior_trees/RSV2.trees}
YULE50=${YULE50:-$HOME/Git/cladeDiscoveryCurve/data/Yule50}
TALK=${TALK:-$HOME/Git/nagoya-talk}
THREADS=${THREADS:-8}
CP="../beast2/build:../BeastFX/build:../beast2/lib/*"   # asm not needed
BUILD=/tmp/ccd-talk-build

mkdir -p "$BUILD"

echo "== compiling experiments =="
javac -cp "$CP" -sourcepath src -d "$BUILD" \
  src/ccd/experiments/regularisation/RSV2ThreeModelTable.java \
  src/ccd/experiments/regularisation/KRegPITSweep.java

echo "== RSV2: three-model table + mu sweep + inflation =="
java -Xmx32g -cp "$BUILD:$CP" \
  ccd.experiments.regularisation.RSV2ThreeModelTable "$RSV2" 0.1 1000 | tee doc/rsv2-threemodel.out
# -> writes doc/rsv2-musweep.tsv and doc/rsv2-inflation.tsv; CV-fitted mu is printed (= 0.00462).

echo "== Yule50 PIT sweep (CCD1, regCCD, KRegCCD; 100 reps, sizes 300/1000/3000) =="
java -Xmx32g -cp "$BUILD:$CP" ccd.experiments.regularisation.KRegPITSweep \
  "$YULE50" yule-n50 1 100 300,1000,3000 ccd1,regccd,kreg 20 100000 42 doc/pit-yule50-3model.tsv "$THREADS"

echo "== plots =="
python3 doc/plot_musweep.py doc/rsv2-musweep.pdf doc/rsv2-musweep.tsv --cvmu=0.00462
python3 doc/plot_pit.py     doc/pit-yule50-3model.pdf doc/pit-yule50-3model.tsv --title=""

echo "== copy figures into the talk =="
cp doc/rsv2-musweep.pdf      "$TALK/figures/kregccd/musweep.pdf"
cp doc/pit-yule50-3model.pdf "$TALK/figures/kregccd/pit-yule50.pdf"

echo
echo "Done."
echo "  Coverage / log-prob table values  -> doc/rsv2-threemodel.out  (RSV2 result slide)"
echo "  Normalisation table values        -> doc/rsv2-inflation.tsv   (Normalisation slide)"
echo "  Figures copied into $TALK/figures/kregccd/"
