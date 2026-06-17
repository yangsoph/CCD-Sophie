# KRegCCD calibration notes

Working notes on assessing how well **KRegCCD** (the full-support regularised clade-conditional
distribution) is calibrated against held-out trees, and how it compares to CCD1 and CCD0. Covers the
PIT histogram experiment done so far, an effective-sample-size caveat on the fitted parameters, the
planned floor-lowering re-run, and the menu of other calibration experiments we have discussed.

Status as of 2026-06-18: PIT histogram experiment complete on Yule50 (99 reps). The other experiments
are proposed, not yet built.

---

## 1. What "calibration" means here

The model `q` is fit to a set of training trees and we ask whether the probabilities it assigns to
*held-out* trees are statistically consistent with how often those trees actually occur. Two
distinct ways a model can fail:

- **Shape miscalibration** — among the trees it *can* score, the assigned probabilities are
  systematically too high or too low (over/under-confident).
- **Coverage failure** — the model assigns probability **zero** to held-out trees (novel clades), so
  those trees can never fall in any credible set. This is the failure full support is meant to cure;
  CCD0/CCD1 suffer from it, KRegCCD does not.

We keep these two axes separate throughout.

---

## 2. The PIT histogram test (done)

### Statistic
For each held-out test tree `T` we compute the **probability-integral-transform** statistic

```
u(T) = P_{S ~ q}[ q(S) >= q(T) ]
```

— the probability that a tree drawn from the model is at least as probable as `T` — estimated by
Monte Carlo (fraction of 100k model samples with log-prob ≥ `log q(T)`). This is exactly the credible
level at which `T` first enters the model's highest-probability credible set, i.e. the statistic
underlying the credible-set coverage curve. **If `q` is calibrated and `T ~ q`, then
`u(T) ~ Uniform(0,1)`.** The coverage curve is the empirical CDF of `u`; the PIT histogram is its
density view, which localises the failure (∪ = over-confident, ∩ = under-confident, slope = bias).

### Randomized PIT
The tree distribution is discrete, so the raw `u(T)` is not exactly uniform (it clumps on achievable
mass levels). We use the **randomized PIT** `V(T) ~ Uniform(u_lo, u_hi)`, with
`u_lo = P[q(S) > q(T)]` and `u_hi = P[q(S) >= q(T)]`, which is exactly `Uniform(0,1)` under `T ~ q`
(Dunn & Smyth 1996). The histograms and tests use `V`.

### Pooling across replicates
Each replicate is an independent simulated dataset with its own posterior. We pool the PIT across
replicates: each tree's `u` is computed under *its own* replicate's fitted model, and a mixture of
uniforms is uniform, so pooling is valid and simply lowers sampling noise.

### Metric: L1, not the p-values
We report the **L1 calibration error** `sum_i | count_i/m - 1/B |` (L1 distance of the histogram from
uniform; the density analogue of the coverage curve's "L1 norm"). We do **not** lead with the
chi-square / KS p-values: pooling makes `m` large (up to ~3×10^5 here), and at that size the tests
reject any infinitesimal deviation — every p-value in the final run is ≈0, including for visibly flat
KRegCCD panels. L1 is an effect size; its H0 noise floor is `≈ 0.798·sqrt((B-1)/m)`, which shrinks
like `1/sqrt(m)`, so it stays meaningful. (Noise floor at B=20: ~0.20 at m=300, ~0.045 at m=6000,
~0.007 at m=3×10^5.) The p-values are kept in the TSV and shown small on the figure for completeness.

### Setup
- Data: Yule, 50 taxa, figshare collection 7102354 (local copy at
  `~/Git/cladeDiscoveryCurve/data/Yule50`). Each replicate has two **independent** MCMC chains:
  `run1` (training) and `run2` (testing), ~35k trees each. Using independent chains for train/test
  means the held-out split is clean (no shared autocorrelation across the split).
- Sample sizes (even-thinned subsamples): n = 300, 1000, 3000, for both training and testing.
- Models: KRegCCD (`withOptimisedParameters`), CCD1, CCD0.
- 99 replicates pooled in this run. rep100 was skipped because its run trees are *misnamed*
  `yule-n50-0.trees` (a 100→0 collision in the data pipeline), not missing — the files are complete
  (35,001 trees, 50 taxa). `KRegPITSweep` now globs the single `<prefix>*.trees` per run dir instead
  of reconstructing the name, so the planned re-run will include all 100.

### Results

| n | KRegCCD L1 (excl.) | CCD1 L1 (excl.) | CCD0 L1 (excl.) |
|------|--------------------|------------------|------------------|
| 300  | **0.102** (0%)     | 0.028 (29%)      | 0.632 (11%)      |
| 1000 | **0.055** (0%)     | 0.022 (14%)      | 0.653 (5%)       |
| 3000 | **0.039** (0%)     | 0.019 (8%)       | 0.672 (2%)       |

Figure: `pit-yule50.pdf` (3×3 grid, rows = n, cols = model), data `pit-yule50.tsv`.

- **KRegCCD** is the only full-support model (0% excluded) and the only one that is both well-shaped
  *and* covers every test tree. Its L1 improves with n but stays above the noise floor — a small
  residual shape miscalibration remains even at n=3000 (see §3 and §4 for a likely cause).
- **CCD1** has the *flattest* histogram on the trees it supports (lowest L1) but **excludes 8–29%**
  of held-out trees — its failure is coverage, not shape. On the figure this is the green L1 next to
  the red "excl." annotation.
- **CCD0** is grossly shape-miscalibrated (L1 ≈ 0.63–0.67, a monotonic ramp toward `u→1`, i.e. it
  assigns the real trees far less probability than its own typical samples), and does not improve
  with n.
- All χ²/KS p-values ≈ 0 — confirms L1 was the right headline.

### Fitted parameters
Figure: `pit-yule50-params.pdf` (log-log scatter of the 297 `(alpha, mu)` pairs coloured by n),
data `pit-yule50.params.tsv`.

- Median `alpha`: 0.46 → 0.40 → 0.36; median `mu`: 0.0043 → 0.0015 → **0.0010** for n = 300/1000/3000.
  Both fall with n; escape mass `mu` falls fastest.
- **72 of 99 fits at n=3000 (and 21 at n=1000) are pinned at the `mu = 1e-3` search floor.** At large
  n the optimiser's lower bound is binding — the data wants *less* escape mass than the grid allows.
  Only one n=300 fit touches the 0.05 ceiling. This both makes the large-n parameter estimates
  unreliable and may cap KRegCCD's large-n calibration (μ is held above its true optimum).

---

## 3. Caveat: ESS ≠ nominal sample size

The "sample size" axis (300/1000/3000 trees) is the **nominal** count of thinned posterior trees, not
the **effective** sample size (ESS) — the number of effectively independent topologies. MCMC trees
are autocorrelated, and thinning 35k → n does not remove that; the topological ESS can be far below n
and may **saturate** (a chain only explores so many effectively independent topologies, so beyond some
point extra correlated trees add little information). Consequences:

1. **The n axis overstates information.** "n = 3000" may carry far less independent signal than 3000;
   the apparent `mu → floor` and `alpha ↓` trends with n could partly reflect ESS rising then
   saturating rather than genuine large-sample behaviour. Parameters should arguably be plotted
   against **ESS**, not nominal n.

2. **The (alpha, mu) selection may be overfit.** `KRegCCDParameterOptimiser.optimise` chooses
   `(alpha, mu)` by maximising **k-fold cross-validated** held-out log-probability, with the default
   **`STRIDED`** fold assignment (tree `i` → fold `i mod b`). On an autocorrelated chain, strided
   folds put each held-out tree's near-neighbours (highly correlated) into the *training* part of the
   same fold, so the CV held-out set is **not independent** of its training set. CV then
   *underestimates* generalisation error and selects **under-regularised** parameters — too-small
   `mu`, too-small `alpha`. This is a plausible driver of both the `mu`-floor pinning at large n and
   the residual KRegCCD shape miscalibration. (The *final* PIT evaluation on `run2` is clean — it is
   an independent chain — so the calibration verdict stands; it is the *parameter fitting* that is
   suspect.)

### Remedies to try
- Measure topological ESS of the chains (e.g. clade-trace / pseudo-ESS) and re-plot fitted parameters
  against ESS instead of nominal n.
- Thin training trees by the autocorrelation time so that nominal n ≈ ESS before fitting.
- **Done:** switched the CV default from `STRIDED` to **`CONTIGUOUS`** (block folds) in
  `KRegCCDParameterOptimiser` — a contiguous held-out block is far less correlated with its training
  set, so the objective is honest for autocorrelated samples. Correct on principle, but a 1-rep check
  (below) shows it barely moves the fitted `mu`, so it is a correctness fix, not the cure for the
  floor-pinning. (Affects all `withOptimisedParameters` callers, incl. the credible-set experiment.)
- Sanity check: select `(alpha, mu)` on the independent `run2` chain rather than internal CV, and see
  how much the parameters move.

---

## 4. Planned re-run: lower the mu floor (+ ESS arm)

Goal: get unconstrained large-n fits and test whether the residual KRegCCD miscalibration is an
artefact of the floor or of overfitting.

**Arm A — lower the search floor (DONE).** In `KRegCCDParameterOptimiser`: `MU_LO` `1e-3` → `1e-4`,
`MU_GRID` `41` → `61` (keep per-decade density over the wider range), and a floor warning mirroring
the existing ceiling warning. Verified: a 2-rep n=3000 check now finds interior optima below the old
floor (rep2 `mu = 0.00031`, impossible under the 1e-3 grid) with no floor warning firing — the lowered
floor is reachable and not itself binding here.

**Arm B — ESS-aware fitting (DONE: CONTIGUOUS is now the default).** The worry was that strided CV
overfits and pins `mu` artificially low. A 1-rep check (STRIDED vs CONTIGUOUS) shows otherwise: fitted
`mu` is essentially unchanged — n=300/1000/3000 `mu` = 0.01046/0.00357/0.00122 (STRIDED) →
0.01046/0.00323/0.00110 (CONTIGUOUS), i.e. slightly *lower*, not higher; `alpha` shifts modestly and
not in a clean direction. So **the floor-pinning is a genuine preference for small `mu` at large n,
not a strided artifact** — Arm A is still required. (Autocorrelation thinning to ESS remains a
possible refinement, but the param effect looks small.)

**Then:** re-run the sweep with both changes (CONTIGUOUS default + lowered floor), now over all 100
reps (~8 h, seeded → reproducible), and re-render `pit-yule50*.pdf`. Expect the floor-pinning to
resolve under Arm A (interior `mu` at large n); whether the n=3000 KRegCCD L1 then improves is the
open question.

Reproduce command (current code; `build/` is stale, so compile against source):
```
javac -cp "../beast2/build:../BeastFX/build:../asm/build:../beast2/lib/*" -sourcepath src -d /tmp/pitbuild \
  src/ccd/experiments/regularisation/KRegPITSweep.java
java -Xmx8g -cp "/tmp/pitbuild:../beast2/build:../BeastFX/build:../asm/build:../beast2/lib/*" \
  ccd.experiments.regularisation.KRegPITSweep \
  ~/Git/cladeDiscoveryCurve/data/Yule50 yule-n50 1 100 300,1000,3000 kreg,ccd1,ccd0 20 100000 42 \
  pit-yule50.tsv
python3 doc/plot_pit.py pit-yule50.pdf pit-yule50.tsv --title="PIT calibration — Yule50"
python3 doc/plot_params.py pit-yule50-params.pdf pit-yule50.params.tsv --title="Fitted KRegCCD parameters — Yule50"
```

---

## 5. Other calibration experiments discussed (not yet built)

Ordered roughly by value-per-effort. PIT histogram (above) was the first.

1. **Clade-marginal reliability diagram** — bin clades by their model marginal probability and compare
   to the empirical frequency of that clade among test trees ("are 90%-supported clades present 90% of
   the time?"). Lower-dimensional than whole-tree calibration, far more data per bin, and localises
   *which* conditional distributions are miscalibrated. High value.
2. **Novel-clade-stratified coverage** — repeat the PIT/coverage check conditioned on the number of
   novel clades in the test tree. This directly stress-tests the full-support claim on the trees
   CCD0/CCD1 zero out — the regime KRegCCD exists for.
3. **Per-tree probability reliability / ECE** — bin test trees by predicted `q(T)`, compare mean
   predicted probability to observed relative frequency; summarise as Expected Calibration Error.
   Calibrates the probability *values*, not just credible-set mass.
4. **KL / cross-entropy fit comparison** — held-out negative log-likelihood is the cross-entropy
   `H(p,q)`; **differences** of held-out NLL between models are exactly differences in `KL(p‖q)`
   (the unknown `H(p)` cancels), so they rank models by KL with no extra assumptions. For an absolute
   KL, estimate `H(p)` with a near-oracle reference model fit on a large sample. Full support is what
   makes `H(p,q)` finite — CCD0/CCD1 give `+∞` on novel-clade trees. Pair with the model self-entropy
   `H(q)` as a sharpness companion (calibration alone is insufficient — maximise sharpness subject to
   calibration). See the entropy writeup (`kregccd-entropy.tex`).
5. **Simulation-Based Calibration (SBC)** — across many simulated datasets, rank the true generating
   tree within the fitted model (via the mass-rank statistic); ranks should be uniform. Stronger than
   single-posterior coverage (tests the whole prior-predictive); the biggest build.
6. **Proper-scoring / Murphy decomposition** — the reliability component of the Brier (or log) score
   *is* a calibration metric, and ties calibration to sharpness in one number.
7. **MMD / energy two-sample test** — kernel two-sample test between model samples and held-out trees
   using a tree kernel/metric (Robinson–Foulds, Kendall–Colijn); an omnibus "are these the same
   distribution?" check needing no binning.

---

## 6. Files

Code:
- `src/ccd/experiments/regularisation/KRegPITExperiment.java` — single (model, test-set) PIT run;
  reusable static helpers (PIT, histogram, uniformity tests, L1).
- `src/ccd/experiments/regularisation/KRegPITSweep.java` — pools the PIT over reps × sizes × models;
  writes the PIT histogram TSV and the `(alpha, mu)` params TSV.
- `src/ccd/model/KRegCCD.java` — added `getAlpha()` / `getMu()` (and an `alpha` field) so the sweep
  can record fitted parameters.

Plots (matplotlib, vector output):
- `doc/plot_pit.py` — PIT histogram grid (L1 headline, 95% band, capped y-axis, off-scale peaks
  labelled).
- `doc/plot_params.py` — `(alpha, mu)` scatter, log-log, coloured by n, with search floor/ceiling
  lines and per-size median stars.

Outputs (Yule50, 99 reps, seed 42):
- `doc/pit-yule50.pdf`, `doc/pit-yule50.tsv` — PIT calibration grid + data.
- `doc/pit-yule50-params.pdf`, `doc/pit-yule50.params.tsv` — fitted-parameter scatter + data.

Reproducibility: the sweep seeds `Randomizer` (and the PIT jitter RNG), so re-running the same command
reproduces the results exactly.
