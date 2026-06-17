#!/usr/bin/env python3
"""Plot PIT (probability-integral-transform) calibration histograms for tree distributions.

Reads one or more tidy TSVs written by
ccd.experiments.regularisation.KRegPITExperiment (columns: model, sampleSize,
binLeft, binRight, count, expected, chiSq, chiP, ksD, ksP, nScored, nZeroed) and
lays them out as a grid of histograms: one row per sample size, one column per
model (KRegCCD, CCD1, CCD0). A calibrated, full-support model gives a flat PIT
(bars on the dashed uniform line, inside the shaded 95% band); a U-shape means the
model is over-confident, an n-shape under-confident, a slope a systematic bias.
Each panel is annotated with the chi-square and KS uniformity p-values, and (for
non-full-support models) the fraction of test trees excluded because q(T)=0.

Usage: plot_pit.py <out.pdf> <hist1.tsv> [hist2.tsv ...] [--title="..."]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# preferred column order and display names
MODEL_ORDER = ["kreg", "ccd1", "ccd0"]
MODEL_LABEL = {"kreg": "KRegCCD", "ccd1": "CCD1", "ccd0": "CCD0"}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    title = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--title=")), "PIT calibration")
    # cap the shared y-axis so the near-uniform panels stay readable; grossly miscalibrated
    # bars (e.g. CCD0's u->1 spike) clip and are labelled with their true height
    ycap = float(next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--ymax=")), 3.0))
    if len(args) < 2:
        sys.exit('usage: plot_pit.py <out.pdf> <hist1.tsv> [hist2.tsv ...] [--title="..."]')
    out, paths = args[0], args[1:]

    rows = np.concatenate([
        np.atleast_1d(np.genfromtxt(p, delimiter="\t", names=True, dtype=None, encoding="utf-8"))
        for p in paths
    ])

    models = [m for m in MODEL_ORDER if m in rows["model"]]
    models += [m for m in np.unique(rows["model"]) if m not in models]
    sizes = sorted(np.unique(rows["sampleSize"]))
    nrows, ncols = len(sizes), len(models)

    fig, axes = plt.subplots(nrows, ncols, figsize=(2.7 * ncols, 2.4 * nrows),
                             sharex=True, sharey=True, squeeze=False)

    for ri, size in enumerate(sizes):
        for ci, model in enumerate(models):
            ax = axes[ri][ci]
            sel = rows[(rows["model"] == model) & (rows["sampleSize"] == size)]
            if sel.size == 0:
                ax.set_axis_off()
                continue
            sel = sel[np.argsort(sel["binLeft"])]
            left, right, count = sel["binLeft"], sel["binRight"], sel["count"].astype(float)
            width = right - left
            m = float(sel["nScored"][0])
            density = count / (m * width)               # area integrates to 1

            # 95% pointwise band under H0: bin count ~ Binomial(m, 1/B), normal approx
            b = len(sel)
            sd = np.sqrt(m * (1.0 / b) * (1.0 - 1.0 / b))
            lo = max(0.0, (m / b - 1.96 * sd)) / (m * width[0])
            hi = (m / b + 1.96 * sd) / (m * width[0])
            ax.axhspan(lo, hi, color="0.85", zorder=0, label="95% band")
            ax.axhline(1.0, ls="--", lw=0.8, color="0.3", zorder=1)

            ax.bar(left, density, width=width, align="edge", color="#3b6fb6",
                   edgecolor="white", linewidth=0.3, zorder=2)

            # label any bar that clips the y-cap with its true density
            if density.max() > ycap:
                bi = int(np.argmax(density))
                ax.text(left[bi] + width[bi] / 2, ycap * 0.97, f"{density.max():.0f}",
                        ha="center", va="top", fontsize=6, color="#b00", zorder=3)

            l1, ksp, nzero = sel["pitL1"][0], sel["ksP"][0], int(sel["nZeroed"][0])
            excl = 100.0 * nzero / (m + nzero)
            # Two independent failure axes, coloured separately:
            #   L1   = shape calibration (sample-size-robust headline; KS p shown small because
            #          pooling inflates its power so it rejects trivially)
            #   excl = coverage (fraction of test trees the model assigns zero probability)
            ax.text(0.04, 0.96, f"L1 = {l1:.3f}", transform=ax.transAxes, va="top", ha="left",
                    fontsize=9, fontweight="bold", color=("#b00" if l1 > 0.15 else "#060"))
            ax.text(0.04, 0.81, f"KS p={ksp:.1g}", transform=ax.transAxes, va="top", ha="left",
                    fontsize=7, color="0.4")
            if nzero > 0:
                ax.text(0.04, 0.71, f"excl. {excl:.0f}%", transform=ax.transAxes, va="top", ha="left",
                        fontsize=7.5, fontweight="bold", color=("#b00" if excl > 1.0 else "#060"))

            if ri == 0:
                ax.set_title(MODEL_LABEL.get(model, model), fontsize=11)
            if ci == 0:
                ax.set_ylabel(f"n = {int(size)}\n\ndensity", fontsize=9)

    for ax in axes[-1]:
        ax.set_xlabel("PIT  $u(T)$", fontsize=9)
    for ax in axes.flat:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, ycap)

    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out)
    print(f"wrote {out}  ({nrows} sizes x {ncols} models)")


if __name__ == "__main__":
    main()
