#!/usr/bin/env python3
"""Plot the fitted KRegCCD (alpha, mu) parameters from the PIT sweep.

Reads the params TSV written by ccd.experiments.regularisation.KRegPITSweep
(columns: rep, sampleSize, alpha, mu) -- one row per KRegCCD fit.

Two modes, chosen automatically from the data:
  * if alpha varies across fits, the (alpha, mu) pairs are scattered on log-log
    axes, coloured by training sample size, with per-size median stars;
  * if alpha is constant (the optimiser pinned it, e.g. alpha = 0.4), the alpha
    axis carries no information, so instead mu is shown as a violin per sample
    size with the individual fits jittered over it.
In both modes a dashed line marks the reliability ceiling mu = 0.05
(KRegCCD.MU_RELIABLE_MAX, the upper bound of the optimiser's mu search) and a
dotted line the mu search floor; points on either are fits whose held-out
likelihood was still improving in mu at that bound.

Usage: plot_params.py <out.pdf> <params.tsv> [--title="..."] [--muceil=0.05] [--mufloor=1e-4]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = ["#3b6fb6", "#d1772f", "#4c9a52", "#9467bd"]


def plot_scatter(ax, rows, sizes, muceil, mufloor):
    ax.axhline(muceil, ls="--", lw=0.9, color="0.4", zorder=1, label=f"$\\mu$ search ceiling = {muceil:g}")
    ax.axhline(mufloor, ls=":", lw=0.9, color="0.4", zorder=1, label=f"$\\mu$ search floor = {mufloor:g}")
    for i, size in enumerate(sizes):
        sel = rows[rows["sampleSize"] == size]
        ax.scatter(sel["alpha"], sel["mu"], s=22, alpha=0.6, color=COLORS[i % len(COLORS)],
                   edgecolor="white", linewidth=0.3, zorder=2,
                   label=f"n = {int(size)}  ({len(sel)} fits)")
        ax.scatter(np.median(sel["alpha"]), np.median(sel["mu"]), s=160, marker="*",
                   color=COLORS[i % len(COLORS)], edgecolor="black", linewidth=0.7, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\alpha$  (backbone smoothing)")
    ax.set_ylabel(r"$\mu$  (escape probability)")
    ax.legend(fontsize=8, framealpha=0.9, loc="best")
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)


def plot_violin(ax, rows, sizes, muceil, mufloor, alpha):
    """Violin of mu per sample size (in log10 space) with the individual fits jittered over it."""
    rng = np.random.RandomState(0)  # fixed seed: jitter is cosmetic, keep the figure reproducible
    positions = np.arange(1, len(sizes) + 1)
    logmu = [np.log10(rows[rows["sampleSize"] == s]["mu"]) for s in sizes]

    parts = ax.violinplot(logmu, positions=positions, widths=0.8, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(COLORS[i % len(COLORS)])
        body.set_alpha(0.25)
        body.set_edgecolor(COLORS[i % len(COLORS)])
    for i, (pos, vals) in enumerate(zip(positions, logmu)):
        x = pos + rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(x, vals, s=14, alpha=0.55, color=COLORS[i % len(COLORS)],
                   edgecolor="white", linewidth=0.3, zorder=3)
        med = np.median(vals)
        ax.hlines(med, pos - 0.4, pos + 0.4, color="black", lw=1.6, zorder=4)
        ax.text(pos + 0.43, med, f"med $\\mu$={10 ** med:.2g}", va="center", fontsize=8)

    ax.axhline(np.log10(muceil), ls="--", lw=0.9, color="0.4", zorder=1,
               label=f"$\\mu$ search ceiling = {muceil:g}")
    ax.axhline(np.log10(mufloor), ls=":", lw=0.9, color="0.4", zorder=1,
               label=f"$\\mu$ search floor = {mufloor:g}")

    # log10 axis labelled with the actual mu decades
    lo = int(np.floor(min(np.log10(mufloor), min(v.min() for v in logmu))))
    hi = int(np.ceil(max(np.log10(muceil), max(v.max() for v in logmu))))
    ax.set_yticks(range(lo, hi + 1))
    ax.set_yticklabels([f"$10^{{{k}}}$" for k in range(lo, hi + 1)])
    ax.set_xticks(positions)
    ax.set_xticklabels([f"n = {int(s)}\n({len(rows[rows['sampleSize'] == s])} fits)" for s in sizes])
    ax.set_ylabel(r"$\mu$  (escape probability)")
    ax.set_xlabel(rf"training sample size   ($\alpha$ fixed at {alpha:g})")
    ax.legend(fontsize=8, framealpha=0.9, loc="best")
    ax.grid(True, axis="y", ls=":", lw=0.4, alpha=0.5)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    title = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--title=")),
                 "Fitted KRegCCD parameters")
    muceil = float(next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--muceil=")), 0.05))
    mufloor = float(next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mufloor=")), 1e-4))
    if len(args) < 2:
        sys.exit('usage: plot_params.py <out.pdf> <params.tsv> [--title="..."] [--muceil=0.05] [--mufloor=1e-4]')
    out, path = args[0], args[1]

    rows = np.atleast_1d(np.genfromtxt(path, delimiter="\t", names=True))
    sizes = sorted(np.unique(rows["sampleSize"]))

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    alpha_varies = np.ptp(rows["alpha"]) > 1e-9
    if alpha_varies:
        plot_scatter(ax, rows, sizes, muceil, mufloor)
        mode = "(alpha, mu) scatter"
    else:
        plot_violin(ax, rows, sizes, muceil, mufloor, float(rows["alpha"][0]))
        mode = f"mu violins (alpha fixed at {float(rows['alpha'][0]):g})"
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out)
    print(f"wrote {out}  ({len(rows)} fits across {len(sizes)} sample sizes; {mode})")


if __name__ == "__main__":
    main()
