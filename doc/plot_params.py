#!/usr/bin/env python3
"""Scatter the fitted KRegCCD (alpha, mu) parameters from the PIT sweep.

Reads the params TSV written by ccd.experiments.regularisation.KRegPITSweep
(columns: rep, sampleSize, alpha, mu) -- one row per KRegCCD fit -- and plots the
(alpha, mu) pairs on log-log axes, coloured by training sample size. The horizontal
dashed line marks the reliability ceiling mu = 0.05 (KRegCCD.MU_RELIABLE_MAX), the
upper bound of the optimiser's mu search; points sitting on it are fits whose held-out
likelihood was still increasing in mu at the cap.

Usage: plot_params.py <out.pdf> <params.tsv> [--title="..."] [--muceil=0.05]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = ["#3b6fb6", "#d1772f", "#4c9a52", "#9467bd"]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    title = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--title=")),
                 "Fitted KRegCCD parameters")
    muceil = float(next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--muceil=")), 0.05))
    mufloor = float(next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mufloor=")), 1e-3))
    if len(args) < 2:
        sys.exit('usage: plot_params.py <out.pdf> <params.tsv> [--title="..."] [--muceil=0.05]')
    out, path = args[0], args[1]

    rows = np.atleast_1d(np.genfromtxt(path, delimiter="\t", names=True))
    sizes = sorted(np.unique(rows["sampleSize"]))

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.axhline(muceil, ls="--", lw=0.9, color="0.4", zorder=1, label=f"$\\mu$ search ceiling = {muceil:g}")
    ax.axhline(mufloor, ls=":", lw=0.9, color="0.4", zorder=1, label=f"$\\mu$ search floor = {mufloor:g}")
    for i, size in enumerate(sizes):
        sel = rows[rows["sampleSize"] == size]
        # median marker per size, to show how the fit shifts with sample size
        ax.scatter(sel["alpha"], sel["mu"], s=22, alpha=0.6, color=COLORS[i % len(COLORS)],
                   edgecolor="white", linewidth=0.3, zorder=2,
                   label=f"n = {int(size)}  ({len(sel)} fits)")
        ax.scatter(np.median(sel["alpha"]), np.median(sel["mu"]), s=160, marker="*",
                   color=COLORS[i % len(COLORS)], edgecolor="black", linewidth=0.7, zorder=3)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\alpha$  (backbone smoothing)")
    ax.set_ylabel(r"$\mu$  (escape probability)")
    ax.set_title(title)
    ax.legend(fontsize=8, framealpha=0.9, loc="best")
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
    fig.tight_layout()
    fig.savefig(out)
    print(f"wrote {out}  ({len(rows)} (alpha, mu) pairs across {len(sizes)} sample sizes; stars = per-size medians)")


if __name__ == "__main__":
    main()
