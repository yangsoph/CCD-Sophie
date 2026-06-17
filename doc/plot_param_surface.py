#!/usr/bin/env python3
"""Plot the KRegCCD (alpha, mu) held-out log-probability surface.

Reads the TSV written by ccd.experiments.KRegCCDParameterSurface (columns:
alpha, mu, heldout_total, heldout_per_tree) and renders a heatmap with overlaid
contours on log-log axes, marking the argmax. Both parameters are log-spaced, so
the axes are log-scaled and the cells are drawn at the grid points.

Usage: plot_param_surface.py <surface.tsv> [out.png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: plot_param_surface.py <surface.tsv> [out.png]")
    path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "kregccd-param-surface.png"

    rows = np.genfromtxt(path, delimiter="\t", names=True)
    alpha = np.unique(rows["alpha"])
    mu = np.unique(rows["mu"])
    # z[i, j] = held-out per-tree logP at alpha[i], mu[j]
    z = np.full((len(alpha), len(mu)), np.nan)
    ai = {a: i for i, a in enumerate(alpha)}
    mi = {m: j for j, m in enumerate(mu)}
    for r in rows:
        z[ai[r["alpha"]], mi[r["mu"]]] = r["heldout_per_tree"]

    # argmax
    k = np.nanargmax(z)
    bi, bj = np.unravel_index(k, z.shape)
    a_star, mu_star, z_star = alpha[bi], mu[bj], z[bi, bj]

    # log-spaced cell edges so pcolormesh places each grid point at a cell centre
    def edges(v):
        lv = np.log(v)
        e = np.empty(len(v) + 1)
        e[1:-1] = (lv[:-1] + lv[1:]) / 2
        e[0] = lv[0] - (lv[1] - lv[0]) / 2
        e[-1] = lv[-1] + (lv[-1] - lv[-2]) / 2
        return np.exp(e)

    A, M = np.meshgrid(edges(alpha), edges(mu), indexing="ij")

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    # clip the colour floor so the deep low-alpha valley does not wash out the ridge
    vmax = np.nanmax(z)
    vmin = np.nanpercentile(z, 5)
    pc = ax.pcolormesh(A, M, z, cmap="viridis", vmin=vmin, vmax=vmax, shading="flat")

    Ac, Mc = np.meshgrid(alpha, mu, indexing="ij")
    levels = np.round(np.linspace(vmin, vmax, 9), 2)
    cs = ax.contour(Ac, Mc, z, levels=levels, colors="white", linewidths=0.6, alpha=0.7)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.2f")

    ax.plot(a_star, mu_star, "*", color="red", markersize=18, markeredgecolor="white",
            markeredgewidth=0.8,
            label=f"argmax: $\\alpha$={a_star:.3f}, $\\mu$={mu_star:.4f}\n{z_star:.3f} nat/tree")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\alpha$ (backbone smoothing)")
    ax.set_ylabel(r"$\mu$ (escape probability)")
    ax.set_title("KRegCCD held-out log-probability surface\n(5-fold CV, RSV2, 1024 thinned trees)")
    cb = fig.colorbar(pc, ax=ax)
    cb.set_label("held-out log-probability (nat/tree)")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"argmax alpha={a_star:.5f} mu={mu_star:.5f} -> {z_star:.4f} nat/tree")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
