#!/usr/bin/env python3
"""Plot a KRegCCD mu-calibration sweep: mean held-out log-probability per tree and the model's
closed-form P(novel-clade tree) against the escape probability mu (log x-axis), with the empirical
held-out novel-clade rate as a reference line.

Reads the TSV written by ccd.experiments.regularisation.RSV2ThreeModelTable
(columns: mu, pnovel, empirical, meanLogP).

Usage: plot_musweep.py <out.pdf> <sweep.tsv> [--cvmu=0.0046] [--title="..."]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE = "#1E5AC5"   # log-prob
RED = "#C51E3A"    # P(novel)
GREY = "0.45"

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    title = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--title=")), None)
    cvmu = next((float(a.split("=", 1)[1]) for a in sys.argv[1:] if a.startswith("--cvmu=")), None)
    if len(args) < 2:
        sys.exit('usage: plot_musweep.py <out.pdf> <sweep.tsv> [--cvmu=0.0046] [--title="..."]')
    out, tsv = args[0], args[1]

    d = np.genfromtxt(tsv, delimiter="\t", names=True)
    mu, pnov, meanLogP = d["mu"], d["pnovel"], d["meanLogP"]
    emp = float(d["empirical"][0])
    optMu = mu[int(np.argmax(meanLogP))]

    fig, axL = plt.subplots(figsize=(6.4, 4.2))
    axR = axL.twinx()

    # left axis: held-out mean log-prob / tree
    lL, = axL.plot(mu, meanLogP, color=BLUE, lw=2.2, label="mean held-out $\\log P$/tree")
    axL.set_xscale("log")
    axL.set_xlabel("escape probability  $\\mu$")
    axL.set_ylabel("mean held-out $\\log P$ / tree", color=BLUE)
    axL.tick_params(axis="y", labelcolor=BLUE)
    axL.axvline(optMu, color=BLUE, ls=":", lw=1.3)
    axL.annotate(f"log $P$ optimum\n$\\mu$ = {optMu:.4f}", xy=(optMu, meanLogP.max()),
                 xytext=(6, -2), textcoords="offset points", fontsize=8, color=BLUE, va="top")

    # right axis: model P(novel) and empirical reference
    lR, = axR.plot(mu, pnov, color=RED, lw=2.2, label="model $P(\\mathrm{novel})$")
    lE = axR.axhline(emp, color=GREY, ls="--", lw=1.5, label=f"empirical rate = {emp:.3f}")
    axR.set_ylabel("$P(\\mathrm{novel\\ clade\\ tree})$", color=RED)
    axR.tick_params(axis="y", labelcolor=RED)
    axR.set_ylim(0, max(pnov) * 1.12)

    handles = [lL, lR, lE]
    if cvmu is not None:
        lC = axL.axvline(cvmu, color="0.2", ls="-.", lw=1.3, label=f"CV-fitted $\\mu$ = {cvmu:.4f}")
        handles.append(lC)

    axL.set_xlim(mu.min(), mu.max())
    axL.legend(handles=handles, loc="lower center", fontsize=8.5, framealpha=0.9, ncol=1)
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97) if title else (0, 0, 1, 1))
    fig.savefig(out)
    print(f"wrote {out}  (optimum mu={optMu:.5f}, empirical={emp:.3f})")


if __name__ == "__main__":
    main()
