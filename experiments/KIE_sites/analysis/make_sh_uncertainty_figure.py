#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Publication figure: true variance attribution of the SH-only alpha13C_OH constraint.

Matches the manuscript figure style in make_manuscript_figures.py (same rcParams,
bar color, axis conventions). This figure decomposes the ACTUAL Phase-6 Monte
Carlo for CGO+SPO by toggling each uncertainty group on/off; it is NOT the
hardcoded guess-sigma diagnostic used by the earlier phase13/figS10 panel.
See analysis/sh_variance_purepy.py for the driver.

Panel (a): isolated Var(alpha) per input group, ranked, as % of summed one-at-a-
           time variance (parts sum to the all-on total; interaction ~ 0).
Panel (b): the same variance rolled up into the three input families
           (observation / wetland / sink), so the family split is legible.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EXPT_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = EXPT_DIR / "figures" / "manuscript"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Values from analysis/sh_variance_purepy.py (120,000-draw SH decomposition).
# Reproduces the reported SH interval: median 1.0049, sd 0.0033.
# Each entry: (label, family, sd_alpha, pct_of_oat_variance).
# ---------------------------------------------------------------------------
GROUPS = [
    ("$\\delta^{13}$C amplitude", "obs",  0.00216, 41.8),
    ("$\\delta$D amplitude",      "obs",  0.00210, 39.7),
    ("Sink $\\rightarrow\\ \\alpha$ (f$_{OH}$, KIEs)", "sink", 0.00091, 7.4),
    ("$\\delta^{13}$C phase",     "obs",  0.00081, 5.9),
    ("Wetland $\\delta^{13}$C signature", "wet", 0.00071, 4.5),
    ("Wetland flux phasor",       "wet",  0.00024, 0.5),
    ("$\\delta$D phase",          "obs",  0.00014, 0.2),
    ("Wetland $\\delta$D signature", "wet", 0.00007, 0.03),
]

# Family palette — kept muted/publication-like, consistent with fig3 hues.
FAM_COLOR = {"obs": "#4c78a8", "wet": "#59a14f", "sink": "#e0a530"}
FAM_LABEL = {"obs": "Observation (harmonic fit)",
             "wet": "Wetland source correction",
             "sink": "Sink $\\rightarrow\\ \\alpha$ conversion"}

# Family rollup of the true decomposition (sums of the GROUPS above):
#   Observation = d13C amp + dD amp + d13C phase + dD phase = 87.6%
#   Sink        = 7.4%
#   Wetland     = d13C sig + flux phasor + dD sig = 5.0%
FAMILY_ORDER = ["obs", "sink", "wet"]
FAMILY_PCT = {
    "obs": 41.8 + 39.7 + 5.9 + 0.2,   # 87.6
    "sink": 7.4,
    "wet": 4.5 + 0.5 + 0.03,          # 5.0
}


def set_common_style():
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    })


def make_figure():
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.3), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.35, 1.0]})

    # ---- Panel (a): ranked per-group attribution -------------------------
    ax = axes[0]
    labels = [g[0] for g in GROUPS]
    fams = [g[1] for g in GROUPS]
    sds = [g[2] for g in GROUPS]
    pcts = [g[3] for g in GROUPS]
    colors = [FAM_COLOR[f] for f in fams]
    y = np.arange(len(GROUPS))

    ax.barh(y, pcts, color=colors, edgecolor="white", linewidth=0.6, zorder=3)
    for yi, pct, sd in zip(y, pcts, sds):
        ax.text(pct + 0.9, yi, f"{pct:.1f}%   sd={sd:.5f}",
                va="center", ha="left", fontsize=8, color="#222222")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of one-at-a-time Var($\\alpha^{13}$C$_{OH}$)  (%)")
    ax.set_title("(a) Variance attribution by input group")
    ax.set_xlim(0, 52)
    ax.grid(axis="x", alpha=0.25)

    handles = [plt.Rectangle((0, 0), 1, 1, color=FAM_COLOR[f]) for f in FAMILY_ORDER]
    ax.legend(handles, [FAM_LABEL[f] for f in FAMILY_ORDER],
              loc="lower right", frameon=False, fontsize=8)

    # ---- Panel (b): family rollup ----------------------------------------
    ax = axes[1]
    yb = np.arange(len(FAMILY_ORDER))
    fpct = [FAMILY_PCT[f] for f in FAMILY_ORDER]
    fcol = [FAM_COLOR[f] for f in FAMILY_ORDER]
    ax.barh(yb, fpct, color=fcol, edgecolor="white", linewidth=0.6, zorder=3)
    for yi, pct in zip(yb, fpct):
        ax.text(pct + 1.5, yi, f"{pct:.1f}%", va="center", ha="left",
                fontsize=9, color="#222222")
    ax.set_yticks(yb)
    ax.set_yticklabels([FAM_LABEL[f].replace(" (harmonic fit)", "\n(harmonic fit)")
                        .replace(" source correction", "\nsource correction")
                        .replace("$\\rightarrow\\ \\alpha$ conversion", "$\\rightarrow\\ \\alpha$\nconversion")
                        for f in FAMILY_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of variance  (%)")
    ax.set_title("(b) Rolled up by input family")
    ax.set_xlim(0, 100)
    ax.grid(axis="x", alpha=0.25)

    out = OUT_DIR / "figS11_sh_true_uncertainty_attribution.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    set_common_style()
    out = make_figure()
    print(out.relative_to(EXPT_DIR))


if __name__ == "__main__":
    main()
