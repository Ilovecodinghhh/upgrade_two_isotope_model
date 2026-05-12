#!/usr/bin/env python3
"""Phase 10 — 3-panel variance decomposition + KIE spread figure."""
import sys
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
ANALYSIS_DIR = Path(__file__).resolve().parent.parent / "analysis"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(ANALYSIS_DIR))
from common import load_data, smooth_5yr
from core import run_2box_flex, OUT_DIR
from variance_decomposition import run_2box

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    data = load_data(REPO_ROOT, two_box=True)
    N = 400; SEED = 42

    # Load bootstrap results
    with open(OUT_DIR / "phase9_bootstrap.json") as f:
        boot = json.load(f)

    # Run models for time series
    print("Running models for time series...")
    FF_c13 = run_2box(data, "d13C_only", N, SEED, use_real_hemi_dD=False)
    FF_off = run_2box(data, "dual", N, SEED, use_real_hemi_dD=False)
    FF_real = run_2box_flex(data, N, SEED)

    years = data.model_years

    # ========== FIGURE ==========
    fig, axes = plt.subplots(1, 3, figsize=(180/25.4, 65/25.4), dpi=300)
    plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9})

    colors = {'d13C': '#2166ac', 'offset': '#b2182b', 'real': '#1b7837'}
    comp_colors = ['#fee08b', '#abdda4', '#66c2a5', '#d9d9d9']  # KIE, Sig, τ, Resid

    # Panel A: Stacked variance bar chart
    ax = axes[0]
    configs = ['d13C_only', 'dual_offset', 'dual_real_hemi']
    labels = ['δ¹³C-only', 'Dual\n(offset)', 'Dual\n(real hemi)']
    x = np.arange(3)

    bottoms = np.zeros(3)
    for i, (comp, clr, cname) in enumerate([
        ('kie_pct', comp_colors[0], 'KIE'),
        ('sig_pct', comp_colors[1], 'Source sig'),
        ('tau_pct', comp_colors[2], 'Lifetime'),
        ('resid_pct', comp_colors[3], 'Residual'),
    ]):
        vals = [boot[c][comp][0] for c in configs]
        ax.bar(x, vals, bottom=bottoms, width=0.6, color=clr, edgecolor='k',
               linewidth=0.5, label=cname)
        bottoms += vals

    for i, c in enumerate(configs):
        sig = boot[c]['sigma_ff'][0]
        ax.text(i, 103, f'σ={sig:.1f}', ha='center', va='bottom', fontsize=7, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel('Variance contribution (%)')
    ax.set_ylim(0, 115)
    ax.legend(loc='upper right', fontsize=6, framealpha=0.9)
    ax.set_title('A. Variance decomposition', fontweight='bold', loc='left')

    # Panel B: FF time series
    ax = axes[1]
    for FF, clr, lbl in [
        (FF_c13, colors['d13C'], 'δ¹³C-only'),
        (FF_off, colors['offset'], 'Dual (offset)'),
        (FF_real, colors['real'], 'Dual (real hemi)'),
    ]:
        FF_s = smooth_5yr(FF)
        med = np.nanmedian(FF_s, axis=1)
        lo = np.nanpercentile(FF_s, 5, axis=1)
        hi = np.nanpercentile(FF_s, 95, axis=1)
        ax.plot(years, med, color=clr, linewidth=1.2, label=lbl)
        ax.fill_between(years, lo, hi, color=clr, alpha=0.15)

    ax.axvline(2007, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.set_xlabel('Year')
    ax.set_ylabel('FF emissions (Tg/yr)')
    ax.legend(loc='upper left', fontsize=6, framealpha=0.9)
    ax.set_title('B. FF time series (5-yr smoothed)', fontweight='bold', loc='left')

    # Panel C: KIE spread comparison
    ax = axes[2]
    spreads = [13.0, 0.7, 0.9, 0.8]
    spread_labels = ['Basu 2022\n(3D)', 'δ¹³C-only', 'Dual\n(offset)', 'Dual\n(real hemi)']
    bar_colors = ['#969696', colors['d13C'], colors['offset'], colors['real']]
    y = np.arange(4)
    ax.barh(y, spreads, color=bar_colors, edgecolor='k', linewidth=0.5, height=0.6)
    for i, v in enumerate(spreads):
        ax.text(v + 0.3, i, f'{v:.1f}', va='center', fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(spread_labels, fontsize=7)
    ax.set_xlabel('KIE spread (Tg/yr)')
    ax.set_xlim(0, 16)
    ax.set_title('C. KIE spread: Basu vs ours', fontweight='bold', loc='left')

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig_variance_v2.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIG_DIR / 'fig_variance_v2.pdf', bbox_inches='tight')
    print(f"Saved: {FIG_DIR / 'fig_variance_v2.png'}")
    print(f"Saved: {FIG_DIR / 'fig_variance_v2.pdf'}")
    plt.close()


if __name__ == "__main__":
    main()
