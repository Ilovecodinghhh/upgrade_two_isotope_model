#!/usr/bin/env python3
"""v2 figures: compare v1 (flawed) vs v2 (improved) results."""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def compute_trends(arr, years, start=2007, end_trim=1):
    end = years[-1] - end_trim
    mask = (years >= start) & (years <= end)
    yrs = years[mask]
    sub = arr[mask, :]
    return np.array([sp_stats.linregress(yrs, sub[:, k]).slope for k in range(sub.shape[1])])


def fig_v1_vs_v2():
    """Compare v1 and v2 model results side by side."""
    # v1 data
    v1_hemi = np.load(RESULTS_DIR / "twobox_hemi" / "all_iterations.npz")
    v1_1box = np.load(RESULTS_DIR / "onebox_reference" / "all_iterations.npz")
    v1_df = pd.read_csv(RESULTS_DIR / "twobox_hemi" / "hemispheric_detail.csv")
    v1_years = v1_df['year'].values

    # v2 data
    v2_2box = np.load(RESULTS_DIR / "v2_improved" / "twobox_v2.npz")
    v2_1box = np.load(RESULTS_DIR / "v2_improved" / "onebox_v2.npz")
    v2_years = v2_2box['years']

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle('Model Improvement: v1 (arbitrary weights) → v2 (uncertainty-based weights)', fontsize=13, y=0.98)

    sources = ['FF', 'Mic', 'BB']
    colors = ['#E64B35', '#4DBBD5', '#7E6148']
    labels = ['Fossil Fuel', 'Microbial', 'Biomass Burning']

    for i, (src, color, label) in enumerate(zip(sources, colors, labels)):
        # v1: top row
        ax_v1 = axes[0, i]
        global_v1 = v1_hemi[f'NH_{src}'] + v1_hemi[f'SH_{src}']
        med_v1 = np.median(global_v1, axis=1)
        p5_v1 = np.percentile(global_v1, 5, axis=1)
        p95_v1 = np.percentile(global_v1, 95, axis=1)
        ax_v1.fill_between(v1_years, p5_v1, p95_v1, alpha=0.2, color=color)
        ax_v1.plot(v1_years, med_v1, color=color, linewidth=2, label='2-box')

        onebox_v1 = v1_1box[src]
        med_1b_v1 = np.median(onebox_v1, axis=1)
        ax_v1.plot(v1_years, med_1b_v1, '--', color='gray', linewidth=2, label='1-box')

        ax_v1.set_title(f'{label}', fontsize=11)
        if i == 0:
            ax_v1.set_ylabel('v1 (W=diag(100,1,0.5))\nTg CH₄/yr', fontsize=9)
        ax_v1.legend(fontsize=7, loc='upper left')
        ax_v1.tick_params(labelsize=8)
        ax_v1.grid(alpha=0.2)

        # v2: bottom row
        ax_v2 = axes[1, i]
        global_v2 = v2_2box[f'NH_{src}'] + v2_2box[f'SH_{src}']
        med_v2 = np.median(global_v2, axis=1)
        p5_v2 = np.percentile(global_v2, 5, axis=1)
        p95_v2 = np.percentile(global_v2, 95, axis=1)
        ax_v2.fill_between(v2_years, p5_v2, p95_v2, alpha=0.2, color=color)
        ax_v2.plot(v2_years, med_v2, color=color, linewidth=2, label='2-box')

        onebox_v2 = v2_1box[src]
        med_1b_v2 = np.median(onebox_v2, axis=1)
        ax_v2.plot(v2_years, med_1b_v2, '--', color='gray', linewidth=2, label='1-box')

        if i == 0:
            ax_v2.set_ylabel('v2 (uncertainty-based W)\nTg CH₄/yr', fontsize=9)
        ax_v2.set_xlabel('Year', fontsize=9)
        ax_v2.legend(fontsize=7, loc='upper left')
        ax_v2.tick_params(labelsize=8)
        ax_v2.grid(alpha=0.2)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGURES_DIR / 'fig_v1_vs_v2_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved fig_v1_vs_v2_comparison")


def fig_v2_hemispheric():
    """v2 hemispheric source attribution."""
    v2 = np.load(RESULTS_DIR / "v2_improved" / "twobox_v2.npz")
    years = v2['years']

    fig, axes = plt.subplots(3, 2, figsize=(10, 10), sharex=True)
    fig.suptitle('v2 Hemispheric Source Attribution\n(uncertainty-based weights, literature IH gradient)', fontsize=12, y=0.99)

    sources = ['FF', 'Mic', 'BB']
    colors = ['#E64B35', '#4DBBD5', '#7E6148']
    labels_full = ['Fossil Fuel', 'Microbial', 'Biomass Burning']

    for i, (src, color, label) in enumerate(zip(sources, colors, labels_full)):
        for j, (hemi, hemi_label) in enumerate(zip(['NH', 'SH'], ['Northern Hemisphere', 'Southern Hemisphere'])):
            ax = axes[i, j]
            arr = v2[f'{hemi}_{src}']
            med = np.median(arr, axis=1)
            p5 = np.percentile(arr, 5, axis=1)
            p95 = np.percentile(arr, 95, axis=1)

            ax.fill_between(years, p5, p95, alpha=0.2, color=color)
            ax.plot(years, med, color=color, linewidth=1.5)

            # Trend (trim last year)
            slopes = compute_trends(arr, years)
            med_slope = np.median(slopes)
            pct_pos = np.mean(slopes > 0) * 100
            ax.text(0.02, 0.95, f'{med_slope:+.2f} Tg/yr²\n({pct_pos:.0f}% pos)',
                    transform=ax.transAxes, fontsize=8, va='top',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

            if i == 0:
                ax.set_title(hemi_label, fontsize=11)
            if j == 0:
                ax.set_ylabel(f'{label}\n(Tg CH₄/yr)', fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=0.2)

    axes[2, 0].set_xlabel('Year', fontsize=10)
    axes[2, 1].set_xlabel('Year', fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGURES_DIR / 'fig_v2_hemispheric.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved fig_v2_hemispheric")


def fig_v2_reconciliation():
    """Updated reconciliation bar chart with v2 results."""
    v2 = np.load(RESULTS_DIR / "v2_improved" / "twobox_v2.npz")
    v2_1box = np.load(RESULTS_DIR / "v2_improved" / "onebox_v2.npz")
    years = v2['years']

    # Compute slopes
    ff_2box = compute_trends(v2['NH_FF'] + v2['SH_FF'], years)
    ff_nh = compute_trends(v2['NH_FF'], years)
    ff_1box = compute_trends(v2_1box['FF'], years)

    fig, ax = plt.subplots(figsize=(8, 5))

    studies = ['Riddell-Young\n(1-box, δ¹³C+δD)', 'This work\n(1-box v2)',
               'This work\n(2-box v2 global)', 'This work\n(2-box v2 NH)',
               'Basu 2022\n(3D, δ¹³C)']
    trends = [-0.3, np.median(ff_1box), np.median(ff_2box), np.median(ff_nh), 1.7]
    errors = [1.5*1.65,
              (np.percentile(ff_1box,95)-np.percentile(ff_1box,5))/2,
              (np.percentile(ff_2box,95)-np.percentile(ff_2box,5))/2,
              (np.percentile(ff_nh,95)-np.percentile(ff_nh,5))/2,
              1.0*1.65]
    colors = ['gray', '#4DBBD5', '#E64B35', '#E64B35', '#00A087']

    ax.barh(studies, trends, xerr=errors, color=colors, alpha=0.7,
            capsize=4, height=0.6, edgecolor='white')
    ax.axvline(0, color='k', linewidth=1)
    ax.set_xlabel('Post-2007 FF Trend (Tg/yr²)', fontsize=11)
    ax.set_title('v2 Reconciliation: FF Emission Trends', fontsize=12)
    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.2, axis='x')
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / 'fig_v2_reconciliation.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved fig_v2_reconciliation")


def main():
    print("=" * 60)
    print("GENERATING v2 FIGURES")
    print("=" * 60)
    fig_v1_vs_v2()
    fig_v2_hemispheric()
    fig_v2_reconciliation()
    print("  Done!")


if __name__ == "__main__":
    main()
