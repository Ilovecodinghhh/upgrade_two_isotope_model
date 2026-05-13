#!/usr/bin/env python3
"""
Phase 5: Main publication figures for the Hemispheric Divergence experiment.

Fig 1: 6-panel hemispheric source time series
Fig 2: 1-box vs 2-box aliasing comparison
Fig 3: Robustness forest plot + exchange rate sensitivity
"""

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
from matplotlib.patches import FancyArrowPatch

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    hemi_npz = np.load(RESULTS_DIR / "twobox_hemi" / "all_iterations.npz")
    hemi_df = pd.read_csv(RESULTS_DIR / "twobox_hemi" / "hemispheric_detail.csv")
    onebox_npz = np.load(RESULTS_DIR / "onebox_reference" / "all_iterations.npz")
    years = hemi_df['year'].values
    return years, hemi_npz, onebox_npz


def add_trend_line(ax, years, y_med, start=2007, color='k'):
    mask = years >= start
    slope, intercept, _, _, se = sp_stats.linregress(years[mask], y_med[mask])
    trend_y = slope * years[mask] + intercept
    ax.plot(years[mask], trend_y, '--', color=color, alpha=0.6, linewidth=1)
    return slope, se


def fig1_hemispheric_sources():
    """6-panel: NH/SH × FF/Mic/BB time series."""
    years, hemi_npz, onebox_npz = load_data()

    fig, axes = plt.subplots(3, 2, figsize=(10, 10), sharex=True)
    fig.suptitle('Hemispheric Source Attribution (2-box dual-isotope model)', fontsize=13, y=0.98)

    sources = ['FF', 'Mic', 'BB']
    colors = ['#E64B35', '#4DBBD5', '#7E6148']
    labels = ['Fossil Fuel', 'Microbial', 'Biomass Burning']

    for i, (src, color, label) in enumerate(zip(sources, colors, labels)):
        for j, (hemi, hemi_label) in enumerate(zip(['NH', 'SH'], ['Northern Hemisphere', 'Southern Hemisphere'])):
            ax = axes[i, j]
            key = f"{hemi}_{src}"
            arr = hemi_npz[key]

            med = np.median(arr, axis=1)
            p5 = np.percentile(arr, 5, axis=1)
            p95 = np.percentile(arr, 95, axis=1)
            p25 = np.percentile(arr, 25, axis=1)
            p75 = np.percentile(arr, 75, axis=1)

            ax.fill_between(years, p5, p95, alpha=0.15, color=color)
            ax.fill_between(years, p25, p75, alpha=0.3, color=color)
            ax.plot(years, med, color=color, linewidth=1.5)

            slope, se = add_trend_line(ax, years, med, color=color)
            ax.text(0.02, 0.95, f'{slope:+.2f} ± {se:.2f} Tg/yr²',
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

    fig.savefig(FIGURES_DIR / 'fig1_hemispheric_sources.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'fig1_hemispheric_sources.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved fig1_hemispheric_sources")


def fig2_aliasing():
    """Key figure: 2-box vs 1-box FF/Mic/BB comparison showing aliasing."""
    years, hemi_npz, onebox_npz = load_data()

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    fig.suptitle('Source Aliasing: 2-box vs 1-box Global Source Attribution', fontsize=13, y=1.02)

    src_names = ['FF', 'Mic', 'BB']
    colors_2box = ['#E64B35', '#4DBBD5', '#7E6148']
    labels_full = ['Fossil Fuel', 'Microbial', 'Biomass Burning']

    for i, (src, color, label) in enumerate(zip(src_names, colors_2box, labels_full)):
        ax = axes[i]

        # 2-box global
        if src in ['FF', 'Mic', 'BB']:
            global_2box = hemi_npz[f'NH_{src}'] + hemi_npz[f'SH_{src}']
        med_2box = np.median(global_2box, axis=1)
        p5_2box = np.percentile(global_2box, 5, axis=1)
        p95_2box = np.percentile(global_2box, 95, axis=1)

        # 1-box
        arr_1box = onebox_npz[src]
        med_1box = np.median(arr_1box, axis=1)
        p5_1box = np.percentile(arr_1box, 5, axis=1)
        p95_1box = np.percentile(arr_1box, 95, axis=1)

        # Plot 2-box
        ax.fill_between(years, p5_2box, p95_2box, alpha=0.2, color=color)
        ax.plot(years, med_2box, color=color, linewidth=2, label='2-box global')

        # Plot 1-box
        ax.fill_between(years, p5_1box, p95_1box, alpha=0.15, color='gray')
        ax.plot(years, med_1box, color='gray', linewidth=2, linestyle='--', label='1-box')

        # Trends
        slope_2b, se_2b = add_trend_line(ax, years, med_2box, color=color)
        slope_1b, se_1b = add_trend_line(ax, years, med_1box, color='gray')

        ax.text(0.02, 0.95, f'2-box: {slope_2b:+.2f} Tg/yr²\n1-box: {slope_1b:+.2f} Tg/yr²',
                transform=ax.transAxes, fontsize=9, va='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        ax.set_title(label, fontsize=11)
        ax.set_xlabel('Year', fontsize=9)
        if i == 0:
            ax.set_ylabel('Emissions (Tg CH₄/yr)', fontsize=9)
        ax.legend(fontsize=8, loc='lower right')
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.2)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'fig2_source_aliasing.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'fig2_source_aliasing.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved fig2_source_aliasing")


def fig3_robustness():
    """Forest plot + exchange rate sensitivity."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Robustness of Hemispheric Divergence', fontsize=13, y=1.02)

    # Panel A: Forest plot from robustness table
    rob = pd.read_csv(RESULTS_DIR / "robustness_table.csv")

    y_pos = np.arange(len(rob))
    ax1.barh(y_pos, rob['NH_FF_slope'], xerr=[rob['NH_FF_slope'] - rob['NH_FF_p5'],
             rob['NH_FF_p95'] - rob['NH_FF_slope']], color='#E64B35', alpha=0.7,
             capsize=3, height=0.6)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(rob['config'], fontsize=9)
    ax1.axvline(0, color='k', linewidth=0.8, linestyle='-')
    ax1.set_xlabel('NH FF Trend (Tg/yr²)', fontsize=10)
    ax1.set_title('(a) NH FF Trend Across Configurations', fontsize=11)

    # Mark significant (>60% positive)
    for i, row in rob.iterrows():
        marker = '✓' if row['pattern_holds'] else '✗'
        ax1.text(row['NH_FF_p95'] + 0.1, i, marker, fontsize=10, va='center')

    ax1.tick_params(labelsize=8)
    ax1.grid(alpha=0.2, axis='x')

    # Panel B: Exchange rate sensitivity
    exch = pd.read_csv(RESULTS_DIR / "exchange_rate_sensitivity.csv")

    ax2.plot(exch['tau_ex'], exch['NH_FF_slope'], 'o-', color='#E64B35',
             linewidth=2, markersize=6, label='NH FF')
    ax2.fill_between(exch['tau_ex'], exch['NH_FF_p5'], exch['NH_FF_p95'],
                     alpha=0.2, color='#E64B35')
    ax2.plot(exch['tau_ex'], exch['SH_Mic_slope'], 's-', color='#4DBBD5',
             linewidth=2, markersize=6, label='SH Mic')

    ax2.axhline(0, color='k', linewidth=0.8, linestyle='-')
    ax2.axvspan(0.9, 1.1, alpha=0.1, color='green', label='Literature τ_ex range')
    ax2.set_xlabel('τ_ex (yr)', fontsize=10)
    ax2.set_ylabel('Trend (Tg/yr²)', fontsize=10)
    ax2.set_title('(b) Exchange Rate Sensitivity', fontsize=11)
    ax2.legend(fontsize=9)
    ax2.tick_params(labelsize=8)
    ax2.grid(alpha=0.2)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'fig3_robustness.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'fig3_robustness.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved fig3_robustness")


def fig4_reconciliation_schematic():
    """Literature reconciliation comparison."""
    fig, ax = plt.subplots(figsize=(8, 5))

    # Bar chart: FF trends across studies
    studies = ['Riddell-Young\n(1-box, δ¹³C+δD)', 'He 2026\n(TROPOMI)',
               'This work\n(1-box)', 'This work\n(2-box global)',
               'This work\n(2-box NH)', 'Basu 2022\n(3D, δ¹³C)']
    trends = [-0.3, -0.5, -0.82, 2.10, 1.04, 1.7]
    errors = [1.5*1.65, 2.0*1.65, (3.16+1.19)/2, (4.22+0.49)/2, (3.14+1.10)/2, 1.0*1.65]
    colors = ['gray', 'gray', '#4DBBD5', '#E64B35', '#E64B35', '#00A087']

    bars = ax.barh(studies, trends, xerr=errors, color=colors, alpha=0.7,
                   capsize=4, height=0.6, edgecolor='white')

    ax.axvline(0, color='k', linewidth=1)
    ax.set_xlabel('Post-2007 FF Trend (Tg/yr²)', fontsize=11)
    ax.set_title('Reconciliation: FF Emission Trends Across Methods', fontsize=12)

    # Add annotation
    ax.annotate('1-box models\n(FF stable/declining)', xy=(-1.5, 1), fontsize=8,
                ha='center', color='gray', style='italic')
    ax.annotate('Spatially-resolved\n(FF increasing)', xy=(2.0, 4.5), fontsize=8,
                ha='center', color='#00A087', style='italic')

    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.2, axis='x')
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / 'fig4_reconciliation.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'fig4_reconciliation.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved fig4_reconciliation")


def main():
    print("=" * 70)
    print("PHASE 5: Generating Figures")
    print("=" * 70)

    fig1_hemispheric_sources()
    fig2_aliasing()
    fig3_robustness()
    fig4_reconciliation_schematic()

    print(f"\n  All figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    main()
