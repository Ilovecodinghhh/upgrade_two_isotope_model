#!/usr/bin/env python3
"""v3 comparison figures: v1 vs v2 vs v3 (fraction vs delta space)."""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"


def main():
    # Load all three versions
    v1 = np.load(RESULTS_DIR / "twobox_hemi" / "all_iterations.npz")
    v1_1b = np.load(RESULTS_DIR / "onebox_reference" / "all_iterations.npz")
    v2 = np.load(RESULTS_DIR / "v2_improved" / "twobox_v2.npz")
    v2_1b = np.load(RESULTS_DIR / "v2_improved" / "onebox_v2.npz")
    v3 = np.load(RESULTS_DIR / "v3_delta_space" / "twobox_v3.npz")
    v3_1b = np.load(RESULTS_DIR / "v3_delta_space" / "onebox_v3.npz")

    import pandas as pd
    v1_years = pd.read_csv(RESULTS_DIR / "twobox_hemi" / "hemispheric_detail.csv")['year'].values
    v2_years = v2['years']
    v3_years = v3['years']

    fig, axes = plt.subplots(3, 3, figsize=(14, 10))
    fig.suptitle('Model Evolution: v1 (fraction, bad W) → v2 (fraction, good W) → v3 (delta-space)',
                 fontsize=12, y=0.99)

    versions = [
        ('v1: W=diag(100,1,0.5)\nfraction space', v1, v1_1b, v1_years),
        ('v2: uncertainty W\nfraction space', v2, v2_1b, v2_years),
        ('v3: delta-space\ncond=13.5', v3, v3_1b, v3_years),
    ]
    sources = ['FF', 'Mic', 'BB']
    colors = ['#E64B35', '#4DBBD5', '#7E6148']

    for row, (vlabel, v2box, v1box, yrs) in enumerate(versions):
        for col, (src, color) in enumerate(zip(sources, colors)):
            ax = axes[row, col]

            # 2-box global
            global_arr = v2box[f'NH_{src}'] + v2box[f'SH_{src}']
            med = np.nanmedian(global_arr, axis=1)
            p5 = np.nanpercentile(global_arr, 5, axis=1)
            p95 = np.nanpercentile(global_arr, 95, axis=1)
            ax.fill_between(yrs, p5, p95, alpha=0.2, color=color)
            ax.plot(yrs, med, color=color, linewidth=2, label='2-box')

            # 1-box
            med_1b = np.nanmedian(v1box[src], axis=1)
            ax.plot(yrs, med_1b, '--', color='gray', linewidth=1.5, label='1-box')

            if row == 0:
                ax.set_title(src, fontsize=11, fontweight='bold')
            if col == 0:
                ax.set_ylabel(f'{vlabel}\nTg/yr', fontsize=8)
            if row == 2:
                ax.set_xlabel('Year', fontsize=9)
            ax.tick_params(labelsize=7)
            ax.grid(alpha=0.15)
            if row == 0 and col == 0:
                ax.legend(fontsize=7)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIGURES_DIR / 'fig_v1_v2_v3_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved fig_v1_v2_v3_comparison.png")

    # Also make a summary bar chart of FF(2010) across versions
    fig2, ax = plt.subplots(figsize=(8, 5))
    
    j = 11  # 2010 index
    versions_data = {
        'v1 2-box': np.nanmedian((v1['NH_FF'] + v1['SH_FF'])[j,:]),
        'v1 1-box': np.nanmedian(v1_1b['FF'][j,:]),
        'v2 2-box': np.nanmedian((v2['NH_FF'] + v2['SH_FF'])[j,:]),
        'v2 1-box': np.nanmedian(v2_1b['FF'][j,:]),
        'v3 2-box': np.nanmedian((v3['NH_FF'] + v3['SH_FF'])[j,:]),
        'v3 1-box': np.nanmedian(v3_1b['FF'][j,:]),
        'EDGAR': 110,
    }
    
    bars = ax.barh(list(versions_data.keys()), list(versions_data.values()),
                   color=['#E64B35','gray','#E64B35','gray','#E64B35','gray','#00A087'],
                   alpha=0.7, height=0.6)
    ax.axvline(110, color='k', linewidth=1, linestyle=':', label='EDGAR ~110')
    ax.set_xlabel('FF Emissions (Tg CH₄/yr)', fontsize=11)
    ax.set_title('FF Emission Level (2010) Across Model Versions', fontsize=12)
    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.2, axis='x')
    plt.tight_layout()
    fig2.savefig(FIGURES_DIR / 'fig_ff_level_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved fig_ff_level_comparison.png")


if __name__ == "__main__":
    main()
