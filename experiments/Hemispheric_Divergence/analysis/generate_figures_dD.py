#!/usr/bin/env python3
"""Figures comparing the 5 δD improvement approaches."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "dD_improvements"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

def main():
    import json
    with open(RESULTS_DIR / "comparison_summary.json") as f:
        data = json.load(f)

    labels = [d['label'].replace('Baseline (v3)', 'v3 Baseline') for d in data]
    short = ['v3\nBaseline', 'A\nSrc-water\nMic δD', 'B\nEDGAR\nFF δD', 
             'C\nC3/C4\nBB δD', 'D\nBayesian\nA+B+C', 'E\nδD gradient\nconstraint']

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle('δD Improvement Approaches: Impact on Source Attribution', fontsize=13, y=0.99)

    x = np.arange(len(data))

    # 1. FF level
    ax = axes[0, 0]
    ff = [d['ff_2010'] for d in data]
    ax.bar(x, ff, color='#E64B35', alpha=0.7)
    ax.axhline(110, color='k', ls=':', label='EDGAR')
    ax.set_ylabel('FF (Tg/yr)')
    ax.set_title('FF Emissions (2010)')
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=7)
    ax.legend(fontsize=8)

    # 2. NH FF share
    ax = axes[0, 1]
    nh = [d['nh_ff_share']*100 for d in data]
    ax.bar(x, nh, color='#3C5488', alpha=0.7)
    ax.axhline(72, color='k', ls=':', label='EDGAR 72%')
    ax.set_ylabel('NH share (%)')
    ax.set_title('NH Share of FF')
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=7)
    ax.legend(fontsize=8)

    # 3. FF trend
    ax = axes[0, 2]
    ft = [d['ff_trend'] for d in data]
    lo = [d['ff_trend_90ci'][0] for d in data]
    hi = [d['ff_trend_90ci'][1] for d in data]
    err = [[f-l for f,l in zip(ft,lo)], [h-f for f,h in zip(ft,hi)]]
    ax.barh(x, ft, xerr=err, color='#E64B35', alpha=0.7, capsize=3, height=0.6)
    ax.axvline(0, color='k', lw=1)
    ax.axvline(2.1, color='green', ls=':', label='EDGAR +2.1')
    ax.set_xlabel('Tg/yr²')
    ax.set_title('FF Trend (2007-2020)')
    ax.set_yticks(x); ax.set_yticklabels(short, fontsize=7)
    ax.legend(fontsize=8)

    # 4. BB level
    ax = axes[1, 0]
    bb = [d['bb_2010'] for d in data]
    ax.bar(x, bb, color='#7E6148', alpha=0.7)
    ax.axhline(30, color='k', ls=':', label='GFED ~30')
    ax.set_ylabel('BB (Tg/yr)')
    ax.set_title('BB Emissions (2010)')
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=7)
    ax.legend(fontsize=8)

    # 5. FF 90% CI width (uncertainty reduction)
    ax = axes[1, 1]
    ci = [d['ff_90ci_width'] for d in data]
    colors = ['gray' if c >= ci[0] else '#00A087' for c in ci]
    ax.bar(x, ci, color=colors, alpha=0.7)
    ax.set_ylabel('90% CI width (Tg/yr)')
    ax.set_title('FF Uncertainty')
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=7)

    # 6. Source fractions
    ax = axes[1, 2]
    ff_pct = [d['ff_pct'] for d in data]
    mic_pct = [100 - d['ff_pct'] - d['bb_pct'] for d in data]
    bb_pct = [d['bb_pct'] for d in data]
    ax.bar(x, ff_pct, color='#E64B35', alpha=0.7, label='FF')
    ax.bar(x, mic_pct, bottom=ff_pct, color='#4DBBD5', alpha=0.7, label='Mic')
    ax.bar(x, bb_pct, bottom=[f+m for f,m in zip(ff_pct, mic_pct)], color='#7E6148', alpha=0.7, label='BB')
    ax.set_ylabel('%')
    ax.set_title('Source Fractions (2010)')
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=7)
    ax.legend(fontsize=7, loc='lower right')

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIGURES_DIR / 'fig_dD_approaches_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved fig_dD_approaches_comparison.png")

if __name__ == "__main__":
    main()
