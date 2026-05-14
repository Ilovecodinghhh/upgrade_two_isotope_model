#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_graphical_abstract.py — TOC / graphical abstract figure
============================================================

Single-panel threshold curve with annotations.
Designed for journal TOC/graphical abstract (85mm × 85mm).
Reads data from phase6 deep dive results.
"""

import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIG_DIR = Path(__file__).resolve().parent


def main():
    # Load phase6 fine-grid sweep
    deep_path = RESULTS_DIR / "phase6_deep_dive" / "deep_dive_results.json"
    with open(deep_path) as f:
        deep = json.load(f)

    fg = deep['A_fine_grid']
    sweep = fg['sweep']
    mults = sorted(sweep.keys(), key=float)
    sigma_vals = [8.25 * float(m) for m in mults]
    ci_vals = [sweep[m]['ci'] for m in mults]

    # d13C-only reference
    d13c_ref = fg['ref_ci']
    thresh_sigma = fg['crossover_0pct']['sigma_permil']

    # --- Figure ---
    plt.rcParams.update({
        'font.size': 9,
        'font.family': 'sans-serif',
        'axes.linewidth': 0.8,
    })

    fig, ax = plt.subplots(figsize=(85/25.4, 85/25.4), dpi=300)

    # Fill regions
    ax.fill_between(sigma_vals, 0, d13c_ref, alpha=0.06, color='#2171b5')
    ax.fill_between(sigma_vals, d13c_ref, max(ci_vals) * 1.1, alpha=0.06, color='#cb181d')

    # Threshold curve
    ax.plot(sigma_vals, ci_vals, 'o-', color='#2171b5', linewidth=2, markersize=4, zorder=3)

    # Reference line
    ax.axhline(d13c_ref, color='#969696', linewidth=1.5, linestyle='--', label='δ¹³C-only')

    # Crossover
    ax.axvline(thresh_sigma, color='#cb181d', linewidth=1, linestyle=':', alpha=0.7)
    ax.annotate(f'Threshold\nσ ≈ {thresh_sigma:.0f}‰', xy=(thresh_sigma, d13c_ref),
                xytext=(thresh_sigma + 25, d13c_ref * 0.65),
                fontsize=8, ha='center', color='#cb181d',
                arrowprops=dict(arrowstyle='->', color='#cb181d', lw=1.2))

    # Current precision
    current_ci = ci_vals[1] if len(ci_vals) > 1 else ci_vals[0]
    ax.annotate('Current\nprecision', xy=(8.25, current_ci),
                xytext=(30, current_ci * 0.7),
                fontsize=8, ha='center', color='#2171b5',
                arrowprops=dict(arrowstyle='->', color='#2171b5', lw=1.2))

    # Labels
    ax.text(15, d13c_ref * 0.85, 'δD helps', fontsize=8, color='#2171b5', fontstyle='italic')
    ax.text(thresh_sigma + 15, d13c_ref * 1.1, 'δD hurts', fontsize=8, color='#cb181d', fontstyle='italic')

    ax.set_xlabel('σ(Mic δD) [‰]', fontsize=9)
    ax.set_ylabel('90% CI on FF emissions [Tg yr⁻¹]', fontsize=9)
    ax.set_xlim(0, max(sigma_vals) * 1.05)
    ax.set_ylim(0, max(ci_vals) * 1.1)
    ax.legend(fontsize=7, loc='upper left')
    ax.tick_params(labelsize=8)

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig_graphical_abstract.png', dpi=300, bbox_inches='tight')
    plt.savefig(FIG_DIR / 'fig_graphical_abstract.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_DIR}/fig_graphical_abstract.png")
    print(f"Saved: {FIG_DIR}/fig_graphical_abstract.pdf")


if __name__ == "__main__":
    main()
