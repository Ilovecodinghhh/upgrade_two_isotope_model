#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_threshold.py — The key figure for Title 1
===============================================

Panel A: "Information gain vs. source-signature uncertainty"
Panel B: "FF emissions CI width — δ¹³C-only vs. dual at different uncertainties"
"""

import sys
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    # Load threshold results
    with open(RESULTS_DIR / "phase3_threshold" / "threshold_results.json") as f:
        data = json.load(f)
    
    multipliers = data['multipliers']
    baseline_ci = data['baseline_d13C_CI']
    improvements = data['improvements']
    
    # Extract data for plotting
    sigma_vals = [improvements[str(m)]['mic_dD_sigma_permil'] for m in multipliers]
    dual_ci_vals = [improvements[str(m)]['dual_CI_width'] for m in multipliers]
    improvement_vals = [improvements[str(m)]['improvement_pct'] for m in multipliers]
    
    # Load DFS results
    with open(RESULTS_DIR / "phase2_dfs" / "dfs_results.json") as f:
        dfs_data = json.load(f)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)
    
    # === Panel A: % Improvement vs. σ(Mic δD) ===
    ax = axes[0]
    ax.plot(sigma_vals, improvement_vals, 'o-', lw=2.5, color='#2171b5', 
            markersize=8, zorder=5)
    ax.axhline(0, color='black', lw=0.8, ls='-')
    ax.axhline(10, color='gray', lw=1.2, ls='--', alpha=0.7, label='10% threshold')
    ax.axvspan(95, 140, alpha=0.15, color='red', label='Thanwerdas (2024) range')
    ax.axvline(25, color='#e34a33', lw=2, ls=':', alpha=0.8, label='Critical threshold (~25‰)')
    
    ax.fill_between(sigma_vals, improvement_vals, 0, 
                     where=[v >= 0 for v in improvement_vals],
                     alpha=0.1, color='green')
    ax.fill_between(sigma_vals, improvement_vals, 0,
                     where=[v < 0 for v in improvement_vals],
                     alpha=0.1, color='red')
    
    ax.set_xlabel('Microbial δD source-signature uncertainty (‰, 1σ)', fontsize=10)
    ax.set_ylabel('Constraint improvement from adding δD (%)', fontsize=10)
    ax.set_title('A. When does δD help?', fontsize=11, fontweight='bold')
    ax.set_xlim(0, 140)
    ax.set_ylim(-180, 60)
    ax.legend(loc='lower left', fontsize=8)
    ax.grid(alpha=0.3)
    
    # Annotate regions
    ax.text(10, 40, 'δD helps\n(current data)', fontsize=9, color='green', 
            ha='center', fontstyle='italic')
    ax.text(80, -100, 'δD hurts\n(Thanwerdas regime)', fontsize=9, color='red',
            ha='center', fontstyle='italic')
    
    # === Panel B: CI width comparison ===
    ax = axes[1]
    ax.plot(sigma_vals, dual_ci_vals, 'o-', lw=2.5, color='#2171b5',
            markersize=8, label='Dual-isotope (δ¹³C + δD)', zorder=5)
    ax.axhline(baseline_ci, color='#e34a33', lw=2.5, ls='--',
               label=f'δ¹³C-only ({baseline_ci:.0f} Tg/yr)')
    ax.axvspan(95, 140, alpha=0.15, color='red')
    ax.axvline(25, color='#e34a33', lw=2, ls=':', alpha=0.8)
    
    ax.set_xlabel('Microbial δD source-signature uncertainty (‰, 1σ)', fontsize=10)
    ax.set_ylabel('FF emissions 90% CI width (Tg CH₄ yr⁻¹)', fontsize=10)
    ax.set_title('B. FF emission constraint degradation', fontsize=11, fontweight='bold')
    ax.set_xlim(0, 140)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(alpha=0.3)
    
    # Annotate crossing
    ax.annotate('Crossover: δD\nbecomes counterproductive',
                xy=(25, baseline_ci), xytext=(50, 60),
                fontsize=8, ha='center',
                arrowprops=dict(arrowstyle='->', color='#e34a33'),
                color='#e34a33')
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig_dD_threshold.png', dpi=300, bbox_inches='tight')
    plt.savefig(FIG_DIR / 'fig_dD_threshold.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_DIR}/fig_dD_threshold.png")
    print(f"       {FIG_DIR}/fig_dD_threshold.pdf")


if __name__ == "__main__":
    main()
