#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_comprehensive.py — 4-panel Nature Communications-style figure
=================================================================

Panel A: DFS gain from δD (1-box vs 2-box)
Panel B: δD constraint improvement vs σ(Mic δD) — THE threshold plot
Panel C: Thanwerdas replication diagnostic
Panel D: Sensitivity analysis (KIE + lifetime) — threshold is robust
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
    # Load all results
    with open(RESULTS_DIR / "phase2_dfs" / "dfs_results.json") as f:
        dfs = json.load(f)
    
    with open(RESULTS_DIR / "phase3_threshold" / "threshold_results.json") as f:
        thresh = json.load(f)
    
    with open(RESULTS_DIR / "phase3b_thanwerdas" / "thanwerdas_comparison.json") as f:
        than = json.load(f)
    
    with open(RESULTS_DIR / "phase5_sensitivity" / "sensitivity_results.json") as f:
        sens = json.load(f)
    
    # ── Figure setup ──
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=150)
    plt.rcParams.update({'font.size': 9})
    
    c_blue = '#2171b5'
    c_red = '#e34a33'
    c_green = '#31a354'
    c_orange = '#fe9929'
    c_purple = '#756bb1'
    
    # ═══ Panel A: DFS ═══
    ax = axes[0, 0]
    models = ['1-box', '2-box']
    dfs_c13 = [dfs['onebox_d13C_only']['DFS'], dfs['twobox_d13C_only']['DFS']]
    dfs_dual = [dfs['onebox_dual']['DFS'], dfs['twobox_dual']['DFS']]
    delta_dfs = [dfs_dual[i] - dfs_c13[i] for i in range(2)]
    
    x = np.arange(2)
    w = 0.3
    bars1 = ax.bar(x - w/2, dfs_c13, w, label='δ¹³C only', color=c_orange, alpha=0.85)
    bars2 = ax.bar(x + w/2, dfs_dual, w, label='δ¹³C + δD', color=c_blue, alpha=0.85)
    
    for i in range(2):
        ax.annotate(f'ΔDFS = +{delta_dfs[i]:.2f}',
                     xy=(x[i] + w/2, dfs_dual[i]),
                     xytext=(x[i] + 0.45, dfs_dual[i] + 0.1),
                     fontsize=7.5, color=c_blue, fontweight='bold',
                     arrowprops=dict(arrowstyle='->', color=c_blue, lw=0.8))
    
    ax.set_xticks(x)
    ax.set_xticklabels(models)
    ax.set_ylabel('Degrees of Freedom for Signal')
    ax.set_title('A. Information content (DFS)', fontsize=10, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_ylim(0, 4.0)
    ax.grid(axis='y', alpha=0.3)
    
    # ═══ Panel B: Threshold (core result) ═══
    ax = axes[0, 1]
    multipliers = thresh['multipliers']
    improvements = thresh['improvements']
    sigma_vals = [improvements[str(m)]['mic_dD_sigma_permil'] for m in multipliers]
    improve_vals = [improvements[str(m)]['improvement_pct'] for m in multipliers]
    
    ax.plot(sigma_vals, improve_vals, 'o-', lw=2.5, color=c_blue, markersize=7, zorder=5)
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(10, color='gray', lw=1, ls='--', alpha=0.6, label='10% criterion')
    ax.axvspan(95, 140, alpha=0.12, color='red', label='Thanwerdas (2024)')
    ax.axvline(25, color=c_red, lw=2, ls=':', alpha=0.8)
    
    ax.fill_between(sigma_vals, improve_vals, 0,
                     where=[v >= 0 for v in improve_vals], alpha=0.08, color='green')
    ax.fill_between(sigma_vals, improve_vals, 0,
                     where=[v < 0 for v in improve_vals], alpha=0.08, color='red')
    
    ax.text(12, 38, 'δD helps', fontsize=8, color=c_green, ha='center', fontstyle='italic')
    ax.text(80, -80, 'δD hurts', fontsize=8, color=c_red, ha='center', fontstyle='italic')
    ax.annotate('Threshold\n~25‰', xy=(25, 0), xytext=(40, 25),
                fontsize=7.5, ha='center', color=c_red,
                arrowprops=dict(arrowstyle='->', color=c_red))
    
    ax.set_xlabel('σ(Mic δD) source signature (‰)', fontsize=9)
    ax.set_ylabel('FF constraint improvement (%)', fontsize=9)
    ax.set_title('B. The δD threshold', fontsize=10, fontweight='bold')
    ax.set_xlim(0, 140)
    ax.set_ylim(-180, 60)
    ax.legend(fontsize=7, loc='lower left')
    ax.grid(alpha=0.2)
    
    # ═══ Panel C: Thanwerdas diagnostic ═══
    ax = axes[1, 0]
    labels = ['δ¹³C only\n(reference)', 'Dual (our σ≈8‰)', 'Dual (Thanwerdas\nσ≈110‰)']
    ci_vals = [
        than['d13C_reference']['CI'],
        than['our_uncertainties']['CI'],
        than['thanwerdas_uncertainties']['CI'],
    ]
    colors = [c_orange, c_blue, c_red]
    bars = ax.bar(labels, ci_vals, color=colors, alpha=0.85, width=0.6)
    
    # Improvement annotations
    ax.annotate(f'−46%', xy=(1, ci_vals[1]), xytext=(1, ci_vals[1]+15),
                fontsize=9, ha='center', fontweight='bold', color=c_green)
    ax.annotate(f'+168%', xy=(2, ci_vals[2]), xytext=(2, ci_vals[2]+15),
                fontsize=9, ha='center', fontweight='bold', color=c_red)
    
    ax.set_ylabel('FF emissions 90% CI width (Tg/yr)', fontsize=9)
    ax.set_title('C. Thanwerdas replication', fontsize=10, fontweight='bold')
    ax.set_ylim(0, 310)
    ax.grid(axis='y', alpha=0.3)
    
    # ═══ Panel D: Sensitivity ═══
    ax = axes[1, 1]
    mults_plot = [1.0, 2.0, 3.0, 5.0, 8.0]
    sigma_plot = [8.25 * m for m in mults_plot]
    ref_ci = 101.3
    
    # KIE
    for label, ls, c in [('Saueressig', '-', c_blue), ('Cantrell', '--', c_purple), ('Sampled', ':', c_green)]:
        key = f'KIE_{label.lower()}'
        vals = [sens[key][str(m)] for m in mults_plot]
        improve = [(ref_ci - v) / ref_ci * 100 for v in vals]
        ax.plot(sigma_plot, improve, ls=ls, lw=2, color=c, marker='o', markersize=5, label=f'KIE: {label}')
    
    # Lifetime
    for label, ls, c in [('τ=9.0', '-', c_orange), ('τ varying', '--', c_red), ('τ=8.5', ':', '#636363')]:
        key_map = {'τ=9.0': 'tau_fixed_9.0', 'τ varying': 'tau_varying', 'τ=8.5': 'tau_fixed_8.5'}
        vals = [sens[key_map[label]][str(m)] for m in mults_plot]
        improve = [(ref_ci - v) / ref_ci * 100 for v in vals]
        ax.plot(sigma_plot, improve, ls=ls, lw=1.5, color=c, marker='s', markersize=4, 
                label=label, alpha=0.7)
    
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(10, color='gray', lw=1, ls='--', alpha=0.5)
    ax.axvline(25, color=c_red, lw=1.5, ls=':', alpha=0.7)
    
    ax.set_xlabel('σ(Mic δD) source signature (‰)', fontsize=9)
    ax.set_ylabel('FF constraint improvement (%)', fontsize=9)
    ax.set_title('D. Threshold robustness', fontsize=10, fontweight='bold')
    ax.legend(fontsize=6.5, loc='lower left', ncol=2)
    ax.set_xlim(0, 70)
    ax.set_ylim(-170, 60)
    ax.grid(alpha=0.2)
    
    ax.text(15, 35, 'Threshold stable\nacross all configs', fontsize=7.5,
            ha='center', fontstyle='italic', color='#636363',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.4))
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig_comprehensive_4panel.png', dpi=300, bbox_inches='tight')
    plt.savefig(FIG_DIR / 'fig_comprehensive_4panel.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_DIR}/fig_comprehensive_4panel.png")


if __name__ == "__main__":
    main()
