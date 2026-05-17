#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_kie_immunity.py — Main figure for KIE immunity paper
=========================================================
3-panel figure showing:
  A: Variance decomposition (where does FF uncertainty come from?)
  B: Total FF variance — δ¹³C-only vs dual isotopes
  C: Comparison with Basu 2022 (KIE spread)
"""

import sys
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS = Path(__file__).resolve().parent.parent / "results"
FIG_DIR = Path(__file__).resolve().parent
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    with open(RESULTS / "variance_decomposition.json") as f:
        var = json.load(f)
    with open(RESULTS / "basu_comparison.json") as f:
        basu = json.load(f)
    
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), dpi=150)
    
    c_kie = '#e34a33'
    c_sig = '#2171b5'
    c_tau = '#31a354'
    c_res = '#bdbdbd'
    c_dual = '#2171b5'
    c_d13c = '#fe9929'
    c_basu = '#636363'
    
    # ═══ Panel A: Variance decomposition ═══
    ax = axes[0]
    modes = ['δ¹³C only', 'δ¹³C + δD']
    keys = ['d13C_only', 'dual']
    kie_p = [var[k]['kie_pct'] for k in keys]
    sig_p = [var[k]['sigs_pct'] for k in keys]
    tau_p = [var[k]['tau_pct'] for k in keys]
    res_p = [var[k]['residual_pct'] for k in keys]
    
    x = np.arange(len(modes))
    w = 0.55
    p1 = ax.bar(x, kie_p, w, label='KIE', color=c_kie)
    p2 = ax.bar(x, sig_p, w, bottom=kie_p, label='Source sigs', color=c_sig)
    p3 = ax.bar(x, tau_p, w, bottom=[a+b for a,b in zip(kie_p,sig_p)], 
                 label='Lifetime', color=c_tau)
    p4 = ax.bar(x, res_p, w, bottom=[a+b+c for a,b,c in zip(kie_p,sig_p,tau_p)], 
                 label='Other', color=c_res)
    
    # Annotate KIE %
    for i, k in enumerate(kie_p):
        if k > 0.5:
            ax.text(i, k/2, f'KIE\n{k:.0f}%', ha='center', va='center',
                    fontsize=9, fontweight='bold', color='white')
    
    ax.set_xticks(x)
    ax.set_xticklabels(modes)
    ax.set_ylabel('Variance contribution (%)', fontsize=10)
    ax.set_title('A. FF variance decomposition', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.set_ylim(0, 105)
    ax.grid(axis='y', alpha=0.3)
    
    # ═══ Panel B: Total variance ═══
    ax = axes[1]
    var_vals = [var['d13C_only']['var_total'], var['dual']['var_total']]
    std_vals = [np.sqrt(v) for v in var_vals]
    bars = ax.bar(modes, std_vals, color=[c_d13c, c_dual], alpha=0.85, width=0.55)
    
    for i, (s, v) in enumerate(zip(std_vals, var_vals)):
        ax.text(i, s + 1, f'σ = {s:.1f}\n(var = {v:.0f})',
                ha='center', fontsize=9)
    
    reduction = (1 - var_vals[1] / var_vals[0]) * 100
    ax.annotate(f'  −{reduction:.0f}% variance  ',
                xy=(0.5, np.mean(std_vals)),
                xytext=(0.5, max(std_vals)*0.85),
                fontsize=10, ha='center', color=c_dual, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', 
                          edgecolor=c_dual, lw=1))
    
    ax.set_ylabel('FF uncertainty σ (Tg/yr)', fontsize=10)
    ax.set_title('B. Total FF uncertainty', fontsize=11, fontweight='bold')
    ax.set_ylim(0, max(std_vals)*1.25)
    ax.grid(axis='y', alpha=0.3)
    
    # ═══ Panel C: Basu 2022 comparison ═══
    ax = axes[2]
    
    labels = ['Basu 2022\n(3D, δ¹³C-only)', 'Our 2-box\n(δ¹³C-only)', 'Our 2-box\n(δ¹³C + δD)']
    spreads = [basu['basu_spread'],
               basu['our_d13C_only']['kie_spread'],
               basu['our_dual']['kie_spread']]
    colors = [c_basu, c_d13c, c_dual]
    
    bars = ax.bar(labels, spreads, color=colors, alpha=0.85, width=0.55)
    for i, s in enumerate(spreads):
        ax.text(i, s + 0.3, f'{s:.1f}', ha='center', fontsize=9, fontweight='bold')
    
    ax.set_ylabel('KIE-driven FF spread (Tg/yr)', fontsize=10)
    ax.set_title('C. KIE spread: 2-box vs 3D inversion', fontsize=11, fontweight='bold')
    ax.set_ylim(0, max(spreads)*1.2)
    ax.grid(axis='y', alpha=0.3)
    
    # Annotate insight
    ax.text(0.5, max(spreads)*0.6,
            'Constrained 2-box (bounded LS\n+ BB-fixed) already eliminates\nKIE ambiguity — even without δD',
            transform=ax.transData, fontsize=7.5, ha='center', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7))
    
    plt.suptitle('KIE immunity: dual isotopes + bounded LS collapse the Saueressig–Cantrell controversy',
                  fontsize=11, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    fig.savefig(FIG_DIR / "fig_kie_immunity.png", dpi=300, bbox_inches='tight')
    fig.savefig(FIG_DIR / "fig_kie_immunity.pdf", bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_DIR}/fig_kie_immunity.png")


if __name__ == "__main__":
    main()
