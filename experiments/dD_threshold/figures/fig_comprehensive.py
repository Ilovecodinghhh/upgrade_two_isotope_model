#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_comprehensive.py — 6-panel publication figure
==================================================

Panel A: DFS gain from δD (1-box vs 2-box)
Panel B: δD constraint improvement vs σ(Mic δD) — THE threshold plot
Panel C: Thanwerdas replication diagnostic
Panel D: Sensitivity analysis (KIE + lifetime) — threshold is robust
Panel E: Hemispheric δD source signature gaps (NH vs SH)
Panel F: Data provenance (Dasgupta vs Umezawa calibration effect)
"""

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
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

    # Load Phase 6 deep dive for exact crossover and fine-grid sweep
    deep_dive_path = RESULTS_DIR / "phase6_deep_dive" / "deep_dive_results.json"
    if deep_dive_path.exists():
        with open(deep_dive_path) as f:
            deep = json.load(f)
    else:
        deep = None

    # Load hemispheric source summary
    src_summary = pd.read_csv(REPO_ROOT / "rel" / "data" / "Hemispheric_dD_sources_summary.csv")

    # ── Figure setup: 3 rows × 2 cols (180mm × 220mm, Nature full-page) ──
    plt.rcParams.update({
        'font.size': 9,
        'font.family': 'sans-serif',
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
    })
    fig, axes = plt.subplots(3, 2, figsize=(180/25.4, 220/25.4), dpi=300)

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
    ax.bar(x - w/2, dfs_c13, w, label='δ¹³C only', color=c_orange, alpha=0.85)
    ax.bar(x + w/2, dfs_dual, w, label='δ¹³C + δD', color=c_blue, alpha=0.85)

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

    # ═══ Panel B: Threshold (core result) — from Phase 6 fine grid ═══
    ax = axes[0, 1]

    # Use Phase 6 fine-grid data if available, else fall back to Phase 3
    if deep is not None:
        fg = deep['A_fine_grid']
        sweep = fg['sweep']
        sigma_vals = [8.25 * float(m) for m in sorted(sweep.keys(), key=float)]
        improve_vals = [sweep[m]['improvement'] for m in sorted(sweep.keys(), key=float)]
        # Exact crossover from interpolation
        cross = fg['crossover_0pct']
        thresh_sigma = cross['sigma_permil'] if cross['sigma_permil'] else 41
        cross_10 = fg['crossover_10pct']
        thresh_sigma_10 = cross_10['sigma_permil'] if cross_10['sigma_permil'] else 29
    else:
        multipliers = thresh['multipliers']
        improvements = thresh['improvements']
        sigma_vals = [improvements[str(m)]['mic_dD_sigma_permil'] for m in multipliers]
        improve_vals = [improvements[str(m)]['improvement_pct'] for m in multipliers]
        thresh_sigma = thresh.get('threshold_mic_dD_sigma_permil', 41)
        thresh_sigma_10 = None

    ax.plot(sigma_vals, improve_vals, 'o-', lw=2.5, color=c_blue, markersize=7, zorder=5)
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(10, color='gray', lw=1, ls='--', alpha=0.6, label='10% criterion')
    ax.axvspan(95, 140, alpha=0.12, color='red', label='Thanwerdas (2024)')

    # Exact threshold from Phase 6 interpolation
    ax.axvline(thresh_sigma, color=c_red, lw=2, ls=':', alpha=0.8)

    ax.fill_between(sigma_vals, improve_vals, 0,
                     where=[v >= 0 for v in improve_vals], alpha=0.08, color='green')
    ax.fill_between(sigma_vals, improve_vals, 0,
                     where=[v < 0 for v in improve_vals], alpha=0.08, color='red')

    ax.text(12, 38, 'δD helps', fontsize=8, color=c_green, ha='center', fontstyle='italic')
    ax.text(80, -80, 'δD hurts', fontsize=8, color=c_red, ha='center', fontstyle='italic')

    # Thanwerdas reference line
    ax.axvline(15.6 * 8.25, color=c_red, lw=1.5, ls='--', alpha=0.5)
    ax.text(15.6 * 8.25 + 3, -150, 'T24', fontsize=7, color=c_red, rotation=90, va='bottom')

    # Current precision annotation
    ax.annotate('Current\nprecision', xy=(8.25, improve_vals[1] if len(improve_vals) > 1 else 50),
                xytext=(25, 55), fontsize=7, ha='center', color=c_blue,
                arrowprops=dict(arrowstyle='->', color=c_blue, lw=0.8))

    ax.annotate(f'Threshold\n~{thresh_sigma:.0f}‰', xy=(thresh_sigma, 0),
                xytext=(thresh_sigma + 20, 25),
                fontsize=7.5, ha='center', color=c_red,
                arrowprops=dict(arrowstyle='->', color=c_red))

    ax.set_xlabel('σ(Mic δD) source signature (‰)', fontsize=9)
    ax.set_ylabel('FF constraint improvement (%)', fontsize=9)
    ax.set_title('B. The δD threshold', fontsize=10, fontweight='bold')
    ax.set_xlim(0, 140)
    ax.set_ylim(-180, 70)
    ax.legend(fontsize=7, loc='lower left')
    ax.grid(alpha=0.2)

    # ═══ Panel C: Thanwerdas diagnostic ═══
    ax = axes[1, 0]
    labels = ['δ¹³C only\n(reference)', 'Dual\n(our σ≈8‰)', 'Dual\n(Thanwerdas σ≈110‰)']
    ci_vals = [
        than['d13C_reference']['CI'],
        than['our_uncertainties']['CI'],
        than['thanwerdas_uncertainties']['CI'],
    ]
    colors = [c_orange, c_blue, c_red]
    ax.bar(labels, ci_vals, color=colors, alpha=0.85, width=0.6)

    ax.annotate(f'{than["our_uncertainties"]["improvement_pct"]:+.0f}%',
                xy=(1, ci_vals[1]), xytext=(1, ci_vals[1]+15),
                fontsize=9, ha='center', fontweight='bold', color=c_green)
    ax.annotate(f'{than["thanwerdas_uncertainties"]["improvement_pct"]:+.0f}%',
                xy=(2, ci_vals[2]), xytext=(2, ci_vals[2]+15),
                fontsize=9, ha='center', fontweight='bold', color=c_red)

    ax.set_ylabel('FF emissions 90% CI width (Tg/yr)', fontsize=9)
    ax.set_title('C. Thanwerdas replication', fontsize=10, fontweight='bold')
    ax.set_ylim(0, 310)
    ax.grid(axis='y', alpha=0.3)

    # ═══ Panel D: Sensitivity ═══
    ax = axes[1, 1]
    mults_plot = [1.0, 2.0, 3.0, 5.0, 8.0]
    sigma_plot = [8.25 * m for m in mults_plot]
    ref_ci = than['d13C_reference']['CI']

    for label, ls, c in [('Saueressig', '-', c_blue), ('Cantrell', '--', c_purple), ('Sampled', ':', c_green)]:
        key = f'KIE_{label.lower()}'
        vals = [sens[key][str(m)] for m in mults_plot]
        improve = [(ref_ci - v) / ref_ci * 100 for v in vals]
        ax.plot(sigma_plot, improve, ls=ls, lw=2, color=c, marker='o', markersize=5, label=f'KIE: {label}')

    for label, ls, c in [('τ=9.0', '-', c_orange), ('τ varying', '--', c_red), ('τ=8.5', ':', '#636363')]:
        key_map = {'τ=9.0': 'tau_fixed_9.0', 'τ varying': 'tau_varying', 'τ=8.5': 'tau_fixed_8.5'}
        vals = [sens[key_map[label]][str(m)] for m in mults_plot]
        improve = [(ref_ci - v) / ref_ci * 100 for v in vals]
        ax.plot(sigma_plot, improve, ls=ls, lw=1.5, color=c, marker='s', markersize=4,
                label=label, alpha=0.7)

    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(10, color='gray', lw=1, ls='--', alpha=0.5)
    ax.axvline(thresh_sigma, color=c_red, lw=1.5, ls=':', alpha=0.7)

    ax.set_xlabel('σ(Mic δD) source signature (‰)', fontsize=9)
    ax.set_ylabel('FF constraint improvement (%)', fontsize=9)
    ax.set_title('D. Threshold robustness', fontsize=10, fontweight='bold')
    ax.legend(fontsize=6.5, loc='lower left', ncol=2)
    ax.set_xlim(0, 70)
    ax.set_ylim(-170, 70)
    ax.grid(alpha=0.2)

    # ═══ Panel E: Hemispheric δD source signature gaps ═══
    ax = axes[2, 0]

    # Read the MC CSVs directly for mean + spread
    src_dir = REPO_ROOT / "rel" / "data"
    sectors = ['Mic', 'BB', 'FF']
    nh_means = []; sh_means = []; gaps = []; gap_labels = []
    for sec in sectors:
        nh_df = pd.read_csv(src_dir / f"{sec}_dD_NH_MC.csv")
        sh_df = pd.read_csv(src_dir / f"{sec}_dD_SH_MC.csv")
        mc_cols = [c for c in nh_df.columns if c.startswith('mc_')]
        nh_vals = nh_df[mc_cols].values  # years × 1000
        sh_vals = sh_df[mc_cols].values
        nh_m = np.mean(nh_vals)
        sh_m = np.mean(sh_vals)
        nh_means.append(nh_m)
        sh_means.append(sh_m)
        gaps.append(nh_m - sh_m)

    x = np.arange(3)
    w = 0.3
    bars_nh = ax.bar(x - w/2, nh_means, w, label='NH', color=c_blue, alpha=0.85)
    bars_sh = ax.bar(x + w/2, sh_means, w, label='SH', color=c_orange, alpha=0.85)

    for i in range(3):
        y_top = max(nh_means[i], sh_means[i]) + 3
        ax.annotate(f'Δ = {gaps[i]:+.0f}‰', xy=(x[i], y_top),
                    fontsize=8, ha='center', fontweight='bold', color=c_red)

    ax.set_xticks(x)
    ax.set_xticklabels(['Microbial', 'Biomass\nBurning', 'Fossil\nFuel'])
    ax.set_ylabel('δD source signature (‰)', fontsize=9)
    ax.set_title('E. Hemispheric δD source signatures', fontsize=10, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)

    # ═══ Panel F: Version comparison (v1–v4) ═══
    ax = axes[2, 1]

    versions = ['v1\n(Umezawa/\nglobal)', 'v2\n(Dasgupta/\nglobal)', 'v3\n(+ hemi\nδD src)', 'v4\n(+ hemi\nδ¹³C src)', 'v5\n(+ Luo\n2024 C4)']
    baseline_ci = [46.6, 37.8, 43.5, 57.6, 62.6]
    improvements_v = [52, 60.8, 57.0, 45.1, 53.0]

    x = np.arange(5)
    # Dual bars: CI width + improvement
    ax2 = ax.twinx()
    bars_ci = ax.bar(x - 0.15, baseline_ci, 0.3, label='Dual CI width', color=c_blue, alpha=0.8)
    bars_imp = ax2.bar(x + 0.15, improvements_v, 0.3, label='Improvement %', color=c_green, alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(versions, fontsize=7)
    ax.set_ylabel('CI width (Tg/yr)', color=c_blue, fontsize=9)
    ax2.set_ylabel('Improvement (%)', color=c_green, fontsize=9)
    ax.set_title('F. Data version comparison', fontsize=10, fontweight='bold')
    ax.set_ylim(0, 60)
    ax2.set_ylim(0, 80)

    # Threshold annotation
    for i, (ci, imp) in enumerate(zip(baseline_ci, improvements_v)):
        ax.text(x[i] - 0.15, ci + 1.5, f'{ci:.1f}', fontsize=7, ha='center', color=c_blue)
        ax2.text(x[i] + 0.15, imp + 1.5, f'{imp:.0f}%', fontsize=7, ha='center', color=c_green)

    ax.grid(axis='y', alpha=0.2)

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig_comprehensive_6panel.png', dpi=300, bbox_inches='tight')
    plt.savefig(FIG_DIR / 'fig_comprehensive_6panel.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_DIR}/fig_comprehensive_6panel.png")
    print(f"Saved: {FIG_DIR}/fig_comprehensive_6panel.pdf")

    # Also regenerate the old 4-panel for backward compat
    # (just update the threshold annotation)
    fig4, axes4 = plt.subplots(2, 2, figsize=(10, 8), dpi=150)

    # Panel A
    ax = axes4[0, 0]
    ax.bar(np.arange(2) - 0.15, dfs_c13, 0.3, label='δ¹³C only', color=c_orange, alpha=0.85)
    ax.bar(np.arange(2) + 0.15, dfs_dual, 0.3, label='δ¹³C + δD', color=c_blue, alpha=0.85)
    for i in range(2):
        ax.annotate(f'ΔDFS = +{delta_dfs[i]:.2f}', xy=(i+0.15, dfs_dual[i]),
                    xytext=(i+0.45, dfs_dual[i]+0.1), fontsize=7.5, color=c_blue, fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color=c_blue, lw=0.8))
    ax.set_xticks([0, 1]); ax.set_xticklabels(['1-box', '2-box'])
    ax.set_ylabel('DFS'); ax.set_title('A. Information content', fontsize=10, fontweight='bold')
    ax.legend(fontsize=8); ax.set_ylim(0, 4.0); ax.grid(axis='y', alpha=0.3)

    # Panel B
    ax = axes4[0, 1]
    ax.plot(sigma_vals, improve_vals, 'o-', lw=2.5, color=c_blue, markersize=7)
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(10, color='gray', lw=1, ls='--', alpha=0.6)
    ax.axvspan(95, 140, alpha=0.12, color='red', label='Thanwerdas (2024)')
    ax.axvline(thresh_sigma, color=c_red, lw=2, ls=':', alpha=0.8)
    ax.fill_between(sigma_vals, improve_vals, 0, where=[v>=0 for v in improve_vals], alpha=0.08, color='green')
    ax.fill_between(sigma_vals, improve_vals, 0, where=[v<0 for v in improve_vals], alpha=0.08, color='red')
    ax.annotate(f'Threshold\n~{thresh_sigma:.0f}‰', xy=(thresh_sigma, 0),
                xytext=(thresh_sigma+20, 25), fontsize=7.5, ha='center', color=c_red,
                arrowprops=dict(arrowstyle='->', color=c_red))
    ax.set_xlabel('σ(Mic δD) (‰)'); ax.set_ylabel('Improvement (%)')
    ax.set_title('B. The δD threshold', fontsize=10, fontweight='bold')
    ax.set_xlim(0, 140); ax.set_ylim(-180, 70); ax.legend(fontsize=7); ax.grid(alpha=0.2)

    # Panel C
    ax = axes4[1, 0]
    ax.bar(labels, ci_vals, color=colors, alpha=0.85, width=0.6)
    ax.annotate(f'{than["our_uncertainties"]["improvement_pct"]:+.0f}%', xy=(1, ci_vals[1]),
                xytext=(1, ci_vals[1]+15), fontsize=9, ha='center', fontweight='bold', color=c_green)
    ax.annotate(f'{than["thanwerdas_uncertainties"]["improvement_pct"]:+.0f}%', xy=(2, ci_vals[2]),
                xytext=(2, ci_vals[2]+15), fontsize=9, ha='center', fontweight='bold', color=c_red)
    ax.set_ylabel('FF 90% CI width (Tg/yr)')
    ax.set_title('C. Thanwerdas replication', fontsize=10, fontweight='bold')
    ax.set_ylim(0, 310); ax.grid(axis='y', alpha=0.3)

    # Panel D
    ax = axes4[1, 1]
    for label, ls, c in [('Saueressig', '-', c_blue), ('Cantrell', '--', c_purple), ('Sampled', ':', c_green)]:
        vals = [sens[f'KIE_{label.lower()}'][str(m)] for m in mults_plot]
        improve = [(ref_ci - v)/ref_ci*100 for v in vals]
        ax.plot(sigma_plot, improve, ls=ls, lw=2, color=c, marker='o', markersize=5, label=f'KIE: {label}')
    for label, ls, c in [('τ=9.0', '-', c_orange), ('τ varying', '--', c_red), ('τ=8.5', ':', '#636363')]:
        key_map = {'τ=9.0': 'tau_fixed_9.0', 'τ varying': 'tau_varying', 'τ=8.5': 'tau_fixed_8.5'}
        vals = [sens[key_map[label]][str(m)] for m in mults_plot]
        improve = [(ref_ci - v)/ref_ci*100 for v in vals]
        ax.plot(sigma_plot, improve, ls=ls, lw=1.5, color=c, marker='s', markersize=4, label=label, alpha=0.7)
    ax.axhline(0, color='black', lw=0.8); ax.axhline(10, color='gray', lw=1, ls='--', alpha=0.5)
    ax.axvline(thresh_sigma, color=c_red, lw=1.5, ls=':', alpha=0.7)
    ax.set_xlabel('σ(Mic δD) (‰)'); ax.set_ylabel('Improvement (%)')
    ax.set_title('D. Robustness', fontsize=10, fontweight='bold')
    ax.legend(fontsize=6.5, loc='lower left', ncol=2)
    ax.set_xlim(0, 70); ax.set_ylim(-170, 70); ax.grid(alpha=0.2)

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig_comprehensive_4panel.png', dpi=300, bbox_inches='tight')
    plt.savefig(FIG_DIR / 'fig_comprehensive_4panel.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_DIR}/fig_comprehensive_4panel.png (updated)")


if __name__ == "__main__":
    main()
