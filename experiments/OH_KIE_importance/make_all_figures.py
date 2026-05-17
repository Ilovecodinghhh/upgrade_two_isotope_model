#!/usr/bin/env python3
"""
OH KIE Importance — Complete Figure Suite
==========================================
Generates all 10 figures from 3×3 and 2×2 experiment results.
Replaces make_figures.py and make_2x2_figures.py.
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

RESULTS = Path(__file__).resolve().parent / "results"
FIGS = Path(__file__).resolve().parent / "figures"
FIGS.mkdir(exist_ok=True)

with open(RESULTS / "summary.json") as f:
    s3 = json.load(f)
with open(RESULTS / "2x2_summary.json") as f:
    s2 = json.load(f)

# ─── Colorblind-safe palette ───
PAL = {
    'base':  '#4575b4',   # steel blue
    'oh13c': '#d73027',   # red
    'ohd':   '#fc8d59',   # orange
    'both':  '#91bfdb',   # light blue
    'allk':  '#1a9850',   # green
    'sau':   '#4575b4',   # blue (Saueressig)
    'can':   '#d73027',   # red  (Cantrell)
    'resid': '#e0e0e0',   # light gray
    'other': '#8073ac',   # purple
    'c13':   '#d73027',   # red for 13C
    'dD':    '#4575b4',   # blue for D
}

DPI = 250

def load_ts(prefix, config):
    return pd.read_csv(RESULTS / f"{prefix}{config}_timeseries.csv")


# ====================================================================
# FIGURE 1 — 3×3 Variance Attribution Bar Chart
# ====================================================================
def fig1():
    fig, ax = plt.subplots(figsize=(10, 6), dpi=DPI)
    s2_base = s3['ALL_SAMPLED']['sigma_ff']**2

    configs = ['ALL_SAMPLED', 'FIX_OH13C', 'FIX_OHD', 'FIX_BOTH_OH', 'ALL_KIE_FIXED']
    labels = ['All sampled\n(baseline)', 'Fix OH-¹³C', 'Fix OH-D', 'Fix both\nOH KIEs', 'Fix all\n8 KIE']
    colors = [PAL['base'], PAL['oh13c'], PAL['ohd'], PAL['both'], PAL['allk']]
    sigmas = [s3[c]['sigma_ff'] for c in configs]
    reds = [0] + [100*(s2_base - s3[c]['sigma_ff']**2)/s2_base for c in configs[1:]]

    bars = ax.bar(range(5), sigmas, color=colors, edgecolor='black', lw=0.8, width=0.65)
    for i, (b, sig, red) in enumerate(zip(bars, sigmas, reds)):
        ax.text(b.get_x()+b.get_width()/2, sig+1.2,
                f'σ = {sig:.1f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
        if i > 0:
            col_txt = 'white' if i != 3 else 'black'
            ax.text(b.get_x()+b.get_width()/2, sig/2,
                    f'−{red:.1f}%' if red > 0 else f'{red:+.1f}%',
                    ha='center', va='center', fontsize=10, color=col_txt, fontweight='bold')

    ax.axhline(s3['ALL_SAMPLED']['sigma_ff'], color='gray', ls='--', alpha=0.4, lw=1)
    ax.set_xticks(range(5)); ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('σ(FF) [Tg yr⁻¹]', fontsize=12)
    ax.set_title('3×3 Dual-Isotope One-Box: FF Uncertainty by KIE Configuration', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 78)
    ax.grid(axis='y', alpha=0.25); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIGS / "fig1_variance_attribution.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 1 ✓  3×3 variance attribution bar chart")


# ====================================================================
# FIGURE 2 — OH-¹³C Level Shift (Saueressig vs Cantrell)
# ====================================================================
def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=DPI, gridspec_kw={'width_ratios': [1.6, 1]})

    # (a) Time series
    ax = axes[0]
    for cfg, lbl, col, ls in [
        ('OH13C_SAUERESSIG', 'Saueressig (α = 1.0039)', PAL['sau'], '-'),
        ('OH13C_CANTRELL',   'Cantrell (α = 1.0054)',    PAL['can'], '-'),
        ('ALL_SAMPLED',      'All sampled (baseline)',    'gray',     '--'),
    ]:
        ts = load_ts('', cfg)
        ax.plot(ts['year'], ts['FF_mean'], color=col, ls=ls, lw=2.2, label=lbl)
        ax.fill_between(ts['year'], ts['FF_mean']-ts['FF_std'],
                        ts['FF_mean']+ts['FF_std'], alpha=0.12, color=col)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xlabel('Year', fontsize=11); ax.set_ylabel('FF Emissions [Tg yr⁻¹]', fontsize=11)
    ax.set_title('(a) OH-¹³C shifts FF level — identical spread', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left'); ax.grid(alpha=0.25)
    ax.set_ylim(-120, 200)

    # (b) Bar chart with shift annotation
    ax = axes[1]
    cfgs = ['OH13C_SAUERESSIG', 'FIX_OH13C', 'OH13C_CANTRELL']
    lbs = ['Saueressig\n(1.0039)', 'Midpoint\n(1.00465)', 'Cantrell\n(1.0054)']
    cols = [PAL['sau'], 'gray', PAL['can']]
    means = [s3[c]['mean_ff'] for c in cfgs]
    sigs = [s3[c]['sigma_ff'] for c in cfgs]
    bars = ax.bar(range(3), means, yerr=sigs, color=cols, edgecolor='black', lw=0.8,
                  width=0.55, capsize=5, error_kw={'lw': 1.5})
    for i, (b, m, s) in enumerate(zip(bars, means, sigs)):
        ax.text(b.get_x()+b.get_width()/2, m+s+2, f'{m:.0f}±{s:.0f}',
                ha='center', fontsize=9.5, fontweight='bold')
    delta = means[2] - means[0]
    mid_y = (means[0] + means[2]) / 2
    ax.annotate('', xy=(2, means[2]+2), xytext=(0, means[0]+2),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text(1, mid_y+8, f'Δ = {delta:+.0f} Tg/yr', ha='center', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.85))
    ax.set_xticks(range(3)); ax.set_xticklabels(lbs, fontsize=9.5)
    ax.set_ylabel('Mean FF [Tg yr⁻¹]', fontsize=11)
    ax.set_title(f'(b) Level shift = {abs(delta):.0f} Tg/yr', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.25); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIGS / "fig2_oh13c_level_shift.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 2 ✓  OH-¹³C level shift (Saueressig vs Cantrell)")


# ====================================================================
# FIGURE 3 — Variance Decomposition Pie + Cumulative Reduction
# ====================================================================
def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=DPI)
    s2v = {c: s3[c]['sigma_ff']**2 for c in s3}
    tv = s2v['ALL_SAMPLED']
    ohd_v = tv - s2v['FIX_OHD']
    oh13c_v = max(tv - s2v['FIX_OH13C'], 0)
    both_v = tv - s2v['FIX_BOTH_OH']
    all_v = tv - s2v['ALL_KIE_FIXED']
    other_v = max(all_v - both_v, 0)
    resid_v = s2v['ALL_KIE_FIXED']

    # (a) Pie
    ax = axes[0]
    sizes = [ohd_v, other_v, oh13c_v if oh13c_v > 5 else 0, resid_v]
    lbls = [f'OH-D\n({100*ohd_v/tv:.1f}%)',
            f'Other KIE\n({100*other_v/tv:.1f}%)',
            f'OH-¹³C\n({100*oh13c_v/tv:.1f}%)' if oh13c_v > 5 else f'OH-¹³C\n(≈0%)',
            f'Source sigs +\ndata noise\n({100*resid_v/tv:.1f}%)']
    cols = [PAL['ohd'], PAL['other'], PAL['oh13c'], PAL['resid']]
    wedges, texts, autotexts = ax.pie(sizes, labels=lbls, colors=cols, explode=(0.05,0,0,0),
                                       startangle=90, autopct=lambda p: f'{p:.0f}%' if p>2 else '',
                                       pctdistance=0.7, textprops={'fontsize': 9.5})
    ax.set_title('(a) σ²(FF) decomposition — 3×3', fontsize=12, fontweight='bold')

    # (b) Cumulative
    ax = axes[1]
    steps = ['Baseline', '+Fix OH-D', '+Fix OH-¹³C', '+Fix all KIE']
    vals = [tv, s2v['FIX_OHD'], s2v['FIX_BOTH_OH'], s2v['ALL_KIE_FIXED']]
    svals = [np.sqrt(v) for v in vals]
    rcum = [0] + [100*(tv-v)/tv for v in vals[1:]]
    bcols = [PAL['base'], PAL['ohd'], PAL['both'], PAL['allk']]
    bars = ax.bar(range(4), svals, color=bcols, edgecolor='black', lw=0.8, width=0.6)
    for i, (b, sv, rc) in enumerate(zip(bars, svals, rcum)):
        ax.text(b.get_x()+b.get_width()/2, sv+0.8, f'σ={sv:.1f}',
                ha='center', fontsize=10, fontweight='bold')
        if i > 0:
            ax.text(b.get_x()+b.get_width()/2, sv/2, f'−{rc:.0f}%',
                    ha='center', fontsize=10, color='white', fontweight='bold')
    ax.set_xticks(range(4)); ax.set_xticklabels(steps, fontsize=9.5)
    ax.set_ylabel('σ(FF) [Tg yr⁻¹]', fontsize=11)
    ax.set_title('(b) Cumulative σ reduction', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 78); ax.grid(axis='y', alpha=0.25)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIGS / "fig3_variance_decomposition.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 3 ✓  Variance decomposition (pie + cumulative)")


# ====================================================================
# FIGURE 4 — 3×3 FF Time Series Comparison (4 configs)
# ====================================================================
def fig4():
    fig, ax = plt.subplots(figsize=(12, 6), dpi=DPI)
    for cfg, lbl, col, ls, lw in [
        ('ALL_SAMPLED',   'All KIE sampled (baseline)', PAL['base'], '-',  2.5),
        ('FIX_OHD',       'Fix OH-D only',              PAL['ohd'],  '-',  2.0),
        ('FIX_BOTH_OH',   'Fix both OH KIEs',           PAL['both'], '--', 1.8),
        ('ALL_KIE_FIXED', 'Fix all 8 KIEs',             PAL['allk'], ':',  2.0),
    ]:
        ts = load_ts('', cfg)
        ax.plot(ts['year'], ts['FF_mean'], color=col, ls=ls, lw=lw, label=lbl)
        ax.fill_between(ts['year'], ts['FF_mean']-ts['FF_std'],
                        ts['FF_mean']+ts['FF_std'], alpha=0.10, color=col)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xlabel('Year', fontsize=12); ax.set_ylabel('FF Emissions [Tg yr⁻¹]', fontsize=12)
    ax.set_title('3×3 One-Box: FF Emissions Under Progressive KIE Freezing', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='upper left'); ax.grid(alpha=0.25); ax.set_ylim(-120, 200)
    plt.tight_layout()
    plt.savefig(FIGS / "fig4_ff_timeseries_comparison.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 4 ✓  3×3 FF time series (4 configs)")


# ====================================================================
# FIGURE 5 — Dual-Role Summary Diagram
# ====================================================================
def fig5():
    fig, ax = plt.subplots(figsize=(12, 5.5), dpi=DPI)
    ax.set_xlim(-0.3, 3.3); ax.set_ylim(-0.6, 2.6); ax.axis('off')

    ax.text(1.5, 2.4, 'OH KIE: Two Distinct Roles in the 3×3 Inversion',
            fontsize=15, fontweight='bold', ha='center')

    # OH-D box
    r1 = FancyBboxPatch((0, 0.35), 1.25, 1.5, boxstyle="round,pad=0.12",
                         facecolor=PAL['ohd'], alpha=0.25, edgecolor='black', lw=2)
    ax.add_patch(r1)
    ax.text(0.625, 1.65, 'OH-D KIE', fontsize=14, fontweight='bold', ha='center')
    ax.text(0.625, 1.30, 'U(1.294, 1.327)', fontsize=9, ha='center', family='monospace')
    ax.text(0.625, 0.95, 'Controls σ(FF)', fontsize=12, ha='center')
    ax.text(0.625, 0.65, '29% of variance', fontsize=13, ha='center',
            fontweight='bold', color='#b35806')

    # OH-13C box
    r2 = FancyBboxPatch((1.75, 0.35), 1.25, 1.5, boxstyle="round,pad=0.12",
                         facecolor=PAL['oh13c'], alpha=0.15, edgecolor='black', lw=2)
    ax.add_patch(r2)
    ax.text(2.375, 1.65, 'OH-¹³C KIE', fontsize=14, fontweight='bold', ha='center')
    ax.text(2.375, 1.30, 'U(1.0039, 1.0054)', fontsize=9, ha='center', family='monospace')
    ax.text(2.375, 0.95, 'Controls FF level', fontsize=12, ha='center')
    ax.text(2.375, 0.65, 'Δ = 17 Tg/yr shift', fontsize=13, ha='center',
            fontweight='bold', color='#b2182b')

    ax.annotate('SPREAD\n(random uncertainty)', xy=(0.625, 0.35), xytext=(0.625, -0.25),
                fontsize=10.5, ha='center', va='top', color='#b35806', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#b35806', lw=2))
    ax.annotate('BIAS\n(systematic shift)', xy=(2.375, 0.35), xytext=(2.375, -0.25),
                fontsize=10.5, ha='center', va='top', color='#b2182b', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#b2182b', lw=2))

    plt.tight_layout()
    plt.savefig(FIGS / "fig5_dual_role_summary.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 5 ✓  Dual-role summary diagram")


# ====================================================================
# FIGURE 6 — 2×2 Variance Attribution (δ¹³C and δD side-by-side)
# ====================================================================
def fig6():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=DPI, sharey=True)
    configs = ['ALL_SAMPLED', 'FIX_OH13C', 'FIX_OHD', 'FIX_BOTH_OH', 'ALL_KIE_FIXED']
    labels = ['All\nsampled', 'Fix\nOH-¹³C', 'Fix\nOH-D', 'Fix both\nOH', 'Fix all\nKIE']

    for panel, key, title_iso, ax in [
        (0, 'sigma_ff_c13', 'δ¹³C', axes[0]),
        (1, 'sigma_ff_dD',  'δD',   axes[1]),
    ]:
        vals = [s2[c][key] for c in configs]
        s2_base = s2['ALL_SAMPLED'][key]**2
        reds = [0] + [100*(s2_base - s2[c][key]**2)/s2_base for c in configs[1:]]

        # Color the active KIE bar differently
        if key == 'sigma_ff_c13':
            colors = [PAL['base'], PAL['oh13c'], '#cccccc', PAL['both'], PAL['allk']]
        else:
            colors = [PAL['base'], '#cccccc', PAL['dD'], PAL['both'], PAL['allk']]

        bars = ax.bar(range(5), vals, color=colors, edgecolor='black', lw=0.8, width=0.6)
        for i, (b, v, r) in enumerate(zip(bars, vals, reds)):
            ax.text(b.get_x()+b.get_width()/2, v+0.5, f'{v:.1f}',
                    ha='center', fontsize=10, fontweight='bold')
            if abs(r) > 0.5:
                tcol = 'white' if colors[i] not in ['#cccccc', PAL['both']] else 'black'
                ax.text(b.get_x()+b.get_width()/2, v/2, f'−{r:.0f}%' if r>0 else '0%',
                        ha='center', fontsize=9, color=tcol, fontweight='bold')

        ax.set_xticks(range(5)); ax.set_xticklabels(labels, fontsize=9)
        side = '(a)' if panel == 0 else '(b)'
        ax.set_title(f'{side} {title_iso}-only inversion', fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.25); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    axes[0].set_ylabel('σ(FF) [Tg yr⁻¹]', fontsize=11)
    axes[0].set_ylim(0, 48)
    fig.suptitle('2×2 One-Box: FF Uncertainty by KIE Configuration (separate isotope inversions)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGS / "fig6_2x2_variance_attribution.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 6 ✓  2×2 variance attribution (δ¹³C + δD)")


# ====================================================================
# FIGURE 7 — 2×2 Saueressig vs Cantrell (δ¹³C and δD panels)
# ====================================================================
def fig7():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=DPI)

    # (a) δ¹³C
    ax = axes[0]
    for cfg, lbl, col, ls in [
        ('OH13C_SAUERESSIG', 'Saueressig (1.0039)', PAL['sau'], '-'),
        ('OH13C_CANTRELL',   'Cantrell (1.0054)',    PAL['can'], '-'),
        ('ALL_SAMPLED',      'All sampled',          'gray',     '--'),
    ]:
        ts = load_ts('2x2_', cfg)
        ax.plot(ts['year'], ts['FF_c13_mean'], color=col, ls=ls, lw=2.2, label=lbl)
        ax.fill_between(ts['year'], ts['FF_c13_mean']-ts['FF_c13_std'],
                        ts['FF_c13_mean']+ts['FF_c13_std'], alpha=0.12, color=col)
    delta = s2['OH13C_CANTRELL']['mean_ff_c13'] - s2['OH13C_SAUERESSIG']['mean_ff_c13']
    ax.set_ylabel('FF from δ¹³C [Tg yr⁻¹]', fontsize=11)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_title(f'(a) δ¹³C inversion: ΔFF = {delta:+.0f} Tg/yr', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.25)

    # (b) δD
    ax = axes[1]
    for cfg, lbl, col, ls in [
        ('OH13C_SAUERESSIG', 'Saueressig', PAL['sau'], '-'),
        ('OH13C_CANTRELL',   'Cantrell',   PAL['can'], '-'),
        ('ALL_SAMPLED',      'All sampled', 'gray',    '--'),
    ]:
        ts = load_ts('2x2_', cfg)
        ax.plot(ts['year'], ts['FF_dD_mean'], color=col, ls=ls, lw=2.2, label=lbl)
        ax.fill_between(ts['year'], ts['FF_dD_mean']-ts['FF_dD_std'],
                        ts['FF_dD_mean']+ts['FF_dD_std'], alpha=0.12, color=col)
    ax.set_ylabel('FF from δD [Tg yr⁻¹]', fontsize=11)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_title('(b) δD inversion: ΔFF = 0 Tg/yr (immune)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.25)

    fig.suptitle('2×2 One-Box: Saueressig vs Cantrell OH-¹³C KIE', fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGS / "fig7_2x2_saueressig_cantrell.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 7 ✓  2×2 Saueressig vs Cantrell")


# ====================================================================
# FIGURE 8 — Grand Comparison: Stacked Bars Across 3 Architectures
# ====================================================================
def fig8():
    fig, ax = plt.subplots(figsize=(14, 6), dpi=DPI)

    def decompose(base_s2, fix13c_s2, fixd_s2, allkie_s2):
        oh13c = max(100*(base_s2 - fix13c_s2)/base_s2, 0)
        ohd   = max(100*(base_s2 - fixd_s2)/base_s2, 0)
        allk  = 100*(base_s2 - allkie_s2)/base_s2
        other = max(allk - oh13c - ohd, 0)
        resid = max(100 - allk, 0)
        return oh13c, ohd, other, resid

    rows = []
    # 2×2 δ¹³C
    b = s2['ALL_SAMPLED']['sigma_ff_c13']**2
    rows.append(('2×2 δ¹³C-only', s2['ALL_SAMPLED']['sigma_ff_c13'],
                 *decompose(b, s2['FIX_OH13C']['sigma_ff_c13']**2,
                            s2['FIX_OHD']['sigma_ff_c13']**2,
                            s2['ALL_KIE_FIXED']['sigma_ff_c13']**2)))
    # 2×2 δD
    b = s2['ALL_SAMPLED']['sigma_ff_dD']**2
    rows.append(('2×2 δD-only', s2['ALL_SAMPLED']['sigma_ff_dD'],
                 *decompose(b, s2['FIX_OH13C']['sigma_ff_dD']**2,
                            s2['FIX_OHD']['sigma_ff_dD']**2,
                            s2['ALL_KIE_FIXED']['sigma_ff_dD']**2)))
    # 3×3
    b = s3['ALL_SAMPLED']['sigma_ff']**2
    rows.append(('3×3 dual-isotope', s3['ALL_SAMPLED']['sigma_ff'],
                 *decompose(b, s3['FIX_OH13C']['sigma_ff']**2,
                            s3['FIX_OHD']['sigma_ff']**2,
                            s3['ALL_KIE_FIXED']['sigma_ff']**2)))

    y = np.arange(len(rows))
    h = 0.5
    for i, (lbl, sig, oh13c, ohd, other, resid) in enumerate(rows):
        left = 0
        for val, col, elbl in [(oh13c, PAL['oh13c'], 'OH-¹³C'),
                                (ohd,   PAL['ohd'],   'OH-D'),
                                (other, PAL['other'],  'Other KIE'),
                                (resid, PAL['resid'],  'Src sigs + noise')]:
            ax.barh(i, val, h, left=left, color=col, edgecolor='black', lw=0.5,
                    label=elbl if i == 0 else None)
            if val > 4:
                tcol = 'white' if col not in [PAL['resid']] else 'black'
                ax.text(left + val/2, i, f'{val:.0f}%', ha='center', va='center',
                        fontsize=9.5, fontweight='bold', color=tcol)
            left += val
        ax.text(102, i, f'σ = {sig:.0f}', va='center', fontsize=11, fontweight='bold')

    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=12)
    ax.set_xlabel('% of σ²(FF)', fontsize=12)
    ax.set_xlim(0, 118)
    ax.set_title('OH KIE Variance Attribution Across Model Architectures', fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9.5, ncol=2, framealpha=0.9)
    ax.grid(axis='x', alpha=0.25); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(FIGS / "fig8_grand_comparison.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 8 ✓  Grand comparison (3 architectures)")


# ====================================================================
# FIGURE 9 — Level Shift Comparison Across Architectures
# ====================================================================
def fig9():
    fig, ax = plt.subplots(figsize=(10, 6), dpi=DPI)
    models = ['2×2 δ¹³C-only', '2×2 δD-only', '3×3 dual-isotope']
    sau = [s2['OH13C_SAUERESSIG']['mean_ff_c13'],
           s2['OH13C_SAUERESSIG']['mean_ff_dD'],
           s3['OH13C_SAUERESSIG']['mean_ff']]
    can = [s2['OH13C_CANTRELL']['mean_ff_c13'],
           s2['OH13C_CANTRELL']['mean_ff_dD'],
           s3['OH13C_CANTRELL']['mean_ff']]
    shifts = [c - s for s, c in zip(sau, can)]

    x = np.arange(3)
    w = 0.28
    bars_s = ax.bar(x - w/2, sau, w, label='Saueressig (1.0039)', color=PAL['sau'],
                    edgecolor='black', lw=0.8)
    bars_c = ax.bar(x + w/2, can, w, label='Cantrell (1.0054)', color=PAL['can'],
                    edgecolor='black', lw=0.8)

    for i in range(3):
        mid = (sau[i] + can[i]) / 2
        offset = max(sau[i], can[i]) + 5
        ax.annotate('', xy=(i+w/2, can[i]), xytext=(i-w/2, sau[i]),
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
        ax.text(i, mid, f'Δ={shifts[i]:+.0f}', ha='center', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='wheat', alpha=0.85))

    ax.set_xticks(x); ax.set_xticklabels(models, fontsize=11)
    ax.set_ylabel('Mean FF Emissions [Tg yr⁻¹]', fontsize=12)
    ax.set_title('OH-¹³C KIE Level Shift: Saueressig vs Cantrell', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10); ax.grid(axis='y', alpha=0.25)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIGS / "fig9_level_shift_comparison.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 9 ✓  Level shift comparison across architectures")


# ====================================================================
# FIGURE 10 — Cross-Isotope Contamination Proof (2×2 clean separation)
# ====================================================================
def fig10():
    """
    Demonstrates that in the decoupled 2×2, each isotope is sensitive ONLY
    to its own OH KIE, with exactly 0% cross-contamination.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), dpi=DPI)

    # Compute percentages
    # δ¹³C side
    b_c = s2['ALL_SAMPLED']['sigma_ff_c13']**2
    oh13c_on_c13 = 100*(b_c - s2['FIX_OH13C']['sigma_ff_c13']**2) / b_c
    ohd_on_c13   = 100*(b_c - s2['FIX_OHD']['sigma_ff_c13']**2) / b_c

    # δD side
    b_d = s2['ALL_SAMPLED']['sigma_ff_dD']**2
    oh13c_on_dD  = 100*(b_d - s2['FIX_OH13C']['sigma_ff_dD']**2) / b_d
    ohd_on_dD    = 100*(b_d - s2['FIX_OHD']['sigma_ff_dD']**2) / b_d

    # (a) Matrix-style heatmap
    ax = axes[0]
    mat = np.array([[oh13c_on_c13, ohd_on_c13],
                    [oh13c_on_dD, ohd_on_dD]])
    im = ax.imshow(mat, cmap='RdYlGn_r', vmin=-2, vmax=30, aspect='auto')
    ax.set_xticks([0, 1]); ax.set_xticklabels(['OH-¹³C KIE', 'OH-D KIE'], fontsize=12)
    ax.set_yticks([0, 1]); ax.set_yticklabels(['δ¹³C → FF', 'δD → FF'], fontsize=12)
    for i in range(2):
        for j in range(2):
            val = mat[i, j]
            txt = f'{val:.1f}%' if abs(val) > 0.05 else '0.0%'
            col = 'white' if val > 15 else 'black'
            fontsize = 16 if abs(val) > 0.05 else 14
            ax.text(j, i, txt, ha='center', va='center', fontsize=fontsize,
                    fontweight='bold', color=col)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('% of σ²(FF) attributable', fontsize=10)
    ax.set_title('(a) Cross-isotope KIE sensitivity matrix\n(2×2 decoupled)', fontsize=12, fontweight='bold')

    # (b) Same for 3×3 for comparison
    ax = axes[1]
    b3 = s3['ALL_SAMPLED']['sigma_ff']**2
    oh13c_3x3 = max(100*(b3 - s3['FIX_OH13C']['sigma_ff']**2)/b3, 0)
    ohd_3x3   = 100*(b3 - s3['FIX_OHD']['sigma_ff']**2)/b3

    # Build a comparison bar chart
    categories = ['OH-¹³C → δ¹³C FF', 'OH-D → δD FF',
                  'OH-¹³C → δ¹³C FF\n(cross: 0%)', 'OH-D → δD FF',
                  'OH-¹³C → 3×3 FF', 'OH-D → 3×3 FF']
    # Simplified: show 3 grouped bars
    x = np.arange(3)
    w = 0.35
    own_vals = [oh13c_on_c13, ohd_on_dD, ohd_3x3]
    cross_vals = [ohd_on_c13, oh13c_on_dD, oh13c_3x3]
    own_labels = ['Own OH KIE', 'Own OH KIE', 'OH-D (dominant)']
    cross_labels = ['Cross OH KIE', 'Cross OH KIE', 'OH-¹³C (minor)']

    bars1 = ax.bar(x - w/2, own_vals, w, label='Primary pathway',
                   color=[PAL['oh13c'], PAL['dD'], PAL['ohd']], edgecolor='black', lw=0.8)
    bars2 = ax.bar(x + w/2, cross_vals, w, label='Cross pathway',
                   color=['#cccccc', '#cccccc', '#cccccc'], edgecolor='black', lw=0.8)

    for b, v in zip(bars1, own_vals):
        ax.text(b.get_x()+b.get_width()/2, v+0.5, f'{v:.1f}%', ha='center',
                fontsize=10, fontweight='bold')
    for b, v in zip(bars2, cross_vals):
        ax.text(b.get_x()+b.get_width()/2, max(v, 0)+0.5, f'{v:.1f}%', ha='center',
                fontsize=10, fontweight='bold', color='gray')

    ax.set_xticks(x); ax.set_xticklabels(['2×2 δ¹³C\ninversion', '2×2 δD\ninversion',
                                            '3×3 dual\ninversion'], fontsize=10)
    ax.set_ylabel('% of σ²(FF)', fontsize=11)
    ax.set_title('(b) Primary vs cross-pathway KIE contribution', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9.5); ax.grid(axis='y', alpha=0.25)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_ylim(0, 35)

    fig.suptitle('Cross-Isotope Contamination: Zero in 2×2, Non-Zero in 3×3',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGS / "fig10_cross_contamination.png", dpi=DPI, bbox_inches='tight')
    plt.close()
    print("Fig 10 ✓  Cross-isotope contamination proof")


# ====================================================================
# RUN ALL
# ====================================================================
if __name__ == "__main__":
    print(f"Generating figures from {RESULTS}/\n")
    fig1(); fig2(); fig3(); fig4(); fig5()
    fig6(); fig7(); fig8(); fig9(); fig10()
    print(f"\n{'='*50}")
    print(f"All 10 figures saved to {FIGS}/")
    print(f"{'='*50}")
