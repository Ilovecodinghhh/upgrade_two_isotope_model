#!/usr/bin/env python3
"""Figures for the 2×2 one-box OH KIE importance experiment."""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS = Path(__file__).resolve().parent / "results"
FIGS = Path(__file__).resolve().parent / "figures"
FIGS.mkdir(exist_ok=True)

with open(RESULTS / "summary.json") as f:
    s3x3 = json.load(f)
with open(RESULTS / "2x2_summary.json") as f:
    s2x2 = json.load(f)

# Colors
C13C = '#d7191c'
CD   = '#2c7bb6'

# ====================================================================
# Figure 6: 2×2 variance attribution — side-by-side for δ¹³C and δD
# ====================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=200, sharey=True)

configs = ['ALL_SAMPLED', 'FIX_OH13C', 'FIX_OHD', 'FIX_BOTH_OH', 'ALL_KIE_FIXED']
labels = ['All\nsampled', 'Fix\nOH-¹³C', 'Fix\nOH-D', 'Fix both\nOH', 'Fix all\nKIE']

# Panel A: δ¹³C-derived FF
ax = axes[0]
vals = [s2x2[c]['sigma_ff_c13'] for c in configs]
s2_base = s2x2['ALL_SAMPLED']['sigma_ff_c13']**2
reds = [0] + [100*(s2_base - s2x2[c]['sigma_ff_c13']**2)/s2_base for c in configs[1:]]
colors_a = ['gray', C13C, '#abd9e9', '#fdae61', '#1a9641']
bars = ax.bar(range(5), vals, color=colors_a, edgecolor='black', lw=0.8, width=0.65)
for i, (b, v, r) in enumerate(zip(bars, vals, reds)):
    ax.text(b.get_x()+b.get_width()/2, v+0.5, f'{v:.1f}', ha='center', fontsize=10, fontweight='bold')
    if i > 0 and abs(r) > 0.5:
        ax.text(b.get_x()+b.get_width()/2, v/2, f'{r:+.0f}%', ha='center', fontsize=9,
                color='white', fontweight='bold')
ax.set_xticks(range(5)); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('σ(FF) [Tg yr⁻¹]', fontsize=11)
ax.set_title('(a) δ¹³C-only inversion → FF uncertainty', fontsize=12, fontweight='bold')
ax.set_ylim(0, 48)
ax.grid(axis='y', alpha=0.3)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# Panel B: δD-derived FF
ax = axes[1]
vals = [s2x2[c]['sigma_ff_dD'] for c in configs]
s2_base = s2x2['ALL_SAMPLED']['sigma_ff_dD']**2
reds = [0] + [100*(s2_base - s2x2[c]['sigma_ff_dD']**2)/s2_base for c in configs[1:]]
colors_b = ['gray', '#abd9e9', CD, '#fdae61', '#1a9641']
bars = ax.bar(range(5), vals, color=colors_b, edgecolor='black', lw=0.8, width=0.65)
for i, (b, v, r) in enumerate(zip(bars, vals, reds)):
    ax.text(b.get_x()+b.get_width()/2, v+0.5, f'{v:.1f}', ha='center', fontsize=10, fontweight='bold')
    if i > 0 and abs(r) > 0.5:
        ax.text(b.get_x()+b.get_width()/2, v/2, f'{r:+.0f}%', ha='center', fontsize=9,
                color='white', fontweight='bold')
ax.set_xticks(range(5)); ax.set_xticklabels(labels, fontsize=9)
ax.set_title('(b) δD-only inversion → FF uncertainty', fontsize=12, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(FIGS / "fig6_2x2_variance_attribution.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 6: 2×2 variance attribution")

# ====================================================================
# Figure 7: Saueressig vs Cantrell level shift in δ¹³C 2×2
# ====================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=200)

ax = axes[0]
for cfg, lbl, col, ls in [
    ('OH13C_SAUERESSIG', 'Saueressig (1.0039)', '#2166ac', '-'),
    ('OH13C_CANTRELL',   'Cantrell (1.0054)',   '#b2182b', '-'),
    ('ALL_SAMPLED',      'All sampled',         'gray',    '--'),
]:
    ts = pd.read_csv(RESULTS / f"2x2_{cfg}_timeseries.csv")
    ax.plot(ts['year'], ts['FF_c13_mean'], color=col, ls=ls, lw=2, label=lbl)
    ax.fill_between(ts['year'], ts['FF_c13_mean']-ts['FF_c13_std'],
                    ts['FF_c13_mean']+ts['FF_c13_std'], alpha=0.15, color=col)
ax.set_ylabel('FF from δ¹³C [Tg yr⁻¹]', fontsize=11)
ax.set_xlabel('Year', fontsize=11)
ax.set_title('(a) δ¹³C inversion: OH-¹³C KIE shifts FF by 38 Tg/yr', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

ax = axes[1]
for cfg, lbl, col, ls in [
    ('OH13C_SAUERESSIG', 'Saueressig', '#2166ac', '-'),
    ('OH13C_CANTRELL',   'Cantrell',   '#b2182b', '-'),
    ('ALL_SAMPLED',      'All sampled','gray',    '--'),
]:
    ts = pd.read_csv(RESULTS / f"2x2_{cfg}_timeseries.csv")
    ax.plot(ts['year'], ts['FF_dD_mean'], color=col, ls=ls, lw=2, label=lbl)
    ax.fill_between(ts['year'], ts['FF_dD_mean']-ts['FF_dD_std'],
                    ts['FF_dD_mean']+ts['FF_dD_std'], alpha=0.15, color=col)
ax.set_ylabel('FF from δD [Tg yr⁻¹]', fontsize=11)
ax.set_xlabel('Year', fontsize=11)
ax.set_title('(b) δD inversion: OH-¹³C KIE has no effect', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(FIGS / "fig7_2x2_saueressig_cantrell.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 7: Saueressig vs Cantrell in 2×2")

# ====================================================================
# Figure 8: Grand comparison — all three model architectures
# ====================================================================
fig, ax = plt.subplots(figsize=(14, 6), dpi=200)

# Data: (label, isotope_system, all_sampled_sigma, oh13c_var%, ohd_var%, other_kie%, residual%)
rows = []

# 2×2 δ¹³C
s2b = s2x2['ALL_SAMPLED']['sigma_ff_c13']**2
oh13c_pct = 100*(s2b - s2x2['FIX_OH13C']['sigma_ff_c13']**2)/s2b
ohd_pct   = 100*(s2b - s2x2['FIX_OHD']['sigma_ff_c13']**2)/s2b
allkie_pct = 100*(s2b - s2x2['ALL_KIE_FIXED']['sigma_ff_c13']**2)/s2b
rows.append(('2×2 δ¹³C-only', s2x2['ALL_SAMPLED']['sigma_ff_c13'],
             oh13c_pct, ohd_pct, allkie_pct - max(oh13c_pct, ohd_pct), 100-allkie_pct))

# 2×2 δD
s2b = s2x2['ALL_SAMPLED']['sigma_ff_dD']**2
oh13c_pct = 100*(s2b - s2x2['FIX_OH13C']['sigma_ff_dD']**2)/s2b
ohd_pct   = 100*(s2b - s2x2['FIX_OHD']['sigma_ff_dD']**2)/s2b
allkie_pct = 100*(s2b - s2x2['ALL_KIE_FIXED']['sigma_ff_dD']**2)/s2b
rows.append(('2×2 δD-only', s2x2['ALL_SAMPLED']['sigma_ff_dD'],
             oh13c_pct, ohd_pct, allkie_pct - max(oh13c_pct, ohd_pct), 100-allkie_pct))

# 3×3 dual
s2b = s3x3['ALL_SAMPLED']['sigma_ff']**2
oh13c_pct = 100*(s2b - s3x3['FIX_OH13C']['sigma_ff']**2)/s2b
ohd_pct   = 100*(s2b - s3x3['FIX_OHD']['sigma_ff']**2)/s2b
allkie_pct = 100*(s2b - s3x3['ALL_KIE_FIXED']['sigma_ff']**2)/s2b
rows.append(('3×3 dual-isotope', s3x3['ALL_SAMPLED']['sigma_ff'],
             oh13c_pct, ohd_pct, allkie_pct - max(oh13c_pct, ohd_pct), 100-allkie_pct))

labels = [r[0] for r in rows]
sigmas = [r[1] for r in rows]
oh13c_vals = [max(r[2], 0) for r in rows]
ohd_vals = [max(r[3], 0) for r in rows]
other_vals = [max(r[4], 0) for r in rows]
resid_vals = [max(r[5], 0) for r in rows]

y = np.arange(len(rows))
h = 0.55

# Stacked horizontal bars (as % of σ²)
ax.barh(y, oh13c_vals, h, label='OH-¹³C KIE', color=C13C, edgecolor='black', lw=0.5)
left = np.array(oh13c_vals)
ax.barh(y, ohd_vals, h, left=left, label='OH-D KIE', color=CD, edgecolor='black', lw=0.5)
left += np.array(ohd_vals)
ax.barh(y, other_vals, h, left=left, label='Other KIE (Cl,Strat,Soil)', color='#7570b3', edgecolor='black', lw=0.5)
left += np.array(other_vals)
ax.barh(y, resid_vals, h, left=left, label='Source sigs + data noise', color='#e7e1ef', edgecolor='black', lw=0.5)

# Annotate σ values
for i, (sig, lab) in enumerate(zip(sigmas, labels)):
    ax.text(102, i, f'σ(FF) = {sig:.0f}', va='center', fontsize=11, fontweight='bold')

# Annotate significant percentages
for i in range(len(rows)):
    x = 0
    for val, col in [(oh13c_vals[i], 'white'), (ohd_vals[i], 'white'),
                      (other_vals[i], 'white'), (resid_vals[i], 'black')]:
        if val > 5:
            ax.text(x + val/2, i, f'{val:.0f}%', ha='center', va='center',
                    fontsize=9, fontweight='bold', color=col)
        x += val

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=12)
ax.set_xlabel('% of σ²(FF)', fontsize=12)
ax.set_xlim(0, 120)
ax.set_title('OH KIE Variance Attribution: 2×2 (separate isotopes) vs 3×3 (dual-isotope)',
             fontsize=13, fontweight='bold')
ax.legend(loc='lower right', fontsize=9, ncol=2)
ax.grid(axis='x', alpha=0.3)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(FIGS / "fig8_grand_comparison.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 8: Grand comparison (3 architectures)")

# ====================================================================
# Figure 9: Level shift comparison — δ¹³C 2×2 vs 3×3
# ====================================================================
fig, ax = plt.subplots(figsize=(10, 5.5), dpi=200)

models = ['2×2 δ¹³C-only', '2×2 δD-only', '3×3 dual-isotope']
sau_vals = [s2x2['OH13C_SAUERESSIG']['mean_ff_c13'],
            s2x2['OH13C_SAUERESSIG']['mean_ff_dD'],
            s3x3['OH13C_SAUERESSIG']['mean_ff']]
can_vals = [s2x2['OH13C_CANTRELL']['mean_ff_c13'],
            s2x2['OH13C_CANTRELL']['mean_ff_dD'],
            s3x3['OH13C_CANTRELL']['mean_ff']]
shifts = [c - s for s, c in zip(sau_vals, can_vals)]

x = np.arange(3)
w = 0.3
bars_s = ax.bar(x - w/2, sau_vals, w, label='Saueressig (1.0039)', color='#2166ac',
                edgecolor='black', lw=0.8)
bars_c = ax.bar(x + w/2, can_vals, w, label='Cantrell (1.0054)', color='#b2182b',
                edgecolor='black', lw=0.8)

for i, (bs, bc, sh) in enumerate(zip(bars_s, bars_c, shifts)):
    mid = (sau_vals[i] + can_vals[i]) / 2
    ax.annotate('', xy=(i+w/2, can_vals[i]), xytext=(i-w/2, sau_vals[i]),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text(i, mid, f'Δ = {abs(sh):.0f}', ha='center', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='wheat', alpha=0.8))

ax.set_xticks(x)
ax.set_xticklabels(models, fontsize=11)
ax.set_ylabel('Mean FF Emissions [Tg yr⁻¹]', fontsize=12)
ax.set_title('OH-¹³C KIE Level Shift Across Model Architectures', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(axis='y', alpha=0.3)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(FIGS / "fig9_level_shift_comparison.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 9: Level shift comparison")

print(f"\nAll figures → {FIGS}/")
