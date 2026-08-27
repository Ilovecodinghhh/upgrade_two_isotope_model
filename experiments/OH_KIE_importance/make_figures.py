#!/usr/bin/env python3
"""Generate all figures for the OH KIE Importance experiment."""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

RESULTS = Path(__file__).resolve().parent / "results"
FIGS = Path(__file__).resolve().parent / "figures"
FIGS.mkdir(exist_ok=True)

with open(RESULTS / "summary.json") as f:
    summary = json.load(f)

# ─── Color palette (colorblind-safe) ───
C_BASE = '#2c7bb6'   # blue
C_13C  = '#d7191c'   # red
C_D    = '#fdae61'    # orange
C_BOTH = '#abd9e9'    # light blue
C_ALL  = '#1a9641'    # green
C_SAU  = '#fee08b'    # yellow
C_CAN  = '#d73027'    # dark red

# ====================================================================
# Figure 1: Variance Attribution Bar Chart
# ====================================================================
fig, ax = plt.subplots(figsize=(10, 6), dpi=200)

s2_base = summary['ALL_SAMPLED']['sigma_ff']**2

configs = ['ALL_SAMPLED', 'FIX_OH13C', 'FIX_OHD', 'FIX_BOTH_OH', 'ALL_KIE_FIXED']
labels = ['All sampled\n(baseline)', 'Fix OH-¹³C\nonly', 'Fix OH-D\nonly', 
          'Fix both\nOH KIEs', 'Fix all\n8 KIE params']
colors = [C_BASE, C_13C, C_D, C_BOTH, C_ALL]
sigmas = [summary[c]['sigma_ff'] for c in configs]
reductions = [0] + [100*(s2_base - summary[c]['sigma_ff']**2)/s2_base for c in configs[1:]]

bars = ax.bar(range(len(configs)), sigmas, color=colors, edgecolor='black', linewidth=0.8, width=0.7)

# Annotate
for i, (bar, sig, red) in enumerate(zip(bars, sigmas, reductions)):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
            f'σ = {sig:.1f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    if i > 0:
        txt = f'Δσ² = {red:+.1f}%'
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                txt, ha='center', va='center', fontsize=9, 
                color='white' if i != 3 else 'black', fontweight='bold')

ax.set_xticks(range(len(configs)))
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel('σ(FF) [Tg yr⁻¹]', fontsize=12)
ax.set_title('One-Box 3×3 Dual-Isotope: FF Uncertainty by KIE Configuration', fontsize=13, fontweight='bold')
ax.set_ylim(0, 75)
ax.axhline(summary['ALL_SAMPLED']['sigma_ff'], color='gray', ls='--', alpha=0.5, lw=1)
ax.grid(axis='y', alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(FIGS / "fig1_variance_attribution.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 1 saved: variance attribution bar chart")

# ====================================================================
# Figure 2: OH-13C Level Shift (Saueressig vs Cantrell)
# ====================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=200)

# Panel A: Time series
ax = axes[0]
for config, label, color, ls in [
    ('OH13C_SAUERESSIG', 'Saueressig (1.0039)', '#2166ac', '-'),
    ('OH13C_CANTRELL',   'Cantrell (1.0054)',   '#b2182b', '-'),
    ('ALL_SAMPLED',      'All sampled',         'gray',    '--'),
]:
    ts = pd.read_csv(RESULTS / f"{config}_timeseries.csv")
    ax.plot(ts['year'], ts['FF_mean'], color=color, ls=ls, lw=2, label=label)
    ax.fill_between(ts['year'], ts['FF_mean']-ts['FF_std'], ts['FF_mean']+ts['FF_std'],
                    alpha=0.15, color=color)

ax.set_xlabel('Year', fontsize=11)
ax.set_ylabel('Fossil Fuel Emissions [Tg yr⁻¹]', fontsize=11)
ax.set_title('(a) OH-¹³C KIE shifts FF level, not spread', fontsize=12, fontweight='bold')
ax.legend(fontsize=9, loc='upper left')
ax.grid(alpha=0.3)
ax.axhline(0, color='black', lw=0.5)
ax.set_ylim(-100, 200)

# Panel B: Level shift summary
ax = axes[1]
configs_b = ['OH13C_SAUERESSIG', 'FIX_OH13C', 'OH13C_CANTRELL']
labels_b = ['Saueressig\n(α=1.0039)', 'Midpoint\n(α=1.00465)', 'Cantrell\n(α=1.0054)']
colors_b = ['#2166ac', 'gray', '#b2182b']
means = [summary[c]['mean_ff'] for c in configs_b]
sigmas_b = [summary[c]['sigma_ff'] for c in configs_b]

bars = ax.bar(range(3), means, yerr=sigmas_b, color=colors_b, 
              edgecolor='black', linewidth=0.8, width=0.6,
              capsize=5, error_kw={'lw': 1.5})
              
for i, (bar, m, s) in enumerate(zip(bars, means, sigmas_b)):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + s + 2,
            f'{m:.0f} ± {s:.0f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

# Annotate the shift
ax.annotate('', xy=(2, means[2]), xytext=(0, means[0]),
            arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
ax.text(1, (means[0]+means[2])/2, f'Δ = {means[2]-means[0]:.0f} Tg/yr',
        ha='center', va='center', fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

ax.set_xticks(range(3))
ax.set_xticklabels(labels_b, fontsize=10)
ax.set_ylabel('Mean FF Emissions [Tg yr⁻¹]', fontsize=11)
ax.set_title('(b) OH-¹³C KIE: level shift = 17 Tg/yr', fontsize=12, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(FIGS / "fig2_oh13c_level_shift.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 2 saved: OH-13C level shift")

# ====================================================================
# Figure 3: Variance Decomposition Pie + Stacked Bar
# ====================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=200)

# Compute individual contributions
s2 = {c: summary[c]['sigma_ff']**2 for c in summary}
total_var = s2['ALL_SAMPLED']

# OH-D contribution
oh_d_var = s2['ALL_SAMPLED'] - s2['FIX_OHD']
# OH-13C contribution  
oh_13c_var = s2['ALL_SAMPLED'] - s2['FIX_OH13C']
# Other KIE (non-OH)
both_oh_var = s2['ALL_SAMPLED'] - s2['FIX_BOTH_OH']
all_kie_var = s2['ALL_SAMPLED'] - s2['ALL_KIE_FIXED']
other_kie_var = all_kie_var - both_oh_var
# Residual (source sigs + data noise + lifetime)
residual_var = s2['ALL_KIE_FIXED']

# Panel A: Pie chart
ax = axes[0]
sizes = [max(oh_d_var, 0), max(other_kie_var, 0), max(oh_13c_var, 0), residual_var]
labels_p = [
    f'OH-D KIE\n({100*oh_d_var/total_var:.1f}%)',
    f'Other KIE\n(Cl, Strat, Soil)\n({100*other_kie_var/total_var:.1f}%)',
    f'OH-¹³C KIE\n({100*oh_13c_var/total_var:.1f}%)',
    f'Source sigs +\ndata noise\n({100*residual_var/total_var:.1f}%)',
]
colors_p = [C_D, '#7570b3', C_13C, '#e7e1ef']
# Suppress tiny/negative slices
if oh_13c_var <= 0:
    sizes[2] = 0
    labels_p[2] = f'OH-¹³C KIE\n(< 1%)'

explode = (0.05, 0, 0, 0)
wedges, texts, autotexts = ax.pie(sizes, labels=labels_p, colors=colors_p,
                                    explode=explode, startangle=90,
                                    autopct=lambda p: f'{p:.0f}%' if p > 2 else '',
                                    pctdistance=0.7, textprops={'fontsize': 9})
ax.set_title('(a) Variance decomposition of σ²(FF)', fontsize=12, fontweight='bold')

# Panel B: Cumulative variance reduction
ax = axes[1]
steps = ['Baseline', '+Fix OH-D', '+Fix OH-¹³C', '+Fix all KIE']
s2_vals = [total_var, s2['FIX_OHD'], s2['FIX_BOTH_OH'], s2['ALL_KIE_FIXED']]
sigma_vals = [np.sqrt(v) for v in s2_vals]
reductions_cum = [0] + [100*(total_var - v)/total_var for v in s2_vals[1:]]

bar_colors = [C_BASE, C_D, C_BOTH, C_ALL]
bars = ax.bar(range(4), sigma_vals, color=bar_colors, edgecolor='black', linewidth=0.8, width=0.65)

for i, (bar, sv, rc) in enumerate(zip(bars, sigma_vals, reductions_cum)):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
            f'σ = {sv:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    if i > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                f'−{rc:.0f}%', ha='center', va='center', fontsize=10,
                color='white', fontweight='bold')

ax.set_xticks(range(4))
ax.set_xticklabels(steps, fontsize=9.5)
ax.set_ylabel('σ(FF) [Tg yr⁻¹]', fontsize=11)
ax.set_title('(b) Cumulative σ(FF) reduction', fontsize=12, fontweight='bold')
ax.set_ylim(0, 75)
ax.grid(axis='y', alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(FIGS / "fig3_variance_decomposition.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 3 saved: variance decomposition")

# ====================================================================
# Figure 4: All-config FF time series comparison
# ====================================================================
fig, ax = plt.subplots(figsize=(12, 6), dpi=200)

plot_configs = [
    ('ALL_SAMPLED',   'All KIE sampled (baseline)', C_BASE, '-',  2.5),
    ('FIX_OHD',       'Fix OH-D only',              C_D,    '-',  2.0),
    ('FIX_BOTH_OH',   'Fix both OH KIEs',           C_BOTH, '--', 1.8),
    ('ALL_KIE_FIXED', 'Fix all 8 KIEs',             C_ALL,  ':',  2.0),
]

for config, label, color, ls, lw in plot_configs:
    ts = pd.read_csv(RESULTS / f"{config}_timeseries.csv")
    ax.plot(ts['year'], ts['FF_mean'], color=color, ls=ls, lw=lw, label=label)
    ax.fill_between(ts['year'], ts['FF_mean']-ts['FF_std'], ts['FF_mean']+ts['FF_std'],
                    alpha=0.10, color=color)

ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('Fossil Fuel Emissions [Tg yr⁻¹]', fontsize=12)
ax.set_title('One-Box 3×3 Dual-Isotope: FF Emissions Under Different KIE Configurations', 
             fontsize=13, fontweight='bold')
ax.legend(fontsize=10, loc='upper left')
ax.grid(alpha=0.3)
ax.axhline(0, color='black', lw=0.5)
ax.set_ylim(-100, 200)

plt.tight_layout()
plt.savefig(FIGS / "fig4_ff_timeseries_comparison.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 4 saved: FF time series comparison")

# ====================================================================
# Figure 5: Dual-role summary diagram
# ====================================================================
fig, ax = plt.subplots(figsize=(11, 5), dpi=200)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.5, 2.5)
ax.axis('off')

# Title
ax.text(1.5, 2.3, 'OH KIE Dual Role in One-Box 3×3 Inversion', 
        fontsize=15, fontweight='bold', ha='center', va='center')

# OH-D box
rect1 = FancyBboxPatch((0, 0.3), 1.2, 1.4, boxstyle="round,pad=0.1",
                        facecolor=C_D, alpha=0.3, edgecolor='black', linewidth=2)
ax.add_patch(rect1)
ax.text(0.6, 1.5, 'OH-D KIE', fontsize=13, fontweight='bold', ha='center')
ax.text(0.6, 1.15, 'U(1.294, 1.327)', fontsize=9, ha='center', family='monospace')
ax.text(0.6, 0.85, f'Controls σ(FF)', fontsize=11, ha='center')
ax.text(0.6, 0.55, f'28% of variance', fontsize=12, ha='center', fontweight='bold', color='#b35806')

# OH-13C box
rect2 = FancyBboxPatch((1.8, 0.3), 1.2, 1.4, boxstyle="round,pad=0.1",
                        facecolor=C_13C, alpha=0.2, edgecolor='black', linewidth=2)
ax.add_patch(rect2)
ax.text(2.4, 1.5, 'OH-¹³C KIE', fontsize=13, fontweight='bold', ha='center')
ax.text(2.4, 1.15, 'U(1.0039, 1.0054)', fontsize=9, ha='center', family='monospace')
ax.text(2.4, 0.85, 'Controls FF level', fontsize=11, ha='center')
ax.text(2.4, 0.55, f'Δ = 17 Tg/yr shift', fontsize=12, ha='center', fontweight='bold', color='#b2182b')

# Arrow annotations
ax.annotate('spread\n(uncertainty)', xy=(0.6, 0.3), xytext=(0.6, -0.2),
            fontsize=10, ha='center', va='top', color='#b35806',
            arrowprops=dict(arrowstyle='->', color='#b35806', lw=1.5))
ax.annotate('bias\n(systematic)', xy=(2.4, 0.3), xytext=(2.4, -0.2),
            fontsize=10, ha='center', va='top', color='#b2182b',
            arrowprops=dict(arrowstyle='->', color='#b2182b', lw=1.5))

plt.tight_layout()
plt.savefig(FIGS / "fig5_dual_role_summary.png", dpi=200, bbox_inches='tight')
plt.close()
print("Fig 5 saved: dual-role summary diagram")

print(f"\nAll figures saved to {FIGS}/")
