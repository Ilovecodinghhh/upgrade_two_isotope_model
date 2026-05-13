#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure improvements from EditorialAssessment_Figure.md
======================================================
1. figM8_before_after_timeseries.png — Before/after filtering FF time series
2. Recolored versions of key figures (Cantrell=teal, Saueressig=coral)
   - figM1_schematic_v3.png
   - figM6_KSR_summary_v3.png
   - figM7_forest_plot_v2.png
   - fig14_temporal_stability_v3.png
   - figM8_before_after_timeseries.png (uses new palette from start)

Requires Phase 9 data. Reruns inversions at N=5000 to get per-year arrays.
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

BASE = Path(__file__).resolve().parent
FIG_DIR = BASE / "figures"
OUT_DIR = BASE / "results" / "phase10_figure_improvements"
FIG_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Import the inversion machinery from phase9
sys.path.insert(0, str(BASE))
from phase9_editorial_fixes import run_inversions, agreement_metrics

# ── Colour palette ──────────────────────────────────────────────────────────
C_CANTRELL = '#2CA6A4'   # teal
C_SAUERESSIG = '#E8655A' # coral
C_REJECT = '#B03A2E'     # dark red for reject paths
C_KEEP = '#1E8449'       # dark green for keep paths
C_NEUTRAL = '#888888'    # gray for neutral elements
C_PURPLE = '#7D3C98'     # purple for discriminant annotations
C_WLS_FAIL = '#E74C3C'   # red for WLS failure bars
C_FILTER_OK = '#27AE60'  # green for filter success bars

# ── Load Phase 9 + Phase 8 results ─────────────────────────────────────────
with open(BASE / "results/phase9_editorial_fixes/high_n_summary.json") as f:
    p9_high_n = json.load(f)
with open(BASE / "results/phase9_editorial_fixes/cl_sensitivity.json") as f:
    p9_cl = json.load(f)
with open(BASE / "results/phase8_fine_thresholds/summary.json") as f:
    p8 = json.load(f)


# ═══════════════════════════════════════════════════════════════════════════
# TASK 1: Before/after FF time series
# ═══════════════════════════════════════════════════════════════════════════
print("[1/5] Running N=5000 inversions for time-series figure...")

CACHE_S = OUT_DIR / "ff_arrays_saueressig.npz"
CACHE_C = OUT_DIR / "ff_arrays_cantrell.npz"

N_ITER = 5000

if CACHE_S.exists() and CACHE_C.exists():
    print("      Using cached arrays...")
    ds = np.load(CACHE_S)
    dc = np.load(CACHE_C)
    FF_c13_S, FF_dD_S, Mic_c13_S = ds['FF_c13'], ds['FF_dD'], ds['Mic_c13']
    FF_c13_C, FF_dD_C, Mic_c13_C = dc['FF_c13'], dc['FF_dD'], dc['Mic_c13']
    years = ds['years']
else:
    print("      Running Saueressig inversions (N=5000)...")
    FF_c13_S, FF_dD_S, Mic_c13_S, years = run_inversions('saueressig', N_ITER, seed=42)
    np.savez_compressed(CACHE_S, FF_c13=FF_c13_S, FF_dD=FF_dD_S, Mic_c13=Mic_c13_S, years=years)

    print("      Running Cantrell inversions (N=5000)...")
    FF_c13_C, FF_dD_C, Mic_c13_C, years = run_inversions('cantrell', N_ITER, seed=42)
    np.savez_compressed(CACHE_C, FF_c13=FF_c13_C, FF_dD=FF_dD_C, Mic_c13=Mic_c13_C, years=years)

print("      Computing filtered ensembles (T=90)...")
T = 90  # headline threshold

def filter_ensemble(FF_c13, FF_dD, threshold, year_frac=0.80):
    """Return mask of iterations passing the agreement filter."""
    n, N = FF_c13.shape
    valid = ~(np.isnan(FF_c13) | np.isnan(FF_dD))
    agree = (np.abs(FF_c13 - FF_dD) < threshold) & valid
    good = agree.sum(axis=0) >= n * year_frac
    return good

mask_S = filter_ensemble(FF_c13_S, FF_dD_S, T)
mask_C = filter_ensemble(FF_c13_C, FF_dD_C, T)

print(f"      Saueressig: {mask_S.sum()}/{N_ITER} pass ({100*mask_S.mean():.1f}%)")
print(f"      Cantrell:   {mask_C.sum()}/{N_ITER} pass ({100*mask_C.mean():.1f}%)")

# Compute percentiles for the time series
def ts_stats(FF, mask=None):
    """Median + 5th/95th for each year."""
    if mask is not None and mask.sum() >= 20:
        arr = FF[:, mask]
    else:
        arr = FF
    arr = np.clip(arr, 0, None)
    med = np.nanmedian(arr, axis=1)
    lo = np.nanpercentile(arr, 5, axis=1)
    hi = np.nanpercentile(arr, 95, axis=1)
    return med, lo, hi

# ── Figure M8: Before/after FF time series ──────────────────────────────
print("      Generating figM8_before_after_timeseries...")

fig, axes = plt.subplots(2, 1, figsize=(12, 9), dpi=300, sharex=True)

for ax_idx, (label, FF_S, FF_C, mS, mC) in enumerate([
    ('Unfiltered', FF_c13_S, FF_c13_C, None, None),
    (f'Filtered (T = {T} Tg/yr)', FF_c13_S, FF_c13_C, mask_S, mask_C),
]):
    ax = axes[ax_idx]

    med_S, lo_S, hi_S = ts_stats(FF_S, mS)
    med_C, lo_C, hi_C = ts_stats(FF_C, mC)

    ax.fill_between(years[:len(med_S)], lo_S, hi_S, alpha=0.15, color=C_SAUERESSIG)
    ax.plot(years[:len(med_S)], med_S, color=C_SAUERESSIG, lw=2, label='Saueressig (α=1.0039)')

    ax.fill_between(years[:len(med_C)], lo_C, hi_C, alpha=0.15, color=C_CANTRELL)
    ax.plot(years[:len(med_C)], med_C, color=C_CANTRELL, lw=2, label='Cantrell (α=1.0054)')

    # KIE offset annotation
    offset = np.nanmean(med_C - med_S)
    ax.annotate(f'Mean offset: {offset:+.1f} Tg/yr',
                xy=(0.98, 0.95), xycoords='axes fraction',
                ha='right', va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    if mS is not None:
        nS = mS.sum() if mS is not None else N_ITER
        nC = mC.sum() if mC is not None else N_ITER
        ax.annotate(f'n(S)={nS}, n(C)={nC}',
                    xy=(0.98, 0.82), xycoords='axes fraction',
                    ha='right', va='top', fontsize=9, color='gray')

    ax.set_ylabel('Fossil Fuel Emissions (Tg/yr)', fontsize=11)
    ax.set_title(f'({"a" if ax_idx==0 else "b"})  {label}', fontsize=12, weight='bold', loc='left')
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(alpha=0.3)

axes[1].set_xlabel('Year', fontsize=11)
fig.suptitle('Fossil Fuel CH₄ Emission Estimates: Effect of Agreement Filter',
             fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig(FIG_DIR / "figM8_before_after_timeseries.png", dpi=300, bbox_inches='tight')
plt.close()
print("      Saved: figM8_before_after_timeseries.png")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 2: Recolored figM1 schematic (v3)
# ═══════════════════════════════════════════════════════════════════════════
print("[2/5] Generating figM1_schematic_v3 (recolored)...")

fig, ax = plt.subplots(figsize=(13, 7.5), dpi=300)
ax.set_xlim(0, 13)
ax.set_ylim(0, 7.5)
ax.axis('off')

def box(x, y, w, h, text, color, fontsize=10, weight='normal', alpha=0.85, ec='black'):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                       linewidth=1.5, edgecolor=ec, facecolor=color, alpha=alpha)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fontsize, weight=weight)

def arrow(x1, y1, x2, y2, label=None, color='black', lw=1.8,
          style='->,head_width=6,head_length=8', label_offset=0.15, fontsize=8):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                         linewidth=lw, color=color)
    ax.add_patch(a)
    if label:
        ax.text((x1+x2)/2, (y1+y2)/2 + label_offset, label,
                ha='center', fontsize=fontsize, style='italic', color=color)

# Inputs
box(0.3, 5.5, 2.2, 0.9, "Atmospheric obs.\nCH₄, δ¹³C, δD", '#cce5ff', fontsize=9)
box(0.3, 4.2, 2.2, 0.9, "Sink fractions\n& KIE draws", '#cce5ff', fontsize=9)
box(0.3, 2.9, 2.2, 0.9, "Source signatures\n(MC sampled)", '#cce5ff', fontsize=9)
box(0.3, 1.6, 2.2, 0.9, "BB emissions\n(GFEDv4.1s)", '#cce5ff', fontsize=9)

# MC gate
box(3.0, 3.3, 1.7, 1.3, "Monte Carlo\n(N = 5000)", '#ffe5b4', fontsize=10, weight='bold')

# Two inversions
box(5.2, 5.2, 2.4, 1.0, "δ¹³C inversion\n(σ ≈ 4 Tg/yr)", '#d4edda', fontsize=10)
box(5.2, 1.5, 2.4, 1.0, "δD inversion\n(σ ≈ 30+ Tg/yr)", '#d4edda', fontsize=10)

# Outputs
box(8.2, 5.2, 1.5, 1.0, "FF¹³ᶜ, Mic¹³ᶜ", '#f8d7da', fontsize=10)
box(8.2, 1.5, 1.5, 1.0, "FFᴰ", '#f8d7da', fontsize=10)

# Agreement test
box(8.2, 3.3, 1.5, 1.3, "Agree?\n|FF¹³ᶜ − FFᴰ|\n< T", '#fff3cd',
    fontsize=9, weight='bold')

# KEEP path
box(10.3, 3.5, 1.5, 1.0, "✓ KEEP", '#d4edda', fontsize=10, weight='bold')
arrow(9.7, 3.95, 10.3, 3.95, label='YES', color=C_KEEP, fontsize=9)

# REJECT path
box(8.45, 0.2, 1.0, 0.8, "✗ REJECT", '#f8d7da', fontsize=9, weight='bold',
    alpha=0.5, ec=C_REJECT)
arrow(8.95, 3.3, 8.95, 1.0, label='NO', color=C_REJECT, lw=1.5, label_offset=-0.35, fontsize=9)

# Input → MC arrows
for y in [5.95, 4.65, 3.35, 2.05]:
    arrow(2.5, y, 3.0, 3.95 + (y - 3.95) * 0.05, color='gray', lw=1.2)

# MC → inversions
arrow(4.7, 4.3, 5.2, 5.7, color='#333')
arrow(4.7, 3.6, 5.2, 2.0, color='#333')

# Inversions → outputs
arrow(7.6, 5.7, 8.2, 5.7, color='#333')
arrow(7.6, 2.0, 8.2, 2.0, color='#333')

# Outputs → agreement test
arrow(8.95, 5.2, 8.95, 4.6, color='#333')
arrow(8.95, 2.5, 8.95, 3.3, color='#333')

# Retention annotation with new colours
ax.text(11.8, 3.0,
        f"Cantrell: 70.8%\nSaueressig: 35.3%\nΔ = 35.5 pp",
        fontsize=9, ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
arrow(10.8, 3.0, 11.1, 3.0, color='gray', lw=1.0,
      style='->,head_width=4,head_length=5')

# Title
ax.text(6.5, 7.2, 'The Dual-Isotope Agreement Filter',
        ha='center', fontsize=14, weight='bold')
ax.text(6.5, 6.8, 'Two independent isotopic budgets → consistency check → KIE discriminant',
        ha='center', fontsize=10, style='italic', color='dimgray')

plt.savefig(FIG_DIR / "figM1_schematic_v3.png", dpi=300, bbox_inches='tight')
plt.close()
print("      Saved: figM1_schematic_v3.png")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 3: Recolored figM6 KSR summary (v3)
# ═══════════════════════════════════════════════════════════════════════════
print("[3/5] Generating figM6_KSR_summary_v3 (recolored)...")

methods = [
    ('1-box WLS\n($w_D=1$)',      0.20, C_WLS_FAIL, None),
    ('1-box WLS\n($w_D=0.01$)',   0.24, C_WLS_FAIL, None),
    ('2-box WLS\n(fixed)',         0.22, C_WLS_FAIL, None),
    ('Filter\nT=200',              1.08, C_NEUTRAL,  None),
    ('Filter\nT=150',              1.08, C_FILTER_OK, (1.06, 1.10)),
    ('Filter\nT=100',              1.08, C_FILTER_OK, (1.03, 1.15)),
    ('Filter\nT=90\n★ headline',   1.12, C_FILTER_OK, (1.02, 1.24)),
    ('Filter\nT=50',               1.59, C_FILTER_OK, (0.93, 4.89)),
]

labels = [m[0] for m in methods]
ksrs = [m[1] for m in methods]
colors = [m[2] for m in methods]

yerr_lo = [m[1] - m[3][0] if m[3] else 0 for m in methods]
yerr_hi = [m[3][1] - m[1] if m[3] else 0 for m in methods]

fig, ax = plt.subplots(figsize=(12, 6), dpi=300)
bars = ax.bar(range(len(methods)), ksrs, color=colors, edgecolor='black', linewidth=1.2,
              yerr=[yerr_lo, yerr_hi], capsize=5, error_kw={'lw': 1.5})

for i, (bar, k) in enumerate(zip(bars, ksrs)):
    y_off = yerr_hi[i] if yerr_hi[i] > 0 else 0.06
    ax.text(bar.get_x() + bar.get_width()/2, k + y_off + 0.08,
            f'{k:.2f}', ha='center', fontsize=10, weight='bold')

ax.axhline(1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
ax.text(7.6, 1.06, 'KSR = 1 (no improvement)', fontsize=8, style='italic')

ax.axhspan(0, 1, alpha=0.07, color='red', zorder=0)
ax.axhspan(1, 6, alpha=0.07, color='green', zorder=0)
ax.text(-0.45, 0.5, 'AMPLIFIES\nKIE error', rotation=90, fontsize=9,
        weight='bold', color='darkred', va='center')
ax.text(-0.45, 2.5, 'REDUCES\nKIE error', rotation=90, fontsize=9,
        weight='bold', color='darkgreen', va='center')

ax.axvline(2.5, color='black', linestyle=':', alpha=0.5)
ax.text(1, 5.2, 'WLS coupling\n(δD as joint constraint)', ha='center', fontsize=10,
        weight='bold', color='darkred')
ax.text(5.5, 5.2, 'Agreement filter (N=5000)\n(δD as consistency check)\n— THIS WORK',
        ha='center', fontsize=10, weight='bold', color='darkgreen')

ax.set_xticks(range(len(methods)))
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('KIE Sensitivity Ratio (KSR)', fontsize=12)
ax.set_ylim(0, 5.8)
ax.set_title('KIE Sensitivity Reduction: All Methods (N = 5000, with bootstrap 95% CIs)',
             fontsize=13, weight='bold', pad=15)
ax.grid(alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(FIG_DIR / "figM6_KSR_summary_v3.png", dpi=300, bbox_inches='tight')
plt.close()
print("      Saved: figM6_KSR_summary_v3.png")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 4: Recolored forest plot (v2)
# ═══════════════════════════════════════════════════════════════════════════
print("[4/5] Generating figM7_forest_plot_v2 (recolored)...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300,
                         gridspec_kw={'width_ratios': [1.3, 1]})

# Panel (a): Full-record
ax = axes[0]
t90 = p9_high_n['thresholds']['90']
s_rate = t90['rate_S'] * 100
c_rate = t90['rate_C'] * 100
disc = t90['discriminant_pp']
disc_lo = t90['discriminant_CI_pp'][0]
disc_hi = t90['discriminant_CI_pp'][1]

y_positions = [2, 1]
labels_fp = ['Cantrell\n(α = 1.0054)', 'Saueressig\n(α = 1.0039)']
rates = [c_rate, s_rate]
ci_half = [0.3, 0.35]
fp_colors = [C_CANTRELL, C_SAUERESSIG]

ax.set_ylim(0, 3.5)
ax.set_xlim(0, 90)

for i, (y, lbl, rate, ch, clr) in enumerate(zip(y_positions, labels_fp, rates, ci_half, fp_colors)):
    ax.errorbar(rate, y, xerr=ch, fmt='o', color=clr, markersize=12,
                capsize=8, capthick=2, elinewidth=2, markeredgecolor='black', markeredgewidth=1)
    ax.text(rate + 1.5, y, f'{rate:.1f}%', va='center', fontsize=12, weight='bold', color=clr)
    ax.text(2, y, lbl, va='center', fontsize=11, ha='left')

ax.annotate('', xy=(c_rate, 2.3), xytext=(s_rate, 2.3),
            arrowprops=dict(arrowstyle='<->', color=C_PURPLE, lw=2))
ax.text((c_rate + s_rate) / 2, 2.55,
        f'Δ = {disc:.1f} pp\n[{disc_lo:.1f}, {disc_hi:.1f}]\np ≪ 0.001',
        ha='center', fontsize=11, weight='bold', color=C_PURPLE,
        bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.8))

ax.set_xlabel('Agreement Rate at T = 90 Tg/yr (%)', fontsize=12)
ax.set_yticks([])
ax.set_title('(a)  Full Record (1999–2022, N = 5000)', fontsize=12, weight='bold')
ax.grid(alpha=0.3, axis='x')
ax.axvline(50, color='gray', ls=':', alpha=0.3)

# Panel (b): Epoch stability
ax = axes[1]
epochs = p8['temporal_stability']
enames = ['epoch1_1999_2006', 'epoch2_2007_2014', 'epoch3_2015_2022']
elabels = ['1999–2006\n(plateau)', '2007–2014\n(growth)', '2015–2022\n(accel.)']

y_pos = [3, 2, 1]
for i, (ename, elbl, yp) in enumerate(zip(enames, elabels, y_pos)):
    ep = epochs[ename]
    s_r = ep['saueressig']['rate'] * 100
    c_r = ep['cantrell']['rate'] * 100
    s_ci_ep = [v * 100 for v in ep['saueressig']['CI']]
    c_ci_ep = [v * 100 for v in ep['cantrell']['CI']]
    disc_ep = ep['discriminant_pp']
    sig = ep['significant']

    ax.errorbar(s_r, yp - 0.12, xerr=[[s_r - s_ci_ep[0]], [s_ci_ep[1] - s_r]],
                fmt='s', color=C_SAUERESSIG, markersize=8, capsize=5, capthick=1.5, elinewidth=1.5)
    ax.errorbar(c_r, yp + 0.12, xerr=[[c_r - c_ci_ep[0]], [c_ci_ep[1] - c_r]],
                fmt='o', color=C_CANTRELL, markersize=8, capsize=5, capthick=1.5, elinewidth=1.5)

    ax.annotate('', xy=(c_r, yp + 0.35), xytext=(s_r, yp + 0.35),
                arrowprops=dict(arrowstyle='<->', color=C_PURPLE, lw=1.5))
    stars = '***' if sig else 'n.s.'
    ax.text((c_r + s_r) / 2, yp + 0.50, f'Δ={disc_ep:.1f}pp {stars}',
            ha='center', fontsize=9, weight='bold', color=C_PURPLE)
    ax.text(2, yp, elbl, va='center', fontsize=9, ha='left')

ax.set_xlim(0, 90)
ax.set_ylim(0.3, 3.8)
ax.set_yticks([])
ax.set_xlabel('Agreement Rate (%)', fontsize=12)
ax.set_title('(b)  By Epoch (T = 100 Tg/yr)', fontsize=12, weight='bold')
ax.grid(alpha=0.3, axis='x')

legend_elements = [
    Line2D([0], [0], marker='o', color=C_CANTRELL, markersize=8, linestyle='None',
           markeredgecolor='black', label='Cantrell (1.0054)'),
    Line2D([0], [0], marker='s', color=C_SAUERESSIG, markersize=8, linestyle='None',
           markeredgecolor='black', label='Saueressig (1.0039)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

fig.suptitle('Agreement-Rate Discriminant: Cantrell vs. Saueressig', fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig(FIG_DIR / "figM7_forest_plot_v2.png", dpi=300, bbox_inches='tight')
plt.close()
print("      Saved: figM7_forest_plot_v2.png")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 5: Recolored temporal stability (v3)
# ═══════════════════════════════════════════════════════════════════════════
print("[5/5] Generating fig14_temporal_stability_v3 (recolored)...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), dpi=300)

elabels_full = ['1999–2006\n(plateau)', '2007–2014\n(renewed growth)', '2015–2022\n(acceleration)']
rates_S_ep = [epochs[e]['saueressig']['rate'] * 100 for e in enames]
rates_C_ep = [epochs[e]['cantrell']['rate'] * 100 for e in enames]
ci_S = [epochs[e]['saueressig']['CI'] for e in enames]
ci_C = [epochs[e]['cantrell']['CI'] for e in enames]
err_S = [[r - c[0]*100 for r, c in zip(rates_S_ep, ci_S)],
         [c[1]*100 - r for r, c in zip(rates_S_ep, ci_S)]]
err_C = [[r - c[0]*100 for r, c in zip(rates_C_ep, ci_C)],
         [c[1]*100 - r for r, c in zip(rates_C_ep, ci_C)]]

# (a) Paired bars
ax = axes[0]
x = np.arange(len(enames))
w = 0.35
ax.bar(x - w/2, rates_S_ep, w, yerr=err_S, color=C_SAUERESSIG, alpha=0.8,
       label='Saueressig (1.0039)', capsize=4, edgecolor='black')
ax.bar(x + w/2, rates_C_ep, w, yerr=err_C, color=C_CANTRELL, alpha=0.8,
       label='Cantrell (1.0054)', capsize=4, edgecolor='black')

for i in range(len(enames)):
    s_top = rates_S_ep[i] + err_S[1][i]
    c_top = rates_C_ep[i] + err_C[1][i]
    bracket_y = max(s_top, c_top) + 3
    ax.plot([i - w/2, i - w/2, i + w/2, i + w/2],
            [bracket_y - 1, bracket_y, bracket_y, bracket_y - 1],
            color='black', lw=1.2)
    ax.text(i, bracket_y + 0.5, '***', ha='center', fontsize=12, weight='bold')

ax.set_xticks(x)
ax.set_xticklabels(elabels_full, fontsize=9)
ax.set_ylabel('Agreement Rate (%)')
ax.set_title('(a) Agreement Rate by Epoch', fontsize=11, weight='bold')
ax.legend(fontsize=9)
ax.grid(alpha=0.3, axis='y')
ax.set_ylim(0, 85)

# (b) Discriminant bars
ax = axes[1]
discs = [epochs[e]['discriminant_pp'] for e in enames]
sigs = [epochs[e]['significant'] for e in enames]
disc_colors = [C_KEEP if s else C_NEUTRAL for s in sigs]
ax.bar(x, discs, color=disc_colors, alpha=0.75, edgecolor='black')
ax.set_xticks(x)
ax.set_xticklabels(elabels_full, fontsize=9)
ax.set_ylabel('Discriminant Δ (pp)')
ax.set_title('(b) Cantrell − Saueressig Discriminant', fontsize=11, weight='bold')
ax.grid(alpha=0.3, axis='y')
ax.axhline(0, color='black', lw=0.5)

for xi, d, s in zip(x, discs, sigs):
    stars = '***' if s else 'n.s.'
    ax.text(xi, d + 0.8, f'+{d:.1f} pp\n{stars}', ha='center', fontsize=10,
            fontweight='bold', color=C_KEEP if s else C_NEUTRAL)

ax.set_ylim(0, 38)

fig.suptitle('Temporal Stability of the KIE Discriminant (T = 100 Tg/yr)',
             fontsize=12, weight='bold')
plt.tight_layout()
plt.savefig(FIG_DIR / "fig14_temporal_stability_v3.png", dpi=300, bbox_inches='tight')
plt.close()
print("      Saved: fig14_temporal_stability_v3.png")


print("\n✅ All figure improvements complete.")
print("   New files:")
print("     figures/figM8_before_after_timeseries.png  (NEW)")
print("     figures/figM1_schematic_v3.png             (recolored)")
print("     figures/figM6_KSR_summary_v3.png           (recolored)")
print("     figures/figM7_forest_plot_v2.png            (recolored)")
print("     figures/fig14_temporal_stability_v3.png     (recolored)")
