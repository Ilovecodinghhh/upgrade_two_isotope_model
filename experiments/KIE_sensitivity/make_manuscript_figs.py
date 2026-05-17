#!/usr/bin/env python3
"""Generate the two missing manuscript figures: schematic + KSR summary."""

import sys
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

BASE = Path(__file__).resolve().parent
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(exist_ok=True)

# =========================================================
# Figure M1 — Conceptual schematic of the Agreement Filter
# =========================================================
fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
ax.set_xlim(0, 12)
ax.set_ylim(0, 7)
ax.axis('off')

def box(x, y, w, h, text, color, fontsize=10, weight='normal'):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                       linewidth=1.5, edgecolor='black', facecolor=color, alpha=0.85)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fontsize, weight=weight, wrap=True)

def arrow(x1, y1, x2, y2, label=None, color='black'):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='->,head_width=6,head_length=8',
                         linewidth=1.8, color=color)
    ax.add_patch(a)
    if label:
        ax.text((x1+x2)/2, (y1+y2)/2 + 0.15, label, ha='center', fontsize=8, style='italic')

# Inputs
box(0.3, 5.5, 2.4, 1.0, "Atmospheric obs.\nCH$_4$, $\\delta^{13}$C, $\\delta$D", '#cce5ff', fontsize=9)
box(0.3, 4.0, 2.4, 1.0, "Sink fractions\n& KIE priors", '#cce5ff', fontsize=9)
box(0.3, 2.5, 2.4, 1.0, "Source signatures\n(Sherwood+, Menoud+)", '#cce5ff', fontsize=9)
box(0.3, 1.0, 2.4, 1.0, "BB emissions\n(GFEDv4.1s)", '#cce5ff', fontsize=9)

# Monte Carlo gate
box(3.2, 3.2, 1.8, 1.4, "Monte Carlo\n(N=1000)", '#ffe5b4', fontsize=10, weight='bold')

# Two parallel inversions
box(5.5, 5.0, 2.6, 1.1, "$\\delta^{13}$C mass balance\n(well-constrained,\n$\\sigma$ ≈ 4 Tg/yr)", '#d4edda', fontsize=9)
box(5.5, 1.7, 2.6, 1.1, "$\\delta$D mass balance\n(loosely constrained,\n$\\sigma$ ≈ 30+ Tg/yr)", '#d4edda', fontsize=9)

# Outputs from each
box(8.7, 5.0, 1.6, 1.1, "$S_{FF}^{13C}(j,k)$\n$S_{Mic}^{13C}(j,k)$", '#f8d7da', fontsize=9)
box(8.7, 1.7, 1.6, 1.1, "$S_{FF}^{D}(j,k)$", '#f8d7da', fontsize=9)

# Agreement test box
box(8.7, 3.2, 1.6, 1.4, "Agree?\n$|S_{FF}^{13C}-S_{FF}^{D}|<T$", '#fff3cd', fontsize=9, weight='bold')

# Final output
box(10.5, 3.2, 1.4, 1.4, "FILTERED\nensemble\n(KSR = 2.5–3.2)", '#d4edda', fontsize=9, weight='bold')

# Arrows from inputs to MC
for y in [6.0, 4.5, 3.0, 1.5]:
    arrow(2.7, y, 3.2, 3.9 + (y-3.75)*0.1, color='gray')

# MC → two inversions
arrow(5.0, 4.2, 5.5, 5.5)
arrow(5.0, 3.6, 5.5, 2.3)

# Inversions → outputs
arrow(8.1, 5.55, 8.7, 5.55)
arrow(8.1, 2.25, 8.7, 2.25)

# Outputs → agreement test
arrow(9.5, 5.0, 9.5, 4.6)
arrow(9.5, 2.8, 9.5, 3.2)

# Agreement → filtered ensemble
arrow(10.3, 3.9, 10.5, 3.9, label='keep if agree')

# Title
ax.text(6, 6.7, 'The Dual-Isotope Agreement Filter',
        ha='center', fontsize=14, weight='bold')
ax.text(6, 6.35, 'Two independent isotopic budgets → consistency check → reduced KIE sensitivity',
        ha='center', fontsize=10, style='italic', color='dimgray')

plt.savefig(FIG_DIR / "figM1_schematic.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {FIG_DIR / 'figM1_schematic.png'}")

# =========================================================
# Figure M6 — Summary KSR across all methods
# =========================================================
methods = [
    ('1-box\nWLS\n($w_D=1$)', 0.20, '#dc3545'),
    ('1-box\nWLS\n($w_D=0.01$)', 0.24, '#dc3545'),
    ('2-box\nWLS\n(fixed)', 0.22, '#dc3545'),
    ('Agreement\nfilter\n(T=200)', 1.09, '#ffc107'),
    ('Agreement\nfilter\n(T=150)', 1.51, '#28a745'),
    ('Agreement\nfilter\n(T=100)', 2.48, '#28a745'),
    ('Agreement\nfilter\n(T=75)', 2.16, '#28a745'),
    ('Agreement\nfilter\n(T=50)', 3.21, '#28a745'),
]
labels = [m[0] for m in methods]
ksrs = [m[1] for m in methods]
colors = [m[2] for m in methods]

fig, ax = plt.subplots(figsize=(11, 6), dpi=300)
bars = ax.bar(range(len(methods)), ksrs, color=colors, edgecolor='black', linewidth=1.2)

# Annotate values
for i, (bar, k) in enumerate(zip(bars, ksrs)):
    ax.text(bar.get_x() + bar.get_width()/2, k + 0.08,
            f'{k:.2f}', ha='center', fontsize=11, weight='bold')

# Reference line
ax.axhline(1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
ax.text(7.5, 1.08, 'KSR = 1\n(no improvement)', fontsize=8, style='italic', color='black')

# Shaded regions
ax.axhspan(0, 1, alpha=0.08, color='red', zorder=0)
ax.axhspan(1, 4, alpha=0.08, color='green', zorder=0)
ax.text(-0.4, 0.5, 'AMPLIFIES\nKIE error', rotation=90, fontsize=9,
        weight='bold', color='darkred', va='center')
ax.text(-0.4, 2.5, 'REDUCES\nKIE error', rotation=90, fontsize=9,
        weight='bold', color='darkgreen', va='center')

# Group separator
ax.axvline(2.5, color='black', linestyle=':', alpha=0.5)
ax.text(1, 3.6, 'WLS coupling\n(δD as joint constraint)', ha='center', fontsize=10,
        weight='bold', color='darkred')
ax.text(5.5, 3.6, 'Agreement filter\n(δD as consistency check) — THIS WORK', ha='center',
        fontsize=10, weight='bold', color='darkgreen')

ax.set_xticks(range(len(methods)))
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('KIE Sensitivity Ratio (KSR)', fontsize=12)
ax.set_ylim(0, 4)
ax.set_title('KIE Sensitivity Reduction Across All Tested Methods',
             fontsize=13, weight='bold', pad=15)
ax.grid(alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(FIG_DIR / "figM6_KSR_summary.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {FIG_DIR / 'figM6_KSR_summary.png'}")

print("\nAll manuscript figures generated.")
