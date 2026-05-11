#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 14:36:20 2025

@author: ryoung
"""

locals().clear()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# load histograms
Mic_Delta_dD_baseline = pd.read_csv('Output/dD_histogram_UmezawaCal_noBUDS.csv').iloc[:,0]
FF_Delta_dD_baseline = pd.read_csv('Output/dD_histogram_UmezawaCal_noBUDS.csv').iloc[:,1]
Mic_Delta_dD_OH_inc = pd.read_csv('Output/dD_histogram_OH_inc.csv').iloc[:,0]
FF_Delta_dD_OH_inc = pd.read_csv('Output/dD_histogram_OH_inc.csv').iloc[:,1]
Mic_Delta_dD_RedCl = pd.read_csv('Output/dD_histogram_RedCl.csv').iloc[:,0]
FF_Delta_dD_RedCl = pd.read_csv('Output/dD_histogram_RedCl.csv').iloc[:,1]
Mic_Delta_dD_RedBB = pd.read_csv('Output/dD_histogram_RedBB.csv').iloc[:,0]
FF_Delta_dD_RedBB = pd.read_csv('Output/dD_histogram_RedBB.csv').iloc[:,1]

Mic_Delta_d13C_baseline = pd.read_csv('Output/d13C_histogram_Cantrell.csv').iloc[:,0]
FF_Delta_d13C_baseline = pd.read_csv('Output/d13C_histogram_Cantrell.csv').iloc[:,1]
Mic_Delta_d13C_OH_inc = pd.read_csv('Output/d13C_histogram_OH_inc.csv').iloc[:,0]
FF_Delta_d13C_OH_inc = pd.read_csv('Output/d13C_histogram_OH_inc.csv').iloc[:,1]
Mic_Delta_d13C_RedCl = pd.read_csv('Output/d13C_histogram_RedCl.csv').iloc[:,0]
FF_Delta_d13C_RedCl = pd.read_csv('Output/d13C_histogram_RedCl.csv').iloc[:,1]
Mic_Delta_d13C_RedBB = pd.read_csv('Output/d13C_histogram_RedBB.csv').iloc[:,0]
FF_Delta_d13C_RedBB = pd.read_csv('Output/d13C_histogram_RedBB.csv').iloc[:,1]

# calculate means and standard deviations
Mic_Delta_dD_baseline_mean = Mic_Delta_dD_baseline.mean()
Mic_Delta_dD_baseline_std = Mic_Delta_dD_baseline.std()
FF_Delta_dD_baseline_mean = FF_Delta_dD_baseline.mean()
FF_Delta_dD_baseline_std = FF_Delta_dD_baseline.std()

Mic_Delta_dD_OH_inc_mean = Mic_Delta_dD_OH_inc.mean()
Mic_Delta_dD_OH_inc_std = Mic_Delta_dD_OH_inc.std()
FF_Delta_dD_OH_inc_mean = FF_Delta_dD_OH_inc.mean()
FF_Delta_dD_OH_inc_std = FF_Delta_dD_OH_inc.std()

Mic_Delta_dD_RedCl_mean = Mic_Delta_dD_RedCl.mean()
Mic_Delta_dD_RedCl_std = Mic_Delta_dD_RedCl.std()
FF_Delta_dD_RedCl_mean = FF_Delta_dD_RedCl.mean()
FF_Delta_dD_RedCl_std = FF_Delta_dD_RedCl.std()

Mic_Delta_dD_RedBB_mean = Mic_Delta_dD_RedBB.mean()
Mic_Delta_dD_RedBB_std = Mic_Delta_dD_RedBB.std()
FF_Delta_dD_RedBB_mean = FF_Delta_dD_RedBB.mean()
FF_Delta_dD_RedBB_std = FF_Delta_dD_RedBB.std()

Mic_Delta_d13C_baseline_mean = Mic_Delta_d13C_baseline.mean()
Mic_Delta_d13C_baseline_std = Mic_Delta_d13C_baseline.std()
FF_Delta_d13C_baseline_mean = FF_Delta_d13C_baseline.mean()
FF_Delta_d13C_baseline_std = FF_Delta_d13C_baseline.std()

Mic_Delta_d13C_OH_inc_mean = Mic_Delta_d13C_OH_inc.mean()
Mic_Delta_d13C_OH_inc_std = Mic_Delta_d13C_OH_inc.std()
FF_Delta_d13C_OH_inc_mean = FF_Delta_d13C_OH_inc.mean()
FF_Delta_d13C_OH_inc_std = FF_Delta_d13C_OH_inc.std()

Mic_Delta_d13C_RedCl_mean = Mic_Delta_d13C_RedCl.mean()
Mic_Delta_d13C_RedCl_std = Mic_Delta_d13C_RedCl.std()
FF_Delta_d13C_RedCl_mean = FF_Delta_d13C_RedCl.mean()
FF_Delta_d13C_RedCl_std = FF_Delta_d13C_RedCl.std()

Mic_Delta_d13C_RedBB_mean = Mic_Delta_d13C_RedBB.mean()
Mic_Delta_d13C_RedBB_std = Mic_Delta_d13C_RedBB.std()
FF_Delta_d13C_RedBB_mean = FF_Delta_d13C_RedBB.mean()
FF_Delta_d13C_RedBB_std = FF_Delta_d13C_RedBB.std()

# Load Carbon tracker methane
data2 = pd.read_excel('data/CarbonTracker_CH4.xlsx')
micCT = data2.iloc[:,7].values
ffCT = data2.iloc[:,3].values
bbCT = data2.iloc[:,11].values

# Calculate carbontracker difference
# 2020 to 2021 average minus 2005 to 2007 average 
Mic_Delta_CT_total = np.mean(micCT[-2:]) - np.mean(micCT[8:11])
FF_Delta_CT_total = np.mean(ffCT[-2:]) - np.mean(ffCT[8:11])

#%% Perform statistical tests and plot supp fig 12 for reference

# Simplified Bayesian Posterior Probaility Test
# Baseline
Mic_Baseline = Mic_Delta_dD_baseline - Mic_Delta_d13C_baseline
Mic_Baseline_Result = np.sum(Mic_Baseline > 0)/len(Mic_Delta_dD_baseline)
FF_Baseline = FF_Delta_dD_baseline - FF_Delta_d13C_baseline
FF_Baseline_Result = np.sum(FF_Baseline > 0)/len(FF_Delta_dD_baseline)
# OH increase (Morgernstern)
Mic_OH_inc = Mic_Delta_dD_OH_inc - Mic_Delta_d13C_OH_inc
Mic_OH_inc_Result = np.sum(Mic_OH_inc > 0)/len(Mic_Delta_dD_OH_inc)
FF_OH_inc = FF_Delta_dD_OH_inc - FF_Delta_d13C_OH_inc
FF_OH_inc_Result = np.sum(FF_OH_inc > 0)/len(FF_Delta_dD_OH_inc)
# Reduced Cl (narrative changing experiment)
Mic_RedCl = Mic_Delta_dD_RedCl - Mic_Delta_d13C_RedCl
Mic_RedCl_Result = np.sum(Mic_RedCl > 0)/len(Mic_Delta_dD_RedCl)
FF_RedCl = FF_Delta_dD_RedCl - FF_Delta_d13C_RedCl
FF_RedCl_Result = np.sum(FF_RedCl > 0)/len(FF_Delta_dD_RedCl)
# Reduced Cl (narrative changing experiment)
Mic_RedBB = Mic_Delta_dD_RedBB - Mic_Delta_d13C_RedBB
Mic_RedBB_Result = np.sum(Mic_RedBB > 0)/len(Mic_Delta_dD_RedBB)
FF_RedBB = FF_Delta_dD_RedBB - FF_Delta_d13C_RedBB
FF_RedBB_Result = np.sum(FF_RedBB > 0)/len(FF_Delta_dD_RedBB)

# Print results: 
results = {
    "Mic dD baseline > Mic d13C baseline": Mic_Baseline_Result,
    "FF dD baseline > FF d13C baseline": FF_Baseline_Result,
    "Mic dD OH_inc > Mic d13C OH_inc": Mic_OH_inc_Result,
    "FF dD OH_inc > FF d13C OH_inc": FF_OH_inc_Result,
    "Mic dD RedCl > Mic d13C RedCl": Mic_RedCl_Result,
    "FF dD RedCl > FF d13C RedCl": FF_RedCl_Result,
    "Mic dD RedBB > Mic d13C RedBB": Mic_RedBB_Result,
    "FF dD RedBB > FF d13C RedBB": FF_RedBB_Result}

# results in subplot order
posterior_probs = [Mic_Baseline_Result, FF_Baseline_Result, Mic_OH_inc_Result, FF_OH_inc_Result, Mic_RedCl_Result, FF_RedCl_Result, Mic_RedBB_Result, FF_RedBB_Result]

for description, value in results.items():
    print(f"Probability that {description}: {value:.3f}")

# Determine bin edges with width of 3
bin_width = 3
binsMic = np.arange(-10, 150 + bin_width, bin_width)
binsFF = np.arange(-60, 60 + bin_width, bin_width)


#%% plot it

import matplotlib.pyplot as plt
import numpy as np

# set font size globally
plt.rcParams.update({'font.size': 10})

# helper to convert cm to inch
def cm2inch(x): return x/2.54

# create figure with 4x2 subplots, large PNAS format
fig, axes = plt.subplots(4, 2, figsize=(cm2inch(18), cm2inch(22)), dpi=200)
axes = axes.flatten()

# --- probability text boxes ---
axes[1].text(0.95, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {Mic_Baseline_Result:.3f}",
             ha='right', va='top', transform=axes[1].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))
axes[0].text(0.05, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {FF_Baseline_Result:.3f}",
             ha='left', va='top', transform=axes[0].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))
axes[3].text(0.95, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {Mic_RedCl_Result:.3f}",
             ha='right', va='top', transform=axes[3].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))
axes[2].text(0.05, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {FF_RedCl_Result:.3f}",
             ha='left', va='top', transform=axes[2].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))
axes[5].text(0.95, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {Mic_RedBB_Result:.3f}",
             ha='right', va='top', transform=axes[5].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))
axes[4].text(0.05, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {FF_RedBB_Result:.3f}",
             ha='left', va='top', transform=axes[4].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))
axes[7].text(0.95, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {Mic_OH_inc_Result:.3f}",
             ha='right', va='top', transform=axes[7].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))
axes[6].text(0.05, 0.8, rf"$P(\delta D > \delta^{{13}}C)$ = {FF_OH_inc_Result:.3f}",
             ha='left', va='top', transform=axes[6].transAxes,
             fontsize=6, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))

# --- microbial baseline ---
closest_bin_mic = min(binsMic, key=lambda x: abs(x - Mic_Delta_CT_total))
axes[1].axvspan(closest_bin_mic - bin_width/2, closest_bin_mic + bin_width/2,
                color='black', alpha=0.3)
axes[1].hist(Mic_Delta_d13C_baseline, bins=binsMic, color='maroon', edgecolor='black', alpha=0.7)
axes[1].hist(Mic_Delta_dD_baseline, bins=binsMic, color='blue', edgecolor='black', alpha=0.7)
axes[1].grid(axis='y', linestyle='--', alpha=0.6)
axes[1].set_title("Mic Emission Change: Baseline", fontsize=10)
axes[1].set_xticklabels([])
axes[1].axvline(0, color='black', linestyle='--', linewidth=1)
axes[1].set_xlim(-20, 140)

# microbial scenarios
axes[3].hist(Mic_Delta_d13C_RedCl, bins=binsMic, color='maroon', edgecolor='black', alpha=0.7)
axes[3].hist(Mic_Delta_dD_RedCl, bins=binsMic, color='blue', edgecolor='black', alpha=0.7)
axes[3].set_title("Mic Emission Change: Cl Decrease", fontsize=10)
axes[3].grid(axis='y', linestyle='--', alpha=0.6)
axes[3].set_xticklabels([])
axes[3].axvline(0, color='black', linestyle='--', linewidth=1)
axes[3].set_xlim(-20, 140)

axes[5].hist(Mic_Delta_d13C_RedBB, bins=binsMic, color='maroon', edgecolor='black', alpha=0.7)
axes[5].hist(Mic_Delta_dD_RedBB, bins=binsMic, color='blue', edgecolor='black', alpha=0.7)
axes[5].set_title("Mic Emission Change: Pyrogenic Decrease", fontsize=10)
axes[5].grid(axis='y', linestyle='--', alpha=0.6)
axes[5].set_xticklabels([])
axes[5].axvline(0, color='black', linestyle='--', linewidth=1)
axes[5].set_xlim(-20, 140)

axes[7].hist(Mic_Delta_d13C_OH_inc, bins=binsMic, color='maroon', edgecolor='black', alpha=0.7)
axes[7].hist(Mic_Delta_dD_OH_inc, bins=binsMic, color='blue', edgecolor='black', alpha=0.7)
axes[7].set_title("Mic Emission Change: OH Increase", fontsize=10)
axes[7].grid(axis='y', linestyle='--', alpha=0.6)
axes[7].axvline(0, color='black', linestyle='--', linewidth=1)
axes[7].set_xlabel("'05-'07 to '20-'22 Mic Emission Change (Tg/yr)", fontsize=10)
axes[7].set_xlim(-20, 140)

# --- ff baseline ---
closest_bin_FF = min(binsFF, key=lambda x: abs(x - FF_Delta_CT_total))
axes[0].axvspan(closest_bin_FF - bin_width/2, closest_bin_FF + bin_width/2,
                color='black', alpha=0.3, label='CT-CH4')
axes[0].hist(FF_Delta_d13C_baseline, bins=binsFF, color='maroon', edgecolor='black', alpha=0.7, label="${\delta}^{13}C$")
axes[0].hist(FF_Delta_dD_baseline, bins=binsFF, color='blue', edgecolor='black', alpha=0.7, label="${\delta}D$")
axes[0].set_title("FF Emission Change: Baseline", fontsize=10)
axes[0].grid(axis='y', linestyle='--', alpha=0.6)
axes[0].legend(loc='lower right', fontsize=9)
axes[0].set_xticklabels([])
axes[0].axvline(0, color='black', linestyle='--', linewidth=1)

# ff scenarios
axes[2].hist(FF_Delta_d13C_RedCl, bins=binsFF, color='maroon', edgecolor='black', alpha=0.7)
axes[2].hist(FF_Delta_dD_RedCl, bins=binsFF, color='blue', edgecolor='black', alpha=0.7)
axes[2].set_title("FF Emission Change: Cl Decrease", fontsize=10)
axes[2].grid(axis='y', linestyle='--', alpha=0.6)
axes[2].set_xticklabels([])
axes[2].axvline(0, color='black', linestyle='--', linewidth=1)

axes[4].hist(FF_Delta_d13C_RedBB, bins=binsFF, color='maroon', edgecolor='black', alpha=0.7)
axes[4].hist(FF_Delta_dD_RedBB, bins=binsFF, color='blue', edgecolor='black', alpha=0.7)
axes[4].set_title("FF Emission Change: Pyrogenic Decrease", fontsize=10)
axes[4].grid(axis='y', linestyle='--', alpha=0.6)
axes[4].set_xticklabels([])
axes[4].axvline(0, color='black', linestyle='--', linewidth=1)

axes[6].hist(FF_Delta_d13C_OH_inc, bins=binsFF, color='maroon', edgecolor='black', alpha=0.7)
axes[6].hist(FF_Delta_dD_OH_inc, bins=binsFF, color='blue', edgecolor='black', alpha=0.7)
axes[6].set_title("FF Emission Change: OH Increase", fontsize=10)
axes[6].grid(axis='y', linestyle='--', alpha=0.6)
axes[6].axvline(0, color='black', linestyle='--', linewidth=1)
axes[6].set_xlabel("'05-'07 to '20-'22 FF Emission Change (Tg/yr)", fontsize=10)

# --- panel labels ---
labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
for ax, label in zip(axes, labels):
    ax.text(0.02, 0.95, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')

plt.tight_layout()

# save in PNAS-compliant format
plt.savefig("Output/Fig3.pdf", dpi=200, bbox_inches="tight")
plt.show()
plt.close()
