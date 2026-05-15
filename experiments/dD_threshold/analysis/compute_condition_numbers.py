#!/usr/bin/env python3
"""
compute_condition_numbers.py — Compute condition numbers for 1-box and 2-box
source-signature matrices, and run W-matrix sensitivity analysis.

Outputs:
  - Condition numbers (1-box global vs 2-box per-hemisphere)
  - W-matrix sensitivity table (real MC simulations)
"""
import sys
from pathlib import Path
import json
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    load_data, sample_source_signatures_hemi,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_KIE, compute_bulk_KIE,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
)

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ANALYSIS_DIR / "analysis"))
from core import run_twobox, ci_width

data = load_data(REPO_ROOT, two_box=True)
rng = np.random.default_rng(42)

print("=" * 60)
print("PART 1: CONDITION NUMBERS")
print("=" * 60)

# Compute mean source signatures over 100 MC samples, year 10
n_samples = 100
j = 10  # representative year

# Collect source sig arrays
bb_d13c_g, ff_d13c_g, mic_d13c_g = [], [], []
bb_dD_g, ff_dD_g, mic_dD_g = [], [], []
bb_d13c_nh, ff_d13c_nh, mic_d13c_nh = [], [], []
bb_dD_nh, ff_dD_nh, mic_dD_nh = [], [], []
bb_d13c_sh, ff_d13c_sh, mic_d13c_sh = [], [], []
bb_dD_sh, ff_dD_sh, mic_dD_sh = [], [], []

for k in range(n_samples):
    sigs = sample_source_signatures_hemi(rng, data, k, data.n_years)
    # Global (mean of NH/SH)
    bb_d13c_g.append(sigs['bb_d13C'][j])
    ff_d13c_g.append(sigs['ff_d13C'][j])
    mic_d13c_g.append(sigs['mic_d13C'][j])
    bb_dD_g.append(sigs['bb_dD'][j])
    ff_dD_g.append(sigs['ff_dD'][j])
    mic_dD_g.append(sigs['mic_dD'][j])
    # NH
    bb_d13c_nh.append(sigs['bb_d13C_NH'][j])
    ff_d13c_nh.append(sigs['ff_d13C_NH'][j])
    mic_d13c_nh.append(sigs['mic_d13C_NH'][j])
    bb_dD_nh.append(sigs['bb_dD_NH'][j])
    ff_dD_nh.append(sigs['ff_dD_NH'][j])
    mic_dD_nh.append(sigs['mic_dD_NH'][j])
    # SH
    bb_d13c_sh.append(sigs['bb_d13C_SH'][j])
    ff_d13c_sh.append(sigs['ff_d13C_SH'][j])
    mic_d13c_sh.append(sigs['mic_d13C_SH'][j])
    bb_dD_sh.append(sigs['bb_dD_SH'][j])
    ff_dD_sh.append(sigs['ff_dD_SH'][j])
    mic_dD_sh.append(sigs['mic_dD_SH'][j])

def mean(arr): return np.mean(arr)

# 1-box dual source-signature matrix [total; d13C; dD] × [BB, FF, Mic]
# In isotope fraction space
f13_bb_g = delta_to_fraction_d13C(mean(bb_d13c_g))
f13_ff_g = delta_to_fraction_d13C(mean(ff_d13c_g))
f13_mic_g = delta_to_fraction_d13C(mean(mic_d13c_g))
fD_bb_g = delta_to_fraction_dD(mean(bb_dD_g))
fD_ff_g = delta_to_fraction_dD(mean(ff_dD_g))
fD_mic_g = delta_to_fraction_dD(mean(mic_dD_g))

A_1box = np.array([
    [1.0, 1.0, 1.0],
    [f13_bb_g, f13_ff_g, f13_mic_g],
    [fD_bb_g, fD_ff_g, fD_mic_g],
])
kappa_1box = np.linalg.cond(A_1box)
print(f"\n1-box dual source-signature matrix:")
print(f"  A = {A_1box}")
print(f"  κ(A) = {kappa_1box:.1f}")

# 2-box: per-hemisphere
f13_bb_nh = delta_to_fraction_d13C(mean(bb_d13c_nh))
f13_ff_nh = delta_to_fraction_d13C(mean(ff_d13c_nh))
f13_mic_nh = delta_to_fraction_d13C(mean(mic_d13c_nh))
fD_bb_nh = delta_to_fraction_dD(mean(bb_dD_nh))
fD_ff_nh = delta_to_fraction_dD(mean(ff_dD_nh))
fD_mic_nh = delta_to_fraction_dD(mean(mic_dD_nh))

A_NH = np.array([
    [1.0, 1.0, 1.0],
    [f13_bb_nh, f13_ff_nh, f13_mic_nh],
    [fD_bb_nh, fD_ff_nh, fD_mic_nh],
])
kappa_NH = np.linalg.cond(A_NH)

f13_bb_sh = delta_to_fraction_d13C(mean(bb_d13c_sh))
f13_ff_sh = delta_to_fraction_d13C(mean(ff_d13c_sh))
f13_mic_sh = delta_to_fraction_d13C(mean(mic_d13c_sh))
fD_bb_sh = delta_to_fraction_dD(mean(bb_dD_sh))
fD_ff_sh = delta_to_fraction_dD(mean(ff_dD_sh))
fD_mic_sh = delta_to_fraction_dD(mean(mic_dD_sh))

A_SH = np.array([
    [1.0, 1.0, 1.0],
    [f13_bb_sh, f13_ff_sh, f13_mic_sh],
    [fD_bb_sh, fD_ff_sh, fD_mic_sh],
])
kappa_SH = np.linalg.cond(A_SH)

kappa_2box_avg = (kappa_NH + kappa_SH) / 2
print(f"\n2-box NH source-signature matrix:")
print(f"  κ(A_NH) = {kappa_NH:.1f}")
print(f"\n2-box SH source-signature matrix:")
print(f"  κ(A_SH) = {kappa_SH:.1f}")
print(f"\n2-box average κ = {kappa_2box_avg:.1f}")
print(f"Reduction: {(1 - kappa_2box_avg/kappa_1box)*100:.0f}%")

# Also compute d13C-only condition numbers (2×2 submatrices)
A_1box_c13 = A_1box[:2, 1:]  # drop dD row, drop BB col (BB prescribed)
A_NH_c13 = A_NH[:2, 1:]
A_SH_c13 = A_SH[:2, 1:]
print(f"\nδ¹³C-only condition numbers:")
print(f"  κ(1-box) = {np.linalg.cond(A_1box_c13):.1f}")
print(f"  κ(NH)    = {np.linalg.cond(A_NH_c13):.1f}")
print(f"  κ(SH)    = {np.linalg.cond(A_SH_c13):.1f}")

print("\n" + "=" * 60)
print("PART 2: W-MATRIX SENSITIVITY (N=200, real MC)")
print("=" * 60)

N_ITER = 200  # enough for stable CI estimates
SEED = 42

# Baseline d13C-only reference
ref = run_twobox(data, 1.0, N_ITER, SEED, mode="d13C_only")
ref_ci = ci_width(ref['FF_G'])
print(f"\nδ¹³C-only reference CI: {ref_ci:.1f} Tg/yr")

# Test W_mass_balance values: 50, 100, 200, 500
# Need to temporarily patch core.W_NH and core.W_SH
import experiments.dD_threshold.analysis.core as core

results = {}
for w_mass in [50, 100, 200, 500]:
    # Set both NH and SH mass balance weight
    core.W_NH = np.diag([float(w_mass), 1.0, 0.5])
    core.W_SH = np.diag([float(w_mass * 2), 1.0, 0.5])  # SH is 2x NH
    
    dual = run_twobox(data, 1.0, N_ITER, SEED, mode="dual")
    dual_ci = ci_width(dual['FF_G'])
    improvement = (1 - dual_ci / ref_ci) * 100
    
    # Also check threshold roughly: run at 4× and 5×
    d4 = run_twobox(data, 4.0, N_ITER, SEED, mode="dual")
    d5 = run_twobox(data, 5.0, N_ITER, SEED, mode="dual")
    ci4 = ci_width(d4['FF_G'])
    ci5 = ci_width(d5['FF_G'])
    # Linear interpolation for threshold
    if ci4 < ref_ci and ci5 > ref_ci:
        thresh = 4.0 + (ref_ci - ci4) / (ci5 - ci4)
    elif ci5 < ref_ci:
        thresh = 5.5  # approximate
    else:
        thresh = 3.5  # approximate
    
    results[w_mass] = {
        'dual_ci': round(dual_ci, 1),
        'improvement': round(improvement, 1),
        'threshold_approx': round(thresh, 1),
    }
    print(f"  W_mass={w_mass}: CI={dual_ci:.1f} Tg/yr, improvement={improvement:.1f}%, threshold≈{thresh:.1f}×")

# Restore defaults
core.W_NH = np.diag([100.0, 1.0, 0.5])
core.W_SH = np.diag([200.0, 1.0, 0.5])

# Test δ¹³C:δD weight ratio
print("\nδ¹³C:δD weight ratio sensitivity:")
for dD_w in [0.25, 0.5, 1.0]:
    core.W_NH = np.diag([100.0, 1.0, dD_w])
    core.W_SH = np.diag([200.0, 1.0, dD_w])
    dual = run_twobox(data, 1.0, N_ITER, SEED, mode="dual")
    dual_ci = ci_width(dual['FF_G'])
    improvement = (1 - dual_ci / ref_ci) * 100
    print(f"  δ¹³C:δD = 1:{dD_w}: CI={dual_ci:.1f}, improvement={improvement:.1f}%")

# Restore defaults
core.W_NH = np.diag([100.0, 1.0, 0.5])
core.W_SH = np.diag([200.0, 1.0, 0.5])

# Save results
output = {
    'condition_numbers': {
        'onebox_dual': round(kappa_1box, 1),
        'twobox_NH': round(kappa_NH, 1),
        'twobox_SH': round(kappa_SH, 1),
        'twobox_avg': round(kappa_2box_avg, 1),
        'reduction_pct': round((1 - kappa_2box_avg/kappa_1box)*100, 0),
    },
    'w_matrix_sensitivity': results,
    'ref_ci': round(ref_ci, 1),
}
out_path = Path(__file__).resolve().parent.parent / "results" / "condition_and_sensitivity.json"
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {out_path}")
