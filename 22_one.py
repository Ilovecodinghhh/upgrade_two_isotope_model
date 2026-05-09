#!/usr/bin/env python3
"""
22_one.py — 2×2 BB-Fixed, One-Box (Global) Model
===================================================
Infra: 2×2 system | 1 global box

Approach:
  Fix BB from CarbonTracker (GFED4 prior).
  Solve δ¹³C and δD SEPARATELY as two independent 2-equation systems:
    S_total = FF + Mic + BB_fixed
    S_total × δ_source = FF × δ_FF + Mic × δ_Mic + BB × δ_BB
  → FF = (S × δ_source − δ_Mic × (S − BB) − δ_BB × BB) / (δ_FF − δ_Mic)
  → Mic = S − BB − FF

  This yields TWO independent estimates (from δ¹³C and from δD) that can
  be cross-validated. Well-conditioned because each is only 1 equation.

References:
  - Riddell-Young et al. (2025, PNAS) — Ben's approach (1-box, separate isotopes)
  - Bao et al. (this work) — KIE sampling, time-varying τ

Lineage: v3.2 (BB-fixed 2×2) stripped to 1-box geometry
"""

import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add parent to path so `models` package is importable
sys.path.insert(0, str(Path(__file__).resolve().parent))
from models import (
    DEFAULT_CONFIG, SENSITIVITY_PRESETS,
    build_kie_sampler, compute_bulk_KIE, compute_lifetime,
    find_data_dirs, load_CH4, load_d13C_hemispheric,
    load_d13C_iterations, load_dD_iterations, load_source_signatures,
    load_BB_emissions, QualityMonitor, smooth_5yr, pad_to_length,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    PT, C13Std, DStd,
)
from models.inputs import SINK_FRACTION_OPTIONS

# ═══════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════
def get_config(preset_name=None):
    if preset_name and preset_name in SENSITIVITY_PRESETS:
        return SENSITIVITY_PRESETS[preset_name]
    return DEFAULT_CONFIG.copy()

config = get_config(sys.argv[1] if len(sys.argv) > 1 else None)
N = config['n_iterations']
SEED = config['seed']

dirs = find_data_dirs()
OUT_DIR = dirs['base'] / "Output_22_one"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("22_one: 2×2 BB-Fixed, One-Box (Global)")
print("=" * 70)

CH4_years, CH4 = load_CH4(dirs['data'])
_, c13_glob, c13_nh, c13_sh = load_d13C_hemispheric(dirs['data'])
d13C_iter = load_d13C_iterations(dirs['data'])
dD_mean, dD_iter = load_dD_iterations(dirs['data'])
sigs = load_source_signatures(dirs['src'], config)

n_years = len(CH4) - 1
model_years = np.arange(1999, 1999 + n_years)
tau = compute_lifetime(model_years, config)

BB_annual = load_BB_emissions(dirs['data'], config, n_years)

# ── Prepare source signature arrays ──────────────────────────────────────
# d13C
ff_d13C_mean = np.array(sigs['FF_d13C_GlobUnc'].iloc[28:, 1]).flatten()
ff_d13C_U = np.array(sigs['FF_d13C_GlobUnc'].iloc[28:, 2]).flatten()
bb_d13C_arr = np.array(sigs['BB_d13C'].iloc[1:, 1]).flatten()
bb_d13C_U_arr = np.array(sigs['BB_d13C'].iloc[1:, 2]).flatten()
bb_d13C_arr = np.concatenate([bb_d13C_arr, [bb_d13C_arr[-1]]])
bb_d13C_U_arr = np.concatenate([bb_d13C_U_arr, [bb_d13C_U_arr[-1]]])
mic_d13C_mean_val = sigs['Mic_d13C_ann'].iloc[:, 1].mean()

# dD
ff_dD_mean = np.array(sigs['FF_dD_GlobUnc'].iloc[34:, 1]).flatten()
ff_dD_U = np.array(sigs['FF_dD_GlobUnc'].iloc[34:, 2]).flatten()
bb_dD_arr = np.array(sigs['BB_dD'].iloc[:, 1]).flatten()
bb_dD_U_arr = np.array(sigs['BB_dD'].iloc[:, 2]).flatten()
bb_dD_arr = np.concatenate([np.full(3, bb_dD_arr[-1]), bb_dD_arr, [bb_dD_arr[-1]]])
bb_dD_U_arr = np.concatenate([np.full(3, bb_dD_U_arr[-1]), bb_dD_arr[:0], bb_dD_U_arr, [bb_dD_U_arr[-1]]])
mic_dD_mean_val = sigs['Mic_dD_ann'].iloc[:, 1].mean()
mic_dD_U_val = config['mic_dD_uncertainty']

# FF MC (EDGAR or CTCH4)
FF_d13C_MC_raw = sigs['FF_d13C_MC']
FF_dD_MC_raw = sigs['FF_dD_MC']
# Trim/pad to align
if 'EDGAR' in str(config.get('FF_signature', 'EDGAR')):
    FF_d13C_MC = FF_d13C_MC_raw.iloc[28:, 1:] if FF_d13C_MC_raw.shape[0] > 28 else FF_d13C_MC_raw.iloc[:, 1:]
    FF_dD_MC = FF_dD_MC_raw.iloc[34:, 1:] if FF_dD_MC_raw.shape[0] > 34 else FF_dD_MC_raw.iloc[:, 1:]
else:
    FF_d13C_MC = FF_d13C_MC_raw.iloc[1:, 1:] if FF_d13C_MC_raw.shape[0] > 1 else FF_d13C_MC_raw.iloc[:, 1:]
    FF_dD_MC = FF_dD_MC_raw.iloc[7:, 1:] if FF_dD_MC_raw.shape[0] > 7 else FF_dD_MC_raw.iloc[:, 1:]

# Pad dD MC iterations
if FF_dD_MC.shape[0] < n_years:
    pad_n = n_years - FF_dD_MC.shape[0]
    FF_dD_MC = pd.concat([pd.concat([FF_dD_MC.iloc[0:1]] * pad_n)] + [FF_dD_MC], ignore_index=True)

Mic_d13C_MC = sigs['Mic_d13C_MC']
Mic_dD_MC = sigs['Mic_dD_MC']
if Mic_dD_MC.shape[0] < n_years:
    pad_n = n_years - Mic_dD_MC.shape[0]
    Mic_dD_MC = pd.concat([pd.concat([Mic_dD_MC.iloc[0:1]] * pad_n)] + [Mic_dD_MC], ignore_index=True)
elif Mic_dD_MC.shape[0] > n_years:
    Mic_dD_MC = Mic_dD_MC.iloc[:n_years]

# Pad dD atmospheric iterations
if dD_iter.shape[0] < n_years + 1:
    pad_n = n_years + 1 - dD_iter.shape[0]
    dD_iter = np.vstack([np.repeat(dD_iter[0:1], pad_n, axis=0), dD_iter])

# ═══════════════════════════════════════════════════════════════════════════
# SINK FRACTIONS (global 1-box)
# ═══════════════════════════════════════════════════════════════════════════
sf_key = config.get('sink_fractions_global', 'global_default')
SINK_FRACS = SINK_FRACTION_OPTIONS[sf_key]

# ═══════════════════════════════════════════════════════════════════════════
# MONTE CARLO
# ═══════════════════════════════════════════════════════════════════════════
print(f"\nRunning MC: {N} iterations, BB fixed = {np.mean(BB_annual):.1f} Tg/yr")
rng = np.random.default_rng(seed=SEED)
kie_sampler = build_kie_sampler(config)

FF_d13C_out = np.zeros((n_years, N))
Mic_d13C_out = np.zeros((n_years, N))
FF_dD_out = np.zeros((n_years, N))
Mic_dD_out = np.zeros((n_years, N))
qm_d13C = QualityMonitor(n_years, N, "d13C")
qm_dD = QualityMonitor(n_years, N, "dD")

for k in range(N):
    if (k + 1) % 250 == 0:
        print(f"  iter {k+1}/{N}")

    kies = kie_sampler(rng)
    Sink_13C, Sink_D = compute_bulk_KIE(kies, SINK_FRACS)
    alpha_13C = 1.0 / Sink_13C
    alpha_D = 1.0 / Sink_D

    # ── Compute total source per year (mass balance) ─────────────────────
    SumSource = np.zeros(n_years)
    for j in range(n_years):
        SumSource[j] = (CH4[j+1] - CH4[j]) * PT + CH4[j] * PT / tau[j]

    # ── Atmospheric isotope observations for this iteration ──────────────
    # δ¹³C
    d13C_mc = d13C_iter[:n_years + 1, min(k, d13C_iter.shape[1] - 1)]
    # δD
    dD_mc = dD_iter[:n_years + 1, min(k, dD_iter.shape[1] - 1)]
    if len(dD_mc) < n_years + 1:
        dD_mc = np.concatenate([np.full(n_years + 1 - len(dD_mc), dD_mc[0]), dD_mc])

    # Convert to fractions
    f13 = delta_to_fraction_d13C(d13C_mc)
    fD = delta_to_fraction_dD(dD_mc)

    # Compute isotopic source signature (what the combined source must be)
    d13C_source = np.zeros(n_years)
    dD_source = np.zeros(n_years)
    for j in range(n_years):
        # ¹³C
        n13_now = f13[j] * CH4[j] * PT
        n13_next = f13[j+1] * CH4[j+1] * PT
        d13C_source[j] = (n13_next - n13_now + n13_now * alpha_13C / tau[j]) / SumSource[j]
        # D
        nD_now = fD[j] * CH4[j] * PT
        nD_next = fD[j+1] * CH4[j+1] * PT
        dD_source[j] = (nD_next - nD_now + nD_now * alpha_D / tau[j]) / SumSource[j]

    # ── Sample source end-member signatures ──────────────────────────────
    g_ff13 = rng.normal()
    g_bb13 = rng.normal()
    g_ffD = rng.normal()
    g_bbD = rng.normal()

    if k < FF_d13C_MC.shape[1]:
        ff_d13C_k = pad_to_length(np.array(FF_d13C_MC.iloc[:, k]).flatten(), n_years)
    else:
        ff_d13C_k = pad_to_length(ff_d13C_mean + g_ff13 * ff_d13C_U, n_years)

    if k < FF_dD_MC.shape[1]:
        ff_dD_k = pad_to_length(np.array(FF_dD_MC.iloc[:, k]).flatten(), n_years)
    else:
        ff_dD_k = pad_to_length(ff_dD_mean + g_ffD * ff_dD_U, n_years)

    bb_d13C_k = pad_to_length(bb_d13C_arr + g_bb13 * bb_d13C_U_arr, n_years)
    bb_dD_k = pad_to_length(bb_dD_arr + g_bbD * bb_dD_U_arr, n_years)

    if k < Mic_d13C_MC.shape[1]:
        mic_d13C_k = pad_to_length(np.array(Mic_d13C_MC.iloc[:n_years, k]).flatten(), n_years)
    else:
        mic_d13C_k = np.full(n_years, mic_d13C_mean_val)

    if k < Mic_dD_MC.shape[1]:
        mic_dD_k = pad_to_length(np.array(Mic_dD_MC.iloc[:n_years, k]).flatten(), n_years)
    else:
        mic_dD_k = np.full(n_years, mic_dD_mean_val)

    # ── 2×2 solve per year ───────────────────────────────────────────────
    for j in range(n_years):
        S = SumSource[j]
        BB_j = BB_annual[j]
        d13C_src_delta = fraction_to_delta_d13C(d13C_source[j])
        dD_src_delta = fraction_to_delta_dD(dD_source[j])

        # δ¹³C inversion
        denom = ff_d13C_k[j] - mic_d13C_k[j]
        if abs(denom) > 0.1:
            FF_j = (S * d13C_src_delta - mic_d13C_k[j] * (S - BB_j) - bb_d13C_k[j] * BB_j) / denom
            Mic_j = S - BB_j - FF_j
            if FF_j < 0 or Mic_j < 0:
                qm_d13C.record_negative()
                FF_j = max(0, FF_j); Mic_j = max(0, S - BB_j - FF_j)
            else:
                qm_d13C.record_ok()
        else:
            FF_j = np.nan; Mic_j = np.nan
        FF_d13C_out[j, k] = FF_j
        Mic_d13C_out[j, k] = Mic_j

        # δD inversion
        denom_D = ff_dD_k[j] - mic_dD_k[j]
        if abs(denom_D) > 1.0:
            FF_jD = (S * dD_src_delta - mic_dD_k[j] * (S - BB_j) - bb_dD_k[j] * BB_j) / denom_D
            Mic_jD = S - BB_j - FF_jD
            if FF_jD < 0 or Mic_jD < 0:
                qm_dD.record_negative()
                FF_jD = max(0, FF_jD); Mic_jD = max(0, S - BB_j - FF_jD)
            else:
                qm_dD.record_ok()
        else:
            FF_jD = np.nan; Mic_jD = np.nan
        FF_dD_out[j, k] = FF_jD
        Mic_dD_out[j, k] = Mic_jD

# ═══════════════════════════════════════════════════════════════════════════
# POST-PROCESSING
# ═══════════════════════════════════════════════════════════════════════════
print("\nSmoothing + saving...")
FF_d13C_s = smooth_5yr(FF_d13C_out)
Mic_d13C_s = smooth_5yr(Mic_d13C_out)
FF_dD_s = smooth_5yr(FF_dD_out)
Mic_dD_s = smooth_5yr(Mic_dD_out)

# Stats
print(f"\n{'='*60}")
print("RESULTS — 22_one (2×2, 1-box, BB-fixed)")
print(f"{'='*60}")
for label, arr in [('FF(δ¹³C)', FF_d13C_s), ('Mic(δ¹³C)', Mic_d13C_s),
                   ('FF(δD)', FF_dD_s), ('Mic(δD)', Mic_dD_s)]:
    print(f"  {label}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr, axis=0)):.1f} Tg/yr")
print(f"  BB (fixed): {np.mean(BB_annual):.1f} Tg/yr")

qm_d13C.summary()
qm_dD.summary()

# Save
results_df = pd.DataFrame({
    'Year': model_years,
    'FF_d13C_mean': np.nanmean(FF_d13C_s, axis=1),
    'FF_d13C_std': np.nanstd(FF_d13C_s, axis=1),
    'Mic_d13C_mean': np.nanmean(Mic_d13C_s, axis=1),
    'Mic_d13C_std': np.nanstd(Mic_d13C_s, axis=1),
    'FF_dD_mean': np.nanmean(FF_dD_s, axis=1),
    'FF_dD_std': np.nanstd(FF_dD_s, axis=1),
    'Mic_dD_mean': np.nanmean(Mic_dD_s, axis=1),
    'Mic_dD_std': np.nanstd(Mic_dD_s, axis=1),
    'BB_fixed': BB_annual,
})
results_df.to_csv(OUT_DIR / 'results_smoothed.csv', index=False)
save_report = {**qm_d13C.summary(), **qm_dD.summary(), 'config': str(config.get('FF_signature', 'EDGAR'))}
with open(OUT_DIR / 'quality_report.json', 'w') as f:
    json.dump(save_report, f, indent=2)

# ── Plot ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150, sharex=True)
fig.suptitle('22_one: 2×2 BB-Fixed, 1-Box — δ¹³C vs δD', fontsize=13)
for ax, name, arr13, arrD in [
    (axes[0], 'Fossil Fuel', FF_d13C_s, FF_dD_s),
    (axes[1], 'Microbial', Mic_d13C_s, Mic_dD_s),
]:
    m13 = np.nanmean(arr13, axis=1); s13 = 2 * np.nanstd(arr13, axis=1)
    mD = np.nanmean(arrD, axis=1); sD = 2 * np.nanstd(arrD, axis=1)
    ax.plot(model_years, m13, 'r-', lw=2, label='δ¹³C-derived')
    ax.fill_between(model_years, m13 - s13, m13 + s13, alpha=0.2, color='red')
    ax.plot(model_years, mD, 'b-', lw=2, label='δD-derived')
    ax.fill_between(model_years, mD - sD, mD + sD, alpha=0.2, color='blue')
    ax.set_ylabel(f'{name} (Tg/yr)')
    ax.set_xlabel('Year')
    ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / 'comparison_d13C_vs_dD.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\nDone. Output → {OUT_DIR}")
