#!/usr/bin/env python3
"""
22_two.py — 2×2 BB-Fixed, Two-Box (NH + SH) Model
====================================================
Infra: 2×2 system | 2 hemispheric boxes

Approach:
  Same as 22_one (fix BB, solve δ¹³C and δD separately) but in a
  two-hemisphere framework with interhemispheric exchange.

  Per hemisphere per year:
    S_hemi = (M_next − M_now) + M_now/τ_hemi − exchange_flux
    Isotopic source fraction from isotopic mass balance with exchange
    Then: FF = (S × δ_src − δ_Mic × (S − BB) − δ_BB × BB) / (δ_FF − δ_Mic)

Hemispheric features:
  - NH/SH CH₄ split via interhemispheric gradient (~80–100 ppb)
  - τ_ex sampled N(1.0, 0.1) yr  [Patra et al. (2011)]
  - δD NH/SH offset ±6‰  [Riddell-Young et al. (2025, PNAS)]
  - NH/SH lifetime ratio 0.95 / 1.05  [Lawrence et al. (2001)]
  - BB NH/SH split 55/45  [GFED4; van der Werf et al. (2017)]

Lineage: v3.2_bb_fixed_2x2.py cleaned up
"""

import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from models import (
    DEFAULT_CONFIG, SENSITIVITY_PRESETS,
    build_kie_sampler, compute_bulk_KIE, compute_lifetime,
    find_data_dirs, load_CH4, load_d13C_hemispheric,
    load_d13C_iterations, load_dD_iterations, load_source_signatures,
    load_BB_emissions, QualityMonitor, smooth_5yr, pad_to_length,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    PT, PT_HEMI,
)
from models.inputs import SINK_FRACTION_OPTIONS

# ═══════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════
config = SENSITIVITY_PRESETS.get(sys.argv[1], DEFAULT_CONFIG) if len(sys.argv) > 1 else DEFAULT_CONFIG
N = config['n_iterations']
SEED = config['seed']

dirs = find_data_dirs()
OUT_DIR = dirs['base'] / "Output_22_two"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("22_two: 2×2 BB-Fixed, Two-Hemisphere")
print("=" * 70)

CH4_years, CH4_global = load_CH4(dirs['data'])
_, c13_glob, c13_NH, c13_SH = load_d13C_hemispheric(dirs['data'])
d13C_iter = load_d13C_iterations(dirs['data'])
dD_mean, dD_iter = load_dD_iterations(dirs['data'])
sigs = load_source_signatures(dirs['src'], config)

n_years = len(CH4_global) - 1
model_years = np.arange(1999, 1999 + n_years)
tau_global = compute_lifetime(model_years, config)

# Hemispheric parameters
DD_OFFSET = config.get('dD_IH_offset', 6.0)
LT_NH = config.get('lifetime_ratio_NH', 0.95)
LT_SH = config.get('lifetime_ratio_SH', 1.05)
tau_NH = tau_global * LT_NH
tau_SH = tau_global * LT_SH
BB_NH_FRAC = config.get('BB_NH_fraction', 0.55)
BB_SH_FRAC = config.get('BB_SH_fraction', 0.45)
IH_GRAD = np.linspace(config.get('IH_gradient_start', 80),
                       config.get('IH_gradient_end', 100), len(CH4_global))
CH4_NH = CH4_global + IH_GRAD / 2.0
CH4_SH = CH4_global - IH_GRAD / 2.0

BB_annual = load_BB_emissions(dirs['data'], config, n_years)
BB_NH = BB_annual * BB_NH_FRAC
BB_SH = BB_annual * BB_SH_FRAC

# Sink fractions
SF_NH = SINK_FRACTION_OPTIONS[config.get('sink_fractions_NH', 'NH_default')]
SF_SH = SINK_FRACTION_OPTIONS[config.get('sink_fractions_SH', 'SH_default')]

# ── Source signatures (same prep as 22_one) ──────────────────────────────
ff_d13C_mean = np.array(sigs['FF_d13C_GlobUnc'].iloc[28:, 1]).flatten()
ff_d13C_U = np.array(sigs['FF_d13C_GlobUnc'].iloc[28:, 2]).flatten()
bb_d13C_arr = np.concatenate([np.array(sigs['BB_d13C'].iloc[1:, 1]).flatten(),
                               [sigs['BB_d13C'].iloc[-1, 1]]])
bb_d13C_U_arr = np.concatenate([np.array(sigs['BB_d13C'].iloc[1:, 2]).flatten(),
                                  [sigs['BB_d13C'].iloc[-1, 2]]])
mic_d13C_mean_val = sigs['Mic_d13C_ann'].iloc[:, 1].mean()

ff_dD_mean = np.array(sigs['FF_dD_GlobUnc'].iloc[34:, 1]).flatten()
ff_dD_U = np.array(sigs['FF_dD_GlobUnc'].iloc[34:, 2]).flatten()
bb_dD_tmp = np.array(sigs['BB_dD'].iloc[:, 1]).flatten()
bb_dD_U_tmp = np.array(sigs['BB_dD'].iloc[:, 2]).flatten()
bb_dD_arr = np.concatenate([np.full(3, bb_dD_tmp[-1]), bb_dD_tmp, [bb_dD_tmp[-1]]])
bb_dD_U_arr = np.concatenate([np.full(3, bb_dD_U_tmp[-1]), bb_dD_U_tmp, [bb_dD_U_tmp[-1]]])
mic_dD_mean_val = sigs['Mic_dD_ann'].iloc[:, 1].mean()
mic_dD_U_val = config['mic_dD_uncertainty']

# FF MC
if 'EDGAR' in str(config.get('FF_signature', 'EDGAR')):
    FF_d13C_MC = sigs['FF_d13C_MC'].iloc[28:, 1:] if sigs['FF_d13C_MC'].shape[0] > 28 else sigs['FF_d13C_MC'].iloc[:, 1:]
    FF_dD_MC = sigs['FF_dD_MC'].iloc[34:, 1:] if sigs['FF_dD_MC'].shape[0] > 34 else sigs['FF_dD_MC'].iloc[:, 1:]
else:
    FF_d13C_MC = sigs['FF_d13C_MC'].iloc[1:, 1:]
    FF_dD_MC = sigs['FF_dD_MC'].iloc[7:, 1:]

if FF_dD_MC.shape[0] < n_years:
    p = n_years - FF_dD_MC.shape[0]
    FF_dD_MC = pd.concat([pd.concat([FF_dD_MC.iloc[0:1]] * p)] + [FF_dD_MC], ignore_index=True)

Mic_d13C_MC = sigs['Mic_d13C_MC']
Mic_dD_MC = sigs['Mic_dD_MC']
if Mic_dD_MC.shape[0] < n_years:
    p = n_years - Mic_dD_MC.shape[0]
    Mic_dD_MC = pd.concat([pd.concat([Mic_dD_MC.iloc[0:1]] * p)] + [Mic_dD_MC], ignore_index=True)
elif Mic_dD_MC.shape[0] > n_years:
    Mic_dD_MC = Mic_dD_MC.iloc[:n_years]

if dD_iter.shape[0] < n_years + 1:
    p = n_years + 1 - dD_iter.shape[0]
    dD_iter = np.vstack([np.repeat(dD_iter[0:1], p, axis=0), dD_iter])

# ═══════════════════════════════════════════════════════════════════════════
# MONTE CARLO
# ═══════════════════════════════════════════════════════════════════════════
print(f"\nMC: {N} iterations | BB NH={np.mean(BB_NH):.1f}, SH={np.mean(BB_SH):.1f} Tg/yr")
rng = np.random.default_rng(seed=SEED)
kie_sampler = build_kie_sampler(config)

result_keys = ['FF_NH_d13C', 'Mic_NH_d13C', 'FF_SH_d13C', 'Mic_SH_d13C',
               'FF_NH_dD', 'Mic_NH_dD', 'FF_SH_dD', 'Mic_SH_dD']
results = {k: np.zeros((n_years, N)) for k in result_keys}
qm_d13C = QualityMonitor(n_years, N, "d13C")
qm_dD = QualityMonitor(n_years, N, "dD")

for k in range(N):
    if (k + 1) % 250 == 0:
        print(f"  iter {k+1}/{N}")

    tau_ex = max(0.5, rng.normal(1.0, 0.1)) if config.get('tau_ex_mode') == 'sampled' else 1.0
    kies = kie_sampler(rng)
    Sink_13C_NH, Sink_D_NH = compute_bulk_KIE(kies, SF_NH)
    Sink_13C_SH, Sink_D_SH = compute_bulk_KIE(kies, SF_SH)
    a13_NH, aD_NH = 1.0 / Sink_13C_NH, 1.0 / Sink_D_NH
    a13_SH, aD_SH = 1.0 / Sink_13C_SH, 1.0 / Sink_D_SH

    # Total source per hemisphere
    S_NH = np.zeros(n_years)
    S_SH = np.zeros(n_years)
    for j in range(n_years):
        M_NH = CH4_NH[j] * PT_HEMI
        M_SH = CH4_SH[j] * PT_HEMI
        ex_NH = (M_SH - M_NH) / tau_ex
        ex_SH = -ex_NH
        S_NH[j] = (CH4_NH[j+1] - CH4_NH[j]) * PT_HEMI + M_NH / tau_NH[j] - ex_NH
        S_SH[j] = (CH4_SH[j+1] - CH4_SH[j]) * PT_HEMI + M_SH / tau_SH[j] - ex_SH

    # Atmospheric δ for this iteration
    d13C_mc = d13C_iter[:n_years + 1, min(k, d13C_iter.shape[1] - 1)]
    n_c13 = min(len(c13_glob), n_years + 1)
    d13C_offset = d13C_mc[:n_c13] - c13_glob[:n_c13]
    d13C_NH_mc = c13_NH[:n_c13] + d13C_offset
    d13C_SH_mc = c13_SH[:n_c13] + d13C_offset

    dD_mc = dD_iter[:n_years + 1, min(k, dD_iter.shape[1] - 1)]
    if len(dD_mc) < n_years + 1:
        dD_mc = np.concatenate([np.full(n_years + 1 - len(dD_mc), dD_mc[0]), dD_mc])
    dD_NH_mc = dD_mc - DD_OFFSET
    dD_SH_mc = dD_mc + DD_OFFSET

    f13_NH = delta_to_fraction_d13C(d13C_NH_mc)
    f13_SH = delta_to_fraction_d13C(d13C_SH_mc)
    fD_NH = delta_to_fraction_dD(dD_NH_mc)
    fD_SH = delta_to_fraction_dD(dD_SH_mc)

    # Isotopic source signatures per hemi
    d13C_src_NH = np.zeros(n_years); d13C_src_SH = np.zeros(n_years)
    dD_src_NH = np.zeros(n_years); dD_src_SH = np.zeros(n_years)
    for j in range(n_years):
        # NH ¹³C
        n13_NH = f13_NH[j] * CH4_NH[j] * PT_HEMI
        n13_NH1 = f13_NH[j+1] * CH4_NH[j+1] * PT_HEMI
        n13_SH = f13_SH[j] * CH4_SH[j] * PT_HEMI
        d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH * a13_NH / tau_NH[j] -
                          (n13_SH - n13_NH) / tau_ex) / S_NH[j]
        # SH ¹³C
        n13_SH1 = f13_SH[j+1] * CH4_SH[j+1] * PT_HEMI
        d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH * a13_SH / tau_SH[j] -
                          (n13_NH - n13_SH) / tau_ex) / S_SH[j]
        # NH D
        nD_NH = fD_NH[j] * CH4_NH[j] * PT_HEMI
        nD_NH1 = fD_NH[j+1] * CH4_NH[j+1] * PT_HEMI
        nD_SH = fD_SH[j] * CH4_SH[j] * PT_HEMI
        dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH * aD_NH / tau_NH[j] -
                        (nD_SH - nD_NH) / tau_ex) / S_NH[j]
        # SH D
        nD_SH1 = fD_SH[j+1] * CH4_SH[j+1] * PT_HEMI
        dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH * aD_SH / tau_SH[j] -
                        (nD_NH - nD_SH) / tau_ex) / S_SH[j]

    # Source end-members (same draw for both hemispheres)
    g13, gD = rng.normal(), rng.normal()
    g_bb13, g_bbD = rng.normal(), rng.normal()
    ff_d13C_k = pad_to_length(np.array(FF_d13C_MC.iloc[:, min(k, FF_d13C_MC.shape[1]-1)]).flatten(), n_years) if k < FF_d13C_MC.shape[1] else pad_to_length(ff_d13C_mean + g13 * ff_d13C_U, n_years)
    ff_dD_k = pad_to_length(np.array(FF_dD_MC.iloc[:, min(k, FF_dD_MC.shape[1]-1)]).flatten(), n_years) if k < FF_dD_MC.shape[1] else pad_to_length(ff_dD_mean + gD * ff_dD_U, n_years)
    bb_d13C_k = pad_to_length(bb_d13C_arr + g_bb13 * bb_d13C_U_arr, n_years)
    bb_dD_k = pad_to_length(bb_dD_arr + g_bbD * bb_dD_U_arr, n_years)
    mic_d13C_k = pad_to_length(np.array(Mic_d13C_MC.iloc[:n_years, min(k, Mic_d13C_MC.shape[1]-1)]).flatten(), n_years) if k < Mic_d13C_MC.shape[1] else np.full(n_years, mic_d13C_mean_val)
    mic_dD_k = pad_to_length(np.array(Mic_dD_MC.iloc[:n_years, min(k, Mic_dD_MC.shape[1]-1)]).flatten(), n_years) if k < Mic_dD_MC.shape[1] else np.full(n_years, mic_dD_mean_val)

    # 2×2 solve per hemisphere per year
    for j in range(n_years):
        d13_NH_d = fraction_to_delta_d13C(d13C_src_NH[j])
        d13_SH_d = fraction_to_delta_d13C(d13C_src_SH[j])
        dD_NH_d = fraction_to_delta_dD(dD_src_NH[j])
        dD_SH_d = fraction_to_delta_dD(dD_src_SH[j])

        denom_13 = ff_d13C_k[j] - mic_d13C_k[j]
        denom_D = ff_dD_k[j] - mic_dD_k[j]

        for hemi, S, BB_j, d13_d, dD_d, suffix in [
            ('NH', S_NH[j], BB_NH[j], d13_NH_d, dD_NH_d, 'NH'),
            ('SH', S_SH[j], BB_SH[j], d13_SH_d, dD_SH_d, 'SH'),
        ]:
            # δ¹³C
            if abs(denom_13) > 0.1:
                FF_j = (S * d13_d - mic_d13C_k[j] * (S - BB_j) - bb_d13C_k[j] * BB_j) / denom_13
                Mic_j = S - BB_j - FF_j
                if FF_j < 0 or Mic_j < 0:
                    qm_d13C.record_negative()
                    FF_j = max(0, FF_j); Mic_j = max(0, S - BB_j - FF_j)
                else:
                    qm_d13C.record_ok()
            else:
                FF_j = Mic_j = np.nan
            results[f'FF_{suffix}_d13C'][j, k] = FF_j
            results[f'Mic_{suffix}_d13C'][j, k] = Mic_j

            # δD
            if abs(denom_D) > 1.0:
                FF_jD = (S * dD_d - mic_dD_k[j] * (S - BB_j) - bb_dD_k[j] * BB_j) / denom_D
                Mic_jD = S - BB_j - FF_jD
                if FF_jD < 0 or Mic_jD < 0:
                    qm_dD.record_negative()
                    FF_jD = max(0, FF_jD); Mic_jD = max(0, S - BB_j - FF_jD)
                else:
                    qm_dD.record_ok()
            else:
                FF_jD = Mic_jD = np.nan
            results[f'FF_{suffix}_dD'][j, k] = FF_jD
            results[f'Mic_{suffix}_dD'][j, k] = Mic_jD

# ═══════════════════════════════════════════════════════════════════════════
# POST-PROCESSING
# ═══════════════════════════════════════════════════════════════════════════
print("\nSmoothing...")
smoothed = {k: smooth_5yr(v) for k, v in results.items()}

# Globals = NH + SH
for iso in ['d13C', 'dD']:
    for src in ['FF', 'Mic']:
        smoothed[f'{src}_Global_{iso}'] = smoothed[f'{src}_NH_{iso}'] + smoothed[f'{src}_SH_{iso}']

print(f"\n{'='*60}")
print("RESULTS — 22_two (2×2, 2-box, BB-fixed)")
print(f"{'='*60}")
for label in ['FF_Global_d13C', 'Mic_Global_d13C', 'FF_Global_dD', 'Mic_Global_dD']:
    arr = smoothed[label]
    print(f"  {label}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr, axis=0)):.1f} Tg/yr")
print(f"  BB (fixed): {np.mean(BB_annual):.1f} Tg/yr")
qm_d13C.summary(); qm_dD.summary()

# Save
df_out = pd.DataFrame({'Year': model_years})
for k, v in smoothed.items():
    if 'Global' in k or 'NH' in k or 'SH' in k:
        df_out[f'{k}_mean'] = np.nanmean(v, axis=1)
        df_out[f'{k}_std'] = np.nanstd(v, axis=1)
df_out['BB_fixed'] = BB_annual
df_out.to_csv(OUT_DIR / 'results_smoothed.csv', index=False)

# Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=150, sharex=True)
fig.suptitle('22_two: 2×2 BB-Fixed, Two-Hemisphere — δ¹³C vs δD', fontsize=13)
for col, (src, label) in enumerate([('FF', 'Fossil Fuel'), ('Mic', 'Microbial')]):
    for row, scope in enumerate(['Global', 'NH']):
        ax = axes[row, col]
        m13 = np.nanmean(smoothed[f'{src}_{scope}_d13C'], axis=1)
        mD = np.nanmean(smoothed[f'{src}_{scope}_dD'], axis=1)
        s13 = 2 * np.nanstd(smoothed[f'{src}_{scope}_d13C'], axis=1)
        sD = 2 * np.nanstd(smoothed[f'{src}_{scope}_dD'], axis=1)
        ax.plot(model_years, m13, 'r-', lw=2, label='δ¹³C')
        ax.fill_between(model_years, m13 - s13, m13 + s13, alpha=0.2, color='red')
        ax.plot(model_years, mD, 'b-', lw=2, label='δD')
        ax.fill_between(model_years, mD - sD, mD + sD, alpha=0.2, color='blue')
        ax.set_ylabel(f'{label} (Tg/yr)')
        ax.set_title(f'{scope} {label}')
        ax.legend(); ax.grid(True, alpha=0.3)
        if row == 1: ax.set_xlabel('Year')
plt.tight_layout()
plt.savefig(OUT_DIR / 'comparison_d13C_vs_dD.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\nDone. Output → {OUT_DIR}")
