#!/usr/bin/env python3
"""
33_two.py — 3×3 Dual-Isotope, Two-Box (NH + SH) Model
========================================================
Infra: 3×3 system | 2 hemispheric boxes

Approach:
  Use δ¹³C AND δD simultaneously in a 3×3 system PER HEMISPHERE:
    [1]  S_hemi = BB_hemi + FF_hemi + Mic_hemi
    [2]  S_hemi × f¹³C_src = BB × f¹³C_BB + FF × f¹³C_FF + Mic × f¹³C_Mic
    [3]  S_hemi × fD_src   = BB × fD_BB   + FF × fD_FF   + Mic × fD_Mic

  With interhemispheric exchange, hemisphere-specific lifetimes,
  and δD NH/SH offset.

  Bounded least-squares enforces non-negativity. Weighting handles
  ill-conditioning from δD (W = diag(100, 1, 0.5) for NH; adjusted for SH).

References:
  - Naus et al. (2019, ACP) — hemispheric box model framework
  - Riddell-Young et al. (2025, PNAS) — δD hemispheric offset
  - Patra et al. (2011, ACP) — interhemispheric exchange time

Lineage: v3.1_optimized_3x3.py cleaned up
"""

import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import lsq_linear

sys.path.insert(0, str(Path(__file__).resolve().parent))
from models import (
    DEFAULT_CONFIG, SENSITIVITY_PRESETS,
    build_kie_sampler, compute_bulk_KIE, compute_lifetime,
    find_data_dirs, load_CH4, load_d13C_hemispheric,
    load_d13C_iterations, load_dD_iterations, load_source_signatures,
    QualityMonitor, smooth_5yr, pad_to_length,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    PT, PT_HEMI,
)
from models.inputs import SINK_FRACTION_OPTIONS

# ═══════════════════════════════════════════════════════════════════════════
config = SENSITIVITY_PRESETS.get(sys.argv[1], DEFAULT_CONFIG) if len(sys.argv) > 1 else DEFAULT_CONFIG
N = config['n_iterations']
SEED = config['seed']

dirs = find_data_dirs()
OUT_DIR = dirs['base'] / "Output_33_two"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("33_two: 3×3 Dual-Isotope, Two-Hemisphere")
print("=" * 70)

CH4_years, CH4_global = load_CH4(dirs['data'])
_, c13_glob, c13_NH, c13_SH = load_d13C_hemispheric(dirs['data'])
d13C_iter = load_d13C_iterations(dirs['data'])
dD_mean, dD_iter = load_dD_iterations(dirs['data'])
sigs = load_source_signatures(dirs['src'], config)

n_years = len(CH4_global) - 1
model_years = np.arange(1999, 1999 + n_years)
tau_global = compute_lifetime(model_years, config)

# Hemispheric config
DD_OFFSET = config.get('dD_IH_offset', 6.0)
LT_NH = config.get('lifetime_ratio_NH', 0.95)
LT_SH = config.get('lifetime_ratio_SH', 1.05)
tau_NH = tau_global * LT_NH
tau_SH = tau_global * LT_SH
IH_GRAD = np.linspace(config.get('IH_gradient_start', 80),
                       config.get('IH_gradient_end', 100), len(CH4_global))
CH4_NH = CH4_global + IH_GRAD / 2.0
CH4_SH = CH4_global - IH_GRAD / 2.0

SF_NH = SINK_FRACTION_OPTIONS[config.get('sink_fractions_NH', 'NH_default')]
SF_SH = SINK_FRACTION_OPTIONS[config.get('sink_fractions_SH', 'SH_default')]

# ── Source signatures ────────────────────────────────────────────────────
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
# MONTE CARLO — 3×3 per hemisphere
# ═══════════════════════════════════════════════════════════════════════════
print(f"\nMC: {N} iterations | 3×3 × 2 hemispheres")
rng = np.random.default_rng(seed=SEED)
kie_sampler = build_kie_sampler(config)

keys = ['BB_NH', 'FF_NH', 'Mic_NH', 'BB_SH', 'FF_SH', 'Mic_SH']
results = {k: np.zeros((n_years, N)) for k in keys}
qm_NH = QualityMonitor(n_years, N, "NH")
qm_SH = QualityMonitor(n_years, N, "SH")

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
    S_NH = np.zeros(n_years); S_SH = np.zeros(n_years)
    for j in range(n_years):
        M_NH = CH4_NH[j] * PT_HEMI; M_SH = CH4_SH[j] * PT_HEMI
        ex = (M_SH - M_NH) / tau_ex
        S_NH[j] = (CH4_NH[j+1] - CH4_NH[j]) * PT_HEMI + M_NH / tau_NH[j] - ex
        S_SH[j] = (CH4_SH[j+1] - CH4_SH[j]) * PT_HEMI + M_SH / tau_SH[j] + ex

    # Atmospheric δ
    d13C_mc = d13C_iter[:n_years+1, min(k, d13C_iter.shape[1]-1)]
    n_c13 = min(len(c13_glob), n_years+1)
    d13C_off = d13C_mc[:n_c13] - c13_glob[:n_c13]
    d13C_NH_mc = c13_NH[:n_c13] + d13C_off
    d13C_SH_mc = c13_SH[:n_c13] + d13C_off

    dD_mc = dD_iter[:n_years+1, min(k, dD_iter.shape[1]-1)]
    if len(dD_mc) < n_years+1:
        dD_mc = np.concatenate([np.full(n_years+1-len(dD_mc), dD_mc[0]), dD_mc])
    dD_NH_mc = dD_mc - DD_OFFSET
    dD_SH_mc = dD_mc + DD_OFFSET

    f13_NH = delta_to_fraction_d13C(d13C_NH_mc)
    f13_SH = delta_to_fraction_d13C(d13C_SH_mc)
    fD_NH = delta_to_fraction_dD(dD_NH_mc)
    fD_SH = delta_to_fraction_dD(dD_SH_mc)

    # Isotopic source fractions per hemisphere
    d13C_src_NH = np.zeros(n_years); d13C_src_SH = np.zeros(n_years)
    dD_src_NH = np.zeros(n_years); dD_src_SH = np.zeros(n_years)
    for j in range(n_years):
        n13_NH = f13_NH[j]*CH4_NH[j]*PT_HEMI; n13_NH1 = f13_NH[j+1]*CH4_NH[j+1]*PT_HEMI
        n13_SH = f13_SH[j]*CH4_SH[j]*PT_HEMI; n13_SH1 = f13_SH[j+1]*CH4_SH[j+1]*PT_HEMI
        nD_NH = fD_NH[j]*CH4_NH[j]*PT_HEMI; nD_NH1 = fD_NH[j+1]*CH4_NH[j+1]*PT_HEMI
        nD_SH = fD_SH[j]*CH4_SH[j]*PT_HEMI; nD_SH1 = fD_SH[j+1]*CH4_SH[j+1]*PT_HEMI

        d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH*a13_NH/tau_NH[j] - (n13_SH-n13_NH)/tau_ex) / S_NH[j]
        d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH*a13_SH/tau_SH[j] - (n13_NH-n13_SH)/tau_ex) / S_SH[j]
        dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH*aD_NH/tau_NH[j] - (nD_SH-nD_NH)/tau_ex) / S_NH[j]
        dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH*aD_SH/tau_SH[j] - (nD_NH-nD_SH)/tau_ex) / S_SH[j]

    # Source end-members
    g13, gD, g_bb13, g_bbD = rng.normal(), rng.normal(), rng.normal(), rng.normal()
    ff_d13C_k = pad_to_length(np.array(FF_d13C_MC.iloc[:, min(k, FF_d13C_MC.shape[1]-1)]).flatten(), n_years) if k < FF_d13C_MC.shape[1] else pad_to_length(ff_d13C_mean + g13*ff_d13C_U, n_years)
    ff_dD_k = pad_to_length(np.array(FF_dD_MC.iloc[:, min(k, FF_dD_MC.shape[1]-1)]).flatten(), n_years) if k < FF_dD_MC.shape[1] else pad_to_length(ff_dD_mean + gD*ff_dD_U, n_years)
    bb_d13C_k = pad_to_length(bb_d13C_arr + g_bb13*bb_d13C_U_arr, n_years)
    bb_dD_k = pad_to_length(bb_dD_arr + g_bbD*bb_dD_U_arr, n_years)
    mic_d13C_k = pad_to_length(np.array(Mic_d13C_MC.iloc[:n_years, min(k, Mic_d13C_MC.shape[1]-1)]).flatten(), n_years) if k < Mic_d13C_MC.shape[1] else np.full(n_years, mic_d13C_mean_val)
    mic_dD_k = pad_to_length(np.array(Mic_dD_MC.iloc[:n_years, min(k, Mic_dD_MC.shape[1]-1)]).flatten(), n_years) if k < Mic_dD_MC.shape[1] else np.full(n_years, mic_dD_mean_val)

    f13_bb = delta_to_fraction_d13C(bb_d13C_k)
    f13_ff = delta_to_fraction_d13C(ff_d13C_k)
    f13_mic = delta_to_fraction_d13C(mic_d13C_k)
    fD_bb = delta_to_fraction_dD(bb_dD_k)
    fD_ff = delta_to_fraction_dD(ff_dD_k)
    fD_mic = delta_to_fraction_dD(mic_dD_k)

    # 3×3 solve per hemisphere per year
    for j in range(n_years):
        for hemi, S_val, d13_src, dD_src, qm, suffix, W in [
            ('NH', S_NH[j], d13C_src_NH[j], dD_src_NH[j], qm_NH, 'NH', np.diag([100., 1., 0.5])),
            ('SH', S_SH[j], d13C_src_SH[j], dD_src_SH[j], qm_SH, 'SH', np.diag([200., 1., 0.5])),
        ]:
            A = np.array([
                [1.0, 1.0, 1.0],
                [f13_bb[j], f13_ff[j], f13_mic[j]],
                [fD_bb[j], fD_ff[j], fD_mic[j]],
            ])
            B = np.array([S_val, S_val * d13_src, S_val * dD_src])
            ub = S_val * 1.5
            try:
                res = lsq_linear(W @ A, W @ B, bounds=(0, ub))
                x = res.x
                qm.record_cond(j, k, A)
                qm.record_ok()
            except Exception:
                x = np.array([np.nan, np.nan, np.nan])
                qm.record_negative()

            results[f'BB_{suffix}'][j, k] = x[0]
            results[f'FF_{suffix}'][j, k] = x[1]
            results[f'Mic_{suffix}'][j, k] = x[2]

# ═══════════════════════════════════════════════════════════════════════════
# POST-PROCESSING
# ═══════════════════════════════════════════════════════════════════════════
print("\nSmoothing...")
smoothed = {k: smooth_5yr(v) for k, v in results.items()}
for src in ['BB', 'FF', 'Mic']:
    smoothed[f'{src}_Global'] = smoothed[f'{src}_NH'] + smoothed[f'{src}_SH']

print(f"\n{'='*60}")
print("RESULTS — 33_two (3×3, 2-box)")
print(f"{'='*60}")
for label in ['BB_Global', 'FF_Global', 'Mic_Global', 'BB_NH', 'FF_NH', 'Mic_NH', 'BB_SH', 'FF_SH', 'Mic_SH']:
    arr = smoothed[label]
    print(f"  {label}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr, axis=0)):.1f} Tg/yr")
qm_NH.summary(); qm_SH.summary()

# Save
df = pd.DataFrame({'Year': model_years})
for k, v in smoothed.items():
    df[f'{k}_mean'] = np.nanmean(v, axis=1)
    df[f'{k}_std'] = np.nanstd(v, axis=1)
df.to_csv(OUT_DIR / 'results_smoothed.csv', index=False)

# Plot
fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=150, sharex=True)
fig.suptitle('33_two: 3×3 Dual-Isotope, Two-Hemisphere', fontsize=13)
for row, (src, c) in enumerate([('BB', 'red'), ('FF', 'blue'), ('Mic', 'green')]):
    for col, scope in enumerate(['NH', 'SH']):
        ax = axes[row, col]
        m = np.nanmean(smoothed[f'{src}_{scope}'], axis=1)
        s = 2 * np.nanstd(smoothed[f'{src}_{scope}'], axis=1)
        ax.plot(model_years, m, '-', lw=2.5, color=c)
        ax.fill_between(model_years, m - s, m + s, alpha=0.3, color=c)
        ax.set_ylabel(f'{src} (Tg/yr)')
        if row == 0: ax.set_title(scope)
        if row == 2: ax.set_xlabel('Year')
        ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / 'sources_hemispheric.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\nDone. Output → {OUT_DIR}")
