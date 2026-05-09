#!/usr/bin/env python3
"""
33_one.py — 3×3 Dual-Isotope, One-Box (Global) Model
======================================================
Infra: 3×3 system | 1 global box

Approach:
  Use δ¹³C AND δD simultaneously in a 3-equation system:
    [1]   S = FF + Mic + BB
    [2]   S × f¹³C_src = FF × f¹³C_FF + Mic × f¹³C_Mic + BB × f¹³C_BB
    [3]   S × fD_src   = FF × fD_FF   + Mic × fD_Mic   + BB × fD_BB

  Solve for all three unknowns (FF, Mic, BB) simultaneously.
  Uses bounded least-squares to enforce non-negativity.

Advantage:  No need to fix BB externally; all 3 sources from isotopes.
Limitation: System can be ill-conditioned when δD end-members overlap.

References:
  - Schaefer et al. (2016, Science) — original 3×3 approach
  - Chandra et al. (2024, Comm. Earth Environ.) — KIE sensitivity
  - Bao et al. (this work) — stochastic KIE, time-varying τ

Lineage: v2.0_upgraded_box_model.py (1-box, 3×3 implicit)
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
    fraction_to_delta_d13C, fraction_to_delta_dD,
    PT, C13Std, DStd,
)
from models.inputs import SINK_FRACTION_OPTIONS

# ═══════════════════════════════════════════════════════════════════════════
config = SENSITIVITY_PRESETS.get(sys.argv[1], DEFAULT_CONFIG) if len(sys.argv) > 1 else DEFAULT_CONFIG
N = config['n_iterations']
SEED = config['seed']

dirs = find_data_dirs()
OUT_DIR = dirs['base'] / "Output_33_one"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("33_one: 3×3 Dual-Isotope, One-Box (Global)")
print("=" * 70)

CH4_years, CH4 = load_CH4(dirs['data'])
_, c13_glob, c13_nh, c13_sh = load_d13C_hemispheric(dirs['data'])
d13C_iter = load_d13C_iterations(dirs['data'])
dD_mean, dD_iter = load_dD_iterations(dirs['data'])
sigs = load_source_signatures(dirs['src'], config)

n_years = len(CH4) - 1
model_years = np.arange(1999, 1999 + n_years)
tau = compute_lifetime(model_years, config)

# Sink fractions (global)
sf_key = config.get('sink_fractions_global', 'global_default')
SINK_FRACS = SINK_FRACTION_OPTIONS[sf_key]

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
# MONTE CARLO — 3×3
# ═══════════════════════════════════════════════════════════════════════════
print(f"\nMC: {N} iterations | 3×3 simultaneous δ¹³C + δD")
rng = np.random.default_rng(seed=SEED)
kie_sampler = build_kie_sampler(config)

BB_out = np.zeros((n_years, N))
FF_out = np.zeros((n_years, N))
Mic_out = np.zeros((n_years, N))
qm = QualityMonitor(n_years, N, "3×3_global")

for k in range(N):
    if (k + 1) % 250 == 0:
        print(f"  iter {k+1}/{N}")

    kies = kie_sampler(rng)
    Sink_13C, Sink_D = compute_bulk_KIE(kies, SINK_FRACS)
    alpha_13C = 1.0 / Sink_13C
    alpha_D = 1.0 / Sink_D

    # Total source
    SumSource = np.zeros(n_years)
    for j in range(n_years):
        SumSource[j] = (CH4[j+1] - CH4[j]) * PT + CH4[j] * PT / tau[j]

    # Atmospheric observations
    d13C_mc = d13C_iter[:n_years + 1, min(k, d13C_iter.shape[1] - 1)]
    dD_mc = dD_iter[:n_years + 1, min(k, dD_iter.shape[1] - 1)]
    if len(dD_mc) < n_years + 1:
        dD_mc = np.concatenate([np.full(n_years + 1 - len(dD_mc), dD_mc[0]), dD_mc])

    f13 = delta_to_fraction_d13C(d13C_mc)
    fD = delta_to_fraction_dD(dD_mc)

    # Isotopic source signatures
    d13C_source = np.zeros(n_years)
    dD_source = np.zeros(n_years)
    for j in range(n_years):
        n13_now = f13[j] * CH4[j] * PT
        n13_next = f13[j+1] * CH4[j+1] * PT
        d13C_source[j] = (n13_next - n13_now + n13_now * alpha_13C / tau[j]) / SumSource[j]
        nD_now = fD[j] * CH4[j] * PT
        nD_next = fD[j+1] * CH4[j+1] * PT
        dD_source[j] = (nD_next - nD_now + nD_now * alpha_D / tau[j]) / SumSource[j]

    # Source end-members
    g_ff13, g_bb13 = rng.normal(), rng.normal()
    g_ffD, g_bbD = rng.normal(), rng.normal()
    ff_d13C_k = pad_to_length(np.array(FF_d13C_MC.iloc[:, min(k, FF_d13C_MC.shape[1]-1)]).flatten(), n_years) if k < FF_d13C_MC.shape[1] else pad_to_length(ff_d13C_mean + g_ff13 * ff_d13C_U, n_years)
    ff_dD_k = pad_to_length(np.array(FF_dD_MC.iloc[:, min(k, FF_dD_MC.shape[1]-1)]).flatten(), n_years) if k < FF_dD_MC.shape[1] else pad_to_length(ff_dD_mean + g_ffD * ff_dD_U, n_years)
    bb_d13C_k = pad_to_length(bb_d13C_arr + g_bb13 * bb_d13C_U_arr, n_years)
    bb_dD_k = pad_to_length(bb_dD_arr + g_bbD * bb_dD_U_arr, n_years)
    mic_d13C_k = pad_to_length(np.array(Mic_d13C_MC.iloc[:n_years, min(k, Mic_d13C_MC.shape[1]-1)]).flatten(), n_years) if k < Mic_d13C_MC.shape[1] else np.full(n_years, mic_d13C_mean_val)
    mic_dD_k = pad_to_length(np.array(Mic_dD_MC.iloc[:n_years, min(k, Mic_dD_MC.shape[1]-1)]).flatten(), n_years) if k < Mic_dD_MC.shape[1] else np.full(n_years, mic_dD_mean_val)

    # Convert to fractions for matrix
    f13_bb = delta_to_fraction_d13C(bb_d13C_k)
    f13_ff = delta_to_fraction_d13C(ff_d13C_k)
    f13_mic = delta_to_fraction_d13C(mic_d13C_k)
    fD_bb = delta_to_fraction_dD(bb_dD_k)
    fD_ff = delta_to_fraction_dD(ff_dD_k)
    fD_mic = delta_to_fraction_dD(mic_dD_k)

    # 3×3 solve per year (bounded least-squares, non-negative)
    for j in range(n_years):
        A = np.array([
            [1.0, 1.0, 1.0],
            [f13_bb[j], f13_ff[j], f13_mic[j]],
            [fD_bb[j], fD_ff[j], fD_mic[j]],
        ])
        B = np.array([
            SumSource[j],
            SumSource[j] * d13C_source[j],
            SumSource[j] * dD_source[j],
        ])
        # Weight: mass balance >> δ¹³C > δD  (δD has larger measurement noise)
        W = np.diag([100.0, 1.0, 0.5])
        ub = SumSource[j] * 1.5
        try:
            result = lsq_linear(W @ A, W @ B, bounds=(0, ub))
            x = result.x
            qm.record_cond(j, k, A)
            if np.any(x < 0):
                qm.record_negative()
            else:
                qm.record_ok()
        except Exception:
            x = np.array([np.nan, np.nan, np.nan])
            qm.record_negative()

        BB_out[j, k] = x[0]
        FF_out[j, k] = x[1]
        Mic_out[j, k] = x[2]

# ═══════════════════════════════════════════════════════════════════════════
# POST-PROCESSING
# ═══════════════════════════════════════════════════════════════════════════
print("\nSmoothing + saving...")
BB_s = smooth_5yr(BB_out)
FF_s = smooth_5yr(FF_out)
Mic_s = smooth_5yr(Mic_out)

print(f"\n{'='*60}")
print("RESULTS — 33_one (3×3, 1-box)")
print(f"{'='*60}")
for label, arr in [('BB', BB_s), ('FF', FF_s), ('Mic', Mic_s)]:
    print(f"  {label}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr, axis=0)):.1f} Tg/yr")
qm.summary()

df = pd.DataFrame({
    'Year': model_years,
    'BB_mean': np.nanmean(BB_s, axis=1), 'BB_std': np.nanstd(BB_s, axis=1),
    'FF_mean': np.nanmean(FF_s, axis=1), 'FF_std': np.nanstd(FF_s, axis=1),
    'Mic_mean': np.nanmean(Mic_s, axis=1), 'Mic_std': np.nanstd(Mic_s, axis=1),
})
df.to_csv(OUT_DIR / 'results_smoothed.csv', index=False)

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=150, sharex=True)
fig.suptitle('33_one: 3×3 Dual-Isotope, 1-Box — All Sources', fontsize=13)
for ax, (name, arr, c) in zip(axes, [('BB', BB_s, 'red'), ('FF', FF_s, 'blue'), ('Mic', Mic_s, 'green')]):
    m = np.nanmean(arr, axis=1); s = 2 * np.nanstd(arr, axis=1)
    ax.plot(model_years, m, '-', lw=2.5, color=c)
    ax.fill_between(model_years, m - s, m + s, alpha=0.3, color=c)
    ax.set_ylabel(f'{name} (Tg/yr)'); ax.set_xlabel('Year')
    ax.set_title(name); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / 'sources_3x3.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\nDone. Output → {OUT_DIR}")
