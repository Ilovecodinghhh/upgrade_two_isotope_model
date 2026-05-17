#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 4b: Fixed Two-Box Model with Correct Exchange Treatment
==============================================================
The original Phase 4 had inflated absolute values because the isotopic
exchange flux was not properly accounted for in the source δ derivation.

Fix: Include exchange isotopic flux in the isotope budget equation.

δ¹³C source isotope budget for NH:
  Σ_sources_NH × δ¹³C_src = d(n13_NH)/dt + Loss_NH × α × R13_NH
                            + F_ex(NH→SH) × R13_NH - F_ex(SH→NH) × R13_SH

Where R13 = n13 / (CH4 × PT_HEMI) is the isotope ratio in each hemisphere.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from common import (
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH, PT, PT_HEMI,
    BB_NH_FRACTION, BB_SH_FRACTION, DD_IH_OFFSET,
    LIFETIME_RATIO_NH, LIFETIME_RATIO_SH, TAU_EX_MEAN, TAU_EX_STD,
    compute_IH_gradient,
)

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "phase4b_two_box_fixed"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 1000
SEED = 42

W_MASS = 1.0
W_C13 = 1.0 / (1.0**2)
W_D = 1.0 / (5.0**2)


def run_two_box_fixed(oh13c_mode: str, label: str, dual: bool = True):
    """Run 2-box model with corrected exchange isotope treatment."""
    mode_str = "dual" if dual else "d13C-only"
    print(f"\n  2-box (fixed) | {mode_str} | OH_13C={oh13c_mode} ({label})")

    data = load_data(REPO_ROOT, two_box=True)
    n = data.n_years
    years = data.model_years
    CH4_NH = data.CH4_NH
    CH4_SH = data.CH4_SH
    rng = np.random.default_rng(SEED)

    tau_global = compute_lifetime(years, 'varying')
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    BB_global = data.BB_annual.copy()

    FF_all = np.zeros((n, N_ITER))
    Mic_all = np.zeros((n, N_ITER))

    for k in range(N_ITER):
        if (k + 1) % 500 == 0:
            print(f"      iter {k+1}/{N_ITER}")

        kies = sample_KIE(rng, 'sampled')
        if oh13c_mode == 'saueressig':
            kies['OH_13C'] = 1.0039
        elif oh13c_mode == 'cantrell':
            kies['OH_13C'] = 1.0054

        KIE_13C_NH, KIE_D_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        KIE_13C_SH, KIE_D_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        alpha_13C_NH = 1.0 / KIE_13C_NH
        alpha_D_NH = 1.0 / KIE_D_NH
        alpha_13C_SH = 1.0 / KIE_13C_SH
        alpha_D_SH = 1.0 / KIE_D_SH

        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))

        # Atmospheric observations
        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm = sample_atm_dD(data, k, n)

        # NH/SH δ¹³C
        c13_NH = data.c13_NH if data.c13_NH is not None else d13C_atm
        c13_SH = data.c13_SH if data.c13_SH is not None else d13C_atm

        # δD hemispheric
        dD_NH = dD_atm - DD_IH_OFFSET
        dD_SH = dD_atm + DD_IH_OFFSET

        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            BB_j = BB_global[j] if j < len(BB_global) else data.BB_global_mean
            BB_NH_j = BB_j * BB_NH_FRACTION
            BB_SH_j = BB_j * BB_SH_FRACTION

            # Exchange mass fluxes (Tg/yr)
            F_ex_NH_to_SH = CH4_NH[j] * PT_HEMI / tau_ex
            F_ex_SH_to_NH = CH4_SH[j] * PT_HEMI / tau_ex

            # Net source needed (mass balance):
            # S_NH = dM_NH/dt + Loss_NH + Exchange_out - Exchange_in
            S_NH = ((CH4_NH[j+1] - CH4_NH[j]) * PT_HEMI
                    + CH4_NH[j] * PT_HEMI / tau_NH[j]
                    + F_ex_NH_to_SH - F_ex_SH_to_NH)
            S_SH = ((CH4_SH[j+1] - CH4_SH[j]) * PT_HEMI
                    + CH4_SH[j] * PT_HEMI / tau_SH[j]
                    + F_ex_SH_to_NH - F_ex_NH_to_SH)

            # === CORRECTED isotope budget ===
            # For NH δ¹³C:
            # Σ_src × f13_src = d(n13_NH)/dt + n13_NH × α/τ + F_ex(→SH) × f13_NH - F_ex(←SH) × f13_SH
            # => f13_src = [d(n13)/dt + n13 × α/τ + F_ex_out × f13_NH - F_ex_in × f13_SH] / S

            # NH δ¹³C
            f13_NH_j = delta_to_fraction_d13C(np.array([c13_NH[j]]))[0]
            f13_NH_j1 = delta_to_fraction_d13C(np.array([c13_NH[min(j+1, len(c13_NH)-1)]]))[0]
            f13_SH_j = delta_to_fraction_d13C(np.array([c13_SH[j]]))[0]

            n13_NH_j = f13_NH_j * CH4_NH[j] * PT_HEMI
            n13_NH_j1 = f13_NH_j1 * CH4_NH[j+1] * PT_HEMI

            d13C_src_NH_f = (n13_NH_j1 - n13_NH_j
                             + n13_NH_j * alpha_13C_NH / tau_NH[j]
                             + F_ex_NH_to_SH * f13_NH_j
                             - F_ex_SH_to_NH * f13_SH_j) / max(S_NH, 1.0)
            d13C_src_NH = fraction_to_delta_d13C(d13C_src_NH_f)

            # SH δ¹³C
            f13_SH_j1 = delta_to_fraction_d13C(np.array([c13_SH[min(j+1, len(c13_SH)-1)]]))[0]
            n13_SH_j = f13_SH_j * CH4_SH[j] * PT_HEMI
            n13_SH_j1 = f13_SH_j1 * CH4_SH[j+1] * PT_HEMI

            d13C_src_SH_f = (n13_SH_j1 - n13_SH_j
                             + n13_SH_j * alpha_13C_SH / tau_SH[j]
                             + F_ex_SH_to_NH * f13_SH_j
                             - F_ex_NH_to_SH * f13_NH_j) / max(S_SH, 1.0)
            d13C_src_SH = fraction_to_delta_d13C(d13C_src_SH_f)

            # NH δD
            fD_NH_j = delta_to_fraction_dD(np.array([dD_NH[j]]))[0]
            fD_NH_j1 = delta_to_fraction_dD(np.array([dD_NH[min(j+1, len(dD_NH)-1)]]))[0]
            fD_SH_j = delta_to_fraction_dD(np.array([dD_SH[j]]))[0]

            nD_NH_j = fD_NH_j * CH4_NH[j] * PT_HEMI
            nD_NH_j1 = fD_NH_j1 * CH4_NH[j+1] * PT_HEMI

            dD_src_NH_f = (nD_NH_j1 - nD_NH_j
                           + nD_NH_j * alpha_D_NH / tau_NH[j]
                           + F_ex_NH_to_SH * fD_NH_j
                           - F_ex_SH_to_NH * fD_SH_j) / max(S_NH, 1.0)
            dD_src_NH = fraction_to_delta_dD(dD_src_NH_f)

            # SH δD
            fD_SH_j1 = delta_to_fraction_dD(np.array([dD_SH[min(j+1, len(dD_SH)-1)]]))[0]
            nD_SH_j = fD_SH_j * CH4_SH[j] * PT_HEMI
            nD_SH_j1 = fD_SH_j1 * CH4_SH[j+1] * PT_HEMI

            dD_src_SH_f = (nD_SH_j1 - nD_SH_j
                           + nD_SH_j * alpha_D_SH / tau_SH[j]
                           + F_ex_SH_to_NH * fD_SH_j
                           - F_ex_NH_to_SH * fD_NH_j) / max(S_SH, 1.0)
            dD_src_SH = fraction_to_delta_dD(dD_src_SH_f)

            # Solve per hemisphere
            ff_nh, mic_nh = _solve_hemi(S_NH, BB_NH_j, d13C_src_NH, dD_src_NH, sigs, j, dual)
            ff_sh, mic_sh = _solve_hemi(S_SH, BB_SH_j, d13C_src_SH, dD_src_SH, sigs, j, dual)

            FF_all[j, k] = ff_nh + ff_sh
            Mic_all[j, k] = mic_nh + mic_sh

    FF_s = smooth_5yr(FF_all)
    Mic_s = smooth_5yr(Mic_all)
    delta_ff, pct_ff = trend_change(FF_s, years)
    delta_mic, pct_mic = trend_change(Mic_s, years)

    suffix = "dual" if dual else "d13C"
    np.savez(OUT_DIR / f"run_{label}_{suffix}.npz",
             FF=FF_s, Mic=Mic_s, years=years,
             delta_ff=delta_ff, delta_mic=delta_mic)

    print(f"    FF: mean={np.median(FF_s, axis=1).mean():.0f} Tg/yr, "
          f"Δ={delta_ff.mean():+.1f}±{delta_ff.std():.1f}")
    print(f"    Mic: mean={np.median(Mic_s, axis=1).mean():.0f} Tg/yr, "
          f"Δ={delta_mic.mean():+.1f}±{delta_mic.std():.1f}")
    return delta_ff, delta_mic, FF_s, Mic_s


def _solve_hemi(S, BB_j, d13C_src, dD_src, sigs, j, dual):
    """Solve FF, Mic in one hemisphere via WLS or analytic."""
    if S <= 0:
        return 0.0, 0.0
    rhs = S - BB_j
    if rhs <= 0:
        return 0.0, 0.0

    if dual:
        A = np.array([
            [1.0, 1.0],
            [sigs['ff_d13C'][j], sigs['mic_d13C'][j]],
            [sigs['ff_dD'][j], sigs['mic_dD'][j]],
        ])
        b = np.array([
            rhs,
            S * d13C_src - BB_j * sigs['bb_d13C'][j],
            S * dD_src - BB_j * sigs['bb_dD'][j],
        ])
        w = np.sqrt(np.array([W_MASS, W_C13, W_D]))
        Aw = A * w[:, None]
        bw = b * w
        try:
            result = lsq_linear(Aw, bw, bounds=(0, np.inf))
            ff, mic = result.x
            # Enforce mass constraint
            if ff + mic > rhs * 1.5:
                scale = rhs / (ff + mic)
                ff *= scale
                mic *= scale
            return float(ff), float(mic)
        except Exception:
            return 0.0, max(0, rhs)
    else:
        denom = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
        if abs(denom) > 0.1:
            ff = (S * d13C_src - sigs['mic_d13C'][j] * rhs -
                  sigs['bb_d13C'][j] * BB_j) / denom
            mic = rhs - ff
            ff = max(0, min(ff, rhs))
            mic = max(0, rhs - ff)
            return ff, mic
        return 0.0, max(0, rhs)


def main():
    print("=" * 60)
    print("  Phase 4b — Two-Box (FIXED exchange isotopes)")
    print("=" * 60)

    results = {}
    for oh13c, lbl in [('saueressig', 'A_saueressig'), ('cantrell', 'B_cantrell'), ('sampled', 'C_sampled')]:
        for d in [True, False]:
            key = f"{lbl}_{'dual' if d else 'd13C'}"
            dff, dmic, FF_s, Mic_s = run_two_box_fixed(oh13c, lbl, dual=d)
            results[key] = {'delta_ff': dff, 'delta_mic': dmic, 'FF': FF_s, 'Mic': Mic_s}

    # KSR
    spread_s_ff = abs(results['B_cantrell_d13C']['delta_ff'].mean() -
                      results['A_saueressig_d13C']['delta_ff'].mean())
    spread_d_ff = abs(results['B_cantrell_dual']['delta_ff'].mean() -
                      results['A_saueressig_dual']['delta_ff'].mean())
    ksr_ff = spread_s_ff / max(spread_d_ff, 1e-6)

    spread_s_mic = abs(results['B_cantrell_d13C']['delta_mic'].mean() -
                       results['A_saueressig_d13C']['delta_mic'].mean())
    spread_d_mic = abs(results['B_cantrell_dual']['delta_mic'].mean() -
                       results['A_saueressig_dual']['delta_mic'].mean())
    ksr_mic = spread_s_mic / max(spread_d_mic, 1e-6)

    print(f"\n  FIXED 2-Box KSR (FF):  {ksr_ff:.2f}  (spreads: single={spread_s_ff:.2f}, dual={spread_d_ff:.2f})")
    print(f"  FIXED 2-Box KSR (Mic): {ksr_mic:.2f}  (spreads: single={spread_s_mic:.2f}, dual={spread_d_mic:.2f})")

    # Sanity check: mean emissions
    FF_c_dual = results['C_sampled_dual']['FF']
    Mic_c_dual = results['C_sampled_dual']['Mic']
    print(f"\n  Sanity: dual-isotope global FF mean = {np.median(FF_c_dual).mean():.0f} Tg/yr")
    print(f"  Sanity: dual-isotope global Mic mean = {np.median(Mic_c_dual).mean():.0f} Tg/yr")

    summary = {
        'KSR_FF_2box_fixed': float(ksr_ff),
        'KSR_Mic_2box_fixed': float(ksr_mic),
        'spread_single_FF': float(spread_s_ff),
        'spread_dual_FF': float(spread_d_ff),
        'spread_single_Mic': float(spread_s_mic),
        'spread_dual_Mic': float(spread_d_mic),
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # Figure
    years = results['C_sampled_dual']['FF'].shape[0]
    yrs = np.arange(1999, 1999 + years)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=300)

    # Top: emission time series (dual isotope, Run C)
    for ax, (key, title) in zip(axes[0], [('FF', 'Fossil Fuel'), ('Mic', 'Microbial')]):
        FF_d = results['C_sampled_dual'][key]
        FF_s = results['C_sampled_d13C'][key]
        med_d = np.median(FF_d, axis=1)
        lo_d, hi_d = np.percentile(FF_d, [2.5, 97.5], axis=1)
        med_s = np.median(FF_s, axis=1)
        lo_s, hi_s = np.percentile(FF_s, [2.5, 97.5], axis=1)

        ax.plot(yrs, med_s, 'r-', lw=2, label='δ¹³C-only')
        ax.fill_between(yrs, lo_s, hi_s, alpha=0.15, color='red')
        ax.plot(yrs, med_d, 'b-', lw=2, label='Dual isotope')
        ax.fill_between(yrs, lo_d, hi_d, alpha=0.15, color='blue')
        ax.set_ylabel('Emissions (Tg/yr)')
        ax.set_title(f'{title} — 2-Box (fixed)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Bottom: trend histograms
    for ax, (key, title) in zip(axes[1], [('delta_ff', 'FF Trend'), ('delta_mic', 'Mic Trend')]):
        for lbl, color, name in [('A_saueressig', 'tab:blue', 'Saueressig'),
                                  ('B_cantrell', 'tab:red', 'Cantrell')]:
            for d_flag, ls in [(True, '-'), (False, '--')]:
                k2 = f"{lbl}_{'dual' if d_flag else 'd13C'}"
                arr = results[k2][key]
                method = 'Dual' if d_flag else 'δ¹³C'
                ax.hist(arr, bins=40, alpha=0.3, color=color,
                        label=f"{name} ({method}): {arr.mean():+.1f}±{arr.std():.1f}",
                        density=True, histtype='stepfilled' if d_flag else 'step',
                        linewidth=2 if not d_flag else 1)
        ax.axvline(0, color='black', ls='--', lw=0.8)
        ax.set_xlabel('Δ Emissions (Tg/yr)')
        ax.set_title(title)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    plt.suptitle(f'2-Box KIE Sensitivity (Fixed Exchange)\nKSR(FF)={ksr_ff:.2f}, KSR(Mic)={ksr_mic:.2f}',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig5_2box_fixed.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure saved: {FIG_DIR / 'fig5_2box_fixed.png'}")
    print("\n✓ Phase 4b complete.")


if __name__ == "__main__":
    main()
