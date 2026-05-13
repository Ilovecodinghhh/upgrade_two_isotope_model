#!/usr/bin/env python3
"""
v3: Delta-space solver — fixes the fundamental numerical problem.

The v1/v2 models solve in FRACTION space:
  A = [[1, 1, 1], [f13_BB, f13_FF, f13_Mic], [fD_BB, fD_FF, fD_Mic]]
  
where f13 ~ 0.0104 and fD ~ 0.00012. This gives effective rank 1 
(cond ~ 170k) because all isotope-fraction values are nearly identical.

v3 solves in DELTA space:
  Step 1: Compute source fractions f_FF, f_Mic, f_BB (summing to 1)
    by solving a well-conditioned (cond ~ 27) 3×3 system in ‰ space
  Step 2: Multiply by total source S to get absolute emissions
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy import stats as sp_stats
import json

from common import (
    ModelConfig, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH, SINK_FRACTIONS_GLOBAL,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    TAU_EX_MEAN, TAU_EX_STD,
)

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results" / "v3_delta_space"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

NI = 1000
TREND_START = 2007
SEED = 42


def realistic_IH_gradient(years):
    anchor_years = np.array([2000, 2010, 2020, 2022])
    anchor_grad = np.array([108.0, 120.0, 140.0, 145.0])
    return np.interp(years, anchor_years, anchor_grad)


def fraction_to_delta_d13C(f):
    """Convert ¹³C fraction to delta (‰ VPDB). PDB R=0.0112372."""
    R_std = 0.0112372
    R_sample = f / (1.0 - f)
    return (R_sample / R_std - 1.0) * 1000.0


def fraction_to_delta_dD(f):
    """Convert D fraction to delta (‰ VSMOW). VSMOW R=0.00015576."""
    R_std = 0.00015576
    R_sample = f / (1.0 - f)
    return (R_sample / R_std - 1.0) * 1000.0


def solve_delta_space(S, d13C_src_delta, dD_src_delta,
                      d13C_FF, d13C_Mic, d13C_BB,
                      dD_FF, dD_Mic, dD_BB):
    """
    Solve for source fractions in delta (‰) space.
    
    System:
      f_BB + f_FF + f_Mic = 1
      f_BB*δ¹³C_BB + f_FF*δ¹³C_FF + f_Mic*δ¹³C_Mic = δ¹³C_src
      f_BB*δD_BB + f_FF*δD_FF + f_Mic*δD_Mic = δD_src
    
    Returns (BB, FF, Mic) in Tg/yr.
    """
    A = np.array([
        [1.0, 1.0, 1.0],
        [d13C_BB, d13C_FF, d13C_Mic],
        [dD_BB, dD_FF, dD_Mic],
    ])
    b = np.array([1.0, d13C_src_delta, dD_src_delta])
    
    # Scale rows for better conditioning
    # Row 0: O(1), Row 1: O(50), Row 2: O(250)
    # Scale all to O(1)
    scale = np.array([1.0, 1.0/50.0, 1.0/250.0])
    A_scaled = A * scale[:, None]
    b_scaled = b * scale
    
    # Solve with non-negativity + upper bound = 1
    try:
        res = lsq_linear(A_scaled, b_scaled, bounds=(0.0, 1.0))
        fracs = res.x
        # Renormalize fractions to sum to 1 (should be close already)
        if fracs.sum() > 0:
            fracs /= fracs.sum()
        return fracs[0] * S, fracs[1] * S, fracs[2] * S  # BB, FF, Mic
    except:
        return np.nan, np.nan, np.nan


def compute_source_delta(f_atm, f_atm_next, M, M_next, alpha, tau, 
                         f_other_hemi=None, M_other=None, tau_ex=None):
    """
    Compute bulk source isotopic composition in FRACTION space,
    then convert to delta.
    
    For 1-box: dn/dt = S*f_src - n*α/τ
    For 2-box: dn/dt = S*f_src - n*α/τ - (n - n_other)/τ_ex
    """
    n = f_atm * M
    n_next = f_atm_next * M_next
    
    if f_other_hemi is not None and M_other is not None and tau_ex is not None:
        n_other = f_other_hemi * M_other
        f_src = (n_next - n + n * alpha / tau - (n_other - n) / tau_ex)
    else:
        f_src = (n_next - n + n * alpha / tau)
    
    return f_src  # This is S * f_src_fraction, divide by S to get f_src


def run_two_box_v3(data, cfg):
    """2-box model with delta-space solver."""
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    all_years = np.arange(years[0], years[-1] + 2)
    IH_grad = realistic_IH_gradient(all_years.astype(float))
    CH4_NH = data.CH4_global + IH_grad / 2.0
    CH4_SH = data.CH4_global - IH_grad / 2.0

    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global

    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    NI = cfg.n_iterations
    FF_NH = np.zeros((n, NI)); FF_SH = np.zeros((n, NI))
    Mic_NH = np.zeros((n, NI)); Mic_SH = np.zeros((n, NI))
    BB_NH = np.zeros((n, NI)); BB_SH = np.zeros((n, NI))

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        kies = sample_KIE(rng, cfg.kie_mode)
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        a13_NH = 1.0 / K13_NH; aD_NH = 1.0 / KD_NH
        a13_SH = 1.0 / K13_SH; aD_SH = 1.0 / KD_SH

        # Total sources
        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i] * PT_HEMI; M_NH1 = CH4_NH[i+1] * PT_HEMI
            M_SH = CH4_SH[i] * PT_HEMI; M_SH1 = CH4_SH[i+1] * PT_HEMI
            S_NH[i] = (M_NH1 - M_NH) + M_NH / tau_NH[i] - (M_SH - M_NH) / tau_ex
            S_SH[i] = (M_SH1 - M_SH) + M_SH / tau_SH[i] - (M_NH - M_SH) / tau_ex

        # Atmospheric observations (MC-sampled)
        d13C_glob_MC = sample_atm_d13C(data, k, n)
        dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off

        # Convert to fractions for mass-balance
        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
        fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
        fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)

        # Compute source isotopic compositions in FRACTION space
        d13C_src_frac_NH = np.zeros(n); d13C_src_frac_SH = np.zeros(n)
        dD_src_frac_NH = np.zeros(n); dD_src_frac_SH = np.zeros(n)

        for j in range(n):
            M_NH = CH4_NH[j] * PT_HEMI; M_NH1 = CH4_NH[j+1] * PT_HEMI
            M_SH = CH4_SH[j] * PT_HEMI; M_SH1 = CH4_SH[j+1] * PT_HEMI

            # 13C NH
            n13 = f13_NH_atm[j] * M_NH; n13_1 = f13_NH_atm[j+1] * M_NH1
            n13_SH = f13_SH_atm[j] * M_SH
            d13C_src_frac_NH[j] = (n13_1 - n13 + n13 * a13_NH / tau_NH[j] - (n13_SH - n13) / tau_ex) / S_NH[j]

            # 13C SH
            n13 = f13_SH_atm[j] * M_SH; n13_1 = f13_SH_atm[j+1] * M_SH1
            n13_NH = f13_NH_atm[j] * M_NH
            d13C_src_frac_SH[j] = (n13_1 - n13 + n13 * a13_SH / tau_SH[j] - (n13_NH - n13) / tau_ex) / S_SH[j]

            # D NH
            nD = fD_NH_atm[j] * M_NH; nD_1 = fD_NH_atm[j+1] * M_NH1
            nD_SH = fD_SH_atm[j] * M_SH
            dD_src_frac_NH[j] = (nD_1 - nD + nD * aD_NH / tau_NH[j] - (nD_SH - nD) / tau_ex) / S_NH[j]

            # D SH
            nD = fD_SH_atm[j] * M_SH; nD_1 = fD_SH_atm[j+1] * M_SH1
            nD_NH = fD_NH_atm[j] * M_NH
            dD_src_frac_SH[j] = (nD_1 - nD + nD * aD_SH / tau_SH[j] - (nD_NH - nD) / tau_ex) / S_SH[j]

        # Convert source fractions to DELTA space
        d13C_src_NH = fraction_to_delta_d13C(d13C_src_frac_NH)
        d13C_src_SH = fraction_to_delta_d13C(d13C_src_frac_SH)
        dD_src_NH = fraction_to_delta_dD(dD_src_frac_NH)
        dD_src_SH = fraction_to_delta_dD(dD_src_frac_SH)

        # Source signatures (already in delta space)
        sigs = sample_source_signatures_hemi(rng, data, k, n)

        for j in range(n):
            # NH solve in delta space
            BB_NH[j,k], FF_NH[j,k], Mic_NH[j,k] = solve_delta_space(
                S_NH[j], d13C_src_NH[j], dD_src_NH[j],
                sigs['ff_d13C_NH'][j], sigs['mic_d13C_NH'][j], sigs['bb_d13C_NH'][j],
                sigs['ff_dD_NH'][j], sigs['mic_dD_NH'][j], sigs['bb_dD_NH'][j])

            # SH solve in delta space
            BB_SH[j,k], FF_SH[j,k], Mic_SH[j,k] = solve_delta_space(
                S_SH[j], d13C_src_SH[j], dD_src_SH[j],
                sigs['ff_d13C_SH'][j], sigs['mic_d13C_SH'][j], sigs['bb_d13C_SH'][j],
                sigs['ff_dD_SH'][j], sigs['mic_dD_SH'][j], sigs['bb_dD_SH'][j])

    return {
        'years': years, 'n': n,
        'FF_NH': FF_NH, 'FF_SH': FF_SH,
        'Mic_NH': Mic_NH, 'Mic_SH': Mic_SH,
        'BB_NH': BB_NH, 'BB_SH': BB_SH,
        'CH4_NH': CH4_NH, 'CH4_SH': CH4_SH,
    }


def run_one_box_v3(data, cfg):
    """1-box model with delta-space solver."""
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    CH4 = data.CH4_global
    c13_glob = data.c13_global
    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)

    NI = cfg.n_iterations
    FF = np.zeros((n, NI)); Mic = np.zeros((n, NI)); BB = np.zeros((n, NI))

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        kies = sample_KIE(rng, cfg.kie_mode)
        K13, KD = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        a13 = 1.0 / K13; aD = 1.0 / KD

        S = np.zeros(n)
        for i in range(n):
            M = CH4[i] * PT; M1 = CH4[i+1] * PT
            S[i] = (M1 - M) + M / tau_global[i]

        d13C_MC = sample_atm_d13C(data, k, n)
        dD_MC = sample_atm_dD(data, k, n)

        f13_atm = delta_to_fraction_d13C(d13C_MC)
        fD_atm = delta_to_fraction_dD(dD_MC)

        # Source isotopic composition in fraction space
        d13C_src_frac = np.zeros(n); dD_src_frac = np.zeros(n)
        for j in range(n):
            M = CH4[j] * PT; M1 = CH4[j+1] * PT
            n13 = f13_atm[j] * M; n13_1 = f13_atm[j+1] * M1
            d13C_src_frac[j] = (n13_1 - n13 + n13 * a13 / tau_global[j]) / S[j]
            nD = fD_atm[j] * M; nD_1 = fD_atm[j+1] * M1
            dD_src_frac[j] = (nD_1 - nD + nD * aD / tau_global[j]) / S[j]

        # Convert to delta space
        d13C_src = fraction_to_delta_d13C(d13C_src_frac)
        dD_src = fraction_to_delta_dD(dD_src_frac)

        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            BB[j,k], FF[j,k], Mic[j,k] = solve_delta_space(
                S[j], d13C_src[j], dD_src[j],
                sigs['ff_d13C'][j], sigs['mic_d13C'][j], sigs['bb_d13C'][j],
                sigs['ff_dD'][j], sigs['mic_dD'][j], sigs['bb_dD'][j])

    return {'years': years, 'n': n, 'FF': FF, 'Mic': Mic, 'BB': BB}


def compute_trends(arr, years, start=TREND_START, end_trim=1):
    end = years[-1] - end_trim
    mask = (years >= start) & (years <= end)
    yrs = years[mask]
    sub = arr[mask, :]
    slopes = np.zeros(sub.shape[1])
    for k in range(sub.shape[1]):
        col = sub[:, k]
        if np.any(np.isnan(col)):
            slopes[k] = np.nan
            continue
        slopes[k] = sp_stats.linregress(yrs, col).slope
    return slopes


def main():
    print("=" * 70)
    print("v3: DELTA-SPACE SOLVER")
    print("Fixes fundamental numerical issue (rank-1 → rank-3 A matrix)")
    print("=" * 70)

    data = load_data(ROOT, two_box=True)
    cfg = ModelConfig(NI, "sampled", "varying", 9.0, SEED)

    # Quick conditioning check
    print("\nConditioning in delta space:")
    A_delta = np.array([
        [1.0, 1.0, 1.0],
        [-44.0, -61.0, -24.8],
        [-194.0, -321.0, -225.0],
    ])
    scale = np.array([1.0, 1.0/50.0, 1.0/250.0])
    A_scaled = A_delta * scale[:, None]
    U, s, Vt = np.linalg.svd(A_scaled)
    print(f"  Singular values: {s}")
    print(f"  Condition number: {s[0]/s[-1]:.1f}")
    print(f"  Effective rank: {np.sum(s > 0.01*s[0])}")

    print("\n--- 2-Box Model (v3, delta-space) ---")
    res_2box = run_two_box_v3(data, cfg)

    print("\n--- 1-Box Model (v3, delta-space) ---")
    res_1box = run_one_box_v3(data, cfg)

    years = res_2box['years']

    # Save
    np.savez(RESULTS_DIR / "twobox_v3.npz",
             NH_FF=res_2box['FF_NH'], NH_Mic=res_2box['Mic_NH'], NH_BB=res_2box['BB_NH'],
             SH_FF=res_2box['FF_SH'], SH_Mic=res_2box['Mic_SH'], SH_BB=res_2box['BB_SH'],
             years=years)
    np.savez(RESULTS_DIR / "onebox_v3.npz",
             FF=res_1box['FF'], Mic=res_1box['Mic'], BB=res_1box['BB'],
             years=years)

    # Trend analysis
    print(f"\n--- Trend Analysis ({TREND_START}-{int(years[-2])}, last year trimmed) ---")

    FF_global = res_2box['FF_NH'] + res_2box['FF_SH']
    Mic_global = res_2box['Mic_NH'] + res_2box['Mic_SH']
    BB_global = res_2box['BB_NH'] + res_2box['BB_SH']

    for name, arr in [('2box_NH_FF', res_2box['FF_NH']), ('2box_SH_FF', res_2box['FF_SH']),
                      ('2box_Global_FF', FF_global), ('2box_NH_Mic', res_2box['Mic_NH']),
                      ('2box_SH_Mic', res_2box['Mic_SH']), ('2box_Global_Mic', Mic_global),
                      ('2box_NH_BB', res_2box['BB_NH']), ('2box_Global_BB', BB_global),
                      ('1box_FF', res_1box['FF']), ('1box_Mic', res_1box['Mic']),
                      ('1box_BB', res_1box['BB'])]:
        slopes = compute_trends(arr, years)
        med = float(np.nanmedian(slopes))
        pct_pos = float(np.nanmean(slopes > 0) * 100)
        p5 = float(np.nanpercentile(slopes, 5))
        p95 = float(np.nanpercentile(slopes, 95))
        sig = p5 > 0 or p95 < 0
        flag = "✓ SIG" if sig else ""
        print(f"  {name:20s}: {med:+.2f} [{p5:+.2f}, {p95:+.2f}] ({pct_pos:.0f}% pos) {flag}")

    # Mean levels
    j2010 = np.where(years == 2010)[0][0]
    ff_2box = np.nanmedian(FF_global[j2010, :])
    mic_2box = np.nanmedian(Mic_global[j2010, :])
    bb_2box = np.nanmedian(BB_global[j2010, :])
    total = ff_2box + mic_2box + bb_2box
    ff_1box = np.nanmedian(res_1box['FF'][j2010, :])
    
    print(f"\n--- Mean Emission Levels (2010) ---")
    print(f"  2-box FF:  {ff_2box:.0f} Tg/yr ({ff_2box/total*100:.0f}%)")
    print(f"  2-box Mic: {mic_2box:.0f} Tg/yr ({mic_2box/total*100:.0f}%)")
    print(f"  2-box BB:  {bb_2box:.0f} Tg/yr ({bb_2box/total*100:.0f}%)")
    print(f"  Total:     {total:.0f} Tg/yr")
    print(f"  1-box FF:  {ff_1box:.0f} Tg/yr")
    print(f"  EDGAR FF:  ~110 Tg/yr")

    # NH share
    ff_nh = np.nanmedian(res_2box['FF_NH'][j2010, :])
    ff_sh = np.nanmedian(res_2box['FF_SH'][j2010, :])
    print(f"  NH FF share: {ff_nh/(ff_nh+ff_sh)*100:.0f}%  (EDGAR: 72%)")
    
    # NaN count
    nan_pct = np.isnan(FF_global).mean() * 100
    print(f"\n  NaN fraction: {nan_pct:.1f}%")

    # Aliasing bias
    slopes_2box = compute_trends(FF_global, years)
    slopes_1box = compute_trends(res_1box['FF'], years)
    bias = np.nanmedian(slopes_2box) - np.nanmedian(slopes_1box)
    print(f"  Aliasing bias (2box-1box FF): {bias:+.2f} Tg/yr²")

    print(f"\n  Results saved to {RESULTS_DIR}/")
    print(f"  ✅ v3 delta-space model complete")


if __name__ == "__main__":
    main()
