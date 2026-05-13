#!/usr/bin/env python3
"""
Phase 1: Run 2-box and 1-box models with full per-iteration output.

Saves all_iterations.npz and hemispheric_detail.csv for downstream analysis.

Two 2-box configs:
  - global_sigs: original 3x3_two.py behavior (global source signatures)
  - hemi_sigs:   hemispheric source signatures (sample_source_signatures_hemi)

Plus 1-box reference (3x3_one equivalent).
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]  # repo root (upgrade_two_isotope_model/)
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

from common import (
    ModelConfig, QualityMonitor, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
    smooth_5yr, pad_to_length,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
)

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def run_two_box(cfg: ModelConfig, use_hemi_sigs: bool = True,
                use_hemi_dD_atm: bool = True, label: str = "hemi"):
    """Run the 2-box model, save per-iteration arrays."""
    print(f"\n{'='*70}")
    print(f"2-BOX MODEL — {label}")
    print(f"  hemi_sigs={use_hemi_sigs}, hemi_dD_atm={use_hemi_dD_atm}")
    print(f"  KIE={cfg.kie_mode}, τ={cfg.lifetime_mode}, N={cfg.n_iterations}")
    print(f"{'='*70}")

    data = load_data(ROOT, two_box=True)
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global

    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    NI = cfg.n_iterations
    BB_NH = np.zeros((n, NI)); FF_NH = np.zeros((n, NI)); Mic_NH = np.zeros((n, NI))
    BB_SH = np.zeros((n, NI)); FF_SH = np.zeros((n, NI)); Mic_SH = np.zeros((n, NI))

    W_NH = np.diag([100.0, 1.0, 0.5])
    W_SH = np.diag([200.0, 1.0, 0.5])

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        kies = sample_KIE(rng, cfg.kie_mode)
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        a13_NH = 1.0 / K13_NH; aD_NH = 1.0 / KD_NH
        a13_SH = 1.0 / K13_SH; aD_SH = 1.0 / KD_SH

        # Total source per hemisphere
        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i] * PT_HEMI; M_NH1 = CH4_NH[i+1] * PT_HEMI
            M_SH = CH4_SH[i] * PT_HEMI; M_SH1 = CH4_SH[i+1] * PT_HEMI
            ex_NH = (M_SH - M_NH) / tau_ex
            ex_SH = (M_NH - M_SH) / tau_ex
            S_NH[i] = (M_NH1 - M_NH) + M_NH / tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1 - M_SH) + M_SH / tau_SH[i] - ex_SH

        # Atmospheric observations
        d13C_glob_MC = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off

        if use_hemi_dD_atm:
            dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
        else:
            dD_glob_MC = sample_atm_dD(data, k, n)
            dD_NH_MC = dD_glob_MC - DD_IH_OFFSET
            dD_SH_MC = dD_glob_MC + DD_IH_OFFSET

        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
        fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
        fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)

        # Isotopic source fractions
        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
        for j in range(n):
            n13_NH = f13_NH_atm[j] * CH4_NH[j] * PT_HEMI
            n13_NH1 = f13_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            n13_SH = f13_SH_atm[j] * CH4_SH[j] * PT_HEMI
            n13_SH1 = f13_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            ex13_NH = (n13_SH - n13_NH) / tau_ex
            ex13_SH = (n13_NH - n13_SH) / tau_ex
            d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH * a13_NH / tau_NH[j] - ex13_NH) / S_NH[j]
            d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH * a13_SH / tau_SH[j] - ex13_SH) / S_SH[j]

            nD_NH = fD_NH_atm[j] * CH4_NH[j] * PT_HEMI
            nD_NH1 = fD_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            nD_SH = fD_SH_atm[j] * CH4_SH[j] * PT_HEMI
            nD_SH1 = fD_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            exD_NH = (nD_SH - nD_NH) / tau_ex
            exD_SH = (nD_NH - nD_SH) / tau_ex
            dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH * aD_NH / tau_NH[j] - exD_NH) / S_NH[j]
            dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH * aD_SH / tau_SH[j] - exD_SH) / S_SH[j]

        # Source signatures
        if use_hemi_sigs:
            sigs = sample_source_signatures_hemi(rng, data, k, n)
        else:
            sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            # NH: use hemispheric sigs if available
            if use_hemi_sigs:
                f13_bb_nh = delta_to_fraction_d13C(sigs['bb_d13C_NH'][j])
                f13_ff_nh = delta_to_fraction_d13C(sigs['ff_d13C_NH'][j])
                f13_mic_nh = delta_to_fraction_d13C(sigs['mic_d13C_NH'][j])
                fD_bb_nh = delta_to_fraction_dD(sigs['bb_dD_NH'][j])
                fD_ff_nh = delta_to_fraction_dD(sigs['ff_dD_NH'][j])
                fD_mic_nh = delta_to_fraction_dD(sigs['mic_dD_NH'][j])
                f13_bb_sh = delta_to_fraction_d13C(sigs['bb_d13C_SH'][j])
                f13_ff_sh = delta_to_fraction_d13C(sigs['ff_d13C_SH'][j])
                f13_mic_sh = delta_to_fraction_d13C(sigs['mic_d13C_SH'][j])
                fD_bb_sh = delta_to_fraction_dD(sigs['bb_dD_SH'][j])
                fD_ff_sh = delta_to_fraction_dD(sigs['ff_dD_SH'][j])
                fD_mic_sh = delta_to_fraction_dD(sigs['mic_dD_SH'][j])
            else:
                f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'][j])
                f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'][j])
                f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'][j])
                fD_bb = delta_to_fraction_dD(sigs['bb_dD'][j])
                fD_ff = delta_to_fraction_dD(sigs['ff_dD'][j])
                fD_mic = delta_to_fraction_dD(sigs['mic_dD'][j])
                f13_bb_nh = f13_bb_sh = f13_bb
                f13_ff_nh = f13_ff_sh = f13_ff
                f13_mic_nh = f13_mic_sh = f13_mic
                fD_bb_nh = fD_bb_sh = fD_bb
                fD_ff_nh = fD_ff_sh = fD_ff
                fD_mic_nh = fD_mic_sh = fD_mic

            # NH solve
            A_nh = np.array([
                [1.0,         1.0,         1.0],
                [f13_bb_nh,   f13_ff_nh,   f13_mic_nh],
                [fD_bb_nh,    fD_ff_nh,    fD_mic_nh],
            ])
            B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
            ub = S_NH[j] * 1.5
            try:
                res = lsq_linear(W_NH @ A_nh, W_NH @ B_nh, bounds=(0, ub))
                BB_NH[j,k], FF_NH[j,k], Mic_NH[j,k] = res.x
            except Exception:
                BB_NH[j,k] = FF_NH[j,k] = Mic_NH[j,k] = np.nan

            # SH solve
            A_sh = np.array([
                [1.0,         1.0,         1.0],
                [f13_bb_sh,   f13_ff_sh,   f13_mic_sh],
                [fD_bb_sh,    fD_ff_sh,    fD_mic_sh],
            ])
            B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
            ub = S_SH[j] * 1.5
            try:
                res = lsq_linear(W_SH @ A_sh, W_SH @ B_sh, bounds=(0, ub))
                BB_SH[j,k], FF_SH[j,k], Mic_SH[j,k] = res.x
            except Exception:
                BB_SH[j,k] = FF_SH[j,k] = Mic_SH[j,k] = np.nan

    print("  MC complete.")
    out = RESULTS_DIR / f"twobox_{label}"
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "all_iterations.npz",
                        NH_FF=FF_NH, NH_Mic=Mic_NH, NH_BB=BB_NH,
                        SH_FF=FF_SH, SH_Mic=Mic_SH, SH_BB=BB_SH,
                        years=years)

    # Summary CSV
    rows = []
    for i, y in enumerate(years):
        row = {'year': y}
        for src, arr_nh, arr_sh in [('FF', FF_NH, FF_SH), ('Mic', Mic_NH, Mic_SH), ('BB', BB_NH, BB_SH)]:
            for hemi, arr in [('NH', arr_nh), ('SH', arr_sh)]:
                v = arr[i]
                row[f'{hemi}_{src}_median'] = np.nanmedian(v)
                row[f'{hemi}_{src}_p5'] = np.nanpercentile(v, 5)
                row[f'{hemi}_{src}_p95'] = np.nanpercentile(v, 95)
                row[f'{hemi}_{src}_mean'] = np.nanmean(v)
                row[f'{hemi}_{src}_std'] = np.nanstd(v)
            glob = arr_nh[i] + arr_sh[i]
            row[f'Global_{src}_median'] = np.nanmedian(glob)
            row[f'Global_{src}_p5'] = np.nanpercentile(glob, 5)
            row[f'Global_{src}_p95'] = np.nanpercentile(glob, 95)
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "hemispheric_detail.csv", index=False)
    print(f"  Saved to {out}/")
    return years, FF_NH, FF_SH, Mic_NH, Mic_SH, BB_NH, BB_SH


def run_one_box(cfg: ModelConfig, label: str = "onebox"):
    """Run the 1-box 3x3 model (Riddell-Young equivalent)."""
    print(f"\n{'='*70}")
    print(f"1-BOX MODEL — {label}")
    print(f"  KIE={cfg.kie_mode}, τ={cfg.lifetime_mode}, N={cfg.n_iterations}")
    print(f"{'='*70}")

    data = load_data(ROOT, two_box=False)
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    CH4 = data.CH4_global
    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)

    NI = cfg.n_iterations
    BB = np.zeros((n, NI)); FF = np.zeros((n, NI)); Mic = np.zeros((n, NI))

    # Global sink fractions (average of NH/SH)
    sf_global = {k: (SINK_FRACTIONS_NH[k] + SINK_FRACTIONS_SH[k]) / 2
                 for k in SINK_FRACTIONS_NH}

    W = np.diag([100.0, 1.0, 0.5])

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        kies = sample_KIE(rng, cfg.kie_mode)
        K13, KD = compute_bulk_KIE(kies, sf_global)
        a13 = 1.0 / K13; aD = 1.0 / KD

        S = np.zeros(n)
        for i in range(n):
            M = CH4[i] * PT; M1 = CH4[i+1] * PT
            S[i] = (M1 - M) + M / tau_global[i]

        d13C_MC = sample_atm_d13C(data, k, n)
        dD_MC = sample_atm_dD(data, k, n)
        f13_atm = delta_to_fraction_d13C(d13C_MC)
        fD_atm = delta_to_fraction_dD(dD_MC)

        d13C_src = np.zeros(n); dD_src = np.zeros(n)
        for j in range(n):
            n13 = f13_atm[j] * CH4[j] * PT
            n13_1 = f13_atm[j+1] * CH4[j+1] * PT
            d13C_src[j] = (n13_1 - n13 + n13 * a13 / tau_global[j]) / S[j]

            nD = fD_atm[j] * CH4[j] * PT
            nD_1 = fD_atm[j+1] * CH4[j+1] * PT
            dD_src[j] = (nD_1 - nD + nD * aD / tau_global[j]) / S[j]

        sigs = sample_source_signatures(rng, data, k, n)
        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
        fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
        fD_mic = delta_to_fraction_dD(sigs['mic_dD'])

        for j in range(n):
            A = np.array([
                [1.0,       1.0,       1.0],
                [f13_bb[j], f13_ff[j], f13_mic[j]],
                [fD_bb[j],  fD_ff[j],  fD_mic[j]],
            ])
            B = np.array([S[j], S[j]*d13C_src[j], S[j]*dD_src[j]])
            ub = S[j] * 1.5
            try:
                res = lsq_linear(W @ A, W @ B, bounds=(0, ub))
                BB[j,k], FF[j,k], Mic[j,k] = res.x
            except Exception:
                BB[j,k] = FF[j,k] = Mic[j,k] = np.nan

    print("  MC complete.")
    out = RESULTS_DIR / f"onebox_{label}"
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "all_iterations.npz",
                        FF=FF, Mic=Mic, BB=BB, years=years)

    rows = []
    for i, y in enumerate(years):
        row = {'year': y}
        for src, arr in [('FF', FF), ('Mic', Mic), ('BB', BB)]:
            row[f'{src}_median'] = np.nanmedian(arr[i])
            row[f'{src}_p5'] = np.nanpercentile(arr[i], 5)
            row[f'{src}_p95'] = np.nanpercentile(arr[i], 95)
            row[f'{src}_mean'] = np.nanmean(arr[i])
            row[f'{src}_std'] = np.nanstd(arr[i])
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "global_detail.csv", index=False)
    print(f"  Saved to {out}/")
    return years, FF, Mic, BB


if __name__ == "__main__":
    cfg = ModelConfig(n_iterations=1000, kie_mode="sampled",
                      lifetime_mode="varying", seed=42)

    # 2-box with hemispheric source signatures (primary)
    run_two_box(cfg, use_hemi_sigs=True, use_hemi_dD_atm=True, label="hemi")

    # 2-box with global source signatures (baseline comparison)
    cfg2 = ModelConfig(n_iterations=1000, kie_mode="sampled",
                       lifetime_mode="varying", seed=42)
    run_two_box(cfg2, use_hemi_sigs=False, use_hemi_dD_atm=False, label="global_sigs")

    # 1-box reference
    cfg3 = ModelConfig(n_iterations=1000, kie_mode="sampled",
                       lifetime_mode="varying", seed=42)
    run_one_box(cfg3, label="reference")

    print("\n✅ All Phase 1 model runs complete.")
