#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core.py — Shared 2-box model runner for the δD threshold experiment.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    smooth_5yr,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    TAU_EX_MEAN, TAU_EX_STD,
    BB_NH_FRACTION, BB_SH_FRACTION,
    sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD_hemi,
)

# Constants
W_NH = np.diag([100.0, 1.0, 0.5])
W_SH = np.diag([200.0, 1.0, 0.5])
DEFAULT_N_ITER = 1000
DEFAULT_SEED = 42


def inflate_dD_uncertainty(sigs, rng, multiplier):
    """Inflate δD signatures by a multiplier around their mean."""
    if multiplier <= 1.0:
        return sigs
    extra_mic = rng.normal() * 8.25 * (multiplier - 1)
    extra_ff = rng.normal() * 0.70 * (multiplier - 1)
    extra_bb = rng.normal() * 7.09 * (multiplier - 1)
    sigs_new = dict(sigs)
    sigs_new['mic_dD'] = sigs['mic_dD'] + extra_mic
    sigs_new['ff_dD'] = sigs['ff_dD'] + extra_ff
    sigs_new['bb_dD'] = sigs['bb_dD'] + extra_bb
    for hemi in ('NH', 'SH'):
        for src, extra in [('mic_dD', extra_mic), ('ff_dD', extra_ff), ('bb_dD', extra_bb)]:
            key = f'{src}_{hemi}'
            if key in sigs:
                sigs_new[key] = sigs[key] + extra
    return sigs_new


def ci_width(arr, start_idx=8):
    """90% CI width over years starting at start_idx."""
    s = smooth_5yr(arr)
    ci5 = np.nanpercentile(s[start_idx:], 5, axis=1)
    ci95 = np.nanpercentile(s[start_idx:], 95, axis=1)
    return float(np.nanmean(ci95 - ci5))


def ci_width_hemi(arr, start_idx=8):
    """90% CI width for a single hemisphere's FF array."""
    s = smooth_5yr(arr)
    ci5 = np.nanpercentile(s[start_idx:], 5, axis=1)
    ci95 = np.nanpercentile(s[start_idx:], 95, axis=1)
    return float(np.nanmean(ci95 - ci5))


def run_twobox(data, multiplier, n_iter, seed, mode="dual",
               kie_mode="sampled", lifetime_mode="varying", tau_fixed=9.0,
               year_start_idx=0, track_bounds=False):
    """
    Unified 2-box (NH/SH) model runner.

    Parameters
    ----------
    data : loaded data object (from load_data with two_box=True)
    multiplier : float, δD uncertainty inflation factor
    n_iter : int, number of MC iterations
    seed : int, random seed
    mode : 'dual' or 'd13C_only'
    kie_mode : 'sampled', 'saueressig', or 'cantrell'
    lifetime_mode : 'varying' or 'fixed'
    tau_fixed : float, fixed lifetime (used when lifetime_mode='fixed')
    year_start_idx : int, start index for year clipping (0 = no clip)
    track_bounds : bool, whether to track bound hits

    Returns
    -------
    dict with keys: FF_NH, FF_SH, FF_G, Mic_NH, Mic_SH, Mic_G,
                    and optionally BB_NH, BB_SH, BB_G, bound_hits
    """
    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years

    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global

    tau_global = compute_lifetime(years, lifetime_mode, tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    BB_hemi_NH = data.BB_global_mean * BB_NH_FRACTION
    BB_hemi_SH = data.BB_global_mean * BB_SH_FRACTION

    FF_NH = np.zeros((n, n_iter))
    FF_SH = np.zeros((n, n_iter))
    Mic_NH = np.zeros((n, n_iter))
    Mic_SH = np.zeros((n, n_iter))
    BB_NH_comp = np.zeros((n, n_iter)) if mode == "dual" else None
    BB_SH_comp = np.zeros((n, n_iter)) if mode == "dual" else None
    bound_hits = np.zeros((n, n_iter), dtype=int) if track_bounds else None

    for k in range(n_iter):
        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        kies = sample_KIE(rng, kie_mode)
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        a13_NH = 1.0 / K13_NH; aD_NH = 1.0 / KD_NH
        a13_SH = 1.0 / K13_SH; aD_SH = 1.0 / KD_SH

        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i] * PT_HEMI; M_NH1 = CH4_NH[i+1] * PT_HEMI
            M_SH = CH4_SH[i] * PT_HEMI; M_SH1 = CH4_SH[i+1] * PT_HEMI
            ex_NH = (M_SH - M_NH) / tau_ex
            ex_SH = (M_NH - M_SH) / tau_ex
            S_NH[i] = (M_NH1 - M_NH) + M_NH / tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1 - M_SH) + M_SH / tau_SH[i] - ex_SH

        d13C_glob_mc = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_mc[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off
        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)

        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        for j in range(n):
            n13_NH = f13_NH_atm[j] * CH4_NH[j] * PT_HEMI
            n13_NH1 = f13_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            n13_SH = f13_SH_atm[j] * CH4_SH[j] * PT_HEMI
            n13_SH1 = f13_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            ex13_NH = (n13_SH - n13_NH) / tau_ex
            ex13_SH = (n13_NH - n13_SH) / tau_ex
            d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH * a13_NH / tau_NH[j] - ex13_NH) / S_NH[j]
            d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH * a13_SH / tau_SH[j] - ex13_SH) / S_SH[j]

        sigs = sample_source_signatures_hemi(rng, data, k, n)
        if multiplier > 1.0 and mode == "dual":
            sigs = inflate_dD_uncertainty(sigs, rng, multiplier)

        f13_bb_NH = delta_to_fraction_d13C(sigs['bb_d13C_NH'])
        f13_ff_NH = delta_to_fraction_d13C(sigs['ff_d13C_NH'])
        f13_mic_NH = delta_to_fraction_d13C(sigs['mic_d13C_NH'])
        f13_bb_SH = delta_to_fraction_d13C(sigs['bb_d13C_SH'])
        f13_ff_SH = delta_to_fraction_d13C(sigs['ff_d13C_SH'])
        f13_mic_SH = delta_to_fraction_d13C(sigs['mic_d13C_SH'])

        if mode == "dual":
            dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
            fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
            fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)

            dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
            for j in range(n):
                nD_NH = fD_NH_atm[j] * CH4_NH[j] * PT_HEMI
                nD_NH1 = fD_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
                nD_SH = fD_SH_atm[j] * CH4_SH[j] * PT_HEMI
                nD_SH1 = fD_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
                exD_NH = (nD_SH - nD_NH) / tau_ex
                exD_SH = (nD_NH - nD_SH) / tau_ex
                dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH * aD_NH / tau_NH[j] - exD_NH) / S_NH[j]
                dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH * aD_SH / tau_SH[j] - exD_SH) / S_SH[j]

            fD_bb_NH = delta_to_fraction_dD(sigs['bb_dD_NH'])
            fD_ff_NH = delta_to_fraction_dD(sigs['ff_dD_NH'])
            fD_mic_NH = delta_to_fraction_dD(sigs['mic_dD_NH'])
            fD_bb_SH = delta_to_fraction_dD(sigs['bb_dD_SH'])
            fD_ff_SH = delta_to_fraction_dD(sigs['ff_dD_SH'])
            fD_mic_SH = delta_to_fraction_dD(sigs['mic_dD_SH'])

            for j in range(n):
                A_nh = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb_NH[j], f13_ff_NH[j], f13_mic_NH[j]],
                    [fD_bb_NH[j], fD_ff_NH[j], fD_mic_NH[j]],
                ])
                B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
                try:
                    res = lsq_linear(W_NH @ A_nh, W_NH @ B_nh, bounds=(0, S_NH[j]*1.5))
                    x = res.x
                except:
                    x = np.array([np.nan, np.nan, np.nan])
                if BB_NH_comp is not None:
                    BB_NH_comp[j, k] = x[0]
                FF_NH[j, k] = x[1]; Mic_NH[j, k] = x[2]
                if track_bounds:
                    if np.any(np.isclose(x, 0)) or np.any(np.isclose(x, S_NH[j]*1.5)):
                        bound_hits[j, k] = 1

                A_sh = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb_SH[j], f13_ff_SH[j], f13_mic_SH[j]],
                    [fD_bb_SH[j], fD_ff_SH[j], fD_mic_SH[j]],
                ])
                B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
                try:
                    res = lsq_linear(W_SH @ A_sh, W_SH @ B_sh, bounds=(0, S_SH[j]*1.5))
                    x = res.x
                except:
                    x = np.array([np.nan, np.nan, np.nan])
                if BB_SH_comp is not None:
                    BB_SH_comp[j, k] = x[0]
                FF_SH[j, k] = x[1]; Mic_SH[j, k] = x[2]
                if track_bounds:
                    if np.any(np.isclose(x, 0)) or np.any(np.isclose(x, S_SH[j]*1.5)):
                        bound_hits[j, k] = 1

        else:  # d13C_only
            for j in range(n):
                # NH
                denom_nh = f13_ff_NH[j] - f13_mic_NH[j]
                if abs(denom_nh) < 1e-15:
                    FF_NH[j, k] = np.nan; Mic_NH[j, k] = np.nan
                else:
                    S_rem = S_NH[j] - BB_hemi_NH
                    rhs = S_NH[j]*d13C_src_NH[j] - BB_hemi_NH*f13_bb_NH[j]
                    FF_NH[j, k] = (rhs - S_rem*f13_mic_NH[j]) / denom_nh
                    Mic_NH[j, k] = S_rem - FF_NH[j, k]
                # SH
                denom_sh = f13_ff_SH[j] - f13_mic_SH[j]
                if abs(denom_sh) < 1e-15:
                    FF_SH[j, k] = np.nan; Mic_SH[j, k] = np.nan
                else:
                    S_rem = S_SH[j] - BB_hemi_SH
                    rhs = S_SH[j]*d13C_src_SH[j] - BB_hemi_SH*f13_bb_SH[j]
                    FF_SH[j, k] = (rhs - S_rem*f13_mic_SH[j]) / denom_sh
                    Mic_SH[j, k] = S_rem - FF_SH[j, k]

    result = {
        'FF_NH': FF_NH, 'FF_SH': FF_SH, 'FF_G': FF_NH + FF_SH,
        'Mic_NH': Mic_NH, 'Mic_SH': Mic_SH, 'Mic_G': Mic_NH + Mic_SH,
    }
    if BB_NH_comp is not None:
        result['BB_NH'] = BB_NH_comp
        result['BB_SH'] = BB_SH_comp
        result['BB_G'] = BB_NH_comp + BB_SH_comp
    if bound_hits is not None:
        result['bound_hits'] = bound_hits
    return result
