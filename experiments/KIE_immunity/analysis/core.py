#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core.py — Shared 2-box engine for KIE_immunity phases 5–13.
=============================================================
v4: Addresses Manuscript_Review_V1.0.md issues:
  - A2: W matrix is now a parameter (default unchanged for reproducibility)
  - A3/B6: fix_kie uses KIE_FIXED from common.py (single source of truth)
  - B4: Added compute_trend_regression() for linear trend + p-value
  - B7: Added solver failure/bound-hit tracking
  - C4: Exception handling improved (no silent pass)
  - BB scaling parameter for sensitivity tests (B2)
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import linregress

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
    smooth_5yr, KIE_FIXED,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
    BB_NH_FRACTION, BB_SH_FRACTION,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Default W matrix — kept for backward compatibility
W_DEFAULT = np.diag([100.0, 1.0, 0.5])


def _make_sink_fractions(cl_frac=None):
    """Return (NH, SH) sink fraction dicts, optionally overriding Cl fraction."""
    if cl_frac is None:
        return SINK_FRACTIONS_NH, SINK_FRACTIONS_SH

    def _adjust(base, new_cl):
        old_cl = base['Cl']
        delta = new_cl - old_cl
        return {
            'OH': base['OH'] - delta,
            'Cl': new_cl,
            'Strat': base['Strat'],
            'Soil': base['Soil'],
        }

    return _adjust(SINK_FRACTIONS_NH, cl_frac), _adjust(SINK_FRACTIONS_SH, cl_frac)


def run_2box_flex(data, n_iter=400, seed=42, *,
                  tau_mode="varying", tau_fixed=9.0,
                  oh_d_fixed=None,
                  cl_frac=None,
                  tau_ex_fixed=None,
                  fix_kie=False, fix_sigs=False,
                  W=None,
                  bb_scale=1.0,
                  track_diagnostics=False):
    """
    Flexible 2-box dual real-hemi-dD model.

    Parameters
    ----------
    W : ndarray (3,3) or None
        Weight matrix for the 3×3 least-squares solver.
        Default: diag(100, 1, 0.5).
    bb_scale : float
        Multiplicative scaling factor for prescribed BB emissions.
        Default: 1.0 (no perturbation).
    track_diagnostics : bool
        If True, return (FF_G, diagnostics_dict) instead of just FF_G.

    Returns
    -------
    FF_G : ndarray (n_years, n_iter)
    diagnostics : dict (only if track_diagnostics=True)
        Keys: 'solver_failures', 'bound_hits', 'total_solves'
    """
    if W is None:
        W = W_DEFAULT

    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years
    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global

    tau_global = compute_lifetime(years, tau_mode, tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    BB_NH = data.BB_global_mean * BB_NH_FRACTION * bb_scale
    BB_SH = data.BB_global_mean * BB_SH_FRACTION * bb_scale

    sf_NH, sf_SH = _make_sink_fractions(cl_frac)

    FF_G = np.zeros((n, n_iter))

    # Diagnostics
    solver_failures = 0
    bound_hits = 0
    total_solves = 0

    if fix_kie:
        kie_base = dict(KIE_FIXED)
        if oh_d_fixed is not None:
            kie_base['OH_D'] = oh_d_fixed
        kie_fixed = kie_base
    if fix_sigs:
        rng_tmp = np.random.default_rng(0)
        sigs_fixed = sample_source_signatures_hemi(rng_tmp, data, 0, n)

    for k in range(n_iter):
        if tau_ex_fixed is not None:
            tau_ex = tau_ex_fixed
        else:
            tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))

        if fix_kie:
            kies = dict(kie_fixed)
        else:
            kies = sample_KIE(rng, "sampled")
            if oh_d_fixed is not None:
                kies['OH_D'] = oh_d_fixed

        K13_NH, KD_NH = compute_bulk_KIE(kies, sf_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, sf_SH)
        a13_NH, aD_NH = 1.0/K13_NH, 1.0/KD_NH
        a13_SH, aD_SH = 1.0/K13_SH, 1.0/KD_SH

        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i]*PT_HEMI; M_NH1 = CH4_NH[i+1]*PT_HEMI
            M_SH = CH4_SH[i]*PT_HEMI; M_SH1 = CH4_SH[i+1]*PT_HEMI
            ex_NH = (M_SH-M_NH)/tau_ex; ex_SH = (M_NH-M_SH)/tau_ex
            S_NH[i] = (M_NH1-M_NH) + M_NH/tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1-M_SH) + M_SH/tau_SH[i] - ex_SH

        d13C_glob = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n+1)
        d13C_off = d13C_glob[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off
        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)

        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        for j in range(n):
            n13_NH = f13_NH_atm[j]*CH4_NH[j]*PT_HEMI
            n13_NH1 = f13_NH_atm[j+1]*CH4_NH[j+1]*PT_HEMI
            n13_SH = f13_SH_atm[j]*CH4_SH[j]*PT_HEMI
            n13_SH1 = f13_SH_atm[j+1]*CH4_SH[j+1]*PT_HEMI
            ex13_NH = (n13_SH-n13_NH)/tau_ex
            ex13_SH = (n13_NH-n13_SH)/tau_ex
            d13C_src_NH[j] = (n13_NH1-n13_NH + n13_NH*a13_NH/tau_NH[j] - ex13_NH)/S_NH[j]
            d13C_src_SH[j] = (n13_SH1-n13_SH + n13_SH*a13_SH/tau_SH[j] - ex13_SH)/S_SH[j]

        if fix_sigs:
            sigs = sigs_fixed
        else:
            sigs = sample_source_signatures_hemi(rng, data, k, n)

        # Hemispheric δ¹³C source signatures
        f13_bb_NH  = delta_to_fraction_d13C(sigs['bb_d13C_NH'])
        f13_ff_NH  = delta_to_fraction_d13C(sigs['ff_d13C_NH'])
        f13_mic_NH = delta_to_fraction_d13C(sigs['mic_d13C_NH'])
        f13_bb_SH  = delta_to_fraction_d13C(sigs['bb_d13C_SH'])
        f13_ff_SH  = delta_to_fraction_d13C(sigs['ff_d13C_SH'])
        f13_mic_SH = delta_to_fraction_d13C(sigs['mic_d13C_SH'])

        dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
        fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
        fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)

        dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
        for j in range(n):
            nD_NH = fD_NH_atm[j]*CH4_NH[j]*PT_HEMI
            nD_NH1 = fD_NH_atm[j+1]*CH4_NH[j+1]*PT_HEMI
            nD_SH = fD_SH_atm[j]*CH4_SH[j]*PT_HEMI
            nD_SH1 = fD_SH_atm[j+1]*CH4_SH[j+1]*PT_HEMI
            exD_NH = (nD_SH-nD_NH)/tau_ex
            exD_SH = (nD_NH-nD_SH)/tau_ex
            dD_src_NH[j] = (nD_NH1-nD_NH + nD_NH*aD_NH/tau_NH[j] - exD_NH)/S_NH[j]
            dD_src_SH[j] = (nD_SH1-nD_SH + nD_SH*aD_SH/tau_SH[j] - exD_SH)/S_SH[j]

        fD_bb_NH = delta_to_fraction_dD(sigs['bb_dD_NH'])
        fD_ff_NH = delta_to_fraction_dD(sigs['ff_dD_NH'])
        fD_mic_NH = delta_to_fraction_dD(sigs['mic_dD_NH'])
        fD_bb_SH = delta_to_fraction_dD(sigs['bb_dD_SH'])
        fD_ff_SH = delta_to_fraction_dD(sigs['ff_dD_SH'])
        fD_mic_SH = delta_to_fraction_dD(sigs['mic_dD_SH'])

        for j in range(n):
            for S, d13C_src, dD_src, f13_bb_h, f13_ff_h, f13_mic_h, fD_bb_h, fD_ff_h, fD_mic_h in [
                (S_NH[j], d13C_src_NH[j], dD_src_NH[j],
                 f13_bb_NH[j], f13_ff_NH[j], f13_mic_NH[j],
                 fD_bb_NH[j], fD_ff_NH[j], fD_mic_NH[j]),
                (S_SH[j], d13C_src_SH[j], dD_src_SH[j],
                 f13_bb_SH[j], f13_ff_SH[j], f13_mic_SH[j],
                 fD_bb_SH[j], fD_ff_SH[j], fD_mic_SH[j]),
            ]:
                total_solves += 1
                A = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb_h, f13_ff_h, f13_mic_h],
                    [fD_bb_h, fD_ff_h, fD_mic_h],
                ])
                B = np.array([S, S*d13C_src, S*dD_src])
                try:
                    res = lsq_linear(W@A, W@B, bounds=(0, S*1.5))
                    FF_G[j, k] += res.x[1]
                    # Track bound hits
                    if any(res.x <= 1e-6) or any(res.x >= S*1.5 - 1e-6):
                        bound_hits += 1
                except Exception:
                    solver_failures += 1

    if track_diagnostics:
        diag = {
            'solver_failures': solver_failures,
            'bound_hits': bound_hits,
            'total_solves': total_solves,
            'failure_rate_pct': solver_failures / max(1, total_solves) * 100,
            'bound_hit_rate_pct': bound_hits / max(1, total_solves) * 100,
        }
        return FF_G, diag

    return FF_G


def compute_trend(FF, years):
    """Post-2007 ΔFF: mean(2010–2018) − mean(2000–2006), per MC iteration.
    
    Note: 2007–2009 excluded as transition years.
    """
    FF_s = smooth_5yr(FF)
    yrs = np.array(years)
    pre = np.where((yrs >= 2000) & (yrs <= 2006))[0]
    post = np.where((yrs >= 2010) & (yrs <= 2018))[0]
    return np.nanmean(FF_s[post], axis=0) - np.nanmean(FF_s[pre], axis=0)


def compute_trend_regression(FF, years):
    """Linear regression slope of FF over 2000–2020, per MC iteration.
    
    Returns dict with slope (Tg/yr²), p-value, and CI for each iteration's
    median statistics.
    """
    FF_s = smooth_5yr(FF)
    yrs = np.array(years)
    mask = (yrs >= 2000) & (yrs <= 2020)
    idx = np.where(mask)[0]
    x = yrs[mask]
    
    n_iter = FF_s.shape[1]
    slopes = np.zeros(n_iter)
    pvalues = np.zeros(n_iter)
    
    for k in range(n_iter):
        y = FF_s[idx, k]
        valid = ~np.isnan(y)
        if valid.sum() < 5:
            slopes[k] = np.nan
            pvalues[k] = np.nan
            continue
        res = linregress(x[valid], y[valid])
        slopes[k] = res.slope
        pvalues[k] = res.pvalue
    
    return {
        'slope_median': float(np.nanmedian(slopes)),
        'slope_5pct': float(np.nanpercentile(slopes, 5)),
        'slope_95pct': float(np.nanpercentile(slopes, 95)),
        'pvalue_median': float(np.nanmedian(pvalues)),
        'pct_significant': float(np.nanmean(pvalues < 0.05) * 100),
    }


def trend_stats(FF, years):
    """Return (median, lo5, hi95) of ΔFF trend (step-change metric)."""
    deltas = compute_trend(FF, years)
    return float(np.median(deltas)), float(np.percentile(deltas, 5)), float(np.percentile(deltas, 95))


def sigma_ff(FF):
    """Post-2007 mean σ(FF)."""
    FF_s = smooth_5yr(FF)
    return float(np.sqrt(np.nanmean(np.nanvar(FF_s[8:], axis=1))))
