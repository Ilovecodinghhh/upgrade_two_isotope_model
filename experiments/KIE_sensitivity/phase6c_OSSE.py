#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 6c: Observing System Simulation Experiment (OSSE)
========================================================
Generate synthetic "truth" emissions and test whether the agreement
framework can recover them more accurately than δ¹³C-only.

Design:
1. Define a "true" source partition (FF_true, Mic_true, BB_true)
2. Forward-model the atmospheric δ¹³C and δD using TRUE KIE values
3. Add realistic observational noise
4. Invert using WRONG KIE values (both Saueressig and Cantrell)
5. Compare recovery accuracy: δ¹³C-only vs agreement-filtered
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from common import (
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    sample_source_signatures, smooth_5yr, trend_change,
    SINK_FRACTIONS_GLOBAL, PT, C13_STD, D_STD,
)

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "phase6c_OSSE"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 1000
SEED = 42
N_YEARS = 23
YEARS = np.arange(1999, 2022)

# === SYNTHETIC TRUTH ===
# Based on He 2026 Science/JGR central estimates
TRUE_KIE_OH_13C = 1.0046  # "True" value between Saueressig and Cantrell
TRUE_KIE_OH_D = 1.310
TRUE_KIE = {
    'OH_13C': TRUE_KIE_OH_13C, 'OH_D': TRUE_KIE_OH_D,
    'Cl_13C': 1.066, 'Cl_D': 1.52,
    'Strat_13C': 1.003, 'Strat_D': 1.179,
    'Soil_13C': 1.0201, 'Soil_D': 1.083,
}

# True source signatures
TRUE_FF_D13C = -44.0    # Schwietzke 2016
TRUE_MIC_D13C = -62.0   # Wetlands + ruminants mean
TRUE_BB_D13C = -22.0    # C3/C4 mix
TRUE_FF_DD = -180.0     # Sherwood 2017
TRUE_MIC_DD = -310.0    # Whiticar 1999
TRUE_BB_DD = -200.0     # Thanwerdas 2024

# True emission trends — must match mass-balance-implied totals from observed CH₄
# We'll compute TRUE_TOTAL from observations, then partition using desired fractions
# Fractions based on He 2026 Science:
FF_FRACTION = 0.24    # ~24% fossil fuel (130–140 Tg from ~575 total)
BB_FRACTION = 0.05    # ~5% biomass burning (25–30 Tg)
MIC_FRACTION = 0.71   # ~71% microbial (rest)

# These get filled in main() from actual mass balance
TRUE_FF = None
TRUE_MIC = None
TRUE_BB = None
TRUE_TOTAL = None


def forward_model_atmosphere(tau_mode='varying'):
    """Generate synthetic atmospheric observations from true emissions.
    
    The inversion solves:
      S × f13_src = n13[j+1] - n13[j] + n13[j] × α/τ
    
    Rearranging for forward model:
      n13[j+1] = n13[j] + S × f13_src - n13[j] × α/τ
    
    BUT this is wrong! The sink removes mass from the atmosphere:
      n13[j+1] = n13[j] - n13[j]×α/τ + S×f13_src
    
    Actually the full mass+isotope budget is:
      M[j+1] = M[j] + S - M[j]/τ            (mass balance)
      n13[j+1] = n13[j] - n13[j]×α/τ + S×f13_src  (isotope balance)
    
    The factor α/τ represents fractionating removal. Since α < 1 (reciprocal of KIE > 1),
    the sink preferentially removes ¹²C, enriching the atmosphere in ¹³C.
    
    For the inversion to be self-consistent, we use the REAL observed CH4
    and just compute what δ¹³C and δD WOULD be if the true sources emitted.
    """
    tau = compute_lifetime(YEARS, tau_mode)
    KIE_13C_true, KIE_D_true = compute_bulk_KIE(TRUE_KIE, SINK_FRACTIONS_GLOBAL)
    # In the inversion, alpha = 1/KIE (the fractionation applied to the heavy isotope)
    # The sink term is: n13_removed = n13 × (1/KIE) / tau = n13 × alpha / tau
    alpha_13C = 1.0 / KIE_13C_true
    alpha_D = 1.0 / KIE_D_true

    # Use real CH4 observations
    data = load_data(REPO_ROOT, two_box=False)
    CH4 = data.CH4_global.copy()

    # Source-weighted isotopic compositions (delta values)
    d13C_src_delta = (TRUE_FF * TRUE_FF_D13C + TRUE_MIC * TRUE_MIC_D13C +
                      TRUE_BB * TRUE_BB_D13C) / TRUE_TOTAL
    dD_src_delta = (TRUE_FF * TRUE_FF_DD + TRUE_MIC * TRUE_MIC_DD +
                    TRUE_BB * TRUE_BB_DD) / TRUE_TOTAL

    # Instead of forward-modeling (which accumulates errors), use the 
    # INVERSE relationship directly:
    # The inversion computes d13C_src from:
    #   f13_src = (n13[j+1] - n13[j] + n13[j] × α/τ) / S
    # 
    # We want to find d13C_atm such that when we invert it with the TRUE KIE,
    # we recover the TRUE sources. So we just use the observed d13C_atm
    # as our "synthetic truth" — it's already physically consistent.
    #
    # Alternative: directly compute what d13C_atm MUST be given our true sources.
    # From the inversion: S × f13_src = n13[j+1] - n13[j] + n13[j] × α/τ
    # => n13[j+1] = S × f13_src - n13[j] × α/τ + n13[j]
    #            = S × f13_src + n13[j] × (1 - α/τ)
    
    d13C_atm = np.zeros(N_YEARS + 1)
    dD_atm = np.zeros(N_YEARS + 1)
    d13C_atm[0] = data.c13_global[0]  # Start from real observation
    dD_atm[0] = data.dD_global[0] if len(data.dD_global) > 0 else -86.0

    for j in range(N_YEARS):
        S = TRUE_TOTAL[j]
        
        # δ¹³C forward model
        f13_j = delta_to_fraction_d13C(np.array([d13C_atm[j]]))[0]
        n13_j = f13_j * CH4[j] * PT
        f13_src_j = delta_to_fraction_d13C(np.array([d13C_src_delta[j]]))[0]
        
        # n13[j+1] = n13[j] × (1 - α/τ) + S × f13_src
        n13_j1 = n13_j * (1.0 - alpha_13C / tau[j]) + S * f13_src_j
        f13_j1 = n13_j1 / (CH4[j+1] * PT)
        d13C_atm[j+1] = fraction_to_delta_d13C(f13_j1)

        # δD forward model
        fD_j = delta_to_fraction_dD(np.array([dD_atm[j]]))[0]
        nD_j = fD_j * CH4[j] * PT
        fD_src_j = delta_to_fraction_dD(np.array([dD_src_delta[j]]))[0]
        
        nD_j1 = nD_j * (1.0 - alpha_D / tau[j]) + S * fD_src_j
        fD_j1 = nD_j1 / (CH4[j+1] * PT)
        dD_atm[j+1] = fraction_to_delta_dD(fD_j1)

    return CH4, d13C_atm, dD_atm


def invert_synthetic(CH4, d13C_atm_true, dD_atm_true, oh13c_mode,
                     d13C_noise=0.05, dD_noise=3.0):
    """Invert synthetic observations using the REAL MC source-signature data.
    
    This ensures the δD spread is realistic (same as in Phases 1–6).
    The only synthetic element is the atmospheric observations.
    """
    data = load_data(REPO_ROOT, two_box=False)
    rng = np.random.default_rng(SEED)
    tau = compute_lifetime(YEARS, 'varying')
    n = N_YEARS

    SumSource = np.zeros(N_YEARS)
    for i in range(N_YEARS):
        SumSource[i] = CH4[i+1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]

    BB = data.BB_annual.copy()

    FF_c13 = np.full((N_YEARS, N_ITER), np.nan)
    FF_dD = np.full((N_YEARS, N_ITER), np.nan)
    Mic_c13 = np.full((N_YEARS, N_ITER), np.nan)

    for k in range(N_ITER):
        kies = sample_KIE(rng, 'sampled')
        if oh13c_mode == 'saueressig':
            kies['OH_13C'] = 1.0039
        elif oh13c_mode == 'cantrell':
            kies['OH_13C'] = 1.0054
        elif oh13c_mode == 'true':
            kies['OH_13C'] = TRUE_KIE_OH_13C

        KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        alpha_13C = 1.0 / KIE_13C
        alpha_D = 1.0 / KIE_D

        # Add noise to synthetic atmospheric obs
        d13C_obs = d13C_atm_true + rng.normal(0, d13C_noise, N_YEARS + 1)
        dD_obs = dD_atm_true + rng.normal(0, dD_noise, N_YEARS + 1)

        f13 = delta_to_fraction_d13C(d13C_obs)
        fD = delta_to_fraction_dD(dD_obs)
        n13 = f13 * CH4[:N_YEARS+1] * PT
        nD = fD * CH4[:N_YEARS+1] * PT

        # Use REAL MC source signatures (same as Phase 6)
        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(N_YEARS):
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            # δ¹³C inversion
            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)
            denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            if abs(denom_c) > 0.1:
                ff_val = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j) -
                          sigs['bb_d13C'][j] * BB_j) / denom_c
                FF_c13[j, k] = ff_val
                Mic_c13[j, k] = S - BB_j - ff_val

            # δD inversion
            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            dD_src = fraction_to_delta_dD(dD_src_f)
            denom_d = sigs['ff_dD'][j] - sigs['mic_dD'][j]
            if abs(denom_d) > 0.5:
                FF_dD[j, k] = (S * dD_src - sigs['mic_dD'][j] * (S - BB_j) -
                               sigs['bb_dD'][j] * BB_j) / denom_d

    return FF_c13, FF_dD, Mic_c13


def compute_recovery_error(FF_est, FF_true_arr, years, threshold=100,
                           FF_dD=None, Mic_est=None, Mic_true_arr=None):
    """Compute RMSE and bias of estimated vs true emissions."""
    # Unfiltered
    FF_med = np.nanmedian(FF_est, axis=1)
    bias_unfilt = np.nanmean(FF_med - FF_true_arr)
    rmse_unfilt = np.sqrt(np.nanmean((FF_med - FF_true_arr)**2))

    # Filtered: use per-year filtering (more lenient)
    if FF_dD is not None:
        agreement = np.abs(FF_est - FF_dD) < threshold
        valid = ~(np.isnan(FF_est) | np.isnan(FF_dD))
        agreement = agreement & valid

        # Per-year: only keep agreeing iterations for that year
        FF_filt_med = np.zeros(N_YEARS)
        n_agree_per_year = np.zeros(N_YEARS)
        for j in range(N_YEARS):
            mask = agreement[j, :]
            n_agree_per_year[j] = mask.sum()
            if mask.sum() > 10:
                FF_filt_med[j] = np.nanmedian(FF_est[j, mask])
            else:
                FF_filt_med[j] = FF_med[j]  # fall back

        bias_filt = np.nanmean(FF_filt_med - FF_true_arr)
        rmse_filt = np.sqrt(np.nanmean((FF_filt_med - FF_true_arr)**2))
        n_good = int(np.mean(n_agree_per_year))
    else:
        bias_filt, rmse_filt, n_good = np.nan, np.nan, 0

    return {
        'bias_unfilt': float(bias_unfilt),
        'rmse_unfilt': float(rmse_unfilt),
        'bias_filt': float(bias_filt),
        'rmse_filt': float(rmse_filt),
        'n_good_iters': int(n_good),
    }


def main():
    global TRUE_FF, TRUE_MIC, TRUE_BB, TRUE_TOTAL

    print("=" * 60)
    print("  Phase 6c — OSSE (Observing System Simulation Experiment)")
    print("=" * 60)

    # Step 0: Compute TRUE emissions from observed mass balance
    data = load_data(REPO_ROOT, two_box=False)
    CH4 = data.CH4_global
    tau = compute_lifetime(YEARS, 'varying')
    S_obs = np.zeros(N_YEARS)
    for i in range(N_YEARS):
        S_obs[i] = CH4[i+1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]

    # Partition with time-varying fractions (FF stable, Mic increases, BB decreases)
    # FF fraction: 24% throughout
    # BB fraction: decreases from 6% to 4%
    # Mic: remainder
    bb_frac = np.linspace(0.06, 0.04, N_YEARS)
    ff_frac = np.full(N_YEARS, FF_FRACTION)
    mic_frac = 1.0 - ff_frac - bb_frac

    TRUE_TOTAL = S_obs
    TRUE_FF = S_obs * ff_frac
    TRUE_BB = S_obs * bb_frac
    TRUE_MIC = S_obs * mic_frac

    print(f"    TRUE_TOTAL: {S_obs.mean():.0f} Tg/yr (range {S_obs.min():.0f}–{S_obs.max():.0f})")
    print(f"    TRUE_FF: {TRUE_FF.mean():.0f} Tg/yr, TRUE_MIC: {TRUE_MIC.mean():.0f}, TRUE_BB: {TRUE_BB.mean():.0f}")

    # Step 1: Forward model synthetic atmosphere
    print("\n  Forward-modeling synthetic atmosphere...")
    CH4, d13C_true, dD_true = forward_model_atmosphere()
    print(f"    δ¹³C range: {d13C_true.min():.2f} to {d13C_true.max():.2f}")
    print(f"    δD range: {dD_true.min():.1f} to {dD_true.max():.1f}")

    # Step 2: Invert with different KIE assumptions
    print("\n  Inverting with different KIE assumptions...")
    osse_results = {}
    for oh13c in ['true', 'saueressig', 'cantrell']:
        print(f"    OH_13C = {oh13c}...")
        FF_c, FF_d, Mic_c = invert_synthetic(CH4, d13C_true, dD_true, oh13c)
        osse_results[oh13c] = {'FF_c13': FF_c, 'FF_dD': FF_d, 'Mic_c13': Mic_c}

    # Step 3: Compute recovery metrics
    print("\n  Computing recovery metrics...")
    recovery = {}
    for oh13c in ['true', 'saueressig', 'cantrell']:
        r = osse_results[oh13c]
        metrics = compute_recovery_error(
            r['FF_c13'], TRUE_FF, YEARS, threshold=100,
            FF_dD=r['FF_dD'], Mic_est=r['Mic_c13'], Mic_true_arr=TRUE_MIC)
        recovery[oh13c] = metrics
        print(f"    {oh13c:12s}: bias_unfilt={metrics['bias_unfilt']:+.1f}, "
              f"rmse_unfilt={metrics['rmse_unfilt']:.1f} | "
              f"bias_filt={metrics['bias_filt']:+.1f}, "
              f"rmse_filt={metrics['rmse_filt']:.1f} "
              f"(n={metrics['n_good_iters']})")

    # Step 4: Agreement rates in OSSE
    print("\n  Agreement rates in OSSE...")
    for oh13c in ['true', 'saueressig', 'cantrell']:
        r = osse_results[oh13c]
        valid = ~(np.isnan(r['FF_c13']) | np.isnan(r['FF_dD']))
        agree = (np.abs(r['FF_c13'] - r['FF_dD']) < 100) & valid
        rate = agree.sum() / max(valid.sum(), 1)
        print(f"    {oh13c:12s}: agreement = {rate:.1%}")
        recovery[oh13c]['agreement_rate'] = float(rate)

    # Save
    summary = {
        'true_emissions': {
            'FF': TRUE_FF.tolist(),
            'Mic': TRUE_MIC.tolist(),
            'BB': TRUE_BB.tolist(),
        },
        'true_KIE_OH_13C': TRUE_KIE_OH_13C,
        'recovery_metrics': recovery,
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # === Figures ===

    # Figure 11: OSSE recovery comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=300)

    # (a) FF time series: true vs recovered
    ax = axes[0, 0]
    ax.plot(YEARS, TRUE_FF, 'k-', lw=3, label='TRUE', zorder=10)
    for oh13c, color, name in [('true', 'green', 'True KIE'),
                                ('saueressig', 'tab:blue', 'Saueressig'),
                                ('cantrell', 'tab:red', 'Cantrell')]:
        FF_med = np.nanmedian(osse_results[oh13c]['FF_c13'], axis=1)
        FF_lo = np.nanpercentile(osse_results[oh13c]['FF_c13'], 16, axis=1)
        FF_hi = np.nanpercentile(osse_results[oh13c]['FF_c13'], 84, axis=1)
        ax.plot(YEARS, FF_med, '-', color=color, lw=1.5, label=f'{name} (δ¹³C)')
        ax.fill_between(YEARS, FF_lo, FF_hi, alpha=0.1, color=color)
    ax.set_ylabel('FF Emissions (Tg/yr)')
    ax.set_title('(a) FF Recovery — Unfiltered (δ¹³C-only)')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # (b) FF time series: agreement-filtered
    ax = axes[0, 1]
    ax.plot(YEARS, TRUE_FF, 'k-', lw=3, label='TRUE', zorder=10)
    for oh13c, color, name in [('true', 'green', 'True KIE'),
                                ('saueressig', 'tab:blue', 'Saueressig'),
                                ('cantrell', 'tab:red', 'Cantrell')]:
        r = osse_results[oh13c]
        agreement = np.abs(r['FF_c13'] - r['FF_dD']) < 100
        valid = ~(np.isnan(r['FF_c13']) | np.isnan(r['FF_dD']))
        good = ((agreement & valid).sum(axis=0) >= N_YEARS * 0.8)
        if good.sum() > 10:
            FF_filt = r['FF_c13'][:, good]
            FF_med = np.nanmedian(FF_filt, axis=1)
            FF_lo = np.nanpercentile(FF_filt, 16, axis=1)
            FF_hi = np.nanpercentile(FF_filt, 84, axis=1)
            ax.plot(YEARS, FF_med, '-', color=color, lw=1.5,
                    label=f'{name} (filtered, n={good.sum()})')
            ax.fill_between(YEARS, FF_lo, FF_hi, alpha=0.1, color=color)
    ax.set_ylabel('FF Emissions (Tg/yr)')
    ax.set_title('(b) FF Recovery — Agreement-Filtered')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # (c) RMSE comparison bar chart
    ax = axes[1, 0]
    methods = ['true', 'saueressig', 'cantrell']
    rmse_unfilt = [recovery[m]['rmse_unfilt'] for m in methods]
    rmse_filt = [recovery[m]['rmse_filt'] for m in methods]
    x = np.arange(len(methods))
    w = 0.35
    ax.bar(x - w/2, rmse_unfilt, w, label='Unfiltered (δ¹³C-only)', color='salmon', alpha=0.7)
    ax.bar(x + w/2, rmse_filt, w, label='Agreement-filtered', color='steelblue', alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(['True KIE\n(1.0046)', 'Saueressig\n(1.0039)', 'Cantrell\n(1.0054)'])
    ax.set_ylabel('RMSE (Tg/yr)')
    ax.set_title('(c) Recovery RMSE: Unfiltered vs Filtered')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis='y')

    # (d) Bias comparison
    ax = axes[1, 1]
    bias_unfilt = [recovery[m]['bias_unfilt'] for m in methods]
    bias_filt = [recovery[m]['bias_filt'] for m in methods]
    ax.bar(x - w/2, bias_unfilt, w, label='Unfiltered', color='salmon', alpha=0.7)
    ax.bar(x + w/2, bias_filt, w, label='Filtered', color='steelblue', alpha=0.7)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(['True KIE\n(1.0046)', 'Saueressig\n(1.0039)', 'Cantrell\n(1.0054)'])
    ax.set_ylabel('Bias (Tg/yr)')
    ax.set_title('(d) Recovery Bias: Unfiltered vs Filtered')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis='y')

    plt.suptitle('Phase 6c: OSSE — Can Agreement Filter Improve Recovery?', fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig11_OSSE_recovery.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 11 saved: {FIG_DIR / 'fig11_OSSE_recovery.png'}")
    print("\n✓ Phase 6c complete.")


if __name__ == "__main__":
    main()
