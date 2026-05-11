#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 4: Two-Box (NH/SH) Extension — KIE Sensitivity
======================================================
Repeats the Phase 1–3 analysis using a hemispheric 2-box model.
Per-hemisphere WLS: 3×2 system (mass + δ¹³C + δD) for each hemisphere.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import ks_2samp
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
OUT_DIR = BASE / "results" / "phase4_two_box"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 1000
SEED = 42

W_MASS = 1.0
W_C13 = 1.0 / (1.0**2)
W_D = 1.0 / (5.0**2)


def run_two_box(oh13c_mode: str, label: str, dual: bool = True):
    """Run 2-box model with either dual-isotope or δ¹³C-only.
    
    dual=True: WLS 3×2 per hemisphere (mass + δ¹³C + δD)
    dual=False: 2×2 per hemisphere (mass + δ¹³C only)
    """
    mode_str = "dual" if dual else "d13C-only"
    print(f"\n  2-box | {mode_str} | OH_13C={oh13c_mode} ({label})")

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

        # Inter-hemispheric exchange
        tau_ex = rng.normal(TAU_EX_MEAN, TAU_EX_STD)
        tau_ex = max(0.5, tau_ex)

        # Atmospheric observations
        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm = sample_atm_dD(data, k, n)

        # NH/SH split for δ¹³C
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

            # Exchange fluxes
            F_ex_NH_to_SH = CH4_NH[j] * PT_HEMI / tau_ex
            F_ex_SH_to_NH = CH4_SH[j] * PT_HEMI / tau_ex

            # Source totals per hemisphere (mass balance)
            S_NH = (CH4_NH[j+1] - CH4_NH[j]) * PT_HEMI + CH4_NH[j] * PT_HEMI / tau_NH[j] + F_ex_NH_to_SH - F_ex_SH_to_NH
            S_SH = (CH4_SH[j+1] - CH4_SH[j]) * PT_HEMI + CH4_SH[j] * PT_HEMI / tau_SH[j] + F_ex_SH_to_NH - F_ex_NH_to_SH

            # Isotopic source compositions — NH
            f13_NH = delta_to_fraction_d13C(np.array([c13_NH[j], c13_NH[j+1] if j+1 < len(c13_NH) else c13_NH[j]]))
            n13_NH = f13_NH * np.array([CH4_NH[j], CH4_NH[j+1]]) * PT_HEMI
            d13C_src_NH_f = (n13_NH[1] - n13_NH[0] + n13_NH[0] * alpha_13C_NH / tau_NH[j]) / max(S_NH, 1e-6) if S_NH > 0 else 0
            d13C_src_NH = fraction_to_delta_d13C(d13C_src_NH_f) if S_NH > 0 else -50.0

            fD_NH = delta_to_fraction_dD(np.array([dD_NH[j], dD_NH[j+1] if j+1 < len(dD_NH) else dD_NH[j]]))
            nD_NH = fD_NH * np.array([CH4_NH[j], CH4_NH[j+1]]) * PT_HEMI
            dD_src_NH_f = (nD_NH[1] - nD_NH[0] + nD_NH[0] * alpha_D_NH / tau_NH[j]) / max(S_NH, 1e-6) if S_NH > 0 else 0
            dD_src_NH = fraction_to_delta_dD(dD_src_NH_f) if S_NH > 0 else -280.0

            # SH
            f13_SH = delta_to_fraction_d13C(np.array([c13_SH[j], c13_SH[j+1] if j+1 < len(c13_SH) else c13_SH[j]]))
            n13_SH = f13_SH * np.array([CH4_SH[j], CH4_SH[j+1]]) * PT_HEMI
            d13C_src_SH_f = (n13_SH[1] - n13_SH[0] + n13_SH[0] * alpha_13C_SH / tau_SH[j]) / max(S_SH, 1e-6) if S_SH > 0 else 0
            d13C_src_SH = fraction_to_delta_d13C(d13C_src_SH_f) if S_SH > 0 else -50.0

            fD_SH = delta_to_fraction_dD(np.array([dD_SH[j], dD_SH[j+1] if j+1 < len(dD_SH) else dD_SH[j]]))
            nD_SH = fD_SH * np.array([CH4_SH[j], CH4_SH[j+1]]) * PT_HEMI
            dD_src_SH_f = (nD_SH[1] - nD_SH[0] + nD_SH[0] * alpha_D_SH / tau_SH[j]) / max(S_SH, 1e-6) if S_SH > 0 else 0
            dD_src_SH = fraction_to_delta_dD(dD_src_SH_f) if S_SH > 0 else -280.0

            # Solve per-hemisphere
            ff_nh, mic_nh = _solve_hemisphere(
                S_NH, BB_NH_j, d13C_src_NH, dD_src_NH, sigs, j, dual)
            ff_sh, mic_sh = _solve_hemisphere(
                S_SH, BB_SH_j, d13C_src_SH, dD_src_SH, sigs, j, dual)

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

    print(f"    FF: Δ={delta_ff.mean():+.1f}±{delta_ff.std():.1f} | Mic: Δ={delta_mic.mean():+.1f}±{delta_mic.std():.1f}")
    return delta_ff, delta_mic


def _solve_hemisphere(S, BB_j, d13C_src, dD_src, sigs, j, dual):
    """Solve for FF, Mic in one hemisphere."""
    if S <= 0:
        return 0.0, 0.0

    if dual:
        A = np.array([
            [1.0, 1.0],
            [sigs['ff_d13C'][j], sigs['mic_d13C'][j]],
            [sigs['ff_dD'][j], sigs['mic_dD'][j]],
        ])
        b = np.array([
            S - BB_j,
            S * d13C_src - BB_j * sigs['bb_d13C'][j],
            S * dD_src - BB_j * sigs['bb_dD'][j],
        ])
        w = np.sqrt(np.array([W_MASS, W_C13, W_D]))
        Aw = A * w[:, None]
        bw = b * w
        try:
            result = lsq_linear(Aw, bw, bounds=(0, np.inf))
            return float(result.x[0]), float(result.x[1])
        except Exception:
            return 0.0, max(0, S - BB_j)
    else:
        # δ¹³C-only 2×2
        denom = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
        if abs(denom) > 0.1:
            ff = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j) -
                  sigs['bb_d13C'][j] * BB_j) / denom
            mic = S - BB_j - ff
            ff = max(0, ff)
            mic = max(0, mic)
            return ff, mic
        return 0.0, max(0, S - BB_j)


def main():
    print("=" * 60)
    print("  Phase 4 — Two-Box (NH/SH) KIE Sensitivity")
    print("=" * 60)

    # Run all 6 combinations: 3 KIE × 2 methods
    results = {}
    for oh13c, lbl in [('saueressig', 'A_saueressig'), ('cantrell', 'B_cantrell'), ('sampled', 'C_sampled')]:
        for dual in [True, False]:
            key = f"{lbl}_{'dual' if dual else 'd13C'}"
            dff, dmic = run_two_box(oh13c, lbl, dual=dual)
            results[key] = {'delta_ff': dff, 'delta_mic': dmic}

    # Compute KSR for 2-box
    spread_s_ff = abs(results['B_cantrell_d13C']['delta_ff'].mean() -
                      results['A_saueressig_d13C']['delta_ff'].mean())
    spread_d_ff = abs(results['B_cantrell_dual']['delta_ff'].mean() -
                      results['A_saueressig_dual']['delta_ff'].mean())
    ksr_ff_2box = spread_s_ff / max(spread_d_ff, 1e-6)

    spread_s_mic = abs(results['B_cantrell_d13C']['delta_mic'].mean() -
                       results['A_saueressig_d13C']['delta_mic'].mean())
    spread_d_mic = abs(results['B_cantrell_dual']['delta_mic'].mean() -
                       results['A_saueressig_dual']['delta_mic'].mean())
    ksr_mic_2box = spread_s_mic / max(spread_d_mic, 1e-6)

    print(f"\n  2-Box KSR (FF):  {ksr_ff_2box:.2f}")
    print(f"  2-Box KSR (Mic): {ksr_mic_2box:.2f}")

    # Load 1-box KSR from phase3
    p3_summary_path = BASE / "results" / "phase3_comparison" / "summary.json"
    if p3_summary_path.exists():
        with open(p3_summary_path) as f:
            p3 = json.load(f)
        ksr_ff_1box = p3['KSR_FF']['value']
        ksr_mic_1box = p3['KSR_Mic']['value']
    else:
        ksr_ff_1box = None
        ksr_mic_1box = None

    # Save
    summary = {
        'KSR_FF_2box': ksr_ff_2box,
        'KSR_Mic_2box': ksr_mic_2box,
        'KSR_FF_1box': ksr_ff_1box,
        'KSR_Mic_1box': ksr_mic_1box,
        'spread_single_FF': float(spread_s_ff),
        'spread_dual_FF': float(spread_d_ff),
        'spread_single_Mic': float(spread_s_mic),
        'spread_dual_Mic': float(spread_d_mic),
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # --- Figure: KSR comparison bar chart ---
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

    labels = ['FF (1-box)', 'FF (2-box)', 'Mic (1-box)', 'Mic (2-box)']
    values = [ksr_ff_1box or 0, ksr_ff_2box, ksr_mic_1box or 0, ksr_mic_2box]
    colors = ['tab:red', 'darkred', 'tab:blue', 'darkblue']

    bars = ax.bar(labels, values, color=colors, alpha=0.7, edgecolor='black')
    ax.axhline(1.0, color='gray', ls='--', lw=1, label='KSR = 1 (no improvement)')
    ax.set_ylabel('KIE Sensitivity Ratio (KSR)')
    ax.set_title('Does Dual-Isotope Reduce KIE Sensitivity?\n1-Box vs 2-Box Comparison')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{val:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig4_KSR_1box_vs_2box.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure saved: {FIG_DIR / 'fig4_KSR_1box_vs_2box.png'}")
    print("\n✓ Phase 4 complete.")


if __name__ == "__main__":
    main()
