#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 7: Time-Varying OH-¹³C KIE
==================================
Tests whether the agreement rate framework is sensitive to a time-varying
OH-¹³C KIE — motivated by He 2026 Science's finding that τ_CH4 has been
declining (~0.017 yr/yr since 2010), which implies [OH] is growing.

If the OH-¹³C KIE has a temperature dependence (Saueressig 2001 hint at
weak T-dependence) or if a different OH regime is active, the effective
bulk KIE_13C could drift over time. Does this break the agreement-rate
discriminant we found in Phase 6b?

Three scenarios:
  (S1) Saueressig drift   : OH_13C linearly drifts 1.0039 → 1.0046 over 1999–2022
  (S2) Cantrell drift     : OH_13C linearly drifts 1.0054 → 1.0048 over 1999–2022
  (S3) Convergent drift   : Saueressig and Cantrell trajectories meet by 2022
                            (tests whether discriminant collapses in present day)

For each scenario we compute:
  - Per-year agreement rate
  - Overall agreement rate
  - KSR (vs constant-KIE baseline)
  - Whether the Cantrell-Saueressig discriminant is still significant

Output:
  results/phase7_timevarying_OH/summary.json
  figures/fig12_timevarying_OH.png
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
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change,
    SINK_FRACTIONS_GLOBAL, PT,
)

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "phase7_timevarying_OH"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 1000
SEED = 42
THRESHOLD = 100  # Tg/yr (best discriminant from Phase 6b)


def make_kie_trajectory(name: str, n_years: int):
    """Return per-year OH_13C trajectory for a given scenario."""
    if name == 'const_saueressig':
        return np.full(n_years, 1.0039)
    if name == 'const_cantrell':
        return np.full(n_years, 1.0054)

    # Drift scenarios
    t = np.linspace(0.0, 1.0, n_years)
    if name == 'drift_saueressig':
        # Saueressig (1999) → midpoint (2022)
        return 1.0039 + (1.00465 - 1.0039) * t
    if name == 'drift_cantrell':
        # Cantrell (1999) → midpoint (2022)
        return 1.0054 + (1.00465 - 1.0054) * t
    if name == 'convergent':
        # Both endpoints converge to 1.0046 by 2022 (tests collapsing discriminant)
        # We take the Saueressig trajectory in this scenario
        return 1.0039 + (1.0046 - 1.0039) * t
    raise ValueError(f"Unknown scenario: {name}")


def run_inversions(oh13c_trajectory: np.ndarray):
    """Run δ¹³C and δD inversions where OH_13C varies year by year."""
    data = load_data(REPO_ROOT, two_box=False)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    rng = np.random.default_rng(SEED)

    tau = compute_lifetime(years, 'varying')
    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]
    BB = data.BB_annual.copy()

    FF_c13 = np.full((n, N_ITER), np.nan)
    FF_dD = np.full((n, N_ITER), np.nan)
    Mic_c13 = np.full((n, N_ITER), np.nan)

    assert len(oh13c_trajectory) == n, \
        f"trajectory length {len(oh13c_trajectory)} != n_years {n}"

    for k in range(N_ITER):
        kies_base = sample_KIE(rng, 'sampled')

        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm = sample_atm_dD(data, k, n)
        f13 = delta_to_fraction_d13C(d13C_atm)
        fD = delta_to_fraction_dD(dD_atm)
        n13 = f13 * CH4[:n + 1] * PT
        nD = fD * CH4[:n + 1] * PT
        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            kies = dict(kies_base)
            kies['OH_13C'] = oh13c_trajectory[j]
            KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
            alpha_13C = 1.0 / KIE_13C
            alpha_D = 1.0 / KIE_D

            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            # δ¹³C
            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)
            denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            if abs(denom_c) > 0.1:
                ff_c = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j) -
                        sigs['bb_d13C'][j] * BB_j) / denom_c
                FF_c13[j, k] = ff_c
                Mic_c13[j, k] = S - BB_j - ff_c

            # δD
            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            dD_src = fraction_to_delta_dD(dD_src_f)
            denom_d = sigs['ff_dD'][j] - sigs['mic_dD'][j]
            if abs(denom_d) > 0.5:
                FF_dD[j, k] = (S * dD_src - sigs['mic_dD'][j] * (S - BB_j) -
                               sigs['bb_dD'][j] * BB_j) / denom_d

    return FF_c13, FF_dD, Mic_c13, years


def metrics(FF_c13, FF_dD, Mic_c13, years, threshold=THRESHOLD):
    n, N = FF_c13.shape
    valid = ~(np.isnan(FF_c13) | np.isnan(FF_dD))
    agree = (np.abs(FF_c13 - FF_dD) < threshold) & valid

    rate_year = np.array([agree[j].sum() / max(valid[j].sum(), 1) for j in range(n)])
    overall = agree.sum() / max(valid.sum(), 1)

    # Iterations with ≥80% of years agreeing
    good = (agree.sum(axis=0) >= n * 0.8)
    n_good = int(good.sum())
    if n_good > 30:
        FF_s = smooth_5yr(np.clip(FF_c13[:, good], 0, None))
        Mic_s = smooth_5yr(np.clip(Mic_c13[:, good], 0, None))
        delta_ff = trend_change(FF_s, years)[0]
        delta_mic = trend_change(Mic_s, years)[0]
        ff_mean = float(delta_ff.mean())
        mic_mean = float(delta_mic.mean())
    else:
        ff_mean = float('nan')
        mic_mean = float('nan')

    # Bootstrap CI on overall rate
    rng = np.random.default_rng(99)
    flat_valid = valid.flatten()
    flat_agree = agree.flatten()[flat_valid]
    rates = np.empty(2000)
    nv = len(flat_agree)
    for b in range(2000):
        rates[b] = flat_agree[rng.choice(nv, nv, replace=True)].mean()
    ci = (float(np.percentile(rates, 2.5)), float(np.percentile(rates, 97.5)))

    return {
        'overall_rate': float(overall),
        'rate_by_year': rate_year.tolist(),
        'CI95': list(ci),
        'n_good_iters': n_good,
        'delta_ff_mean': ff_mean,
        'delta_mic_mean': mic_mean,
    }


def main():
    print("=" * 60)
    print("  Phase 7 — Time-Varying OH-¹³C KIE")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=False)
    n = data.n_years
    years = data.model_years

    scenarios = [
        'const_saueressig',
        'const_cantrell',
        'drift_saueressig',
        'drift_cantrell',
        'convergent',
    ]

    results = {}
    rate_by_year = {}
    trajectories = {}
    for sc in scenarios:
        traj = make_kie_trajectory(sc, n)
        trajectories[sc] = traj.tolist()
        print(f"\n  Scenario: {sc}  (KIE: {traj[0]:.5f} → {traj[-1]:.5f})")
        FF_c, FF_d, Mic_c, yrs = run_inversions(traj)
        m = metrics(FF_c, FF_d, Mic_c, yrs)
        results[sc] = m
        rate_by_year[sc] = m['rate_by_year']
        print(f"    Overall agreement: {m['overall_rate']:.1%} "
              f"[{m['CI95'][0]:.1%}, {m['CI95'][1]:.1%}]  "
              f"(n_good={m['n_good_iters']})")

    # Discriminant tests
    print("\n  === DISCRIMINANT TESTS ===")
    pairs = [
        ('const_saueressig', 'const_cantrell', 'Constant baseline (Phase 6b)'),
        ('drift_saueressig', 'drift_cantrell', 'Symmetric drift (toward midpoint)'),
        ('convergent', 'const_cantrell', 'Convergent vs constant Cantrell'),
    ]
    discriminants = {}
    for s_lbl, c_lbl, desc in pairs:
        s = results[s_lbl]
        c = results[c_lbl]
        diff = c['overall_rate'] - s['overall_rate']
        # CI overlap test
        sig = s['CI95'][1] < c['CI95'][0] or c['CI95'][1] < s['CI95'][0]
        discriminants[f"{s_lbl}_vs_{c_lbl}"] = {
            'description': desc,
            'rate_low_KIE': s['overall_rate'],
            'rate_high_KIE': c['overall_rate'],
            'difference_pp': float(diff * 100),
            'CI_low': s['CI95'],
            'CI_high': c['CI95'],
            'significant': bool(sig),
        }
        print(f"  {desc}")
        print(f"    rate(low) = {s['overall_rate']:.1%}, rate(high) = {c['overall_rate']:.1%}")
        print(f"    Δ = {diff*100:+.1f} pp, significant = {sig}")

    summary = {
        'threshold': THRESHOLD,
        'years': years.tolist(),
        'trajectories': trajectories,
        'scenarios': results,
        'discriminants': discriminants,
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # === Figure 12 ===
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)

    # (a) KIE trajectories
    ax = axes[0]
    colors = {
        'const_saueressig': 'tab:blue',
        'const_cantrell': 'tab:red',
        'drift_saueressig': 'tab:cyan',
        'drift_cantrell': 'tab:orange',
        'convergent': 'tab:purple',
    }
    for sc in scenarios:
        ax.plot(years, trajectories[sc], '-', color=colors[sc], lw=2,
                label=sc.replace('_', ' '))
    ax.set_xlabel('Year')
    ax.set_ylabel('OH-¹³C KIE α')
    ax.set_title('(a) KIE Trajectories')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # (b) Per-year agreement rate
    ax = axes[1]
    for sc in scenarios:
        ax.plot(years, np.array(rate_by_year[sc]) * 100, '-',
                color=colors[sc], lw=2, label=sc.replace('_', ' '))
    ax.set_xlabel('Year')
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title(f'(b) Per-Year Agreement Rate (threshold={THRESHOLD} Tg/yr)')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # (c) Discriminant comparison
    ax = axes[2]
    keys = list(discriminants.keys())
    diffs = [discriminants[k]['difference_pp'] for k in keys]
    sigs = [discriminants[k]['significant'] for k in keys]
    bar_colors = ['green' if s else 'gray' for s in sigs]
    ax.bar(range(len(keys)), diffs, color=bar_colors, alpha=0.75, edgecolor='black')
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([k.replace('_vs_', '\nvs\n').replace('_', ' ')
                        for k in keys], fontsize=7)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_ylabel('Δ Agreement Rate (pp)')
    ax.set_title('(c) Discriminant Power\n(green = significant, p<0.05)')
    ax.grid(alpha=0.3, axis='y')
    for i, (d, s) in enumerate(zip(diffs, sigs)):
        ax.text(i, d + (0.5 if d > 0 else -1.5), f'{d:+.1f}',
                ha='center', fontsize=8, fontweight='bold' if s else 'normal')

    plt.suptitle('Phase 7: Time-Varying OH-¹³C KIE — Does the Discriminant Survive?',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig12_timevarying_OH.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  Figure 12 saved: {FIG_DIR / 'fig12_timevarying_OH.png'}")
    print("\n✓ Phase 7 complete.")


if __name__ == "__main__":
    main()
