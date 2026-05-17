#!/usr/bin/env python3
"""
phase14_W_sensitivity.py — Weight matrix sensitivity analysis (Review A2)
=========================================================================
Tests impact of W on all headline numbers.
"""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, compute_trend_regression, load_data, REPO_ROOT, OUT_DIR
from common import smooth_5yr

def main():
    print("=" * 70)
    print("PHASE 14: W matrix sensitivity (Review issue A2)")
    print("=" * 70)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400; SEED = 42

    w_configs = [
        ("Identity",               np.diag([1.0, 1.0, 1.0])),
        ("Equal isotopes",         np.diag([100.0, 1.0, 1.0])),
        ("Default (100,1,0.5)",    np.diag([100.0, 1.0, 0.5])),
        ("δD upweighted",         np.diag([100.0, 1.0, 2.0])),
        ("δD dominant",           np.diag([100.0, 0.5, 2.0])),
        ("Inverse-variance est",   np.diag([100.0, 1.0/0.05, 1.0/1.0])),  # ~inv measurement uncertainty
    ]

    results = {}
    print(f"\n{'Config':<25} {'σ(FF)':>7} {'ΔFF step':>10} {'slope':>8} {'KIE%':>6} {'Sig%':>6}")
    print("-" * 72)

    for label, W in w_configs:
        FF = run_2box_flex(data, N, SEED, W=W)
        sig = sigma_ff(FF)
        med, lo, hi = trend_stats(FF, data.model_years)
        reg = compute_trend_regression(FF, data.model_years)

        # Variance decomposition
        var_total = np.nanmean(np.nanvar(smooth_5yr(FF)[8:], axis=1))

        FF_fk = run_2box_flex(data, N, SEED, W=W, fix_kie=True)
        var_fk = np.nanmean(np.nanvar(smooth_5yr(FF_fk)[8:], axis=1))
        kie_pct = max(0, (var_total - var_fk) / var_total * 100) if var_total > 0 else 0

        FF_fs = run_2box_flex(data, N, SEED, W=W, fix_sigs=True)
        var_fs = np.nanmean(np.nanvar(smooth_5yr(FF_fs)[8:], axis=1))
        sig_pct = max(0, (var_total - var_fs) / var_total * 100) if var_total > 0 else 0

        print(f"{label:<25} {sig:>6.1f}  {med:>+7.1f}     {reg['slope_median']:>+6.2f}  {kie_pct:>5.1f} {sig_pct:>6.1f}")

        results[label] = {
            'W_diag': W.diagonal().tolist(),
            'sigma_ff': sig,
            'trend_step_median': med, 'trend_step_5pct': lo, 'trend_step_95pct': hi,
            'trend_regression': reg,
            'kie_pct': kie_pct, 'sig_pct': sig_pct,
        }

    # Summary
    sigmas = [r['sigma_ff'] for r in results.values()]
    print(f"\n  σ(FF) range across W: {min(sigmas):.1f} – {max(sigmas):.1f} Tg/yr")
    signs = [r['trend_step_median'] for r in results.values()]
    if all(s > 0 for s in signs):
        print("  ✓ All W configs give POSITIVE FF trend")
    elif all(s < 0 for s in signs):
        print("  ✓ All W configs give NEGATIVE FF trend")
    else:
        print("  ⚠ FF trend SIGN DEPENDS on W choice!")

    with open(OUT_DIR / "phase14_W_sensitivity.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase14_W_sensitivity.json'}")


if __name__ == "__main__":
    main()
