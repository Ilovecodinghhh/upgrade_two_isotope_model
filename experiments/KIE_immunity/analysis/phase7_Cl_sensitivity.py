#!/usr/bin/env python3
"""Phase 7 — Cl fraction sensitivity on FF trend + variance decomposition."""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR
from common import smooth_5yr

def main():
    print("=" * 60)
    print("PHASE 7: Cl fraction sensitivity (dual real-hemi-dD)")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400

    cl_values = [
        ("Thanwerdas (0.6%)", 0.006),
        ("Low (2.0%)",        0.020),
        ("Default (3.5%)",    0.035),
        ("Medium (5.0%)",     0.050),
        ("High (6.5%)",       0.065),
        ("Allan upper (10%)", 0.100),
    ]

    results = {}
    print(f"\n{'Config':<22} {'ΔFF med':>8} {'[5%':>7} {'95%]':>7} {'σ(FF)':>7} {'KIE%':>6} {'Sig%':>6} {'τ%':>6}")
    print("-" * 80)

    for label, cl in cl_values:
        FF = run_2box_flex(data, N, 42, cl_frac=cl)
        med, lo, hi = trend_stats(FF, data.model_years)
        sig = sigma_ff(FF)

        # Variance decomposition
        var_total = np.nanmean(np.nanvar(smooth_5yr(FF)[8:], axis=1))

        FF_fk = run_2box_flex(data, N, 42, cl_frac=cl, fix_kie=True)
        var_fk = np.nanmean(np.nanvar(smooth_5yr(FF_fk)[8:], axis=1))
        kie_pct = max(0, (var_total - var_fk) / var_total * 100) if var_total > 0 else 0

        FF_fs = run_2box_flex(data, N, 42, cl_frac=cl, fix_sigs=True)
        var_fs = np.nanmean(np.nanvar(smooth_5yr(FF_fs)[8:], axis=1))
        sig_pct = max(0, (var_total - var_fs) / var_total * 100) if var_total > 0 else 0

        FF_ft = run_2box_flex(data, N, 42, cl_frac=cl, tau_mode="fixed", tau_fixed=9.0)
        var_ft = np.nanmean(np.nanvar(smooth_5yr(FF_ft)[8:], axis=1))
        tau_pct = max(0, (var_total - var_ft) / var_total * 100) if var_total > 0 else 0

        print(f"{label:<22} {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}  {sig:>6.1f} {kie_pct:>5.1f} {sig_pct:>6.1f} {tau_pct:>5.1f}")

        results[label] = {
            'cl_frac': cl,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
            'sigma_ff': sig,
            'kie_pct': kie_pct, 'sig_pct': sig_pct, 'tau_pct': tau_pct,
        }

    print(f"\n{'=' * 60}")
    print("SIGN ANALYSIS")
    print(f"{'=' * 60}")
    for label, r in results.items():
        sign = "NEGATIVE" if r['trend_median'] < 0 else "POSITIVE"
        print(f"  {label:<22} ΔFF = {r['trend_median']:+.1f}  {sign}")

    with open(OUT_DIR / "phase7_Cl_sensitivity.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase7_Cl_sensitivity.json'}")


if __name__ == "__main__":
    main()
