#!/usr/bin/env python3
"""Phase 5 — Lifetime sensitivity on FF trend."""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR

def main():
    print("=" * 60)
    print("PHASE 5: Lifetime sensitivity (dual real-hemi-dD)")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400

    configs = [
        ("τ=8.0 fixed",   "fixed", 8.0),
        ("τ=8.5 fixed",   "fixed", 8.5),
        ("τ=9.0 fixed",   "fixed", 9.0),
        ("τ=9.5 fixed",   "fixed", 9.5),
        ("τ=10.0 fixed",  "fixed", 10.0),
        ("He 2026 varying", "varying", 9.0),
    ]

    results = {}
    print(f"\n{'Config':<22} {'ΔFF med':>8} {'[5%':>7} {'95%]':>7} {'σ(FF)':>7}")
    print("-" * 55)

    for label, mode, tau_val in configs:
        FF = run_2box_flex(data, N, 42, tau_mode=mode, tau_fixed=tau_val)
        med, lo, hi = trend_stats(FF, data.model_years)
        sig = sigma_ff(FF)

        # Also compute KIE contribution
        FF_fix_kie = run_2box_flex(data, N, 42, tau_mode=mode, tau_fixed=tau_val,
                                   fix_kie=True)
        var_total = np.nanmean(np.nanvar(FF[8:], axis=1))
        var_no_kie = np.nanmean(np.nanvar(FF_fix_kie[8:], axis=1))
        kie_pct = max(0, (var_total - var_no_kie) / var_total * 100) if var_total > 0 else 0

        print(f"{label:<22} {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}  {sig:>6.1f}")

        results[label] = {
            'tau_mode': mode, 'tau_value': tau_val,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
            'sigma_ff': sig, 'kie_pct': kie_pct,
        }

    # Identify if trend flips sign
    print(f"\n{'=' * 60}")
    print("SIGN ANALYSIS")
    print(f"{'=' * 60}")
    for label, r in results.items():
        sign = "NEGATIVE" if r['trend_median'] < 0 else "POSITIVE"
        ci_contains_zero = r['trend_5pct'] < 0 < r['trend_95pct']
        robust = "✓ sign robust" if not ci_contains_zero else "⚠ CI includes zero"
        print(f"  {label:<22} ΔFF = {r['trend_median']:+.1f}  {sign:<9} {robust}")

    with open(OUT_DIR / "phase5_tau_sensitivity.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase5_tau_sensitivity.json'}")


if __name__ == "__main__":
    main()
