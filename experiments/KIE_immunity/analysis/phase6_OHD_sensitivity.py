#!/usr/bin/env python3
"""Phase 6 — OH-D KIE sensitivity on FF trend."""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR

def main():
    print("=" * 60)
    print("PHASE 6: OH-D KIE sensitivity (dual real-hemi-dD)")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400

    oh_d_values = [
        ("Saueressig (1.294)", 1.294),
        ("Midpoint (1.310)",   1.310),
        ("Cantrell (1.327)",   1.327),
        ("He 2026 upper (1.35)", 1.350),
        ("Sampled (default)",  None),
    ]

    results = {}
    print(f"\n{'Config':<25} {'ΔFF med':>8} {'[5%':>7} {'95%]':>7} {'σ(FF)':>7}")
    print("-" * 58)

    for label, oh_d in oh_d_values:
        FF = run_2box_flex(data, N, 42, oh_d_fixed=oh_d)
        med, lo, hi = trend_stats(FF, data.model_years)
        sig = sigma_ff(FF)
        print(f"{label:<25} {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}  {sig:>6.1f}")

        results[label] = {
            'oh_d': oh_d,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
            'sigma_ff': sig,
        }

    print(f"\n{'=' * 60}")
    print("SIGN ANALYSIS")
    print(f"{'=' * 60}")
    all_negative = True
    for label, r in results.items():
        sign = "NEGATIVE" if r['trend_median'] < 0 else "POSITIVE"
        if r['trend_median'] >= 0:
            all_negative = False
        ci_zero = r['trend_5pct'] < 0 < r['trend_95pct']
        robust = "✓ robust" if not ci_zero else "⚠ CI includes zero"
        print(f"  {label:<25} ΔFF = {r['trend_median']:+.1f}  {sign:<9} {robust}")

    if all_negative:
        print("\n  ✓ FF trend reversal is ROBUST to OH-D KIE choice")
    else:
        print("\n  ⚠ FF trend depends on OH-D KIE — threshold exists")

    with open(OUT_DIR / "phase6_OHD_sensitivity.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase6_OHD_sensitivity.json'}")


if __name__ == "__main__":
    main()
