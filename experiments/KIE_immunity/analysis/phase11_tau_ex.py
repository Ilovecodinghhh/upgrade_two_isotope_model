#!/usr/bin/env python3
"""Phase 11 — Interhemispheric exchange time sensitivity."""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR
from common import smooth_5yr

def main():
    print("=" * 60)
    print("PHASE 11: τ_ex sensitivity (dual real-hemi-dD)")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400

    tau_ex_values = [
        ("Fast (0.5 yr)", 0.5),
        ("Default ~N(1.0, 0.1)", None),  # None = use normal sampling
        ("Fixed (1.0 yr)", 1.0),
        ("Slow (1.5 yr)", 1.5),
        ("Very slow (2.0 yr)", 2.0),
    ]

    results = {}
    print(f"\n{'Config':<25} {'ΔFF med':>8} {'[5%':>7} {'95%]':>7} {'σ(FF)':>7}")
    print("-" * 58)

    for label, tau_ex in tau_ex_values:
        FF = run_2box_flex(data, N, 42, tau_ex_fixed=tau_ex)
        med, lo, hi = trend_stats(FF, data.model_years)
        sig = sigma_ff(FF)

        # Also variance decomposition
        var_total = np.nanmean(np.nanvar(smooth_5yr(FF)[8:], axis=1))
        FF_fk = run_2box_flex(data, N, 42, tau_ex_fixed=tau_ex, fix_kie=True)
        var_fk = np.nanmean(np.nanvar(smooth_5yr(FF_fk)[8:], axis=1))
        kie_pct = max(0, (var_total - var_fk) / var_total * 100) if var_total > 0 else 0

        FF_fs = run_2box_flex(data, N, 42, tau_ex_fixed=tau_ex, fix_sigs=True)
        var_fs = np.nanmean(np.nanvar(smooth_5yr(FF_fs)[8:], axis=1))
        sig_pct = max(0, (var_total - var_fs) / var_total * 100) if var_total > 0 else 0

        print(f"{label:<25} {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}  {sig:>6.1f}")

        results[label] = {
            'tau_ex': tau_ex,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
            'sigma_ff': sig, 'kie_pct': kie_pct, 'sig_pct': sig_pct,
        }

    # Analysis
    print(f"\n{'=' * 60}")
    print("ANALYSIS: Does τ_ex affect σ(FF)?")
    print(f"{'=' * 60}")
    sigs = [(r['tau_ex'] or 1.0, r['sigma_ff']) for r in results.values()]
    sig_range = max(s[1] for s in sigs) - min(s[1] for s in sigs)
    sig_mean = np.mean([s[1] for s in sigs])
    print(f"  σ(FF) range: {min(s[1] for s in sigs):.1f} – {max(s[1] for s in sigs):.1f} Tg/yr")
    print(f"  σ(FF) spread: {sig_range:.1f} Tg/yr ({sig_range/sig_mean*100:.0f}% of mean)")

    if sig_range > 3:
        print(f"  → τ_ex MATTERS: hemispheric transport separation provides real constraint")
    else:
        print(f"  → τ_ex has MODEST effect: 2-box advantage is mainly from source-sig separation")

    with open(OUT_DIR / "phase11_tau_ex.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase11_tau_ex.json'}")


if __name__ == "__main__":
    main()
