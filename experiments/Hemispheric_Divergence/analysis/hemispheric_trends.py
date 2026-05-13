#!/usr/bin/env python3
"""
Phase 2: Hemispheric trend analysis and significance testing.

Loads Phase 1 output and computes:
  - Linear trends (2007-2022) per iteration → slope distributions
  - Significance tests (fraction positive/negative)
  - 1-box vs 2-box comparison
  - Spatial aliasing quantification
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]  # repo root
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
import json

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results"


def load_twobox(label="hemi"):
    d = np.load(RESULTS_DIR / f"twobox_{label}" / "all_iterations.npz")
    return d


def load_onebox(label="reference"):
    d = np.load(RESULTS_DIR / f"onebox_{label}" / "all_iterations.npz")
    return d


def compute_trends(arr, years, start=2007, end=None):
    """
    For arr[n_years, n_iter], fit linear trend over [start, end] for each iteration.
    Returns slopes[n_iter] in Tg/yr/yr.
    """
    if end is None:
        end = years[-1]
    mask = (years >= start) & (years <= end)
    yrs = years[mask]
    sub = arr[mask, :]
    n_iter = sub.shape[1]
    slopes = np.zeros(n_iter)
    for k in range(n_iter):
        col = sub[:, k]
        if np.any(np.isnan(col)):
            slopes[k] = np.nan
            continue
        slope, _, _, _, _ = sp_stats.linregress(yrs, col)
        slopes[k] = slope
    return slopes


def compute_trend_change(arr, years, start=2007, end=None):
    """Cumulative change = slope × (end - start), in Tg/yr."""
    if end is None:
        end = years[-1]
    slopes = compute_trends(arr, years, start, end)
    return slopes * (end - start)


def trend_summary(slopes, label=""):
    valid = slopes[~np.isnan(slopes)]
    if len(valid) == 0:
        return {"label": label, "n_valid": 0}
    return {
        "label": label,
        "slope_median": float(np.median(valid)),
        "slope_mean": float(np.mean(valid)),
        "slope_p5": float(np.percentile(valid, 5)),
        "slope_p95": float(np.percentile(valid, 95)),
        "slope_std": float(np.std(valid)),
        "p_positive": float(np.mean(valid > 0) * 100),
        "p_negative": float(np.mean(valid < 0) * 100),
        "significant_pos": bool(np.percentile(valid, 5) > 0),
        "significant_neg": bool(np.percentile(valid, 95) < 0),
        "n_valid": int(len(valid)),
    }


def main():
    print("=" * 70)
    print("PHASE 2: Hemispheric Trend Analysis")
    print("=" * 70)

    # Load 2-box hemispheric
    d2 = load_twobox("hemi")
    years = d2['years']

    # Load 1-box
    d1 = load_onebox("reference")
    years_1box = d1['years']

    # Load 2-box with global sigs
    d2g = load_twobox("global_sigs")

    results = []

    print("\n--- 2-box hemispheric source sigs ---")
    for src in ['FF', 'Mic', 'BB']:
        for hemi in ['NH', 'SH']:
            key = f'{hemi}_{src}'
            arr = d2[key]
            slopes = compute_trends(arr, years, start=2007)
            s = trend_summary(slopes, f"2box_hemi_{key}")
            results.append(s)
            sig = "✓ SIG" if s.get('significant_pos') or s.get('significant_neg') else ""
            print(f"  {key:12s}: slope = {s['slope_median']:+.2f} [{s['slope_p5']:+.2f}, {s['slope_p95']:+.2f}] Tg/yr²  "
                  f"({s['p_positive']:.0f}% pos) {sig}")

        # Global aggregate from 2-box
        glob = d2[f'NH_{src}'] + d2[f'SH_{src}']
        slopes = compute_trends(glob, years, start=2007)
        s = trend_summary(slopes, f"2box_hemi_Global_{src}")
        results.append(s)
        sig = "✓ SIG" if s.get('significant_pos') or s.get('significant_neg') else ""
        print(f"  {'Global_'+src:12s}: slope = {s['slope_median']:+.2f} [{s['slope_p5']:+.2f}, {s['slope_p95']:+.2f}] Tg/yr²  "
              f"({s['p_positive']:.0f}% pos) {sig}")

    print("\n--- 2-box global source sigs ---")
    for src in ['FF', 'Mic', 'BB']:
        for hemi in ['NH', 'SH']:
            key = f'{hemi}_{src}'
            arr = d2g[key]
            slopes = compute_trends(arr, years, start=2007)
            s = trend_summary(slopes, f"2box_global_{key}")
            results.append(s)
            sig = "✓ SIG" if s.get('significant_pos') or s.get('significant_neg') else ""
            print(f"  {key:12s}: slope = {s['slope_median']:+.2f} [{s['slope_p5']:+.2f}, {s['slope_p95']:+.2f}] Tg/yr²  "
                  f"({s['p_positive']:.0f}% pos) {sig}")

    print("\n--- 1-box reference ---")
    for src in ['FF', 'Mic', 'BB']:
        arr = d1[src]
        n_yrs = min(len(years_1box), arr.shape[0])
        slopes = compute_trends(arr[:n_yrs], years_1box[:n_yrs], start=2007)
        s = trend_summary(slopes, f"1box_{src}")
        results.append(s)
        sig = "✓ SIG" if s.get('significant_pos') or s.get('significant_neg') else ""
        print(f"  {src:12s}: slope = {s['slope_median']:+.2f} [{s['slope_p5']:+.2f}, {s['slope_p95']:+.2f}] Tg/yr²  "
              f"({s['p_positive']:.0f}% pos) {sig}")

    # ---- Spatial aliasing test ----
    print("\n" + "=" * 70)
    print("SPATIAL ALIASING TEST")
    print("=" * 70)
    NH_FF = d2['NH_FF']; SH_FF = d2['SH_FF']
    glob_FF = NH_FF + SH_FF

    nh_slopes = compute_trends(NH_FF, years, start=2007)
    sh_slopes = compute_trends(SH_FF, years, start=2007)
    glob_slopes = compute_trends(glob_FF, years, start=2007)

    print(f"  NH_FF trend: {np.median(nh_slopes):+.3f} Tg/yr² ({np.mean(nh_slopes > 0)*100:.0f}% positive)")
    print(f"  SH_FF trend: {np.median(sh_slopes):+.3f} Tg/yr² ({np.mean(sh_slopes > 0)*100:.0f}% positive)")
    print(f"  Global FF trend (NH+SH): {np.median(glob_slopes):+.3f} Tg/yr²")

    # 1-box FF
    ff_1box = d1['FF']
    n_yrs = min(len(years_1box), ff_1box.shape[0])
    onebox_slopes = compute_trends(ff_1box[:n_yrs], years_1box[:n_yrs], start=2007)
    print(f"  1-box FF trend: {np.median(onebox_slopes):+.3f} Tg/yr²")

    nh_pos = np.mean(nh_slopes > 0) * 100
    glob_neg = np.mean(glob_slopes < 0) * 100
    aliasing = nh_pos > 60 and glob_neg > 40
    print(f"\n  Aliasing detected: NH positive ({nh_pos:.0f}%) but Global ~zero/negative ({glob_neg:.0f}% neg)")
    print(f"  → {'YES — spatial aliasing confirmed!' if aliasing else 'Inconclusive'}")

    # ---- Microbial divergence ----
    NH_Mic = d2['NH_Mic']; SH_Mic = d2['SH_Mic']
    nh_mic_slopes = compute_trends(NH_Mic, years, start=2007)
    sh_mic_slopes = compute_trends(SH_Mic, years, start=2007)
    print(f"\n  NH_Mic trend: {np.median(nh_mic_slopes):+.3f} Tg/yr²")
    print(f"  SH_Mic trend: {np.median(sh_mic_slopes):+.3f} Tg/yr²")
    print(f"  SH_Mic positive: {np.mean(sh_mic_slopes > 0)*100:.0f}%")

    # Save
    out_file = RESULTS_DIR / "trend_analysis.json"
    with open(out_file, 'w') as f:
        json.dump(results, f, indent=2)

    # CSV summary
    df = pd.DataFrame(results)
    df.to_csv(RESULTS_DIR / "trend_analysis.csv", index=False)
    print(f"\n  Saved to {RESULTS_DIR}/trend_analysis.csv")

    # ---- Also save aliasing summary ----
    aliasing_summary = {
        "NH_FF_slope_median": float(np.median(nh_slopes)),
        "SH_FF_slope_median": float(np.median(sh_slopes)),
        "Global_FF_slope_median": float(np.median(glob_slopes)),
        "Onebox_FF_slope_median": float(np.median(onebox_slopes)),
        "NH_FF_pct_positive": float(nh_pos),
        "Global_FF_pct_negative": float(glob_neg),
        "SH_Mic_slope_median": float(np.median(sh_mic_slopes)),
        "SH_Mic_pct_positive": float(np.mean(sh_mic_slopes > 0)*100),
        "aliasing_detected": bool(aliasing),
        "hypothesis_supported": bool(aliasing and np.mean(sh_mic_slopes > 0)*100 > 60),
        # Additional: 1-box vs 2-box divergence
        "onebox_vs_twobox_FF_divergence": float(np.median(glob_slopes) - np.median(onebox_slopes)),
    }
    with open(RESULTS_DIR / "aliasing_test.json", 'w') as f:
        json.dump(aliasing_summary, f, indent=2)
    print(f"  Saved aliasing test to {RESULTS_DIR}/aliasing_test.json")


if __name__ == "__main__":
    main()
