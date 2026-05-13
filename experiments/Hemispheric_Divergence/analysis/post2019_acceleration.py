#!/usr/bin/env python3
"""
Phase 6.1: Post-2019 acceleration — hemispheric decomposition.

Tests: Is the 2020-2022 CH₄ growth surge dominated by SH-Mic
(tropical wetlands, La Niña)? Does NH-FF show a COVID-19 dip?
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import json

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


def main():
    print("=" * 70)
    print("PHASE 6.1: Post-2019 Hemispheric Decomposition")
    print("=" * 70)

    hemi_df = pd.read_csv(RESULTS_DIR / "twobox_hemi" / "hemispheric_detail.csv")
    hemi_npz = np.load(RESULTS_DIR / "twobox_hemi" / "all_iterations.npz")
    onebox_npz = np.load(RESULTS_DIR / "onebox_reference" / "all_iterations.npz")

    years = hemi_df['year'].values
    NH_FF = hemi_npz['NH_FF']
    NH_Mic = hemi_npz['NH_Mic']
    NH_BB = hemi_npz['NH_BB']
    SH_FF = hemi_npz['SH_FF']
    SH_Mic = hemi_npz['SH_Mic']
    SH_BB = hemi_npz['SH_BB']

    FF_1box = onebox_npz['FF']
    Mic_1box = onebox_npz['Mic']

    # Reference period: 2010-2018 mean
    ref_mask = (years >= 2010) & (years <= 2018)
    # Post-2019 period
    post_mask = years >= 2019

    if not np.any(post_mask):
        print("  ⚠ No post-2019 data in model output!")
        print(f"  Year range: {years[0]}–{years[-1]}")
        print("  Cannot perform post-2019 analysis.")
        return

    ref_years = years[ref_mask]
    post_years = years[post_mask]

    print(f"\n  Reference period: {ref_years[0]:.0f}–{ref_years[-1]:.0f}")
    print(f"  Post-2019 years: {post_years}")

    # Anomalies relative to 2010-2018 mean
    sources = {
        'NH_FF': NH_FF, 'NH_Mic': NH_Mic, 'NH_BB': NH_BB,
        'SH_FF': SH_FF, 'SH_Mic': SH_Mic, 'SH_BB': SH_BB,
    }

    print(f"\n  ── Anomalies vs 2010–2018 Mean (Tg/yr, median [90% CI]) ──")
    anomaly_results = {}

    for name, arr in sources.items():
        ref_mean = np.median(np.mean(arr[ref_mask, :], axis=0))
        for yi, yr in enumerate(post_years):
            idx = np.where(years == yr)[0][0]
            anom = arr[idx, :] - np.mean(arr[ref_mask, :], axis=0)
            med = np.median(anom)
            p5, p95 = np.percentile(anom, [5, 95])
            key = f"{name}_{int(yr)}"
            anomaly_results[key] = {"median": float(med), "p5": float(p5), "p95": float(p95)}

    # Print formatted
    for yr in post_years:
        yr_int = int(yr)
        print(f"\n  Year {yr_int}:")
        for name in sources:
            r = anomaly_results[f"{name}_{yr_int}"]
            flag = ""
            if abs(r['median']) > 2.0:
                flag = " ★"
            print(f"    {name:8s}: {r['median']:+6.1f} [{r['p5']:+6.1f}, {r['p95']:+6.1f}]{flag}")

    # COVID test: is 2020 NH_FF lower than 2019?
    if 2020 in post_years and 2019 in post_years:
        idx_2019 = np.where(years == 2019)[0][0]
        idx_2020 = np.where(years == 2020)[0][0]
        covid_dip = NH_FF[idx_2020, :] - NH_FF[idx_2019, :]
        pct_neg = np.mean(covid_dip < 0) * 100
        print(f"\n  ── COVID-19 Test ──")
        print(f"  NH_FF 2020 vs 2019: {np.median(covid_dip):+.1f} Tg/yr ({pct_neg:.0f}% negative)")
        if pct_neg > 60:
            print(f"  ✓ Suggestive COVID dip in NH FF emissions")
        else:
            print(f"  ✗ No clear COVID signal (may be within noise)")

    # La Niña test: is 2020-2022 SH_Mic elevated?
    if np.any(post_mask):
        sh_mic_post = np.mean(SH_Mic[post_mask, :], axis=0)
        sh_mic_ref = np.mean(SH_Mic[ref_mask, :], axis=0)
        enhancement = sh_mic_post - sh_mic_ref
        print(f"\n  ── La Niña Wetland Test ──")
        print(f"  SH_Mic post-2019 enhancement: {np.median(enhancement):+.1f} Tg/yr [{np.percentile(enhancement,5):+.1f}, {np.percentile(enhancement,95):+.1f}]")
        print(f"  Chandra 2024 tropical wetland surge: ~10-15 Tg/yr above baseline")
        if np.median(enhancement) > 3:
            print(f"  ✓ SH_Mic shows significant post-2019 enhancement")
        else:
            print(f"  ✗ SH_Mic enhancement is modest")

    # Total growth attribution
    print(f"\n  ── Growth Attribution (post-2019 vs 2010-2018) ──")
    total_growth = {}
    for name, arr in sources.items():
        ref_mean = np.mean(np.median(arr[ref_mask, :], axis=1))
        post_mean = np.mean(np.median(arr[post_mask, :], axis=1))
        total_growth[name] = post_mean - ref_mean

    total = sum(total_growth.values())
    for name, val in total_growth.items():
        pct = val / total * 100 if total != 0 else 0
        print(f"    {name:8s}: {val:+6.1f} Tg/yr ({pct:+5.1f}%)")
    print(f"    {'TOTAL':8s}: {total:+6.1f} Tg/yr")

    # Save
    with open(RESULTS_DIR / "post2019_analysis.json", 'w') as f:
        json.dump({
            "anomalies": anomaly_results,
            "growth_attribution": {k: float(v) for k, v in total_growth.items()},
            "total_growth": float(total),
        }, f, indent=2)
    print(f"\n  Saved to {RESULTS_DIR}/post2019_analysis.json")


if __name__ == "__main__":
    main()
