#!/usr/bin/env python3
"""
Phase 4.2: Predicted NH-SH δD gradient — a novel testable prediction.

Uses the 2-box model's source partitioning + sink KIE to predict
what the interhemispheric δD-CH₄ gradient should be, given the
hemispheric source mix. This gradient has never been observed because
δD monitoring is too sparse for hemispheric resolution.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import json

from common import (
    load_data, KIE_FIXED,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


def predict_dD_gradient(data):
    """
    Estimate the expected NH-SH δD gradient from source partitioning.
    
    Logic: At steady state for a 2-box model, the atmospheric δD in each
    hemisphere reflects the balance of source δD and sink fractionation.
    
    δD_atm ≈ δD_source - ε_sink  (simplified)
    
    where ε_sink = (α_sink - 1) × 1000
    
    The NH-SH gradient comes from:
    1. Different source mixes (more FF in NH → higher δD_source_NH)
    2. Different sink fractionation (more Cl in NH → larger ε_D_NH)
    """
    
    # Source δD values (typical, from common.py source sig ranges)
    dD_FF = -180.0   # ‰ (thermogenic: -150 to -200)
    dD_Mic = -310.0  # ‰ (microbial: -280 to -340)
    dD_BB = -215.0   # ‰ (biomass burning: -200 to -230)
    
    # Hemispheric source fractions from our 2-box model (median values)
    # From the trend_analysis results, we know approximate source mix
    hemi_npz = np.load(RESULTS_DIR / "twobox_hemi" / "all_iterations.npz")
    NH_FF = hemi_npz['NH_FF']
    NH_Mic = hemi_npz['NH_Mic']
    NH_BB = hemi_npz['NH_BB']
    SH_FF = hemi_npz['SH_FF']
    SH_Mic = hemi_npz['SH_Mic']
    SH_BB = hemi_npz['SH_BB']
    
    n_years = NH_FF.shape[0]
    
    # Source-weighted δD for each hemisphere
    NH_total = NH_FF + NH_Mic + NH_BB
    SH_total = SH_FF + SH_Mic + SH_BB
    
    dD_source_NH = (NH_FF * dD_FF + NH_Mic * dD_Mic + NH_BB * dD_BB) / NH_total
    dD_source_SH = (SH_FF * dD_FF + SH_Mic * dD_Mic + SH_BB * dD_BB) / SH_total
    
    # Sink fractionation (ε = (α-1) × 1000)
    # Bulk KIE for δD
    OH_D = KIE_FIXED['OH_D']      # 1.3105
    Cl_D = KIE_FIXED['Cl_D']      # 1.52
    Strat_D = KIE_FIXED['Strat_D']  # 1.179
    Soil_D = KIE_FIXED['Soil_D']    # 1.083
    
    # NH sink ε_D
    sf_nh = SINK_FRACTIONS_NH
    bulk_alpha_NH = (sf_nh['OH']*OH_D + sf_nh['Cl']*Cl_D + 
                     sf_nh['Strat']*Strat_D + sf_nh['Soil']*Soil_D)
    eps_D_NH = (bulk_alpha_NH - 1) * 1000  # ‰
    
    sf_sh = SINK_FRACTIONS_SH
    bulk_alpha_SH = (sf_sh['OH']*OH_D + sf_sh['Cl']*Cl_D + 
                     sf_sh['Strat']*Strat_D + sf_sh['Soil']*Soil_D)
    eps_D_SH = (bulk_alpha_SH - 1) * 1000  # ‰
    
    print(f"  Sink fractionation ε_D:")
    print(f"    NH: {eps_D_NH:.1f}‰  (α = {bulk_alpha_NH:.4f})")
    print(f"    SH: {eps_D_SH:.1f}‰  (α = {bulk_alpha_SH:.4f})")
    print(f"    Difference (NH-SH): {eps_D_NH - eps_D_SH:.1f}‰")
    
    # Source δD gradient (NH - SH) per iteration, per year
    dD_source_gradient = dD_source_NH - dD_source_SH  # shape (n_years, n_iter)
    
    # Expected atmospheric gradient:
    # δD_atm ≈ δD_source - ε_sink
    # ΔδD_atm = (δD_source_NH - ε_NH) - (δD_source_SH - ε_SH)
    #         = (δD_source_NH - δD_source_SH) - (ε_NH - ε_SH)
    dD_atm_gradient = dD_source_gradient - (eps_D_NH - eps_D_SH)
    
    # Summary stats
    med_gradient = np.median(np.median(dD_atm_gradient, axis=1))
    p5 = np.percentile(np.median(dD_atm_gradient, axis=1), 5)
    p95 = np.percentile(np.median(dD_atm_gradient, axis=1), 95)
    
    print(f"\n  ── Predicted NH-SH δD Gradient ──")
    print(f"  Source δD gradient (NH-SH): {np.median(np.median(dD_source_gradient, axis=1)):+.1f}‰")
    print(f"  Sink ε correction:          {-(eps_D_NH - eps_D_SH):+.1f}‰")
    print(f"  Net atmospheric gradient:   {med_gradient:+.1f}‰ [{p5:+.1f}, {p95:+.1f}]")
    
    # Time series
    med_by_year = np.median(dD_atm_gradient, axis=1)
    p5_by_year = np.percentile(dD_atm_gradient, 5, axis=1)
    p95_by_year = np.percentile(dD_atm_gradient, 95, axis=1)
    
    # Load actual dD data used in model
    dD_NH = data.dD_NH if hasattr(data, 'dD_NH') else None
    dD_SH = data.dD_SH if hasattr(data, 'dD_SH') else None
    
    # Our model uses dD_offset = ±3‰ (from inputs.py)
    print(f"\n  Model input assumption: ±3‰ NH-SH offset (6‰ gradient)")
    print(f"  Model prediction:       {med_gradient:+.1f}‰ gradient")
    
    if abs(med_gradient - 6.0) > 3.0:
        print(f"  ⚠ Predicted gradient differs substantially from input assumption!")
        print(f"    This self-inconsistency could be important for future work.")
    else:
        print(f"  ✓ Predicted gradient is broadly consistent with input assumption.")
    
    # Trend in gradient
    from scipy import stats as sp_stats
    hemi_df = pd.read_csv(RESULTS_DIR / "twobox_hemi" / "hemispheric_detail.csv")
    years = hemi_df['year'].values
    grad_trends = np.array([sp_stats.linregress(years, dD_atm_gradient[:, k]).slope 
                            for k in range(dD_atm_gradient.shape[1])])
    
    print(f"\n  Gradient trend: {np.median(grad_trends):+.3f} ‰/yr [{np.percentile(grad_trends,5):+.3f}, {np.percentile(grad_trends,95):+.3f}]")
    if np.median(grad_trends) > 0:
        print(f"  → NH-SH δD gradient is WIDENING (NH becoming more FF-enriched)")
    else:
        print(f"  → NH-SH δD gradient is narrowing")
    
    # Save
    results = {
        "predicted_dD_gradient_median": float(med_gradient),
        "predicted_dD_gradient_p5": float(p5),
        "predicted_dD_gradient_p95": float(p95),
        "source_dD_gradient": float(np.median(np.median(dD_source_gradient, axis=1))),
        "sink_epsilon_NH": float(eps_D_NH),
        "sink_epsilon_SH": float(eps_D_SH),
        "gradient_trend_per_yr": float(np.median(grad_trends)),
        "model_input_offset": 6.0,
        "testable_prediction": (
            f"The NH-SH δD-CH₄ gradient should be {med_gradient:+.0f}‰ "
            f"[{p5:+.0f}, {p95:+.0f}‰, 90% CI]. "
            f"This can be tested as the IRMS/TILDAS monitoring network expands."
        ),
    }
    
    with open(RESULTS_DIR / "dD_gradient_prediction.json", 'w') as f:
        json.dump(results, f, indent=2)
    
    ts_df = pd.DataFrame({
        'year': years,
        'dD_gradient_median': med_by_year,
        'dD_gradient_p5': p5_by_year,
        'dD_gradient_p95': p95_by_year,
    })
    ts_df.to_csv(RESULTS_DIR / "dD_gradient_timeseries.csv", index=False)
    
    print(f"\n  Saved to {RESULTS_DIR}/dD_gradient_prediction.json")


def main():
    print("=" * 70)
    print("PHASE 4.2: δD Gradient Prediction")
    print("=" * 70)
    
    data = load_data(ROOT, two_box=True)
    predict_dD_gradient(data)


if __name__ == "__main__":
    main()
