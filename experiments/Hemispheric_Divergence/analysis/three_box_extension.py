#!/usr/bin/env python3
"""
Phase 6.2: Three-box extension — mathematical framework and feasibility.

This is conceptual/theoretical: defines the 3-box system, identifies data
requirements, and assesses feasibility for future work.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import json

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


def main():
    print("=" * 70)
    print("PHASE 6.2: Three-Box Extension Framework")
    print("=" * 70)

    framework = {
        "title": "3-Box Dual-Isotope Methane Source Attribution",
        "boxes": {
            "NHext": {"lat_range": "30°N–90°N", "mass_fraction": 0.25,
                      "dominant_sources": "FF (industry, gas), boreal wetlands, ruminants",
                      "monitoring": "NOAA flask network, AGAGE, Greenland stations"},
            "Trop": {"lat_range": "30°S–30°N", "mass_fraction": 0.50,
                     "dominant_sources": "tropical wetlands, rice, livestock, tropical BB",
                     "monitoring": "Sparse! Mauna Loa, Cape Grim edge, Ascension Island"},
            "SHext": {"lat_range": "90°S–30°S", "mass_fraction": 0.25,
                      "dominant_sources": "minimal local sources, clean background",
                      "monitoring": "South Pole, Cape Grim"},
        },
        "exchange_rates": {
            "NHext_Trop": {"tau_yr": 0.8, "note": "Faster (Hadley cell + eddies)"},
            "Trop_SHext": {"tau_yr": 1.2, "note": "Slower (ITCZ barrier)"},
            "NHext_SHext": {"tau_yr": "∞", "note": "No direct exchange (via Trop only)"},
        },
        "equations": """
        Mass balance for each box i ∈ {NHext, Trop, SHext}:

        dM_i/dt = S_i - M_i/τ_i + Σ_j (M_j - M_i)/τ_{ij}

        where:
          M_i = CH₄ mass in box i
          S_i = total source in box i (= BB_i + FF_i + Mic_i)
          τ_i = chemical lifetime in box i
          τ_{ij} = exchange time between boxes i and j

        Isotope mass balance (for ¹³CH₄):
          d(R_i·M_i)/dt = R_src_i·S_i - α_i·R_i·M_i/τ_i + Σ_j (R_j·M_j - R_i·M_i)/τ_{ij}

        Same for CH₃D (δD).

        System: 9 unknowns (3 sources × 3 boxes)
                9 equations (3 mass + 3 δ¹³C + 3 δD)
        
        → Exactly determined! (vs. 2-box: 6 unknowns, 6 equations)
        → Each box needs: CH₄ concentration, δ¹³C, δD (3 observables)
        """,
        "data_requirements": {
            "CH4_zonal": {
                "status": "AVAILABLE",
                "source": "NOAA GML flask + in-situ, zonal means available",
                "uncertainty": "±2 ppb for annual zonal means"
            },
            "d13C_zonal": {
                "status": "PARTIALLY AVAILABLE",
                "source": "NOAA INSTAAR, sparse in tropics",
                "issue": "Tropical δ¹³C is interpolated from NH/SH stations, not independent",
                "uncertainty": "±0.05‰ at well-sampled sites, ±0.2‰ in tropics"
            },
            "dD_zonal": {
                "status": "NOT AVAILABLE for 3 zones",
                "source": "Only a few sites globally (ALT, MLO, CGO, SPO)",
                "issue": "Cannot independently constrain 3 zonal δD values",
                "uncertainty": "±3‰ per station, but too few stations for tropical mean"
            },
            "source_signatures_zonal": {
                "status": "PARTIALLY AVAILABLE",
                "source": "Regional inventories, process studies",
                "issue": "Tropical source sigs less well-characterized than NH"
            }
        },
        "feasibility_assessment": {
            "verdict": "NOT YET FEASIBLE for full implementation",
            "limiting_factor": "Tropical δD observations",
            "timeline": "Feasible within 5-10 years if IRMS/TILDAS network expands",
            "workaround": "Could use satellite-derived CH₄ + δ¹³C (GOSAT) with modeled δD priors",
            "advantage_over_2box": (
                "Separates tropical wetland signal from NH industrial + SH background. "
                "Current 2-box lumps tropical into NH+SH, losing the dominant signal."
            ),
        },
        "key_question": (
            "Would the 3-box resolve the tropical wetland surge that drives post-2019 "
            "acceleration? The 2-box lumps tropical sources into NH and SH equally, "
            "but most tropical wetland emissions are near-equatorial. A 3-box could "
            "independently estimate the tropical microbial trend."
        ),
    }

    # Save
    with open(RESULTS_DIR / "three_box_framework.json", 'w') as f:
        json.dump(framework, f, indent=2)

    # Print summary
    print(f"\n  ── 3-Box System ──")
    for box, info in framework['boxes'].items():
        print(f"  {box} ({info['lat_range']}): {info['dominant_sources']}")

    print(f"\n  ── Exchange Rates ──")
    for pair, info in framework['exchange_rates'].items():
        print(f"  {pair}: τ = {info['tau_yr']} yr — {info['note']}")

    print(f"\n  ── Data Availability ──")
    for obs, info in framework['data_requirements'].items():
        print(f"  {obs}: {info['status']}")

    print(f"\n  ── Feasibility ──")
    fa = framework['feasibility_assessment']
    print(f"  Verdict: {fa['verdict']}")
    print(f"  Limiting factor: {fa['limiting_factor']}")
    print(f"  Timeline: {fa['timeline']}")
    print(f"  Advantage: {fa['advantage_over_2box']}")

    print(f"\n  Saved to {RESULTS_DIR}/three_box_framework.json")


if __name__ == "__main__":
    main()
