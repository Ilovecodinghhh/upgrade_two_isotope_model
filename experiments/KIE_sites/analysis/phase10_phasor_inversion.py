#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 10: direct phasor-space inversion diagnostics for KIE_sites."""

from pathlib import Path
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPT_DIR = Path(__file__).resolve().parent.parent
PHASE6_JSON = EXPT_DIR / "results" / "phase6_phasor" / "phasor_results.json"
PHASE8_JSON = EXPT_DIR / "results" / "phase8_biomass_burning" / "bb_correction_results.json"
OUT_DIR = EXPT_DIR / "results" / "phase10_phasor_inversion"
FIG_DIR = EXPT_DIR / "figures"
OUT_JSON = OUT_DIR / "phasor_inversion_results.json"
OUT_FIG = FIG_DIR / "fig19_phasor_inversion_diagnostics.png"

ALPHA_D_OH = 1.294


def ratio_to_alpha_13c(R_sink, alpha_d_oh=ALPHA_D_OH):
    """Convert sink phasor amplitude ratio to pure-OH alpha."""
    return 1.0 + R_sink * (alpha_d_oh - 1.0)


def invert_sink_phasors(sink_13c, sink_dD):
    """Invert complex sink phasors into R_sink and alpha diagnostics."""
    A13 = abs(sink_13c)
    AD = abs(sink_dD)
    R_sink = A13 / AD if AD > 0 else np.nan
    alpha = ratio_to_alpha_13c(R_sink) if np.isfinite(R_sink) else np.nan
    phase_diff = np.angle(sink_13c) - np.angle(sink_dD)
    phase_diff = (phase_diff + np.pi) % (2 * np.pi) - np.pi
    return {
        "A_sink_13c": float(A13),
        "A_sink_dD": float(AD),
        "R_sink": float(R_sink),
        "alpha_13c_oh": float(alpha),
        "phase_diff_rad": float(phase_diff),
        "phase_diff_months": float(abs(phase_diff) / (2 * np.pi) * 12.0),
    }


def phasor_residual(observed, source, sink):
    """Return observed - source - sink as a complex residual."""
    return observed - source - sink


def _complex_from_pair(pair):
    return complex(float(pair[0]), float(pair[1]))


def analyze_phase6_site(site):
    """Compute direct phasor inversion diagnostics for one Phase 6 site."""
    observed_13c = _complex_from_pair(site["Z_obs_13c"])
    observed_dD = _complex_from_pair(site["Z_obs_dD"])
    source_13c = _complex_from_pair(site["Z_src_13c"])
    source_dD = _complex_from_pair(site["Z_src_dD"])
    sink_13c = _complex_from_pair(site["Z_sink_13c"])
    sink_dD = _complex_from_pair(site["Z_sink_dD"])
    inv = invert_sink_phasors(sink_13c, sink_dD)
    res_13c = phasor_residual(observed_13c, source_13c, sink_13c)
    res_dD = phasor_residual(observed_dD, source_dD, sink_dD)
    inv.update(
        {
            "residual_13c_abs": float(abs(res_13c)),
            "residual_dD_abs": float(abs(res_dD)),
            "R_observed": float(site["R_obs"]),
            "R_phase6_corrected": float(site["R_corrected"]),
            "source_band": site["source_band"],
        }
    )
    return inv


def analyze_phase8_site(site):
    """Compute direct inversion after wetland+BB source subtraction."""
    sink_13c = _complex_from_pair(site["Z_sink_wetland_plus_bb_13c"])
    sink_dD = _complex_from_pair(site["Z_sink_wetland_plus_bb_dD"])
    inv = invert_sink_phasors(sink_13c, sink_dD)
    inv.update(
        {
            "R_observed": float(site["R_obs"]),
            "R_wetland_plus_bb": float(site["R_wetland_plus_bb"]),
            "source_band": site["band"],
        }
    )
    return inv


def plot_alpha_by_site(results):
    codes = list(results["phase6_wetland_only"].keys())
    x = np.arange(len(codes))
    phase6 = [results["phase6_wetland_only"][c]["alpha_13c_oh"] for c in codes]
    phase8 = [
        results["phase8_wetland_plus_bb"][c]["alpha_13c_oh"]
        for c in codes
        if c in results["phase8_wetland_plus_bb"]
    ]

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, phase6, "o-", label="Wetland-only correction")
    if len(phase8) == len(codes):
        ax.plot(x, phase8, "s-", label="Wetland + biomass burning")
    ax.axhline(1.0039, color="C2", ls="--", lw=1, label="Saueressig")
    ax.axhline(1.0054, color="C3", ls="--", lw=1, label="Cantrell")
    ax.set_xticks(x)
    ax.set_xticklabels(codes)
    ax.set_ylabel("alpha13C_OH inferred from corrected phasors")
    ax.set_title("Direct phasor inversion diagnostics")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=180)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    phase6 = json.loads(PHASE6_JSON.read_text(encoding="utf-8"))
    phase8 = json.loads(PHASE8_JSON.read_text(encoding="utf-8")) if PHASE8_JSON.exists() else {"sites": {}}
    results = {
        "metadata": {
            "method": "Direct complex-phasor inversion from existing corrected sink phasors",
            "alpha_d_oh": ALPHA_D_OH,
        },
        "phase6_wetland_only": {
            code: analyze_phase6_site(site) for code, site in phase6["sites"].items()
        },
        "phase8_wetland_plus_bb": {
            code: analyze_phase8_site(site) for code, site in phase8.get("sites", {}).items()
        },
    }
    OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    plot_alpha_by_site(results)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
