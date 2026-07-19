#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 13: uncertainty attribution and NH residual diagnostics.

DEPRECATED (superseded 2026-07-04) — do not use the attribution output for the
manuscript. The `attribution` block here is a prioritization-only diagnostic
built from HARDCODED, ASSUMED per-group sigma scales (see
build_default_perturbations: observation 0.0040, wetland_phasor 0.0030, etc.),
NOT from the actual Monte Carlo. It overstated the wetland contribution and
understated the observational contribution.

The real variance budget comes from analysis/sh_variance_purepy.py, which
decomposes the actual Phase-6 SH (CGO+SPO) Monte Carlo by toggling each
uncertainty group on/off. That true decomposition (observations ~88%, sink ~7%,
wetland ~5%) is what feeds the manuscript figures:
  - fig3 panel (c)                        (make_manuscript_figures.py)
  - figS10_uncertainty_attribution.png    (make_manuscript_figures.py)
  - figS11_sh_true_uncertainty_attribution.png (make_sh_uncertainty_figure.py)
Neither figure loads this script's JSON any longer.

The `nh_residual_summary` helper below is unaffected and may still be useful as
a standalone NH corrected-ratio diagnostic.
"""

from pathlib import Path
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPT_DIR = Path(__file__).resolve().parent.parent
PHASE6_JSON = EXPT_DIR / "results" / "phase6_phasor" / "phasor_results.json"
PHASE8_JSON = EXPT_DIR / "results" / "phase8_biomass_burning" / "bb_correction_results.json"
OUT_DIR = EXPT_DIR / "results" / "phase13_uncertainty_attribution"
FIG_DIR = EXPT_DIR / "figures"
OUT_JSON = OUT_DIR / "uncertainty_attribution_results.json"
OUT_FIG = FIG_DIR / "fig22_uncertainty_attribution.png"

ALPHA_D_OH = 1.294
R_CANTRELL = (1.0054 - 1.0) / (ALPHA_D_OH - 1.0)


def one_at_a_time_attribution(model, perturbations):
    """Estimate output variance from each perturbation group independently."""
    details = {}
    for name, values in perturbations.items():
        draws = {key: np.zeros_like(np.asarray(values, dtype=float)) for key, values in perturbations.items()}
        draws[name] = np.asarray(values, dtype=float)
        output = np.asarray(model(draws), dtype=float)
        variance = float(np.nanvar(output))
        details[name] = {
            "variance": variance,
            "std": float(np.sqrt(variance)),
        }
    total = sum(v["variance"] for v in details.values())
    for value in details.values():
        value["fraction_of_oat_variance"] = float(value["variance"] / total) if total > 0 else 0.0
    ranking = sorted(details, key=lambda k: details[k]["variance"], reverse=True)
    details["ranking"] = ranking
    details["total_oat_variance"] = float(total)
    return details


def nh_residual_summary(sites, oh_ratio_high=R_CANTRELL):
    """Summarize NH corrected-ratio excess above the Cantrell OH band edge."""
    nh_values = []
    sh_values = []
    for site in sites.values():
        band = site.get("source_band", site.get("band", ""))
        value = site.get("R_corrected", site.get("R_wetland_plus_bb"))
        if value is None:
            continue
        if str(band).startswith("NH") or band == "Tropics":
            nh_values.append(float(value))
        elif str(band).startswith("SH"):
            sh_values.append(float(value))
    nh = np.asarray(nh_values, dtype=float)
    sh = np.asarray(sh_values, dtype=float)
    return {
        "n_nh_sites": int(len(nh)),
        "n_sh_sites": int(len(sh)),
        "mean_nh_R": float(np.nanmean(nh)) if len(nh) else np.nan,
        "mean_sh_R": float(np.nanmean(sh)) if len(sh) else np.nan,
        "mean_excess_above_cantrell": float(np.nanmean(nh - oh_ratio_high)) if len(nh) else np.nan,
        "oh_ratio_high": float(oh_ratio_high),
    }


def _alpha_from_ratio(R, alpha_d_oh=ALPHA_D_OH):
    return 1.0 + R * (alpha_d_oh - 1.0)


def build_default_perturbations(n=2000, seed=123):
    """Return deterministic perturbation draws for attribution groups."""
    rng = np.random.default_rng(seed)
    return {
        "observation": rng.normal(0.0, 0.0040, n),
        "wetland_phasor": rng.normal(0.0, 0.0030, n),
        "wetland_isotopes": rng.normal(0.0, 0.0025, n),
        "bb_correction": rng.normal(0.0, 0.0015, n),
        "sink_fractions": rng.normal(0.0, 0.0020, n),
        "alpha_D_OH": rng.normal(0.0, 0.0018, n),
        "non_oh_kie": rng.normal(0.0, 0.0012, n),
    }


def run_attribution(base_alpha=1.0044, n=2000, seed=123):
    """Run one-at-a-time alpha perturbation attribution."""
    perturbations = build_default_perturbations(n=n, seed=seed)

    def model(draws):
        total = np.zeros(n, dtype=float)
        for values in draws.values():
            total = total + np.asarray(values, dtype=float)
        return base_alpha + total

    return one_at_a_time_attribution(model, perturbations)


def plot_attribution(attribution):
    ranking = attribution["ranking"]
    values = [attribution[name]["fraction_of_oat_variance"] for name in ranking]
    x = np.arange(len(ranking))
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(x, values)
    ax.set_xticks(x)
    ax.set_xticklabels(ranking, rotation=25, ha="right")
    ax.set_ylabel("Fraction of one-at-a-time variance")
    ax.set_title("Grouped alpha uncertainty attribution")
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=180)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    phase6 = json.loads(PHASE6_JSON.read_text(encoding="utf-8"))
    phase8 = json.loads(PHASE8_JSON.read_text(encoding="utf-8")) if PHASE8_JSON.exists() else {"sites": {}}
    attribution = run_attribution()
    output = {
        "metadata": {
            "method": "One-at-a-time grouped alpha uncertainty attribution and NH corrected-ratio residual summary",
            "note": "Perturbation scales are diagnostic defaults for prioritization, not a replacement for Phase 6 MC.",
        },
        "attribution": attribution,
        "nh_residual_phase6": nh_residual_summary(phase6["sites"]),
        "nh_residual_phase8": nh_residual_summary(phase8.get("sites", {})) if phase8.get("sites") else None,
        "alpha_from_cantrell_ratio": _alpha_from_ratio(R_CANTRELL),
    }
    OUT_JSON.write_text(json.dumps(output, indent=2), encoding="utf-8")
    plot_attribution(attribution)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
