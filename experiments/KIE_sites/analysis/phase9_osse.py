#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 9: synthetic OSSE checks for KIE_sites alpha recovery.

This phase tests the inversion logic on controlled monthly cycles with known
alpha_13C_OH. It is intentionally compact: the purpose is to expose when the
scalar amplitude-ratio workflow is unbiased, and when source seasonality pushes
the retrieved alpha high.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPT_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = EXPT_DIR / "results" / "phase9_osse"
FIG_DIR = EXPT_DIR / "figures"
OUT_JSON = OUT_DIR / "osse_results.json"
OUT_FIG = FIG_DIR / "fig18_osse_recovery.png"

ALPHA_D_OH = 1.294
D13C_SOURCE = -62.0
DD_SOURCE = -310.0
D13C_ATM = -47.3
DD_ATM = -86.0


@dataclass(frozen=True)
class Scenario:
    """Configuration for one synthetic monthly OSSE case."""

    name: str
    alpha_13c_oh: float
    sink_dD_amplitude: float
    source_fraction_amplitude: float
    noise_13c: float
    noise_dD: float
    n_years: int = 6
    source_phase_month: float = 7.0
    sink_phase_month: float = 7.0
    seed: int = 123


def _cos_cycle(months, peak_month, amplitude):
    """Return a cosine cycle with the requested 1-indexed peak month."""
    return amplitude * np.cos(2 * np.pi * (months - peak_month) / 12.0)


def invert_ratio_to_alpha(R_obs, alpha_d_oh=ALPHA_D_OH):
    """Invert pure-OH ratio R = (alpha_13C - 1) / (alpha_D - 1)."""
    return 1.0 + R_obs * (alpha_d_oh - 1.0)


def fit_annual_phasor(months, values):
    """Fit y = c0 + B sin(wm) + C cos(wm) and return amplitude and phasor."""
    x = 2 * np.pi * (np.asarray(months, dtype=float) - 0.5) / 12.0
    design = np.column_stack([np.ones_like(x), np.sin(x), np.cos(x)])
    coeffs, _, _, _ = np.linalg.lstsq(design, np.asarray(values, dtype=float), rcond=None)
    _, B, C = coeffs
    return {"B": float(B), "C": float(C), "amplitude": float(np.hypot(B, C))}


def synthesize_monthly_series(scenario):
    """Build synthetic paired monthly d13C and dD cycles for one scenario."""
    rng = np.random.default_rng(scenario.seed)
    months_one_year = np.arange(1, 13, dtype=float)
    months = np.tile(months_one_year, scenario.n_years)

    R_true = (scenario.alpha_13c_oh - 1.0) / (ALPHA_D_OH - 1.0)
    sink_dD = _cos_cycle(months, scenario.sink_phase_month, scenario.sink_dD_amplitude)
    sink_13c = _cos_cycle(months, scenario.sink_phase_month, scenario.sink_dD_amplitude * R_true)

    source_frac = _cos_cycle(months, scenario.source_phase_month, scenario.source_fraction_amplitude)
    source_13c = (D13C_SOURCE - D13C_ATM) * source_frac
    source_dD = (DD_SOURCE - DD_ATM) * source_frac

    d13c = sink_13c + source_13c
    dD = sink_dD + source_dD
    if scenario.noise_13c:
        d13c = d13c + rng.normal(0.0, scenario.noise_13c, size=d13c.shape)
    if scenario.noise_dD:
        dD = dD + rng.normal(0.0, scenario.noise_dD, size=dD.shape)

    return {
        "month": months,
        "d13C": d13c,
        "dD": dD,
        "sink_13c": sink_13c,
        "sink_dD": sink_dD,
        "source_13c": source_13c,
        "source_dD": source_dD,
    }


def run_osse_scenario(scenario):
    """Run one synthetic case and return true/retrieved alpha diagnostics."""
    data = synthesize_monthly_series(scenario)
    fit_13c = fit_annual_phasor(data["month"], data["d13C"])
    fit_dD = fit_annual_phasor(data["month"], data["dD"])
    R_observed = fit_13c["amplitude"] / fit_dD["amplitude"]
    alpha_retrieved = invert_ratio_to_alpha(R_observed)

    source_fit_13c = fit_annual_phasor(data["month"], data["source_13c"])
    source_fit_dD = fit_annual_phasor(data["month"], data["source_dD"])
    sink_fit_13c = fit_annual_phasor(data["month"], data["sink_13c"])
    sink_fit_dD = fit_annual_phasor(data["month"], data["sink_dD"])

    return {
        "scenario": asdict(scenario),
        "R_true": float((scenario.alpha_13c_oh - 1.0) / (ALPHA_D_OH - 1.0)),
        "R_observed": float(R_observed),
        "alpha_true": float(scenario.alpha_13c_oh),
        "alpha_retrieved": float(alpha_retrieved),
        "bias": float(alpha_retrieved - scenario.alpha_13c_oh),
        "A_obs_13c": fit_13c["amplitude"],
        "A_obs_dD": fit_dD["amplitude"],
        "A_source_13c": source_fit_13c["amplitude"],
        "A_source_dD": source_fit_dD["amplitude"],
        "A_sink_13c": sink_fit_13c["amplitude"],
        "A_sink_dD": sink_fit_dD["amplitude"],
    }


def default_scenarios():
    """Return the standard OSSE suite for this phase."""
    return [
        Scenario("pure_oh_saueressig", 1.0039, 5.0, 0.0, 0.0, 0.0),
        Scenario("pure_oh_cantrell", 1.0054, 5.0, 0.0, 0.0, 0.0),
        Scenario("wetland_contaminated", 1.0039, 5.0, 0.010, 0.0, 0.0),
        Scenario("sparse_like_noise", 1.0039, 5.0, 0.006, 0.03, 0.8, seed=456),
    ]


def plot_results(results):
    """Create a compact true-vs-retrieved alpha diagnostic."""
    labels = {
        "pure_oh_saueressig": "Pure OH\nSaueressig input",
        "pure_oh_cantrell": "Pure OH\nCantrell input",
        "wetland_contaminated": "Wetland source\nadded",
        "sparse_like_noise": "Wetland source\n+ sparse noise",
    }
    names = [labels.get(r["scenario"]["name"], r["scenario"]["name"]) for r in results]
    true_alpha = np.array([r["alpha_true"] for r in results])
    got_alpha = np.array([r["alpha_retrieved"] for r in results])
    x = np.arange(len(results))

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, true_alpha, "o", label="Input alpha")
    ax.plot(x, got_alpha, "s", label="Retrieved from amplitude ratio")
    ax.axhline(1.0039, color="C2", ls="--", lw=1, label="Saueressig")
    ax.axhline(1.0054, color="C3", ls="--", lw=1, label="Cantrell")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=18, ha="right")
    ax.set_xlabel("Synthetic scenario")
    ax.set_ylabel("alpha13C_OH")
    ax.set_title("OSSE recovery from scalar amplitude ratios")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=180)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    results = [run_osse_scenario(s) for s in default_scenarios()]
    output = {
        "metadata": {
            "method": "Synthetic monthly OSSE for scalar amplitude-ratio alpha recovery",
            "alpha_d_oh": ALPHA_D_OH,
        },
        "scenarios": results,
    }
    OUT_JSON.write_text(json.dumps(output, indent=2), encoding="utf-8")
    plot_results(results)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
