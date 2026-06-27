import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase9_osse import (  # noqa: E402
    Scenario,
    invert_ratio_to_alpha,
    run_osse_scenario,
)


def test_pure_oh_osse_recovers_input_alpha():
    scenario = Scenario(
        name="pure",
        alpha_13c_oh=1.0039,
        sink_dD_amplitude=5.0,
        source_fraction_amplitude=0.0,
        noise_13c=0.0,
        noise_dD=0.0,
        n_years=6,
    )

    result = run_osse_scenario(scenario)

    assert abs(result["alpha_retrieved"] - 1.0039) < 1e-6
    assert result["bias"] == result["alpha_retrieved"] - scenario.alpha_13c_oh


def test_wetland_contamination_biases_uncorrected_ratio_high():
    clean = Scenario(
        name="clean",
        alpha_13c_oh=1.0039,
        sink_dD_amplitude=5.0,
        source_fraction_amplitude=0.0,
        noise_13c=0.0,
        noise_dD=0.0,
        n_years=6,
    )
    wetland = Scenario(
        name="wetland",
        alpha_13c_oh=1.0039,
        sink_dD_amplitude=5.0,
        source_fraction_amplitude=0.01,
        noise_13c=0.0,
        noise_dD=0.0,
        n_years=6,
    )

    clean_result = run_osse_scenario(clean)
    wetland_result = run_osse_scenario(wetland)

    assert wetland_result["R_observed"] > clean_result["R_observed"]
    assert wetland_result["alpha_retrieved"] > clean_result["alpha_retrieved"]


def test_ratio_inversion_round_trips_alpha_values():
    for alpha in np.array([1.0039, 1.0054, 1.011]):
        ratio = (alpha - 1.0) / (1.294 - 1.0)

        assert np.isclose(invert_ratio_to_alpha(ratio, alpha_d_oh=1.294), alpha)
