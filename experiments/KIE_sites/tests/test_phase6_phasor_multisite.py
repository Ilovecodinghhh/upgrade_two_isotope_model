import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase6_phasor import (  # noqa: E402
    combine_ratio_samples_by_site,
    ratio_to_alpha_13c,
    summarize_multisite_ratio_constraint,
)


def test_combine_ratio_samples_averages_sites_samplewise():
    ratio_samples = {
        "CGO": np.array([0.030, 0.034, 0.032]),
        "SPO": np.array([0.020, 0.024, 0.022]),
    }
    sigmas = {"CGO": 1.0, "SPO": 1.0}

    result = combine_ratio_samples_by_site(["CGO", "SPO"], ratio_samples, sigmas)

    assert np.allclose(result["R_weighted_samples"], [0.025, 0.029, 0.027])
    assert np.isclose(result["R_weighted_mean"], 0.027)


def test_multisite_alpha_uses_combined_ratio_not_pooled_site_samples():
    ratio_samples = {
        "CGO": np.full(10, 0.033),
        "SPO": np.full(10, 0.023),
    }
    sigmas = {"CGO": 1.0, "SPO": 1.0}

    result = summarize_multisite_ratio_constraint(
        ["CGO", "SPO"],
        ratio_samples,
        sigmas,
        include_sink_uncertainty=False,
    )

    expected = ratio_to_alpha_13c(0.028)

    assert np.isclose(result["R_weighted_mean"], 0.028)
    assert np.isclose(result["alpha_13c_oh_median"], expected)
