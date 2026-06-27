import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase13_uncertainty_attribution import (  # noqa: E402
    nh_residual_summary,
    one_at_a_time_attribution,
)


def test_one_at_a_time_attribution_identifies_only_active_source():
    def model(draws):
        return 1.0 + draws["active"]

    perturbations = {
        "inactive": np.zeros(100),
        "active": np.linspace(-1.0, 1.0, 100),
    }

    result = one_at_a_time_attribution(model, perturbations)

    assert result["active"]["variance"] > 0
    assert result["inactive"]["variance"] == 0
    assert result["ranking"][0] == "active"


def test_nh_residual_summary_reports_positive_excess():
    sites = {
        "ALT": {"source_band": "NH_high", "R_corrected": 0.05},
        "CGO": {"source_band": "SH_extra", "R_corrected": 0.03},
    }

    result = nh_residual_summary(sites, oh_ratio_high=0.0174)

    assert result["n_nh_sites"] == 1
    assert result["mean_excess_above_cantrell"] > 0
