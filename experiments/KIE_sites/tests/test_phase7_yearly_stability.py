import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase7_yearly_stability import (  # noqa: E402
    MIN_MONTHS,
    classify_stability,
    fit_site_year,
)


def _monthly_frame(n_months):
    months = np.arange(1, n_months + 1)
    phase = 2 * np.pi * (months - 0.5) / 12.0
    return pd.DataFrame(
        {
            "year": np.full(n_months, 2006),
            "month": months,
            "decimal_year": 2006 + (months - 0.5) / 12.0,
            "d13C_mean": 0.1 * np.sin(phase),
            "dD_mean": 2.0 * np.sin(phase),
        }
    )


def test_site_year_requires_at_least_eight_paired_months():
    assert MIN_MONTHS == 8
    assert fit_site_year("TST", 2006, _monthly_frame(7)) is None
    assert fit_site_year("TST", 2006, _monthly_frame(8)) is not None


def test_stability_classification_requires_two_usable_years():
    assert classify_stability(1, np.nan, np.nan) == "insufficient_years"
    assert classify_stability(2, 0.5, 1.0) == "yearly_ratios_broadly_consistent"
