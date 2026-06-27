import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase11_block_bootstrap import (  # noqa: E402
    block_bootstrap_years,
    fit_resampled_site,
)


def test_block_bootstrap_years_is_reproducible_with_seed():
    years = np.array([2005, 2006, 2007, 2008])

    first = block_bootstrap_years(years, n_boot=3, seed=123)
    second = block_bootstrap_years(years, n_boot=3, seed=123)

    assert first == second
    assert all(len(draw) == len(years) for draw in first)


def test_fit_resampled_site_preserves_known_ratio():
    rows = []
    for year in [2005, 2006, 2007]:
        for month in range(1, 13):
            phase = 2 * np.pi * (month - 0.5) / 12.0
            rows.append(
                {
                    "year": year,
                    "month": month,
                    "decimal_year": year + (month - 0.5) / 12.0,
                    "d13C_mean": 0.2 * np.sin(phase),
                    "dD_mean": 4.0 * np.sin(phase),
                }
            )
    df = pd.DataFrame(rows)

    result = fit_resampled_site(df, [2005, 2006, 2007])

    assert np.isclose(result["ratio"], 0.05)
