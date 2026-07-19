import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase12_harmonic_sensitivity import (  # noqa: E402
    fit_annual,
    fit_annual_plus_semiannual,
    fit_monthly_fixed_effect,
    json_safe,
    leave_one_year_out_ratios,
)


def test_annual_separates_linear_trend_from_seasonal_amplitude():
    t = []
    y = []
    for year in range(2005, 2011):
        for month in range(1, 13):
            decimal = year + (month - 0.5) / 12.0
            t.append(decimal)
            y.append(
                1.7 * (decimal - 2008.0)
                + 0.3 * np.sin(2 * np.pi * decimal)
                - 0.4 * np.cos(2 * np.pi * decimal)
            )

    result = fit_annual(np.array(t), np.array(y))

    assert np.isclose(result["trend"], 1.7)
    assert np.isclose(result["amplitude"], 0.5)


def test_annual_plus_semiannual_preserves_pure_annual_amplitude():
    t = []
    y = []
    for year in [2005, 2006, 2007]:
        for month in range(1, 13):
            decimal = year + (month - 0.5) / 12.0
            t.append(decimal)
            y.append(
                1.7 * (decimal - 2006.0)
                + 0.3 * np.sin(2 * np.pi * decimal)
                - 0.4 * np.cos(2 * np.pi * decimal)
            )

    result = fit_annual_plus_semiannual(np.array(t), np.array(y))

    assert np.isclose(result["trend"], 1.7)
    assert np.isclose(result["annual_amplitude"], 0.5)
    assert result["semiannual_amplitude"] < 1e-10


def test_monthly_fixed_effect_returns_half_peak_to_trough_amplitude():
    months = np.tile(np.arange(1, 13), 3)
    decimal_years = np.concatenate(
        [year + (np.arange(1, 13) - 0.5) / 12.0 for year in [2005, 2006, 2007]]
    )
    seasonal = np.tile(
        np.array([1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2], dtype=float),
        3,
    )
    values = seasonal + 4.0 * (decimal_years - decimal_years.mean())

    result = fit_monthly_fixed_effect(decimal_years, months, values)

    assert np.isclose(result["amplitude"], 3.0)


def test_monthly_fixed_effect_is_not_estimable_when_a_calendar_month_is_missing():
    months = np.tile(np.arange(1, 12), 3)
    decimal_years = np.concatenate(
        [year + (np.arange(1, 12) - 0.5) / 12.0 for year in [2005, 2006, 2007]]
    )
    values = np.sin(2 * np.pi * decimal_years)

    result = fit_monthly_fixed_effect(decimal_years, months, values)

    assert result["estimable"] is False
    assert np.isnan(result["amplitude"])


def test_leave_one_year_out_ratios_returns_one_entry_per_year():
    records = []
    for year in [2005, 2006, 2007]:
        for month in range(1, 13):
            phase = 2 * np.pi * (month - 0.5) / 12.0
            decimal = year + (month - 0.5) / 12.0
            records.append(
                (
                    year,
                    decimal,
                    0.1 * np.sin(phase) + 0.4 * (decimal - 2006.0),
                    2.0 * np.sin(phase) - 0.3 * (decimal - 2006.0),
                )
            )

    result = leave_one_year_out_ratios(records)

    assert sorted(result) == [2005, 2006, 2007]
    assert all(np.isclose(v, 0.05) for v in result.values())


def test_json_safe_converts_nan_to_none():
    result = json_safe({"outer": {"missing": np.nan, "valid": 1.2}})

    assert result == {"outer": {"missing": None, "valid": 1.2}}
