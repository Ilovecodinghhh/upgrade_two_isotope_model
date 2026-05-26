import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


EXPT_DIR = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = EXPT_DIR / "analysis"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


phase7 = load_module("phase7_yearly_stability", ANALYSIS_DIR / "phase7_yearly_stability.py")


class YearlyStabilityTests(unittest.TestCase):
    def test_fit_site_year_requires_minimum_months(self):
        months = np.arange(1, 8)
        df = pd.DataFrame({
            "year": 2006,
            "month": months,
            "decimal_year": 2006 + (months - 0.5) / 12.0,
            "d13C_mean": np.sin(2 * np.pi * (months - 0.5) / 12.0),
            "dD_mean": 10 * np.sin(2 * np.pi * (months - 0.5) / 12.0),
        })

        self.assertIsNone(phase7.fit_site_year("TST", 2006, df, min_months=8))

    def test_fit_site_year_returns_ratio_for_synthetic_cycle(self):
        months = np.arange(1, 13)
        phase = 2 * np.pi * (months - 0.5) / 12.0
        df = pd.DataFrame({
            "year": 2006,
            "month": months,
            "decimal_year": 2006 + (months - 0.5) / 12.0,
            "d13C_mean": 0.2 * np.sin(phase),
            "dD_mean": 4.0 * np.sin(phase),
        })

        result = phase7.fit_site_year("TST", 2006, df, min_months=8)

        self.assertIsNotNone(result)
        self.assertEqual(result["n_months"], 12)
        self.assertAlmostEqual(result["ratio"], 0.05, places=10)
        self.assertAlmostEqual(result["phase_diff_months"], 0.0, places=10)


if __name__ == "__main__":
    unittest.main()
