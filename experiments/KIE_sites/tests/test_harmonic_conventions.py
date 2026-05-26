import importlib.util
from pathlib import Path
import unittest

import numpy as np


EXPT_DIR = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = EXPT_DIR / "analysis"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


phase2 = load_module("phase2_harmonics", ANALYSIS_DIR / "phase2_harmonics.py")
phase6 = load_module("phase6_phasor", ANALYSIS_DIR / "phase6_phasor.py")


class HarmonicConventionTests(unittest.TestCase):
    def test_wetland_month_index_coefficients_align_with_phase2_month_midpoints(self):
        """A July peak in wetland month-index space must align with Phase 2 B,C."""
        months = np.arange(12)
        years = np.arange(2005, 2011)

        # Wetland convention: m=0..11, Jan=0. B=0, C=-1 peaks at July (m=6).
        b_wetland = 0.0
        c_wetland = -1.0
        monthly_cycle = (
            b_wetland * np.sin(2 * np.pi * months / 12.0)
            + c_wetland * np.cos(2 * np.pi * months / 12.0)
        )

        t = []
        y = []
        for year in years:
            for month_index, value in zip(months, monthly_cycle):
                # Phase 1 stores monthly means at month midpoints.
                t.append(year + (month_index + 0.5) / 12.0)
                y.append(value)

        fit = phase2.fit_harmonic(np.array(t), np.array(y))
        b_expected, c_expected = phase6.convert_wetland_to_phase2_phasor(
            b_wetland, c_wetland
        )

        self.assertAlmostEqual(fit["B"], b_expected, places=12)
        self.assertAlmostEqual(fit["C"], c_expected, places=12)

    def test_mc_phasor_draws_vectorized_random_samples(self):
        """MC uncertainty should draw arrays, not loop one scalar draw per sample."""

        class RecordingRng:
            def __init__(self):
                self.scalar_calls = 0
                self.vector_calls = 0
                self._rng = np.random.default_rng(123)

            def normal(self, loc=0.0, scale=1.0, size=None):
                if size is None:
                    self.scalar_calls += 1
                else:
                    self.vector_calls += 1
                return self._rng.normal(loc, scale, size)

        rng = RecordingRng()
        n_mc = 37
        R, alpha = phase6.mc_phasor(
            B_obs_13c=0.2,
            C_obs_13c=-0.1,
            amp_ci_13c=[0.15, 0.25],
            peak_ci_13c=[4.5, 5.5],
            B_obs_dD=1.0,
            C_obs_dD=-2.0,
            amp_ci_dD=[1.5, 3.0],
            peak_ci_dD=[5.5, 6.5],
            B_Q=-0.2,
            C_Q=-1.2,
            Q_mean=0.9,
            dD_wetland=-374.0,
            dD_sigma=10.0,
            n_mc=n_mc,
            rng=rng,
        )

        self.assertEqual(R.shape, (n_mc,))
        self.assertEqual(alpha.shape, (n_mc,))
        self.assertGreaterEqual(rng.vector_calls, 12)
        self.assertEqual(rng.scalar_calls, 0)

    def test_alpha_filter_summary_reports_narrow_and_wide_windows(self):
        alpha = np.array([0.979, 0.985, 0.995, 1.000, 1.010, 1.019, 1.026, 1.049, 1.051])

        summary = phase6.summarize_alpha_filters(alpha)

        self.assertEqual(summary["narrow"]["filter"], [0.99, 1.02])
        self.assertEqual(summary["wide"]["filter"], [0.98, 1.05])
        self.assertEqual(summary["narrow"]["n_samples"], 4)
        self.assertEqual(summary["wide"]["n_samples"], 7)
        self.assertGreater(summary["impact_wide_minus_narrow"]["ci95"][1], 0.0)


if __name__ == "__main__":
    unittest.main()
