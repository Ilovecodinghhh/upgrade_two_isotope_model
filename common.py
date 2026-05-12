#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
common.py — Shared utilities for the methane isotope box model suite
=====================================================================

Provides:
  - Isotope conversion functions (delta ↔ fraction ↔ ratio)
  - KIE configuration and sampling (fixed or literature distributions)
  - Time-varying methane lifetime
  - Data loading (CH₄, δ¹³C, δD, source signatures, CarbonTracker BB)
  - 5-year smoothing
  - Solution quality monitoring
  - Configuration dataclass

Used by: 2x2_one.py, 2x2_two.py, 3x3_one.py, 3x3_two.py
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# ============================================================================
# CONFIGURATION
# ============================================================================

@dataclass
class ModelConfig:
    """Configuration for all model variants.

    Attributes
    ----------
    n_iterations : int
        Number of Monte Carlo iterations.
    kie_mode : str
        'fixed' — use literature central values (no sampling).
        'sampled' — draw KIE from literature distributions each iteration.
    lifetime_mode : str
        'fixed' — constant τ for all years (uses ``tau_fixed``).
        'varying' — τ(t) = 9.0 − 0.017·(t − 2010), following He et al. (2026).
    tau_fixed : float
        Fixed lifetime in years (only used when ``lifetime_mode='fixed'``).
    seed : int
        RNG seed for reproducibility.
    """
    n_iterations: int = 1000
    kie_mode: str = "sampled"       # 'fixed' | 'sampled'
    lifetime_mode: str = "varying"  # 'fixed' | 'varying'
    tau_fixed: float = 9.0
    seed: int = 42


# ============================================================================
# CONSTANTS
# ============================================================================

C13_STD = 0.011113       # ¹³C/¹²C VPDB standard
D_STD = 0.00015576       # D/H VSMOW standard
PT = 2.815               # ppb → Tg conversion factor
PT_HEMI = PT / 2.0       # per hemisphere


# ============================================================================
# ISOTOPE MATH
# ============================================================================

def delta_to_R(delta_permil: np.ndarray, std: float) -> np.ndarray:
    """δ (‰) → isotope ratio R."""
    return (delta_permil / 1000.0 + 1.0) * std

def R_to_fraction(R: np.ndarray) -> np.ndarray:
    """R → heavy-isotope mole fraction f = R / (1 + R)."""
    return R / (1.0 + R)

def fraction_to_R(f: np.ndarray) -> np.ndarray:
    """f → R = f / (1 − f)."""
    return f / (1.0 - f)

def delta_to_fraction_d13C(delta: np.ndarray) -> np.ndarray:
    return R_to_fraction(delta_to_R(delta, C13_STD))

def delta_to_fraction_dD(delta: np.ndarray) -> np.ndarray:
    return R_to_fraction(delta_to_R(delta, D_STD))

def fraction_to_delta_d13C(f: np.ndarray) -> np.ndarray:
    R = fraction_to_R(f)
    return ((R - C13_STD) / C13_STD) * 1000.0

def fraction_to_delta_dD(f: np.ndarray) -> np.ndarray:
    R = fraction_to_R(f)
    return ((R - D_STD) / D_STD) * 1000.0


# ============================================================================
# KIE CONFIGURATION & SAMPLING
# ============================================================================

# Literature distributions for Kinetic Isotope Effects
KIE_DISTRIBUTIONS = {
    'OH_13C':    {'dist': 'uniform', 'low': 1.0039, 'high': 1.0054},
    'OH_D':      {'dist': 'uniform', 'low': 1.294,  'high': 1.327},
    'Cl_13C':    {'dist': 'normal',  'mean': 1.066,  'std': 0.002},
    'Cl_D':      {'dist': 'normal',  'mean': 1.52,   'std': 0.02},
    'Strat_13C': {'dist': 'normal',  'mean': 1.003,  'std': 0.001},
    'Strat_D':   {'dist': 'normal',  'mean': 1.179,  'std': 0.01},
    'Soil_13C':  {'dist': 'normal',  'mean': 1.0201, 'std': 0.003},
    'Soil_D':    {'dist': 'normal',  'mean': 1.083,  'std': 0.01},
}

# Central values (used when kie_mode='fixed')
KIE_FIXED = {
    'OH_13C':    (1.0039 + 1.0054) / 2,  # 1.00465
    'OH_D':      (1.294 + 1.327) / 2,     # 1.3105
    'Cl_13C':    1.066,
    'Cl_D':      1.52,
    'Strat_13C': 1.003,
    'Strat_D':   1.179,
    'Soil_13C':  1.0201,
    'Soil_D':    1.083,
}


def sample_KIE(rng: np.random.Generator, mode: str = "sampled") -> dict:
    """Draw KIE values — either fixed central values or from literature distributions.

    Parameters
    ----------
    rng : numpy Generator
    mode : 'fixed' or 'sampled'

    Returns
    -------
    dict with keys OH_13C, OH_D, Cl_13C, Cl_D, Strat_13C, Strat_D, Soil_13C, Soil_D
    """
    if mode == "fixed":
        return dict(KIE_FIXED)

    kies = {}
    for key, cfg in KIE_DISTRIBUTIONS.items():
        if cfg['dist'] == 'uniform':
            kies[key] = rng.uniform(cfg['low'], cfg['high'])
        elif cfg['dist'] == 'normal':
            kies[key] = rng.normal(cfg['mean'], cfg['std'])
    return kies


# Sink fractions by hemisphere (also used for global one-box)
SINK_FRACTIONS_GLOBAL = {'OH': 0.835, 'Cl': 0.035, 'Strat': 0.070, 'Soil': 0.060}
SINK_FRACTIONS_NH = {'OH': 0.825, 'Cl': 0.040, 'Strat': 0.070, 'Soil': 0.065}
SINK_FRACTIONS_SH = {'OH': 0.850, 'Cl': 0.028, 'Strat': 0.070, 'Soil': 0.052}


def compute_bulk_KIE(kies: dict, sink_fracs: dict):
    """Compute weighted-mean KIE for ¹³C and D from individual sink KIEs.

    Returns (KIE_13C, KIE_D).
    """
    kie_13C = sum(kies[f'{s}_13C'] * sink_fracs[s] for s in ('OH', 'Cl', 'Strat', 'Soil'))
    kie_D   = sum(kies[f'{s}_D']   * sink_fracs[s] for s in ('OH', 'Cl', 'Strat', 'Soil'))
    return kie_13C, kie_D


# ============================================================================
# TIME-VARYING LIFETIME
# ============================================================================

def compute_lifetime(years: np.ndarray, mode: str = "varying",
                     tau_fixed: float = 9.0) -> np.ndarray:
    """Compute CH₄ lifetime array.

    Parameters
    ----------
    years : array of calendar years
    mode : 'fixed' or 'varying'
    tau_fixed : constant lifetime (used when mode='fixed')

    Returns
    -------
    τ : array of lifetimes in years
    """
    years = np.asarray(years, dtype=float)
    if mode == "fixed":
        return np.full_like(years, tau_fixed)
    # He et al. (2026): τ(t) = 9.0 − 0.017·(t − 2010)
    return 9.0 - 0.017 * (years - 2010)


LIFETIME_RATIO_NH = 0.95
LIFETIME_RATIO_SH = 1.05


# ============================================================================
# TWO-HEMISPHERE PARAMETERS
# ============================================================================

# NH/SH BB split from GFED4
BB_NH_FRACTION = 0.55
BB_SH_FRACTION = 0.45

# Interhemispheric exchange
TAU_EX_MEAN = 1.0   # years
TAU_EX_STD = 0.1

# δD hemispheric offset (Riddell-Young 2025: NH ~12‰ lower than SH)
DD_IH_OFFSET = 6.0  # ‰ — LEGACY: NH = global − 6, SH = global + 6 (unused when real hemi data loaded)

# IH gradient in CH₄ concentration
def compute_IH_gradient(n_points: int) -> np.ndarray:
    return np.linspace(80, 100, n_points)


# ============================================================================
# 5-YEAR SMOOTHING
# ============================================================================

def smooth_5yr(arr_2d: np.ndarray) -> np.ndarray:
    """Apply 5-year moving average along axis 0. Shape: [years, iterations]."""
    n_years = arr_2d.shape[0]
    if n_years < 5:
        return arr_2d.copy()
    result = np.zeros_like(arr_2d)
    result[0] = np.mean(arr_2d[0:3], axis=0)
    result[1] = np.mean(arr_2d[0:4], axis=0)
    for i in range(2, n_years - 2):
        result[i] = np.mean(arr_2d[i-2:i+3], axis=0)
    result[-2] = np.mean(arr_2d[-4:], axis=0)
    result[-1] = np.mean(arr_2d[-3:], axis=0)
    return result


# ============================================================================
# SOLUTION QUALITY MONITOR
# ============================================================================

class QualityMonitor:
    """Track condition numbers and non-physical solutions."""

    def __init__(self, n_years: int, n_iter: int, label: str = ""):
        self.label = label
        self.n_years = n_years
        self.n_iter = n_iter
        self.condition_numbers = np.zeros((n_years, n_iter))
        self.is_nonphysical = np.zeros((n_years, n_iter), dtype=bool)
        self.is_nan = np.zeros((n_years, n_iter), dtype=bool)

    def record(self, yr: int, it: int, A: np.ndarray, x: np.ndarray):
        cond = np.linalg.cond(A)
        self.condition_numbers[yr, it] = cond
        if np.any(~np.isfinite(x)):
            self.is_nan[yr, it] = True
            self.is_nonphysical[yr, it] = True
        elif np.any(x < 0):
            self.is_nonphysical[yr, it] = True

    def summary(self) -> dict:
        total = self.n_years * self.n_iter
        n_nonphys = int(np.sum(self.is_nonphysical))
        pct = 100.0 * n_nonphys / total
        mean_c = float(np.mean(self.condition_numbers))
        max_c = float(np.max(self.condition_numbers))
        report = {
            'label': self.label,
            'total_solves': total,
            'nonphysical_count': n_nonphys,
            'nonphysical_pct': round(pct, 2),
            'nan_count': int(np.sum(self.is_nan)),
            'mean_condition': round(mean_c, 1),
            'max_condition': round(max_c, 1),
        }
        print(f"  [{self.label}] {pct:.1f}% non-physical, "
              f"mean cond = {mean_c:.1f} (max = {max_c:.1f})")
        return report


# ============================================================================
# DATA LOADING
# ============================================================================

def _resolve_data_dirs(base_dir: Path):
    """Resolve data and source-signature directories."""
    # Data lives in sibling repo ../TwoIsotopeBoxModel/rel/
    rel_dir = base_dir.parent / "TwoIsotopeBoxModel" / "rel"
    # Fallback: local rel/ copy
    if not rel_dir.exists():
        rel_dir = base_dir / "rel"
    data_dir = rel_dir / "data"
    src_dir = rel_dir / "output"
    return data_dir, src_dir


def pad_to_length(arr: np.ndarray, length: int) -> np.ndarray:
    """Pad or truncate a 1-D array to exact ``length``."""
    arr = np.asarray(arr).flatten()
    if len(arr) >= length:
        return arr[:length]
    return np.concatenate([np.full(length - len(arr), arr[0]), arr])[:length]


def _annual_average(dates, values):
    """Group sub-annual data into annual means (≥6 points required)."""
    years_floor = np.floor(dates).astype(int)
    unique_years = np.unique(years_floor)
    ann_years, ann_means = [], []
    for yr in unique_years:
        mask = years_floor == yr
        if np.sum(mask) >= 6:
            ann_years.append(yr)
            ann_means.append(np.nanmean(values[mask]))
    return np.array(ann_years), np.array(ann_means)


@dataclass
class LoadedData:
    """Container for all loaded observational / source-signature data."""

    # CH₄ concentrations (1999–2022, 24 values)
    CH4_global: np.ndarray = field(default_factory=lambda: np.array([]))
    CH4_years: np.ndarray = field(default_factory=lambda: np.array([]))
    # For two-box models
    CH4_NH: Optional[np.ndarray] = None
    CH4_SH: Optional[np.ndarray] = None

    # δ¹³C annual means (1999–2022)
    c13_global: np.ndarray = field(default_factory=lambda: np.array([]))
    c13_NH: Optional[np.ndarray] = None
    c13_SH: Optional[np.ndarray] = None

    # δD annual global mean
    dD_global: np.ndarray = field(default_factory=lambda: np.array([]))
    # δD hemispheric (NH/SH) — actual observations when available
    dD_NH: Optional[np.ndarray] = None
    dD_SH: Optional[np.ndarray] = None

    # MC iteration matrices (rows=years, cols=iterations)
    d13C_MC: np.ndarray = field(default_factory=lambda: np.array([]))
    dD_MC: np.ndarray = field(default_factory=lambda: np.array([]))
    dD_NH_MC: Optional[np.ndarray] = None
    dD_SH_MC: Optional[np.ndarray] = None
    # δD data year range (may differ from CH₄/δ¹³C range of 1999–2022)
    dD_start_year: int = 2005

    # Source signatures — central values + uncertainties
    ff_d13C: np.ndarray = field(default_factory=lambda: np.array([]))
    ff_d13C_U: np.ndarray = field(default_factory=lambda: np.array([]))
    bb_d13C: np.ndarray = field(default_factory=lambda: np.array([]))
    bb_d13C_U: np.ndarray = field(default_factory=lambda: np.array([]))
    mic_d13C_mean: float = 0.0

    ff_dD: np.ndarray = field(default_factory=lambda: np.array([]))
    ff_dD_U: np.ndarray = field(default_factory=lambda: np.array([]))
    bb_dD: np.ndarray = field(default_factory=lambda: np.array([]))
    bb_dD_U: np.ndarray = field(default_factory=lambda: np.array([]))
    mic_dD_mean: float = 0.0

    # MC-trend matrices for source signatures
    FF_d13C_MC_EDGAR: Optional[pd.DataFrame] = None
    FF_dD_MC_EDGAR: Optional[pd.DataFrame] = None
    Mic_d13C_MC: Optional[pd.DataFrame] = None
    Mic_dD_MC: Optional[pd.DataFrame] = None

    # Hemispheric δD source-signature MC matrices (rows=years, cols=1000 MC)
    FF_dD_NH_MC: Optional[np.ndarray] = None
    FF_dD_SH_MC: Optional[np.ndarray] = None
    Mic_dD_NH_MC: Optional[np.ndarray] = None
    Mic_dD_SH_MC: Optional[np.ndarray] = None
    BB_dD_NH_MC: Optional[np.ndarray] = None
    BB_dD_SH_MC: Optional[np.ndarray] = None

    # CarbonTracker BB (for 2×2 models)
    BB_annual: np.ndarray = field(default_factory=lambda: np.array([]))
    BB_global_mean: float = 0.0

    # Model dimensions
    n_years: int = 0
    model_years: np.ndarray = field(default_factory=lambda: np.array([]))


def load_data(base_dir: Path, two_box: bool = False) -> LoadedData:
    """Load all observational and source-signature data.

    Parameters
    ----------
    base_dir : Path to the repository root
    two_box : if True, also compute NH/SH CH₄ and δ¹³C splits

    Returns
    -------
    LoadedData instance
    """
    data_dir, src_dir = _resolve_data_dirs(base_dir)
    d = LoadedData()

    # === CH₄ concentrations ===
    CH4raw = pd.read_excel(data_dir / "GML_CH4_AnnualMean.xlsx").to_numpy()
    d.CH4_global = CH4raw[15:39, 1].astype(float)  # 1999–2022 (24 values)
    d.CH4_years = CH4raw[15:39, 0].astype(float)

    if two_box:
        IH_grad = compute_IH_gradient(len(d.CH4_global))
        d.CH4_NH = d.CH4_global + IH_grad / 2.0
        d.CH4_SH = d.CH4_global - IH_grad / 2.0

    # === δ¹³C ===
    C13raw = pd.read_excel(data_dir / "ch4c13_nh_sh_mean.xlsx", header=None).to_numpy()
    c13_dates = C13raw[:, 0]
    c13_ann_years, c13_ann_global = _annual_average(c13_dates, C13raw[:, 1])

    si = np.where(c13_ann_years == 1999)[0][0]
    ei = np.where(c13_ann_years == 2022)[0][0] + 1
    d.c13_global = c13_ann_global[si:ei]

    if two_box:
        _, c13_ann_NH = _annual_average(c13_dates, C13raw[:, 2])
        _, c13_ann_SH = _annual_average(c13_dates, C13raw[:, 3])
        d.c13_NH = c13_ann_NH[si:ei]
        d.c13_SH = c13_ann_SH[si:ei]

    # === δ¹³C MC iterations ===
    d13C_raw = np.loadtxt(str(data_dir / "d13C_dei_compiled.txt"))
    d.d13C_MC = d13C_raw[1:, 1:]  # rows=years, cols=iterations

    # === δD (Dasgupta calibration) ===
    dD_df = pd.read_excel(data_dir / "GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx")
    # Drop any non-numeric columns (e.g. 'Unnamed: ...')
    dD_df = dD_df.loc[:, [isinstance(c, (int, float)) for c in dD_df.columns]]
    dD_num = dD_df.apply(pd.to_numeric, errors="coerce")
    # Reconstruct full annual MC: header row (year0 = 2005) + data rows (2006..2023)
    d.dD_start_year = int(float(dD_df.columns[0]))  # 2005
    _dD_first = dD_num.columns[1:].to_numpy(dtype=np.float64).reshape(1, -1)
    _dD_rest = dD_num.iloc[:, 1:].to_numpy(dtype=np.float64)
    d.dD_MC = np.vstack([_dD_first, _dD_rest])  # (19, ~998) — years × iterations
    d.dD_global = np.nanmean(d.dD_MC, axis=1)    # annual global mean from MC

    # === δD hemispheric MC iterations (Dasgupta calibration) ===
    if two_box:
        dD_NH_file = data_dir / "NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx"
        dD_SH_file = data_dir / "SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx"
        if dD_NH_file.exists() and dD_SH_file.exists():
            dD_NH_df = pd.read_excel(dD_NH_file)
            dD_SH_df = pd.read_excel(dD_SH_file)
            dD_NH_df = dD_NH_df.loc[:, [isinstance(c, (int, float)) for c in dD_NH_df.columns]]
            dD_SH_df = dD_SH_df.loc[:, [isinstance(c, (int, float)) for c in dD_SH_df.columns]]
            dD_NH_num = dD_NH_df.apply(pd.to_numeric, errors="coerce")
            dD_SH_num = dD_SH_df.apply(pd.to_numeric, errors="coerce")
            # Reconstruct: header year + data rows → full annual series
            nh_first = dD_NH_num.columns[1:].to_numpy(dtype=np.float64).reshape(1, -1)
            nh_rest = dD_NH_num.iloc[:, 1:].to_numpy(dtype=np.float64)
            d.dD_NH_MC = np.vstack([nh_first, nh_rest])  # (19, 1000)
            sh_first = dD_SH_num.columns[1:].to_numpy(dtype=np.float64).reshape(1, -1)
            sh_rest = dD_SH_num.iloc[:, 1:].to_numpy(dtype=np.float64)
            d.dD_SH_MC = np.vstack([sh_first, sh_rest])  # (19, 1000)
            # Annual means for reference
            d.dD_NH = np.nanmean(d.dD_NH_MC, axis=1)
            d.dD_SH = np.nanmean(d.dD_SH_MC, axis=1)

            # Fill NaN rows (e.g. 2020-2023 station gaps) using global MC +
            # mean hemispheric offset from years with good data
            _nan_rows_nh = np.where(np.all(np.isnan(d.dD_NH_MC), axis=1))[0]
            if len(_nan_rows_nh) > 0:
                # Compute mean offset from valid rows that exist in both global and hemi
                _n_common = min(d.dD_NH_MC.shape[1], d.dD_MC.shape[1])
                _valid = np.where(~np.all(np.isnan(d.dD_NH_MC[:, :_n_common]), axis=1))[0]
                _nh_offset = np.nanmean(
                    d.dD_NH_MC[np.ix_(_valid, range(_n_common))] -
                    d.dD_MC[np.ix_(_valid, range(_n_common))]
                )
                _sh_offset = np.nanmean(
                    d.dD_SH_MC[np.ix_(_valid, range(_n_common))] -
                    d.dD_MC[np.ix_(_valid, range(_n_common))]
                )
                for _r in _nan_rows_nh:
                    if _r < d.dD_MC.shape[0]:
                        # Use global MC + offset; broadcast to match column count
                        _g = d.dD_MC[_r, :_n_common]
                        d.dD_NH_MC[_r, :_n_common] = _g + _nh_offset
                        d.dD_SH_MC[_r, :_n_common] = _g + _sh_offset
                        # Fill any remaining columns with repeat
                        if _n_common < d.dD_NH_MC.shape[1]:
                            d.dD_NH_MC[_r, _n_common:] = d.dD_NH_MC[_r, _n_common - 1]
                            d.dD_SH_MC[_r, _n_common:] = d.dD_SH_MC[_r, _n_common - 1]
                # Recompute means
                d.dD_NH = np.nanmean(d.dD_NH_MC, axis=1)
                d.dD_SH = np.nanmean(d.dD_SH_MC, axis=1)

    # === Hemispheric δD MC iterations (from Riddell-Young station-level pipeline) ===
    if two_box:
        dD_NH_file = data_dir / "NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx"
        dD_SH_file = data_dir / "SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx"
        if dD_NH_file.exists() and dD_SH_file.exists():
            dD_NH_df = pd.read_excel(dD_NH_file, header=None)
            dD_SH_df = pd.read_excel(dD_SH_file, header=None)
            # Col 0 = year (2005–2023), cols 1–1000 = MC iterations
            dD_NH_years = dD_NH_df.iloc[:, 0].to_numpy(dtype=int)
            dD_NH_mc = dD_NH_df.iloc[:, 1:].to_numpy(dtype=np.float64)
            dD_SH_mc = dD_SH_df.iloc[:, 1:].to_numpy(dtype=np.float64)
            # Forward-fill NaN rows (2020+ may be all-NaN)
            for arr in [dD_NH_mc, dD_SH_mc]:
                for i in range(1, arr.shape[0]):
                    mask = np.isnan(arr[i])
                    if mask.all():
                        arr[i] = arr[i - 1]
            d.dD_NH_MC = dD_NH_mc
            d.dD_SH_MC = dD_SH_mc
            d.dD_NH = np.nanmean(d.dD_NH_MC, axis=1)
            d.dD_SH = np.nanmean(d.dD_SH_MC, axis=1)
            d._dD_hemi_years = dD_NH_years

    # === Hemispheric δD source signatures MC ===
    if two_box:
        for attr, fname in [
            ('FF_dD_NH_MC', 'FF_dD_NH_MC.csv'),
            ('FF_dD_SH_MC', 'FF_dD_SH_MC.csv'),
            ('Mic_dD_NH_MC', 'Mic_dD_NH_MC.csv'),
            ('Mic_dD_SH_MC', 'Mic_dD_SH_MC.csv'),
            ('BB_dD_NH_MC', 'BB_dD_NH_MC.csv'),
            ('BB_dD_SH_MC', 'BB_dD_SH_MC.csv'),
        ]:
            fpath = data_dir / fname
            if fpath.exists():
                mc_df = pd.read_csv(fpath)
                # Col 0 = year (1998–2021), cols 1–1000 = MC iterations
                setattr(d, attr, mc_df.iloc[:, 1:].to_numpy(dtype=np.float64))

    # === Source signatures — d13C ===
    BB_d13C_data = pd.read_csv(src_dir / "BB_d13C_annual.csv", header=None)
    Mic_d13C_data = pd.read_csv(src_dir / "Mic_d13C_annual.csv", header=None)
    FF_d13C_data = pd.read_csv(src_dir / "FF_d13C_GlobUnc.csv")
    FF_d13C_MC_EDGAR = pd.read_csv(src_dir / "FF_d13C_GlobMC_EDGAR.csv")
    Mic_d13C_MC = pd.read_csv(src_dir / "Mic_d13C_MC.csv", header=None)

    d.mic_d13C_mean = Mic_d13C_data.iloc[:, 1].mean()
    d.ff_d13C = np.array(FF_d13C_data.iloc[28:, 1]).flatten()
    d.ff_d13C_U = np.array(FF_d13C_data.iloc[28:, 2]).flatten()

    bb_c = np.array(BB_d13C_data.iloc[1:, 1]).flatten()
    bb_cU = np.array(BB_d13C_data.iloc[1:, 2]).flatten()
    d.bb_d13C = np.concatenate([bb_c, [bb_c[-1]]])
    d.bb_d13C_U = np.concatenate([bb_cU, [bb_cU[-1]]])

    d.FF_d13C_MC_EDGAR = FF_d13C_MC_EDGAR.iloc[28:, 1:]
    d.Mic_d13C_MC = Mic_d13C_MC.iloc[:, 1:]

    # === Source signatures — dD ===
    BB_dD_data = pd.read_csv(src_dir / "BB_dD_annual.csv", header=None)
    Mic_dD_data = pd.read_csv(src_dir / "Mic_dD_AnnGlob.csv", header=None)
    FF_dD_data = pd.read_csv(src_dir / "FF_dD_GlobUnc.csv")
    FF_dD_MC_EDGAR = pd.read_csv(src_dir / "FF_dD_GlobMC_EDGAR.csv")
    Mic_dD_MC_trends = pd.read_csv(src_dir / "Mic_dD_MC.csv", header=None)

    d.mic_dD_mean = Mic_dD_data.iloc[:, 1].mean()
    d.ff_dD = np.array(FF_dD_data.iloc[34:, 1]).flatten()
    d.ff_dD_U = np.array(FF_dD_data.iloc[34:, 2]).flatten()

    bb_d = np.array(BB_dD_data.iloc[:, 1]).flatten()
    bb_dU = np.array(BB_dD_data.iloc[:, 2]).flatten()
    d.bb_dD = np.concatenate([np.full(3, bb_d[-1]), bb_d, [bb_d[-1]]])
    d.bb_dD_U = np.concatenate([np.full(3, bb_dU[-1]), bb_dU, [bb_dU[-1]]])

    d.FF_dD_MC_EDGAR = FF_dD_MC_EDGAR.iloc[34:, 1:]
    # Pad FF_dD_MC_EDGAR to ≥24 rows
    if d.FF_dD_MC_EDGAR.shape[0] < 24:
        pad_n = 24 - d.FF_dD_MC_EDGAR.shape[0]
        pad = pd.concat([d.FF_dD_MC_EDGAR.iloc[0:1]] * pad_n, ignore_index=True)
        d.FF_dD_MC_EDGAR = pd.concat([pad, d.FF_dD_MC_EDGAR], ignore_index=True)

    d.Mic_dD_MC = Mic_dD_MC_trends.iloc[6:, 1:]

    # === Hemispheric δD source signatures (MC iterations) ===
    if two_box:
        _hemi_src_files = {
            'FF_dD_NH_MC': 'FF_dD_NH_MC.csv',
            'FF_dD_SH_MC': 'FF_dD_SH_MC.csv',
            'Mic_dD_NH_MC': 'Mic_dD_NH_MC.csv',
            'Mic_dD_SH_MC': 'Mic_dD_SH_MC.csv',
            'BB_dD_NH_MC': 'BB_dD_NH_MC.csv',
            'BB_dD_SH_MC': 'BB_dD_SH_MC.csv',
        }
        for attr, fname in _hemi_src_files.items():
            fpath = data_dir / fname
            if fpath.exists():
                _df = pd.read_csv(fpath)
                # col 0 = year, cols 1..1000 = MC iterations
                _mat = _df.iloc[:, 1:].to_numpy(dtype=np.float64)
                setattr(d, attr, _mat)

    # === CarbonTracker BB ===
    data_CT = pd.read_excel(data_dir / "CarbonTracker_CH4.xlsx")
    bbCT = data_CT.iloc[:, 9].values
    d.BB_global_mean = float(np.mean(bbCT))
    d.BB_annual = bbCT

    # === Model dimensions ===
    d.n_years = len(d.CH4_global) - 1  # 23
    d.model_years = np.arange(1999, 1999 + d.n_years)

    # Pad source-signature arrays to n_years
    tl = d.n_years
    d.ff_d13C = pad_to_length(d.ff_d13C, tl)
    d.ff_d13C_U = pad_to_length(d.ff_d13C_U, tl)
    d.bb_d13C = pad_to_length(d.bb_d13C, tl)
    d.bb_d13C_U = pad_to_length(d.bb_d13C_U, tl)
    d.ff_dD = pad_to_length(d.ff_dD, tl)
    d.ff_dD_U = pad_to_length(d.ff_dD_U, tl)
    d.bb_dD = pad_to_length(d.bb_dD, tl)
    d.bb_dD_U = pad_to_length(d.bb_dD_U, tl)

    # Pad MC matrices
    # δD data starts at dD_start_year (2005), model starts at 1999 → front-pad
    pad_dD = max(0, tl + 1 - d.dD_MC.shape[0])
    if pad_dD > 0:
        d.dD_MC = np.vstack([np.repeat(d.dD_MC[0:1], pad_dD, axis=0), d.dD_MC])
    # Also pad hemispheric δD MC matrices if available
    if d.dD_NH_MC is not None:
        pad_hemi = max(0, tl + 1 - d.dD_NH_MC.shape[0])
        if pad_hemi > 0:
            d.dD_NH_MC = np.vstack([np.repeat(d.dD_NH_MC[0:1], pad_hemi, axis=0), d.dD_NH_MC])
            d.dD_SH_MC = np.vstack([np.repeat(d.dD_SH_MC[0:1], pad_hemi, axis=0), d.dD_SH_MC])

    if d.Mic_dD_MC.shape[0] < tl:
        pn = tl - d.Mic_dD_MC.shape[0]
        d.Mic_dD_MC = pd.concat(
            [pd.concat([d.Mic_dD_MC.iloc[0:1]] * pn, ignore_index=True), d.Mic_dD_MC],
            ignore_index=True)
    elif d.Mic_dD_MC.shape[0] > tl:
        d.Mic_dD_MC = d.Mic_dD_MC.iloc[:tl]

    # BB annual → align to model years
    if len(d.BB_annual) >= d.n_years:
        d.BB_annual = d.BB_annual[:d.n_years].astype(float)
    else:
        d.BB_annual = np.full(d.n_years, d.BB_global_mean)

    return d


# ============================================================================
# MC SAMPLING HELPERS
# ============================================================================

def sample_source_signatures(rng, data: LoadedData, k: int, target_length: int):
    """Sample source end-member signatures for MC iteration k.

    Returns dict with keys:
        ff_d13C, bb_d13C, mic_d13C, ff_dD, bb_dD, mic_dD
    (each is a 1-D array of length target_length)
    """
    tl = target_length
    g_ff_c = rng.normal()
    g_bb_c = rng.normal()
    g_ff_d = rng.normal()
    g_bb_d = rng.normal()

    # FF δ¹³C
    if k < data.FF_d13C_MC_EDGAR.shape[1]:
        ff_c = np.array(data.FF_d13C_MC_EDGAR.iloc[:, k]).flatten()
    else:
        ff_c = data.ff_d13C + g_ff_c * data.ff_d13C_U
    # FF δD
    if k < data.FF_dD_MC_EDGAR.shape[1]:
        ff_d = np.array(data.FF_dD_MC_EDGAR.iloc[:, k]).flatten()
    else:
        ff_d = data.ff_dD + g_ff_d * data.ff_dD_U
    # BB
    bb_c = data.bb_d13C + g_bb_c * data.bb_d13C_U
    bb_d = data.bb_dD + g_bb_d * data.bb_dD_U
    # Mic
    if k < data.Mic_d13C_MC.shape[1]:
        mic_c = np.array(data.Mic_d13C_MC.iloc[:tl, k]).flatten()
    else:
        mic_c = np.full(tl, data.mic_d13C_mean)
    if k < data.Mic_dD_MC.shape[1]:
        mic_d = np.array(data.Mic_dD_MC.iloc[:tl, k]).flatten()
    else:
        mic_d = np.full(tl, data.mic_dD_mean)

    return {
        'ff_d13C': pad_to_length(ff_c, tl),
        'bb_d13C': pad_to_length(bb_c, tl),
        'mic_d13C': pad_to_length(mic_c, tl),
        'ff_dD': pad_to_length(ff_d, tl),
        'bb_dD': pad_to_length(bb_d, tl),
        'mic_dD': pad_to_length(mic_d, tl),
    }


def sample_source_signatures_hemi(rng, data: LoadedData, k: int, target_length: int):
    """Sample hemispheric δD source signatures for MC iteration k.

    Uses actual NH/SH MC iterations when available; falls back to global.

    Returns dict with keys:
        ff_dD_NH, ff_dD_SH, bb_dD_NH, bb_dD_SH, mic_dD_NH, mic_dD_SH
        ff_d13C, bb_d13C, mic_d13C  (global — δ¹³C unchanged)
    """
    tl = target_length
    # Get global signatures first (includes δ¹³C)
    global_sigs = sample_source_signatures(rng, data, k, tl)

    # δD: use hemispheric MC if available, else duplicate global
    def _pick_hemi(mc_mat, global_arr, k, tl):
        if mc_mat is not None:
            col = min(k, mc_mat.shape[1] - 1)
            arr = mc_mat[:tl, col].copy()
            return pad_to_length(arr, tl)
        return global_arr

    result = {
        'ff_d13C': global_sigs['ff_d13C'],
        'bb_d13C': global_sigs['bb_d13C'],
        'mic_d13C': global_sigs['mic_d13C'],
        'ff_dD_NH': _pick_hemi(data.FF_dD_NH_MC, global_sigs['ff_dD'], k, tl),
        'ff_dD_SH': _pick_hemi(data.FF_dD_SH_MC, global_sigs['ff_dD'], k, tl),
        'bb_dD_NH': _pick_hemi(data.BB_dD_NH_MC, global_sigs['bb_dD'], k, tl),
        'bb_dD_SH': _pick_hemi(data.BB_dD_SH_MC, global_sigs['bb_dD'], k, tl),
        'mic_dD_NH': _pick_hemi(data.Mic_dD_NH_MC, global_sigs['mic_dD'], k, tl),
        'mic_dD_SH': _pick_hemi(data.Mic_dD_SH_MC, global_sigs['mic_dD'], k, tl),
        # Also keep global for backward compat
        'ff_dD': global_sigs['ff_dD'],
        'bb_dD': global_sigs['bb_dD'],
        'mic_dD': global_sigs['mic_dD'],
    }
    return result


def sample_atm_d13C(data: LoadedData, k: int, target_length: int) -> np.ndarray:
    """Return sampled global δ¹³C time series (length target_length+1)."""
    tl1 = target_length + 1
    if k < data.d13C_MC.shape[1]:
        return data.d13C_MC[:tl1, k]
    return data.d13C_MC[:tl1, -1]


def sample_atm_dD(data: LoadedData, k: int, target_length: int) -> np.ndarray:
    """Return sampled global δD time series (length target_length+1)."""
    tl1 = target_length + 1
    if k < data.dD_MC.shape[1]:
        arr = data.dD_MC[:tl1, k]
    else:
        arr = data.dD_MC[:tl1, -1]
    if len(arr) < tl1:
        arr = np.concatenate([np.full(tl1 - len(arr), arr[0]), arr])
    return arr


def sample_atm_dD_hemi(data: LoadedData, k: int, target_length: int):
    """Return sampled NH and SH δD time series (each length target_length+1).

    Uses actual hemispheric MC iterations when available; falls back to
    global +/- DD_IH_OFFSET for backward compatibility.

    Returns (dD_NH, dD_SH) as numpy arrays.
    """
    tl1 = target_length + 1
    if data.dD_NH_MC is not None and data.dD_SH_MC is not None:
        # Use real hemispheric data (already padded during load_data)
        col = min(k, data.dD_NH_MC.shape[1] - 1)
        nh = data.dD_NH_MC[:tl1, col].copy()
        sh = data.dD_SH_MC[:tl1, col].copy()
        if len(nh) < tl1:
            nh = np.concatenate([np.full(tl1 - len(nh), nh[0]), nh])
            sh = np.concatenate([np.full(tl1 - len(sh), sh[0]), sh])
        return nh, sh
    else:
        # Legacy fallback: global ± fixed offset
        dD_glob = sample_atm_dD(data, k, target_length)
        return dD_glob - DD_IH_OFFSET, dD_glob + DD_IH_OFFSET


def sample_source_signatures_hemi(rng, data: LoadedData, k: int, target_length: int):
    """Sample hemispheric δD source signatures for MC iteration k.

    Uses actual NH/SH MC iterations when available; falls back to global.

    Returns dict with keys:
        ff_dD_NH, ff_dD_SH, bb_dD_NH, bb_dD_SH, mic_dD_NH, mic_dD_SH
        ff_d13C, bb_d13C, mic_d13C  (global — δ¹³C unchanged)
    """
    tl = target_length
    # Get global signatures first (includes δ¹³C)
    global_sigs = sample_source_signatures(rng, data, k, tl)

    # δD: use hemispheric MC if available, else duplicate global
    def _pick_hemi(mc_mat, global_arr, k, tl):
        if mc_mat is not None:
            col = min(k, mc_mat.shape[1] - 1)
            # Source sig MC covers 1998–2021 (24 rows); model starts 1999 → skip row 0
            arr = mc_mat[1:tl+1, col].copy()
            return pad_to_length(arr, tl)
        return global_arr

    result = {
        'ff_d13C': global_sigs['ff_d13C'],
        'bb_d13C': global_sigs['bb_d13C'],
        'mic_d13C': global_sigs['mic_d13C'],
        'ff_dD_NH': _pick_hemi(data.FF_dD_NH_MC, global_sigs['ff_dD'], k, tl),
        'ff_dD_SH': _pick_hemi(data.FF_dD_SH_MC, global_sigs['ff_dD'], k, tl),
        'bb_dD_NH': _pick_hemi(data.BB_dD_NH_MC, global_sigs['bb_dD'], k, tl),
        'bb_dD_SH': _pick_hemi(data.BB_dD_SH_MC, global_sigs['bb_dD'], k, tl),
        'mic_dD_NH': _pick_hemi(data.Mic_dD_NH_MC, global_sigs['mic_dD'], k, tl),
        'mic_dD_SH': _pick_hemi(data.Mic_dD_SH_MC, global_sigs['mic_dD'], k, tl),
        # Also keep global for backward compat
        'ff_dD': global_sigs['ff_dD'],
        'bb_dD': global_sigs['bb_dD'],
        'mic_dD': global_sigs['mic_dD'],
    }
    return result
# ============================================================================
# TREND ANALYSIS HELPERS
# ============================================================================

def trend_change(compiled: np.ndarray, model_years: np.ndarray,
                 base=(2005, 2007), recent_n=3):
    """Compute Δ(recent vs base) across MC iterations.

    Returns (delta_array, pct_positive).
    """
    yr0 = int(model_years[0])
    i0 = base[0] - yr0
    i1 = base[1] - yr0 + 1
    base_mean = compiled[i0:i1].mean(axis=0)
    recent_mean = compiled[-recent_n:].mean(axis=0)
    delta = recent_mean - base_mean
    pct_pos = (delta > 0).sum() / len(delta) * 100
    return delta, pct_pos
