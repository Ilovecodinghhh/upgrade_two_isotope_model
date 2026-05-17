"""
Shared infrastructure for all 4 model variants.
=================================================
Contains: data loading, KIE sampling, isotope math, lifetime calculation,
solution quality monitor, smoothing, and plotting utilities.
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import lsq_linear

from .inputs import (
    KIE_OPTIONS, SINK_FRACTION_OPTIONS, LIFETIME_OPTIONS,
    FF_SIGNATURE_OPTIONS, DEFAULT_CONFIG
)

# ═══════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════
C13Std = 0.011113        # ¹³C/¹²C reference ratio (IUPAC 2024)
DStd = 0.00015576        # D/H reference ratio
PT = 2.815               # ppb → Tg conversion factor
PT_HEMI = PT / 2.0       # per hemisphere

# ═══════════════════════════════════════════════════════════════════════════
# ISOTOPE MATH
# ═══════════════════════════════════════════════════════════════════════════
def delta_to_fraction_d13C(delta_permil):
    R = (delta_permil / 1000.0 + 1.0) * C13Std
    return R / (1.0 + R)

def delta_to_fraction_dD(delta_permil):
    R = (delta_permil / 1000.0 + 1.0) * DStd
    return R / (1.0 + R)

def fraction_to_delta_d13C(f):
    R = f / (1.0 - f)
    return ((R - C13Std) / C13Std) * 1000.0

def fraction_to_delta_dD(f):
    R = f / (1.0 - f)
    return ((R - DStd) / DStd) * 1000.0


# ═══════════════════════════════════════════════════════════════════════════
# KIE SAMPLING
# ═══════════════════════════════════════════════════════════════════════════
def build_kie_sampler(config):
    """Build a sampler function from config dict (KIE key → option name).
    
    Returns a function: sampler(rng) → dict of KIE values.
    """
    specs = {}
    for key, option_name in config['KIE'].items():
        spec = KIE_OPTIONS[key][option_name]
        specs[key] = spec

    def sampler(rng):
        kies = {}
        for key, spec in specs.items():
            if spec['dist'] == 'fixed':
                kies[key] = spec['value']
            elif spec['dist'] == 'uniform':
                kies[key] = rng.uniform(spec['low'], spec['high'])
            elif spec['dist'] == 'normal':
                kies[key] = rng.normal(spec['mean'], spec['std'])
        return kies
    return sampler


def compute_bulk_KIE(kies, sink_fracs):
    """Compute flux-weighted bulk KIE for ¹³C and D from individual sink KIEs.
    
    kies: dict with keys OH_13C, OH_D, Cl_13C, Cl_D, etc.
    sink_fracs: dict with keys OH, Cl, Strat, Soil (fractions summing to ~1)
    """
    kie_13C = (kies['OH_13C'] * sink_fracs['OH'] +
               kies['Cl_13C'] * sink_fracs['Cl'] +
               kies['Strat_13C'] * sink_fracs['Strat'] +
               kies['Soil_13C'] * sink_fracs['Soil'])
    kie_D = (kies['OH_D'] * sink_fracs['OH'] +
             kies['Cl_D'] * sink_fracs['Cl'] +
             kies['Strat_D'] * sink_fracs['Strat'] +
             kies['Soil_D'] * sink_fracs['Soil'])
    return kie_13C, kie_D


# ═══════════════════════════════════════════════════════════════════════════
# LIFETIME
# ═══════════════════════════════════════════════════════════════════════════
def compute_lifetime(years, config):
    """Compute CH₄ lifetime array from config.

    Returns np.ndarray of shape (len(years),).
    """
    lt_key = config['lifetime']
    lt = LIFETIME_OPTIONS[lt_key]
    years = np.asarray(years, dtype=float)
    if lt['mode'] == 'fixed':
        return np.full_like(years, lt['value'])
    elif lt['mode'] == 'linear':
        return lt['tau0'] + lt['slope'] * (years - lt['ref_year'])
    else:
        raise ValueError(f"Unknown lifetime mode: {lt['mode']}")


# ═══════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════
def find_data_dirs():
    """Locate data and source-signature directories."""
    base = Path(__file__).resolve().parent.parent
    original = base.parent / "TwoIsotopeBoxModel"
    rel_dir = original / "rel"
    return {
        'base': base,
        'data': rel_dir / "data",
        'src': rel_dir / "output",
    }


def load_CH4(data_dir, year_start=1999, year_end=2022):
    """Load global annual mean CH₄ and return years + values (ppb)."""
    CH4_raw = pd.read_excel(data_dir / "GML_CH4_AnnualMean.xlsx").to_numpy()
    # Find the rows spanning year_start to year_end+1 (we need n+1 for differencing)
    all_years = CH4_raw[:, 0].astype(float)
    idx_start = np.where(all_years == year_start)[0][0]
    idx_end = np.where(all_years == year_end)[0][0] + 1  # inclusive
    CH4_years = all_years[idx_start:idx_end].astype(float)
    CH4_values = CH4_raw[idx_start:idx_end, 1].astype(float)
    return CH4_years, CH4_values


def load_d13C_hemispheric(data_dir, year_start=1999, year_end=2022):
    """Load NH/SH δ¹³C annual averages."""
    C13data = pd.read_excel(data_dir / "ch4c13_nh_sh_mean.xlsx", header=None).to_numpy()
    dates = C13data[:, 0]
    glob_vals = C13data[:, 1]
    nh_vals = C13data[:, 2]
    sh_vals = C13data[:, 3]

    def annual_avg(dates, values):
        yrs = np.floor(dates).astype(int)
        unique = np.unique(yrs)
        out_y, out_v = [], []
        for y in unique:
            mask = yrs == y
            if mask.sum() >= 6:
                out_y.append(y)
                out_v.append(np.nanmean(values[mask]))
        return np.array(out_y), np.array(out_v)

    y_g, v_g = annual_avg(dates, glob_vals)
    _, v_nh = annual_avg(dates, nh_vals)
    _, v_sh = annual_avg(dates, sh_vals)

    i0 = np.where(y_g == year_start)[0][0]
    i1 = np.where(y_g == year_end)[0][0] + 1
    return y_g[i0:i1], v_g[i0:i1], v_nh[i0:i1], v_sh[i0:i1]


def load_d13C_iterations(data_dir):
    """Load δ¹³C MC iteration matrix [years × iterations]."""
    raw = np.loadtxt(str(data_dir / "d13C_dei_compiled.txt"))
    return raw[1:, 1:]  # skip header row and year column


def load_dD_iterations(data_dir):
    """Load δD MC iteration matrix + annual mean."""
    df = pd.read_excel(data_dir / "GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx")
    num = df.apply(pd.to_numeric, errors="coerce")
    ann_mean = num.iloc[:, 1].to_numpy(dtype=np.float64) - 0.5
    mc_matrix = num.iloc[:, 5:].to_numpy(dtype=np.float64)
    return ann_mean, mc_matrix


def load_source_signatures(src_dir, config):
    """Load all source-signature files and return a dict of arrays/dataframes."""
    sigs = {}

    # BB (always same files)
    sigs['BB_d13C'] = pd.read_csv(src_dir / "BB_d13C_annual.csv", header=None)
    sigs['BB_dD'] = pd.read_csv(src_dir / "BB_dD_annual.csv", header=None)

    # Microbial
    sigs['Mic_d13C_ann'] = pd.read_csv(src_dir / "Mic_d13C_annual.csv", header=None)
    sigs['Mic_dD_ann'] = pd.read_csv(src_dir / "Mic_dD_AnnGlob.csv", header=None)
    sigs['Mic_d13C_MC'] = pd.read_csv(src_dir / "Mic_d13C_MC.csv", header=None).iloc[:, 1:]
    sigs['Mic_dD_MC'] = pd.read_csv(src_dir / "Mic_dD_MC.csv", header=None).iloc[6:, 1:]

    # FF — select from catalog
    ff_key = config.get('FF_signature', 'EDGAR')
    ff_opt = FF_SIGNATURE_OPTIONS[ff_key]
    sigs['FF_d13C_MC'] = pd.read_csv(src_dir / ff_opt['d13C_file'], header=None if 'MC' in ff_opt['d13C_file'] else 0)
    sigs['FF_dD_MC'] = pd.read_csv(src_dir / ff_opt['dD_file'], header=None if 'MC' in ff_opt['dD_file'] else 0)

    # Also load the GlobUnc for fallback means
    sigs['FF_d13C_GlobUnc'] = pd.read_csv(src_dir / "FF_d13C_GlobUnc.csv")
    sigs['FF_dD_GlobUnc'] = pd.read_csv(src_dir / "FF_dD_GlobUnc.csv")

    return sigs


def load_BB_emissions(data_dir, config, n_years):
    """Load/compute BB emissions array for 2×2 models."""
    data_CT = pd.read_excel(data_dir / "CarbonTracker_CH4.xlsx")
    bbCT = data_CT.iloc[:, 9].values  # GFED4 prior

    mode = config.get('BB_mode', 'CT_GFED4_mean')
    if mode == 'CT_GFED4_mean':
        return np.full(n_years, np.mean(bbCT))
    elif mode == 'CT_GFED4_annual':
        if len(bbCT) >= n_years:
            return bbCT[:n_years]
        else:
            return np.full(n_years, np.mean(bbCT))
    elif mode == 'declining':
        bb_mean = np.mean(bbCT)
        end_frac = config.get('BB_end_fraction', 0.09)
        return np.concatenate([
            np.full(6, bb_mean),
            np.linspace(bb_mean, bb_mean * end_frac, n_years - 6)
        ])
    else:
        return np.full(n_years, np.mean(bbCT))


# ═══════════════════════════════════════════════════════════════════════════
# SOLUTION QUALITY MONITORING
# ═══════════════════════════════════════════════════════════════════════════
class QualityMonitor:
    """Track non-physical solutions and condition numbers."""

    def __init__(self, n_years, n_iter, label=""):
        self.label = label
        self.n_years = n_years
        self.n_iter = n_iter
        self.cond = np.zeros((n_years, n_iter))
        self.nonphys = np.zeros((n_years, n_iter), dtype=bool)
        self.n_negative = 0
        self.n_total = 0

    def record_cond(self, j, k, A):
        self.cond[j, k] = np.linalg.cond(A)

    def record_negative(self):
        self.n_negative += 1
        self.n_total += 1

    def record_ok(self):
        self.n_total += 1

    def summary(self):
        pct = 100.0 * self.n_negative / max(1, self.n_total)
        mean_c = np.mean(self.cond[self.cond > 0]) if np.any(self.cond > 0) else 0
        report = {
            'label': self.label,
            'total_solves': self.n_total,
            'negative_pct': round(pct, 2),
            'mean_cond': round(float(mean_c), 1),
        }
        print(f"  [{self.label}] {pct:.1f}% negative | mean cond = {mean_c:.1f}")
        return report


# ═══════════════════════════════════════════════════════════════════════════
# SMOOTHING
# ═══════════════════════════════════════════════════════════════════════════
def smooth_5yr(arr_2d):
    """5-year running mean along axis 0. Input: [years × iterations]."""
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


# ═══════════════════════════════════════════════════════════════════════════
# UTILITIES
# ═══════════════════════════════════════════════════════════════════════════
def pad_to_length(arr, length):
    """Pad array to `length` by repeating first element at front."""
    arr = np.asarray(arr, dtype=float).flatten()
    if len(arr) >= length:
        return arr[:length]
    pad = np.full(length - len(arr), arr[0])
    return np.concatenate([pad, arr])[:length]


def save_quality_report(report_dict, out_dir, filename='quality_report.json'):
    with open(Path(out_dir) / filename, 'w') as f:
        json.dump(report_dict, f, indent=2)


def trend_analysis(compiled, model_years, base_slice=slice(6, 9), recent_slice=slice(-3, None)):
    """Compute emission trends. Returns dict with stats."""
    delta = compiled[recent_slice, :].mean(axis=0) - compiled[base_slice, :].mean(axis=0)
    return {
        'mean': float(np.nanmean(delta)),
        'std': float(np.nanstd(delta)),
        'pct_positive': float((delta > 0).sum() / len(delta) * 100),
    }
