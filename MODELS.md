# Model Version Comparison

This document summarises every model version in the repository, how they differ, and traces the data lineage for output files.

---

## Version Overview

| Version | File | Core Approach | Isotopes | Geometry | Solver | Key Innovation |
|---------|------|--------------|----------|----------|--------|----------------|
| **v1.0** | `v1.0_d13C_dD_MassBalance.py` | Deterministic 2×2 mass-balance | δ¹³C + δD | Global 1-box | `np.linalg.solve` (2×2) | Original baseline — fixed KIE values |
| **v2.0** | `v2.0_upgraded_box_model.py` | Stochastic 2×2 mass-balance | δ¹³C + δD | Global 1-box | `np.linalg.solve` (2×2) | KIE sampling, time-varying lifetime, solution-quality monitoring |
| **v3.0** (two files) | `v3.0_bayesian_box_model.py` + `v3.0_two_hemisphere_box_model.py` | Bayesian MCMC (PyMC) **or** 2-hemisphere 3×3 | δ¹³C + δD | 1-box (Bayesian) or NH/SH 2-box | PyMC NUTS / `np.linalg.solve` (3×3) | Prior-based inference; hemispheric split with inter-hemispheric exchange |
| **v3.1** (two files) | `v3.1_bayesian_negative.py` + `v3.1_optimized_3x3.py` | Bayesian with sign constraints **or** optimised 3×3 | δ¹³C + δD | 1-box / NH+SH | PyMC with `pm.Potential` / `np.linalg.solve` (3×3) | Penalise negative (non-physical) emissions; condition-number guard |
| **v3.2** | `v3.2_bb_fixed_2x2.py` | Fixed-BB diagnostic 2×2 | δ¹³C + δD | Global 1-box | `np.linalg.solve` (2×2) | Fixes biomass-burning to external estimate; solves only for FF + Mic |
| **v3.3** | `v3.3_dD_comparison.py` | Multi-version δD overlay | δD only | Global 1-box | (reads pre-computed results) | Comparison/visualisation of δD from v1.0–v3.2 |
| **v4.0** | `v4.0_mic_vs_nonmic.py` | 2-source partitioning (microbial vs non-microbial) | δ¹³C + δD | Global 1-box | `np.linalg.solve` (2×2) | Lumps FF+BB→"non-microbial"; solves Mic vs NonMic |

---

## Detailed Comparison

### v1.0 — Deterministic Two-Isotope Mass Balance

**File:** `v1.0_d13C_dD_MassBalance.py`

**What it does:**
- Classic 2×2 isotopic mass-balance to partition three CH₄ sources (fossil-fuel, microbial, biomass-burning) using δ¹³C and δD observations.
- Monte Carlo (N = 1 000) propagates observational uncertainty in atmospheric δ¹³C/δD and source signatures.
- KIE (Kinetic Isotope Effect) values for OH, Cl, and stratospheric sinks are **fixed constants** (e.g., OH-¹³C KIE = 1.0039).
- Atmospheric lifetime τ = **9 yr** (fixed).
- The 2×2 system solves for two unknowns (FF and Mic fractions) after subtracting the BB contribution; BB is drawn from its own uncertainty distribution.

**Outputs:** Per-year median ± uncertainty for FF, Mic, BB emissions; δ¹³C and δD source signatures.

---

### v2.0 — Upgraded Box Model

**File:** `v2.0_upgraded_box_model.py`

**Differences from v1.0:**

| Feature | v1.0 | v2.0 |
|---------|------|------|
| KIE values | Fixed constants | **Sampled from literature distributions** each MC iteration (OH-¹³C ~ Uniform(1.0039, 1.0054); OH-D ~ Uniform(1.294, 1.327); Cl-¹³C ~ N(1.066, 0.002); Cl-D ~ N(1.52, 0.02)) |
| CH₄ lifetime | τ = 9 yr (fixed) | **Time-varying:** τ(t) = 9.0 − 0.017·(t − 2010), following He et al. (2026) |
| Solution monitoring | None | **Condition-number tracking**, non-physical solution counting, rejection statistics |
| Debug mode | No | `--debug` flag writes per-iteration diagnostics to JSON |

**Why it matters:** The dominant systematic uncertainty in isotope-based source attribution is the KIE of OH oxidation. Saueressig et al. (2001) and Cantrell et al. (1990) differ by ~1.2‰ for ¹³C. v2.0 propagates this by sampling across the literature range, typically widening uncertainty bands by 20–40%.

---

### v3.0 — Two Variants

#### v3.0a — Bayesian Box Model (PyMC)

**File:** `v3.0_bayesian_box_model.py`

**Differences from v2.0:**

| Feature | v2.0 | v3.0a (Bayesian) |
|---------|------|-------------------|
| Statistical framework | Frequentist MC | **Bayesian MCMC** (PyMC + NUTS sampler) |
| Priors | Implicit (uniform draws) | Explicit priors: FF ~ N(100, 30), Mic ~ N(230, 50), BB ~ N(30, 15) Tg yr⁻¹; source δ values ~ Normal |
| KIE | Sampled in MC loop | Sampled as PyMC random variables |
| Output | Point estimates + percentiles | **Full posterior distributions**, trace plots, HDI intervals |
| Negative-emission handling | Post-hoc filter | Implicit (but not constrained — see v3.1) |

#### v3.0b — Two-Hemisphere Box Model

**File:** `v3.0_two_hemisphere_box_model.py`

**Differences from v2.0:**

| Feature | v2.0 | v3.0b (2-hemisphere) |
|---------|------|----------------------|
| Spatial geometry | 1 global box | **2 boxes** (Northern Hemisphere + Southern Hemisphere) |
| Linear system | 2×2 | **3×3** (NH-FF, NH-Mic, SH-Mic) with inter-hemispheric exchange flux |
| Observations used | Global mean δ¹³C, δD | **Separate NH and SH** δ¹³C, δD, plus total CH₄ |
| New parameter | — | `k_ex` = inter-hemispheric exchange rate (~1.1 yr⁻¹) |
| BB treatment | Random draw | **Split into NH-BB and SH-BB** from separate distributions |

---

### v3.1 — Two Variants

#### v3.1a — Bayesian with Negativity Penalty

**File:** `v3.1_bayesian_negative.py`

**Differences from v3.0a:**

| Feature | v3.0a | v3.1a |
|---------|-------|-------|
| Non-physical solutions | Not constrained | **`pm.Potential` penalty** applies a large negative log-likelihood when any source emission goes below zero |
| Prior widths | Fixed | **Tighter priors** to reduce posterior mass in non-physical region |
| Diagnostics | Basic traces | Enhanced reporting of R̂, ESS, divergence counts |

#### v3.1b — Optimised 3×3 System

**File:** `v3.1_optimized_3x3.py`

**Differences from v3.0b:**

| Feature | v3.0b | v3.1b |
|---------|-------|-------|
| Condition-number guard | None | **Skips iterations** where cond(A) > threshold (default 100) |
| NH/SH source signatures | Identical | **Separate NH and SH** δ¹³C and δD signature distributions for FF and Mic |
| Exchange rate | Fixed `k_ex` | `k_ex` **sampled** from N(1.1, 0.15) |
| Output | Emissions only | Also outputs **per-hemisphere δ¹³C/δD reconstructions** |

---

### v3.2 — Fixed-BB Diagnostic

**File:** `v3.2_bb_fixed_2x2.py`

**Differences from v2.0:**

| Feature | v2.0 | v3.2 |
|---------|------|------|
| BB emissions | Sampled from uncertainty distribution each MC iteration | **Fixed to external inventory** (e.g., GFED-based estimate per year) |
| Unknown variables | 2 (FF-fraction, Mic-fraction) | 2 (same), but BB is **prescribed**, not random |
| Purpose | Full 3-source partition | **Sensitivity test**: how much does BB uncertainty drive total uncertainty? |

This version quantifies the uncertainty reduction when BB is independently constrained (e.g., from satellite fire data).

---

### v3.3 — δD Comparison Across Versions

**File:** `v3.3_dD_comparison.py`

This is **not a model** per se — it is a **plotting/comparison script** that:
1. Runs or loads δD results from v1.0, v2.0, v3.0b, v3.1b, and v3.2.
2. Overlays the δD time series from each version on a single figure.
3. Highlights where versions diverge most (typically in the 2006–2020 period where lifetime changes and KIE sampling have the largest effect).

There is also `compare_dD_all_versions.py` which performs a similar cross-version overlay with additional statistics (RMSE, bias, spread).

---

### v4.0 — Microbial vs Non-Microbial

**File:** `v4.0_mic_vs_nonmic.py`

**Differences from v2.0:**

| Feature | v2.0 | v4.0 |
|---------|------|------|
| Source categories | 3 (FF, Mic, BB) | **2** (Microbial, Non-Microbial = FF + BB) |
| Linear system | 2×2 (solve FF, Mic after removing BB) | **2×2** (solve Mic, NonMic directly) |
| BB handling | Separate random variable | **Lumped into non-microbial** |
| Source signatures | 3 sets (FF, Mic, BB) | 2 sets (Mic, NonMic where δ_NonMic is a weighted average of FF + BB signatures) |
| Purpose | Full partition | **Robust 2-source split** when BB/FF separation is unreliable |

This version is useful when the goal is simply to assess the microbial contribution (e.g., for wetland/agriculture policy) without needing to separately resolve fossil and fire sources.

---

## Data Lineage — `rel/output/` Files

All files in `rel/output/` are **pre-computed Monte Carlo simulation results** used as input data by the Python model scripts. They are **not generated** by scripts in this repository. Each CSV contains 1 000 Monte Carlo realisations (columns) for each year (rows), with the first column being the year.

| File | Contents | Data Source / How Calculated |
|------|----------|------------------------------|
| `BB_d13C_annual.csv` | Annual biomass-burning δ¹³C source signature (1000 MC draws × years) | Derived from GFED (Global Fire Emissions Database) fire-type fractions weighted by isotopic signatures of C3/C4 vegetation burning. MC draws sample vegetation-type uncertainty. |
| `BB_dD_annual.csv` | Annual biomass-burning δD source signature (1000 MC draws × years) | Same approach as above but for deuterium; based on pyrolysis δD measurements from literature (e.g., Whiticar 1999, Sherwood et al. 2017). |
| `FF_d13C_GlobUnc.csv` | Annual fossil-fuel δ¹³C with global uncertainty (1000 MC draws × years) | Based on thermogenic methane δ¹³C from EDGAR emission categories (coal, oil, gas) weighted by their relative contributions. MC samples category-fraction and per-category isotopic uncertainty. |
| `FF_dD_GlobMC_CTCH4.csv` | Annual fossil-fuel δD with global uncertainty (1000 MC draws × years) | Same structure for δD; "CTCH4" indicates the source inventory is CT-CH4 or a compatible product. Based on thermogenic δD measurements. |
| `Mic_d13C_MC.csv` | Annual microbial δ¹³C source signature (1000 MC draws × years) | Weighted average of wetland, ruminant, rice, and waste δ¹³C. MC propagates uncertainty in sub-source fractions and their isotopic values (e.g., Dlugokencky et al., Schaefer et al.). |
| `Mic_dD_MC.csv` | Annual microbial δD source signature (1000 MC draws × years) | Same as above for δD; microbial (biogenic) methane δD is strongly influenced by local water δD and methanogenic pathway fractionation. |
| `Mic_dD_AnnGlob.csv` | Annual global-mean microbial δD (1000 MC draws × years) | Similar to `Mic_dD_MC.csv` — likely a variant with different aggregation (annual global mean vs. sub-annual). Used by specific model versions for comparison. |

**How they are used in the models:** Each Python script reads these CSV files via `SRC_DIR = REL_DIR / "output"`. Within the Monte Carlo loop, iteration `i` reads column `i` from each CSV to obtain that iteration's source signatures, ensuring correlated sampling across all source types within each realisation.

---

## Summary of Evolution

```
v1.0  Fixed KIE, fixed τ, deterministic 2×2
  │
  ▼
v2.0  + KIE sampling + time-varying τ + diagnostics
  │
  ├──► v3.0a  Bayesian (PyMC MCMC) — full posteriors
  │      │
  │      ▼
  │    v3.1a  + negativity penalty + tighter priors
  │
  ├──► v3.0b  Two-hemisphere 3×3
  │      │
  │      ▼
  │    v3.1b  + condition guard + separate NH/SH signatures
  │
  ├──► v3.2   Fixed-BB diagnostic (sensitivity test)
  │
  ├──► v3.3   δD comparison overlay (visualisation)
  │
  └──► v4.0   2-source (Mic vs Non-Mic) simplification
```

Each branch addresses a different modelling concern: **v3.0a/3.1a** tackle statistical rigour, **v3.0b/3.1b** tackle spatial resolution, **v3.2** isolates BB uncertainty, and **v4.0** simplifies the source partition for policy-relevant questions.

---

## Three-Box Models (NHext / Trop / SHext)

*Branch: `three-box`*

Two new models extend the NH/SH 2-box framework to three latitude bands:
- **NHext** (>30°N): boreal / high-latitude sources
- **Trop** (30°S–30°N): dominant tropical sources
- **SHext** (<30°S): Southern Hemisphere extratropics

### 2×2_three.py — BB-Fixed, Separate Isotopes, Three-Box

| Feature | 2×2_two.py (2-box) | **2×2_three.py (3-box)** |
|---------|-------------------|--------------------------|
| Geometry | NH + SH | **NHext + Trop + SHext** |
| Exchange | 1 path (NH↔SH, τ_ex~N(1.0,0.1)) | **2 paths** (NHext↔Trop τ_NT~N(0.8,0.1); Trop↔SHext τ_TS~N(1.2,0.1)) |
| δ¹³C isotopes | Real NH/SH obs | Interpolated to 3 boxes from NH/SH |
| δD atmospheric | Real NH/SH MC | **Real 3-box obs** (ThreeBox_atm_dD_annual.csv, 2005–2024) |
| δD source sigs | NH/SH MC (1000 iter) | **Per-box MC** (NHext/Trop/SHext, 1000 iter each) |
| BB split | NH 55%, SH 45% | **NHext 30%, Trop 55%, SHext 15%** (GFED4) |
| Solver | Analytic 2×2 per box | Analytic 2×2 per box per isotope |
| Sink fractions | NH/SH | **Box-specific** (tropics: more OH; extratropics: more strat) |
| Lifetime ratio | NH 0.95×, SH 1.05× | **NHext 1.05×, Trop 0.90×, SHext 1.08×** |

### 3×3_three.py — Coupled Dual-Isotope, Three-Box

| Feature | 3×3_two.py (2-box) | **3×3_three.py (3-box)** |
|---------|-------------------|--------------------------|
| Geometry | NH + SH | **NHext + Trop + SHext** |
| Exchange | 1 path | **2 paths** (NHext↔Trop; Trop↔SHext) |
| System | 3×3 coupled (BB, FF, Mic) | **3×3 coupled per box** |
| Solver | `lsq_linear` (bounded LS) | `lsq_linear` per box with **box-specific weights** |
| δD source sigs | NH/SH MC | **Per-box MC** |
| δ¹³C source sigs | Global (same for NH/SH) | Global (same for all boxes) |
| Condition numbers | NH/SH | **Per-box** (NHext, Trop, SHext reported separately) |

### New Data Files Used

| File | Description | Coverage |
|------|-------------|----------|
| `ThreeBox_atm_dD_annual.csv` | Station-derived atmospheric δD annual means | 2005–2024, 3 boxes |
| `FF_dD_{NHext,Trop,SHext}_MC.csv` | Fossil-fuel δD source signatures | 1998–2021, 1000 MC × 24 yr |
| `BB_dD_{NHext,Trop,SHext}_MC.csv` | Biomass-burning δD source signatures | 1998–2021, 1000 MC × 24 yr |
| `Mic_dD_{NHext,Trop,SHext}_MC.csv` | Microbial δD source signatures | 1998–2021, 1000 MC × 24 yr |
| `ThreeBox_dD_sources_summary.csv` | Summary stats (mean ± std per box per year) | 1998–2021 |

### New Parameters in `common.py`

```python
SINK_FRACTIONS_NHEXT = {'OH': 0.810, 'Cl': 0.040, 'Strat': 0.080, 'Soil': 0.070}
SINK_FRACTIONS_TROP  = {'OH': 0.860, 'Cl': 0.035, 'Strat': 0.055, 'Soil': 0.050}
SINK_FRACTIONS_SHEXT = {'OH': 0.840, 'Cl': 0.025, 'Strat': 0.080, 'Soil': 0.055}

LIFETIME_RATIO_NHEXT = 1.05   # slower destruction in extratropics
LIFETIME_RATIO_TROP  = 0.90   # faster destruction in tropics (more OH)
LIFETIME_RATIO_SHEXT = 1.08

BB_NHEXT_FRACTION = 0.30; BB_TROP_FRACTION = 0.55; BB_SHEXT_FRACTION = 0.15

TAU_EX_NT_MEAN = 0.8; TAU_EX_TS_MEAN = 1.2  # years
```

### Updated Model Evolution Tree

```
v1.0  Fixed KIE, fixed τ, deterministic 2×2
  │
  ▼
v2.0  + KIE sampling + time-varying τ + diagnostics
  │
  ├──► v3.0b  Two-hemisphere 3×3 (NH/SH)
  │      │
  │      ▼
  │    v3.1b  + condition guard + separate NH/SH signatures
  │      │
  │      ├──► 2×2_two.py  BB-fixed, separate isotopes, 2-box
  │      │       │
  │      │       ▼
  │      │     **2×2_three.py  BB-fixed, separate isotopes, 3-box** ← NEW
  │      │
  │      └──► 3×3_two.py  Coupled dual-isotope, 2-box
  │              │
  │              ▼
  │            **3×3_three.py  Coupled dual-isotope, 3-box** ← NEW
  │
  └──► v4.0   2-source (Mic vs Non-Mic)
```
