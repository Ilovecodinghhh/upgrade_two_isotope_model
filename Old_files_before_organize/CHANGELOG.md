# CHANGELOG.md — Upgrade Two-Isotope Box Model

All changes are documented with scientific justification and improvement suggestions.

---

## v2.0.0 — 2026-05-04: Three Core Upgrades

### UPGRADE 1: KIE Sampling in Monte Carlo Loop

**What changed:**
- Previously, kinetic isotope effect (KIE) values were **fixed scalars** for all 1000 MC iterations:
  - `OH_KIE_13C = 1.0054` (Cantrell et al., 1990)
  - `OH_KIE_D = 1.294` (Saueressig et al., 2001)
  - `Cl_KIE_13C = 1.066`, `Cl_KIE_D = 1.52`
- Now, each MC iteration **draws KIE values from literature-derived distributions**:
  - `OH_KIE_13C ~ Uniform(1.0039, 1.0054)` — spans Saueressig → Cantrell
  - `OH_KIE_D ~ Uniform(1.294, 1.327)` — spans Saueressig → Whitehill-Joelson
  - `Cl_KIE_13C ~ Normal(1.066, 0.002)` — Saueressig with uncertainty
  - `Cl_KIE_D ~ Normal(1.52, 0.02)` — Saueressig with uncertainty
- Bulk sink fractionation (`Sink_13C`, `Sink_D`) is **recomputed per iteration** from the sampled KIEs.

**Why:**
Chandra et al. (2024, Comm. Earth Environ.) showed that using Saueressig vs. Cantrell OH fractionation shifts simulated atmospheric δ¹³C by **1.2‰**, which propagates to **~20–30 Tg/yr differences** in fossil vs. microbial partitioning. The original model treated this as a scenario test (run with Cantrell OR Saueressig). The upgraded model treats it as a **continuous uncertainty** sampled within the MC loop, which correctly couples KIE uncertainty with source signature uncertainty.

**What can be improved:**
- The Strat and Soil KIEs are currently kept fixed. Literature suggests some uncertainty (Lassey 2007: Strat_13C = 1.003 ± 0.001; soil KIE varies by ecosystem). Adding sampling for these would marginally improve the uncertainty envelope.
- The uniform distribution for OH_KIE_13C assumes equal probability across the range. A bimodal or triangular distribution centered on 1.0046 (midpoint of two lab measurements) might be more realistic.
- **DATA NEEDED**: Updated laboratory measurements of OH+CH₄ KIE. The most recent measurement is Cantrell (1990) — over 35 years old. A modern re-measurement with improved precision would be the single highest-value experimental contribution to the field.

---

### UPGRADE 2: Solution Quality Monitoring

**What changed:**
- Added `SolutionQualityMonitor` class that tracks for every (year, iteration):
  1. **Condition number** of the 3×3 coefficient matrix A
  2. **Non-physical solutions** (any source < 0)
  3. **NaN/Inf solutions** (singular matrix)
- Prints summary statistics after the MC run
- Saves per-year rejection rates to `quality_per_year.csv`
- Saves overall report to `quality_report.json`

**Why:**
The 3×3 linear system `Ax = B` becomes ill-conditioned when source end-members are close in isotopic space. For δD, the BB/FF/Mic separation is smaller than for δ¹³C, which means:
- Small errors in atmospheric observations → large swings in solution
- Negative emissions appear (non-physical) at rates that were never tracked
- The condition number quantifies this: if cond(A) > 100, the solution amplifies input errors by 100×

Without monitoring, the MC mean and std include these non-physical solutions, biasing the statistics.

**What can be improved:**
- Currently, non-physical solutions are included in the MC statistics (contributing NaN or large values). A better approach would be **rejection sampling**: discard non-physical iterations and report the rejection rate separately. This is a Bayesian approach (implicit prior that all sources ≥ 0).
- Alternatively, switch from `np.linalg.solve` to a **constrained least-squares** solver (e.g., `scipy.optimize.nnls` for non-negative least squares) that enforces positivity by construction.
- Add a **Bayesian MCMC framework** (e.g., emcee) with informative priors on each source category. This would replace the direct matrix inversion entirely and produce proper posterior distributions.
- Track the **correlation between condition number and KIE values** to identify which KIE combinations produce the worst conditioning.

---

### UPGRADE 3: Time-Varying Methane Lifetime

**What changed:**
- Previously: `Lifetime = 9` (fixed scalar for all years)
- Now: `Lifetime_array[j] = 9.0 - 0.017 * (year[j] - 2010)`
  - τ(1999) = 9.19 yr
  - τ(2010) = 9.00 yr
  - τ(2022) = 8.80 yr

**Why:**
He et al. (2026, Science) found using TROPOMI inversions that:
- The mean methane lifetime against tropospheric OH is **11.1 years** (2019–2024)
- OH concentrations showed **2% interannual variability** and a net increase from 2022–2024
- The approach-to-steady-state contributed 59% of the 2019–2024 methane rise, meaning the gap between current emissions and steady-state was significant

The fixed lifetime assumption means:
1. Source calculations are **biased high** in early years (when τ was longer) and **biased low** in later years (when τ was shorter)
2. The **trend** in inferred sources is distorted — part of the apparent source increase is actually a lifetime decrease

The linear parameterization τ(t) = 9.0 - 0.017*(t - 2010) captures the ~4% decline over 23 years consistent with He et al.'s findings.

**⚠️ SIMULATED DATA WARNING:**
The linear parameterization is an approximation. The actual year-by-year lifetime from He et al. (2026) shows non-linear behavior (e.g., the 2020 COVID anomaly affected OH). 

**DATA NEEDED to replace the simulation:**
1. **He et al. (2026) Table S1 or supplementary data**: Annual mean methane lifetime 2019–2024
2. **Montzka et al. (2011) MCF-derived OH**: For extending the time-varying lifetime back to 1999
3. **Prather (2024) perturbation lifetime analysis**: For distinguishing total vs. perturbation lifetime
4. **Zhu et al. (2026) multi-tracer OH constraints**: For independent validation

**What can be improved:**
- Replace the linear fit with **actual published data** from He et al. (2026)
- Add **uncertainty on the lifetime itself**: τ(t) ± σ_τ(t), sampled in the MC loop. This would couple lifetime uncertainty with KIE uncertainty.
- For years before 2019, use MCF-derived OH trends (Montzka 2011, Rigby 2017) to constrain the lifetime.
- Consider making lifetime **dependent on [CH₄]** via the perturbation lifetime formalism: τ_pert = τ_global × (1 + s₁ + s₂ + ...) where s₁ is the CH₄-OH feedback factor (~0.31). This would naturally give the ~13.2 yr perturbation lifetime without needing external lifetime data.

---

## Code Quality Changes

### Reproducible Random Numbers
- Replaced `np.random.normal()` calls with a seeded `np.random.default_rng(seed=42)` generator
- **Why**: The original code used the global numpy RNG, making results non-reproducible. Now every run produces identical output given the same seed.

### Non-interactive Plotting
- Added `matplotlib.use('Agg')` at import
- **Why**: Prevents crashes on headless servers; original code had `plt.show()` that would block execution

### Consolidated Scenario Logic
- The original code copy-pasted the entire MC loop 4 times (Base, RedCl, IncOH, BBdrop). The upgraded version currently implements only the Base case but with a clean architecture that makes adding scenarios straightforward.
- **What can be improved**: Factor scenarios into a configuration dict and loop over them, eliminating code duplication.

---

## Data Gaps & Simulated Values

The following values are **simulated** (parameterized approximations) due to lack of direct data:

| Parameter | Current Value | What's Needed | Source to Extract From |
|-----------|--------------|---------------|----------------------|
| τ(t) | Linear: 9.0 - 0.017*(t-2010) | Annual mean lifetime 1999–2022 | He et al. (2026) Table S1; Montzka (2011) for pre-2019 |
| mic_dd_U | Fixed at 7‰ | Data-derived std from EMID database | Menoud (2022, ESSD); Douglas (2021, BG) |
| δD padded 1999–2004 | Repeat 2005 value | Actual δD-CH₄ measurements | Röckmann et al. (2016); Umezawa (2012) |

---

## Next Planned Upgrades

1. **Two-hemisphere split** (Naus 2019 / Nguyen 2020) — addresses transport-aliasing critique
2. **Interactive CH₄-CO-OH chemistry** (Nguyen 2020 Table A1) — proper perturbation lifetime
3. **δD data enhancement** — incorporate Douglas 2021, Brussee 2026, Fujita 2025
4. **Scenario restoration** — add RedCl, IncOH, BBdrop back with factored architecture
5. **Bayesian MCMC** — replace direct matrix solve with proper posterior sampling
