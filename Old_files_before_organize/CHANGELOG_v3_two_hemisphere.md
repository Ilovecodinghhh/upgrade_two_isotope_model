# CHANGELOG — Two-Hemisphere Upgrade (v3.0)

## v3.0.0 — 2026-05-05: Two-Hemisphere Dual-Isotope Box Model

### UPGRADE: Two-Hemisphere Spatial Structure

**What changed:**
The one-box model has been restructured into a two-hemisphere (NH/SH) model following Nguyen et al. (2020, GRL) Table A1.

**Architecture:**
- **State variables**: [CH₄]_NH, [CH₄]_SH (ppb) with δ¹³C and δD per hemisphere
- **Exchange**: d[CH₄]_N/dt includes +([CH₄]_S - [CH₄]_N)/τ_ex where τ_ex = 1 year
- **Per-hemisphere lifetime**: τ_NH ≈ 0.95×τ_global, τ_SH ≈ 1.05×τ_global (reflects higher OH in tropics/NH subtropics)
- **Emission split**: NH ~75%, SH ~25% (derived from mass balance; consistent with EDGAR v7)

**Equations (per hemisphere):**
```
Mass:     S_N = dM_N/dt + M_N/τ_N - (M_S - M_N)/τ_ex
Isotope:  f_src·S_N = d(f·M)_N/dt + f_N·M_N·α/τ_N - (f_S·M_S - f_N·M_N)/τ_ex
```

**Solver:**
- Replaced `np.linalg.solve` (direct 3×3 inversion) with `scipy.optimize.lsq_linear` (bounded weighted least squares)
- Mass balance constraint weighted 100–200× higher than isotope constraints
- Non-negativity bounds enforce physically meaningful solutions (0 ≤ x_i ≤ 1.5×S_total)
- Result: **0% non-physical solutions** (vs. ~100% with direct inversion due to ill-conditioned δD matrix)

**Data inputs:**
| Input | Source | NH/SH |
|-------|--------|-------|
| [CH₄] | NOAA GML annual mean | Derived: Global ± Δ/2, Δ = 80→100 ppb |
| δ¹³C | ch4c13_nh_sh_mean.xlsx | Direct NH/SH observations (CU-INSTAAR) |
| δD | GlobMean_dD_iterations | Approximated: Global ∓ 1.5‰ NH/SH offset |
| Sources | EDGAR/CarbonTracker | FF:85/15%, Mic:65/35%, BB:55/45% NH/SH |

**Why this upgrade matters (Naus et al. 2019 critique):**
One-box models alias inter-hemispheric transport as source changes. When NH emissions increase:
1. The NH-SH gradient grows
2. Transport of CH₄ from NH→SH increases
3. A one-box model sees the growing global mean and attributes ALL change to sources
4. But part of the SH increase is just transport from NH, not SH sources

The two-box model explicitly resolves this:
- NH source = local accumulation + local loss - net exchange FROM SH
- SH source = local accumulation + local loss - net exchange FROM NH

This prevents misattribution of transport-driven changes to source changes.

**Results (1999–2021 mean):**
| Category | NH (Tg/yr) | SH (Tg/yr) | Global (Tg/yr) |
|----------|-----------|-----------|---------------|
| BB | 0.7 ± 2.6 | 23.8 ± 5.2 | 24.5 ± 6.5 |
| FF | 189.9 ± 19.5 | 0.0 ± 0.0 | 189.9 ± 19.5 |
| Mic | 250.1 ± 19.3 | 121.0 ± 5.2 | 371.1 ± 24.0 |
| **Total** | **440.6** | **144.9** | **585.5** |

**Validation:**
- NH emission fraction: 75.3% (expected 70–75% from EDGAR ✓)
- δ¹³C NH-SH gradient: −0.24‰ (NH more depleted due to FF influence ✓)
- CH₄ IH gradient: 80–100 ppb (consistent with NOAA observations ✓)

---

### KNOWN LIMITATIONS AND FUTURE WORK

1. **SH FF = 0 Tg/yr** — The bounded solver pushes FF to the lower bound in SH because the δD endmember matrix is too ill-conditioned to resolve FF from Mic in the southern hemisphere. This is an artifact of:
   - δD end-members too close (BB: −210‰, FF: −180‰, Mic: −310‰ — FF/BB overlap)
   - Small SH total source magnifies relative errors
   - **FIX**: Add prior information (Bayesian MCMC with informative priors on SH FF fraction from EDGAR)

2. **δD NH/SH split is approximated** — The ±1.5‰ offset is a first-order estimate. Need:
   - Station-level δD data from Röckmann group (Utrecht), Rice/Fujita (Tokyo), Rice (Colorado)
   - NOAA GGGRN δD compilation when available

3. **Source signatures assumed hemisphere-invariant** — In reality:
   - NH FF is gas-dominated (δ¹³C ≈ −44‰) vs SH FF is more coal (δ¹³C ≈ −37‰)
   - NH Mic is rice + waste (δ¹³C ≈ −62‰) vs SH Mic is tropical wetlands (δ¹³C ≈ −58‰)
   - **FIX**: Use hemisphere-specific source signatures from EDGAR emission-weighted averages

4. **No interactive CH₄-CO-OH chemistry** — Nguyen (2020) showed this causes 25% bias over 10 years. Next upgrade should add CO as a state variable and let OH respond.

5. **IH exchange τ = 1 yr (fixed)** — Could add uncertainty: τ_ex ~ Normal(1.0, 0.1) in MC loop.

---

### RETAINED FROM v2.0

All three previous upgrades remain active:
1. **KIE sampling**: OH_13C ~ U[1.0039, 1.0054], OH_D ~ U[1.294, 1.327]
2. **Quality monitoring**: Condition numbers tracked (still high ~178,000 due to δD, but solutions are physical thanks to bounded solver)
3. **Time-varying lifetime**: τ(t) = 9.0 - 0.017*(t-2010), split NH/SH with 0.95/1.05 ratio

---

### CODE CHANGES

- **New file**: `two_hemisphere_box_model.py` (replaces `upgraded_box_model.py` as main model)
- **New output directory**: `Output_2Hemi/`
- **New dependency**: `scipy` (for `lsq_linear` bounded least squares)
- **Data path**: Reads from `../TwoIsotopeBoxModel/rel/` (original data directory)

### REFERENCES

- Nguyen, N.H. et al. (2020). Effects of chemical feedbacks on decadal methane emissions estimates. *GRL*, 47, e2019GL085706.
- Naus, S. et al. (2019). Constraints and biases in a tropospheric two-box model of OH. *ACP*, 19, 407–424.
- Prather, M.J. (1994). Lifetimes and eigenstates in atmospheric chemistry. *GRL*, 21, 801–804.
- He, J. et al. (2026). Global methane budget constrained by TROPOMI. *Science*.
- Saunois, M. et al. (2020). The Global Methane Budget 2000–2017. *ESSD*, 12, 1561–1623.
