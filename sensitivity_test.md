# Sensitivity Test Matrix

A complete inventory of every tunable input in the four-model suite (`2x2_one`, `2x2_two`, `3x3_one`, `3x3_two`), with possible variations, physical justification, and literature references.

---

## 1. Model Architecture

The choice of model structure is itself the first-order sensitivity test.

| Model | System | Geometry | BB Treatment | Solver |
|-------|--------|----------|--------------|--------|
| `2x2_one` | Separate 2×2 (δ¹³C and δD solved independently) | Global 1-box | Fixed (CarbonTracker) | Analytic formula |
| `2x2_two` | Separate 2×2 per hemisphere | NH/SH 2-box | Fixed (CarbonTracker × GFED NH/SH split) | Analytic formula |
| `3x3_one` | Simultaneous 3×3 (δ¹³C + δD coupled) | Global 1-box | Free variable (solved) | `np.linalg.solve` |
| `3x3_two` | Simultaneous 3×3 per hemisphere | NH/SH 2-box | Free variable (solved) | Bounded LS (`scipy.optimize.lsq_linear`) |

**Key structural trade-offs:**
- 2×2 models are well-conditioned but require external BB → cannot independently test BB trends.
- 3×3 models solve BB freely but the δD row is ~100× smaller in absolute scale than the δ¹³C row, causing ill-conditioning (κ ~ 10⁴–10⁵).
- 2×2 models yield *independent* FF/Mic estimates from each isotope — their disagreement quantifies sensitivity to KIE and sink assumptions.

---

## 2. Kinetic Isotope Effect (KIE)

**CLI flag:** `--kie fixed | sampled`

The KIE of each sink reaction determines how the atmospheric isotope ratio relates to the source mixture. It is the single largest systematic uncertainty in isotope-based source partitioning.

### 2.1 Sink-Specific KIE Values

| Sink | Isotope | Distribution (sampled mode) | Fixed Value | References |
|------|---------|----------------------------|-------------|------------|
| OH | ¹³C | Uniform(1.0039, 1.0054) | 1.00465 | Saueressig et al. (2001); Cantrell et al. (1990) |
| OH | D | Uniform(1.294, 1.327) | 1.3105 | Saueressig et al. (2001); Gierczak et al. (1997) |
| Cl | ¹³C | N(1.066, 0.002) | 1.066 | Saueressig et al. (1995); Tyler et al. (2000) |
| Cl | D | N(1.52, 0.02) | 1.52 | Saueressig et al. (1996) |
| Stratospheric | ¹³C | N(1.003, 0.001) | 1.003 | Rice et al. (2003) |
| Stratospheric | D | N(1.179, 0.01) | 1.179 | Rice et al. (2003); Röckmann et al. (2003) |
| Soil uptake | ¹³C | N(1.0201, 0.003) | 1.0201 | Tyler et al. (1994); King et al. (1989) |
| Soil uptake | D | N(1.083, 0.01) | 1.083 | Estimated from limited field data |

### 2.2 Sensitivity Tests

| Test | KIE Mode | Purpose |
|------|----------|---------|
| **A** | `fixed` | Baseline with literature central values — isolates other uncertainties |
| **B** | `sampled` | Full uncertainty propagation — realistic error bars |
| **C** (manual) | Fix OH-¹³C at low end (1.0039) | Test Cantrell et al. (1990) value specifically |
| **D** (manual) | Fix OH-¹³C at high end (1.0054) | Test Saueressig et al. (2001) value specifically |

**Expected impact:** KIE sampling typically widens FF/Mic uncertainty bands by 20–40% (documented in v2.0 → v1.0 comparison). The OH-¹³C KIE is the dominant contributor because OH accounts for ~84% of the total sink.

### 2.3 Sink Fractions (Weights)

The bulk KIE is a weighted sum over individual sink KIEs. Weights differ by geometry:

| Sink | Global | NH | SH | References |
|------|--------|----|----|------------|
| OH | 0.835 | 0.825 | 0.850 | Prather et al. (2012); Zhao et al. (2020) |
| Cl | 0.035 | 0.040 | 0.028 | Allan et al. (2007); Hossaini et al. (2016) |
| Stratospheric | 0.070 | 0.070 | 0.070 | Prather et al. (2012) |
| Soil | 0.060 | 0.065 | 0.052 | Dutaur & Verchot (2007) |

**Potential additional test:** Vary Cl fraction within ±0.01 to test marine boundary layer Cl sensitivity. Recent studies (Wang et al. 2021) suggest Cl may be higher than the commonly assumed 3.5%.

---

## 3. CH₄ Atmospheric Lifetime

**CLI flag:** `--lifetime fixed | varying` and `--tau <value>`

| Mode | Formula | Reference |
|------|---------|-----------|
| `fixed` | τ = 9.0 yr (constant for all years) | IPCC AR6, Prather et al. (2012) |
| `varying` | τ(t) = 9.0 − 0.017 × (t − 2010) | He et al. (2026) — trend driven by rising OH and tropospheric ozone |

For two-box models, the global lifetime is split:
- τ_NH = τ_global × 0.95 (more OH in tropics/NH)
- τ_SH = τ_global × 1.05

### 3.1 Sensitivity Tests

| Test | Lifetime | τ Value | Purpose |
|------|----------|---------|---------|
| **A** | `fixed` | 9.0 yr | Standard assumption (IPCC AR6) |
| **B** | `fixed` | 8.5 yr | Lower bound — tests short-lifetime scenario |
| **C** | `fixed` | 9.5 yr | Upper bound — tests long-lifetime scenario |
| **D** | `fixed` | 10.0 yr | Extreme — Rigby et al. (2017) lower OH estimate |
| **E** | `varying` | N/A | Time-varying He et al. (2026) trend |
| **F** (manual) | `varying` | Modified slope | Change slope from −0.017 to −0.010 or −0.025 |

**Expected impact:** A ±0.5 yr change in τ shifts total source by ~60 Tg/yr. This propagates roughly equally into FF and Mic in the 2×2 models.

---

## 4. Source Isotopic Signatures

### 4.1 Fossil Fuel (FF) Signatures

Two alternative inventories are available:

| Inventory | δ¹³C File | δD File | Rows | Description |
|-----------|-----------|---------|------|-------------|
| **EDGAR** | `FF_d13C_GlobMC_EDGAR.csv` | `FF_dD_GlobMC_EDGAR.csv` | 53 (1970–2022) | EDGAR v8 emission categories (coal, oil, gas) weighted by relative contribution. 1000 MC columns. |
| **CT-CH4** | `FF_d13C_GlobMC_CTCH4.csv` | `FF_dD_GlobMC_CTCH4.csv` | 24 (1999–2022) | CarbonTracker-CH₄ emission categories. 1000 MC columns. |

Central values with uncertainty are also available in `FF_d13C_GlobUnc.csv` and `FF_dD_GlobUnc.csv` (mean ± 1σ per year).

**Sensitivity tests:**

| Test | FF Signature Source | Purpose |
|------|---------------------|---------|
| **A** | EDGAR MC (default) | Standard — largest temporal coverage |
| **B** | CT-CH4 MC | Test inventory dependence (coal/gas ratio differs between EDGAR and CT-CH4) |
| **C** | Gaussian perturbation around GlobUnc central values | Simpler uncertainty model — no correlated sub-source trends |

**Expected impact:** The FF δ¹³C signature has trended less negative over 2000–2022 (shift from coal to gas). This trend drives a ~10–15 Tg/yr difference in inferred FF depending on which inventory is used.

**References:** Schwietzke et al. (2016); Sherwood et al. (2017); EDGAR v8 documentation.

### 4.2 Microbial (Mic) Signatures

| File | Isotope | Rows | Description |
|------|---------|------|-------------|
| `Mic_d13C_MC.csv` | δ¹³C | 24 | Weighted average of wetland, ruminant, rice, waste sub-sources. 1000 MC columns. |
| `Mic_d13C_annual.csv` | δ¹³C | 24 | Central value + uncertainty per year. |
| `Mic_dD_MC.csv` | δD | 24 | Same sub-source weighting for δD. 1000 MC columns. |
| `Mic_dD_AnnGlob.csv` | δD | 24 | Central value + uncertainty per year. |

Central values: δ¹³C_Mic ≈ −61.7 ‰, δD_Mic ≈ −306.5 ‰

**Sensitivity tests:**

| Test | Variation | Purpose |
|------|-----------|---------|
| **A** | Full MC (default) | Propagates sub-source fraction + isotopic uncertainty |
| **B** | Fixed to annual mean (no MC) | Isolate effect of source-signature uncertainty on results |
| **C** (manual) | Shift Mic δD by ±10 ‰ | Test sensitivity to poorly constrained biogenic δD |

**References:** Dlugokencky et al. (2011); Schaefer et al. (2016); Ganesan et al. (2018).

### 4.3 Biomass Burning (BB) Signatures

| File | Isotope | Rows | Description |
|------|---------|------|-------------|
| `BB_d13C_annual.csv` | δ¹³C | 25 | GFED fire-type fractions × C3/C4 vegetation δ¹³C. Central + 1σ. |
| `BB_dD_annual.csv` | δD | 25 | Same for δD from pyrolysis measurements. Central + 1σ. |

Central values: δ¹³C_BB ≈ −22 to −25 ‰, δD_BB ≈ −210 to −230 ‰

**References:** van der Werf et al. (2017) [GFED4]; Whiticar (1999); Sherwood et al. (2017).

---

## 5. Biomass Burning Emissions (Mass)

| Model Type | BB Treatment | Source | Reference |
|------------|-------------|--------|-----------|
| 2×2 models | **Fixed** from CarbonTracker | `CarbonTracker_CH4.xlsx` column 10 | Bruhwiler et al. (2014); CT-CH4 |
| 3×3 models | **Free variable** (solved from dual-isotope system) | N/A (determined by model) | — |

Global mean BB ≈ 29 Tg/yr (CarbonTracker).

For two-box models, BB is split:
- NH: 55% (BB × 0.55)
- SH: 45% (BB × 0.45)

Based on GFED4 hemispheric fire distribution (van der Werf et al. 2017).

**Sensitivity tests (2×2 models):**

| Test | Variation | Purpose |
|------|-----------|---------|
| **A** | CarbonTracker annual (default) | Year-specific BB from CT-CH4 |
| **B** | Fixed global mean (29 Tg/yr, all years) | Remove interannual BB variability |
| **C** | BB ± 30% | Test sensitivity to BB magnitude uncertainty |
| **D** | NH/SH split: 60/40 or 50/50 | Test sensitivity to hemispheric fire distribution |

---

## 6. Atmospheric Observations

### 6.1 CH₄ Concentration

**Source:** `GML_CH4_AnnualMean.xlsx` — NOAA GML global annual mean (rows 16–39 → 1999–2022, 24 values).

For two-box models, NH/SH concentrations are derived from global + interhemispheric (IH) gradient:
- CH₄_NH = CH₄_global + IH_gradient / 2
- CH₄_SH = CH₄_global − IH_gradient / 2
- IH_gradient: linear 80 → 100 ppb over 1999–2022

**Sensitivity test:** Vary IH gradient (e.g., constant 90 ppb, or 60→120 ppb) to test NH/SH partitioning.

**References:** Dlugokencky et al. (NOAA GML); Nisbet et al. (2019).

### 6.2 δ¹³C Observations

**Source:** `ch4c13_nh_sh_mean.xlsx` — sub-annual δ¹³C-CH₄ with global, NH, SH columns. Annualised with ≥6 observations/year filter.

MC uncertainty: `d13C_dei_compiled.txt` — 1000 MC iterations (columns) × years (rows).

For two-box models, NH/SH MC is constructed as:
```
offset = d13C_global_MC − d13C_global_observed
d13C_NH_MC = d13C_NH_observed + offset
d13C_SH_MC = d13C_SH_observed + offset
```

**References:** White & Vaughn (2011); Nisbet et al. (2019).

### 6.3 δD Observations

**Source:** `GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx` — Global mean δD-CH₄ with 1000 MC iterations (Ben Riddell-Young's construction using Umezawa calibration, excluding BUDS stations).

An improved pipeline exists (`improved_dD_pipeline.py`) with:
- Area-weighted hemispheric averaging (sin(lat) weighting)
- Cosine-latitude weighting within latitude bands
- Better handling of station coverage gaps

For two-box models, NH/SH δD is estimated as:
- δD_NH = δD_global − 6 ‰
- δD_SH = δD_global + 6 ‰

Based on Riddell-Young (2025) observation that NH δD is ~12 ‰ lower than SH.

**Sensitivity tests:**

| Test | δD Source | Purpose |
|------|-----------|---------|
| **A** | Ben's original (default) | Standard — used in Riddell-Young (2025) |
| **B** | Improved pipeline output | Test area-weighting impact |
| **C** | Vary IH offset: ±4 ‰ or ±8 ‰ instead of ±6 ‰ | Test δD hemispheric gradient sensitivity |

**References:** Riddell-Young (2025); Umezawa et al. (2012); Rice et al. (2016).

---

## 7. Interhemispheric Exchange (Two-Box Only)

| Parameter | Distribution | Default | Reference |
|-----------|------------|---------|-----------|
| τ_exchange | N(1.0, 0.1) yr | 1.0 yr | Patra et al. (2011); Liang et al. (2017) |

The exchange time controls how fast air (and CH₄) mixes between NH and SH. It appears in the mass and isotope balance as:

```
Exchange_NH = (M_SH − M_NH) / τ_ex
```

**Sensitivity tests:**

| Test | τ_ex | Purpose |
|------|------|---------|
| **A** | N(1.0, 0.1) (default) | Standard with uncertainty |
| **B** | Fixed 1.0 yr | Remove exchange uncertainty |
| **C** | Fixed 0.8 yr | Fast mixing — more homogeneous hemispheres |
| **D** | Fixed 1.3 yr | Slow mixing — larger hemispheric differences |

---

## 8. Solver Configuration (3×3 Models Only)

### 8.1 Direct vs Bounded Least Squares

| Model | Solver | Non-Negativity | Weighting |
|-------|--------|----------------|-----------|
| `3x3_one` | `np.linalg.solve` (direct) | None — allows negative solutions | None |
| `3x3_two` | `scipy.optimize.lsq_linear` (bounded) | Enforced: 0 ≤ x ≤ 1.5·S | Diagonal weighting matrix |

### 8.2 Weighting Matrix (3×3 Two-Box)

The 3×3 system rows have very different magnitudes:
- Row 1 (mass): ~580 Tg
- Row 2 (¹³C): ~6.4 (fraction units)
- Row 3 (δD): ~0.09 (fraction units)

Current weights:

| Hemisphere | W₁₁ (mass) | W₂₂ (¹³C) | W₃₃ (δD) |
|------------|------------|------------|-----------|
| NH | 100 | 1 | 0.5 |
| SH | 200 | 1 | 0.5 |

**Sensitivity tests:**

| Test | Weights | Purpose |
|------|---------|---------|
| **A** | Default (100/200, 1, 0.5) | Standard |
| **B** | Equal (1, 1, 1) — unweighted | Test whether weighting matters |
| **C** | Heavy δD (100, 1, 5) | Force δD to have more influence |
| **D** | δ¹³C only (100, 1, 0.01) | Effectively ignore δD (approaches 2×2 behaviour) |

---

## 9. Monte Carlo Configuration

| Parameter | Default | Range | Purpose |
|-----------|---------|-------|---------|
| `--iterations` | 1000 | 100–10000 | Convergence testing |
| `--seed` | 42 | Any integer | Reproducibility / ensemble comparison |

**Sensitivity tests:**

| Test | N | Purpose |
|------|---|---------|
| **A** | 100 | Quick screening runs (~10× faster) |
| **B** | 1000 | Standard (default) |
| **C** | 5000 | Convergence check — do medians/percentiles stabilise? |
| **D** | 10000 | Publication-grade uncertainty quantification |

---

## 10. Recommended Sensitivity Matrix

The full combinatorial space is large (4 models × 2 KIE × 6 lifetime × ...). We recommend a structured approach: hold all inputs at default except one, then vary that one input.

### 10.1 Core Matrix (24 runs)

Run all 4 models × 6 configurations:

| Run | KIE | Lifetime | τ | Label |
|-----|-----|----------|---|-------|
| 1 | sampled | varying | — | **Default** |
| 2 | fixed | varying | — | KIE sensitivity |
| 3 | sampled | fixed | 9.0 | Lifetime sensitivity |
| 4 | sampled | fixed | 8.5 | Short lifetime |
| 5 | sampled | fixed | 9.5 | Long lifetime |
| 6 | fixed | fixed | 9.0 | Minimal uncertainty (lower bound) |

```bash
# Example: run entire core matrix
for model in 2x2_one 2x2_two 3x3_one 3x3_two; do
  python3 ${model}.py --kie sampled --lifetime varying              # Run 1
  python3 ${model}.py --kie fixed   --lifetime varying              # Run 2
  python3 ${model}.py --kie sampled --lifetime fixed --tau 9.0      # Run 3
  python3 ${model}.py --kie sampled --lifetime fixed --tau 8.5      # Run 4
  python3 ${model}.py --kie sampled --lifetime fixed --tau 9.5      # Run 5
  python3 ${model}.py --kie fixed   --lifetime fixed --tau 9.0      # Run 6
done
```

### 10.2 Extended Tests (Require Code Modifications)

These tests require editing `common.py` or the model files:

| Test ID | Input Changed | Variation | Files to Edit |
|---------|---------------|-----------|---------------|
| E1 | FF signatures | Switch EDGAR → CT-CH4 | `common.py` → `load_data()` |
| E2 | δD observations | Use improved pipeline output | `common.py` → `load_data()` |
| E3 | BB mass | ±30% scaling | Model files → BB array |
| E4 | NH/SH BB split | 60/40 or 50/50 | `common.py` → `BB_NH_FRACTION` |
| E5 | IH gradient | Constant 90 ppb | `common.py` → `compute_IH_gradient()` |
| E6 | τ_exchange | Fixed 0.8 or 1.3 yr | `common.py` → `TAU_EX_MEAN` |
| E7 | δD IH offset | ±4 or ±8 ‰ | `common.py` → `DD_IH_OFFSET` |
| E8 | Cl sink fraction | 0.035 → 0.050 | `common.py` → `SINK_FRACTIONS_*` |
| E9 | Solver weights | Various W matrices | `3x3_two.py` → `W_NH`, `W_SH` |
| E10 | Mic δD shift | ±10 ‰ offset | `common.py` → `load_data()` |

---

## 11. Expected Sensitivity Rankings

Based on prior runs and literature analysis, the approximate sensitivity ranking (largest → smallest impact on FF/Mic partition):

1. **CH₄ lifetime (τ)** — ±0.5 yr → ~60 Tg/yr shift in total source
2. **OH-¹³C KIE** — Cantrell vs Saueressig → ~30–40 Tg/yr shift in FF (via δ¹³C)
3. **Model architecture** (2×2 vs 3×3) — ~50–90 Tg/yr difference in FF
4. **FF source signature inventory** (EDGAR vs CT-CH4) — ~10–15 Tg/yr
5. **δD observations** (Ben vs improved) — ~5–10 Tg/yr (mainly affects δD-based FF)
6. **BB magnitude** (2×2 models) — ~10 Tg/yr per 10 Tg BB change (directly subtracted)
7. **Interhemispheric exchange** — ~5–10 Tg/yr (two-box only)
8. **Cl sink fraction** — ~5 Tg/yr (small sink, but high KIE amplifies effect)
9. **Solver weighting** (3×3 two) — ~5–15 Tg/yr (mainly affects BB/FF split)
10. **δD IH offset** — ~3–5 Tg/yr (two-box δD only)

---

## 12. References

- Allan, W., et al. (2007). Methane carbon isotope effects caused by atomic chlorine in the marine boundary layer. *JGR*, 112, D04306.
- Bruhwiler, L., et al. (2014). CarbonTracker-CH4: an assimilation system for estimating emissions of atmospheric methane. *ACP*, 14, 8269–8293.
- Cantrell, C. A., et al. (1990). Carbon kinetic isotope effect in the oxidation of methane by hydroxyl radicals. *JGR*, 95, 22455–22462.
- Dlugokencky, E. J., et al. (2011). Global atmospheric methane. *Phil. Trans. R. Soc. A*, 369, 2058–2072.
- Dutaur, L. & Verchot, L. V. (2007). A global inventory of the soil CH₄ sink. *Global Biogeochem. Cycles*, 21, GB4013.
- Ganesan, A. L., et al. (2018). Advancing scientific understanding of the global methane budget. *Environ. Res. Lett.*, 14, 063004.
- Gierczak, T., et al. (1997). Rate coefficients for the reactions of hydroxyl radicals with methane isotopologues. *J. Phys. Chem. A*, 101, 3125–3134.
- He, J., et al. (2026). Declining methane lifetime implies rising hydroxyl concentrations. *Nature*, (in press).
- Hossaini, R., et al. (2016). A global model of tropospheric chlorine chemistry. *JGR Atmos.*, 121, 13312–13338.
- King, G. M. (1989). Oxidation of methane in soil. *FEMS Microbiol. Rev.*, 63, 431–442.
- Liang, Q., et al. (2017). Deriving global OH abundance and atmospheric lifetimes for long-lived gases. *JGR Atmos.*, 122, 10462–10487.
- Nisbet, E. G., et al. (2019). Very strong atmospheric methane growth in the 4 years 2014–2017. *Global Biogeochem. Cycles*, 33, 318–342.
- Patra, P. K., et al. (2011). TransCom model simulations of CH₄ and related species. *ACP*, 11, 12813–12837.
- Prather, M. J., et al. (2012). Reactive greenhouse gas scenarios. *GRL*, 39, L09803.
- Rice, A. L., et al. (2003). Atmospheric methane isotopic record favors fossil sources. *PNAS*, 100, 10168–10173.
- Rice, A. L., et al. (2016). Isotopic constraints on the atmospheric methane sink. *Science*, 352, 428–431.
- Riddell-Young, B. (2025). Dual-isotope constraints on the global methane budget. PhD thesis, Royal Holloway, University of London.
- Rigby, M., et al. (2017). Role of atmospheric oxidation in recent methane growth. *PNAS*, 114, 5373–5377.
- Röckmann, T., et al. (2003). The isotopic composition of methane in the stratosphere. *ACP*, 3, 2003–2024.
- Saueressig, G., et al. (1995). Carbon kinetic isotope effect in the reaction of CH₄ with Cl atoms. *GRL*, 22, 1225–1228.
- Saueressig, G., et al. (1996). The reaction of CH₃D + Cl. *GRL*, 23, 3581–3584.
- Saueressig, G., et al. (2001). Carbon 13 and D kinetic isotope effects in the reactions of CH₄ with O(¹D) and OH. *Chem. Phys. Lett.*, 345, 321–328.
- Schaefer, H., et al. (2016). A 21st-century shift from fossil-fuel to biogenic methane emissions indicated by ¹³CH₄. *Science*, 352, 80–84.
- Schwietzke, S., et al. (2016). Upward revision of global fossil fuel methane emissions based on isotope database. *Nature*, 538, 88–91.
- Sherwood, O. A., et al. (2017). Global inventory of gas geochemistry data from fossil fuel, microbial, and biomass burning sources. *ESSD*, 9, 639–656.
- Tyler, S. C., et al. (1994). Carbon and hydrogen isotope fractionation during oxidation of methane by soil. *GCA*, 58, 1625–1633.
- Tyler, S. C., et al. (2000). Stable carbon isotope ratios and ¹⁴C as tracers for fossil fuel influence. *GRL*, 27, 1111–1114.
- Umezawa, T., et al. (2012). Interlaboratory comparison of δ¹³C and δD measurements of atmospheric CH₄. *AMT*, 5, 2807–2820.
- van der Werf, G. R., et al. (2017). Global fire emissions estimates during 1997–2016. *ESSD*, 9, 697–720.
- Wang, X., et al. (2021). Global tropospheric halogen (Cl, Br, I) chemistry and its impact on oxidants. *ACP*, 21, 13973–13996.
- White, J. W. C. & Vaughn, B. H. (2011). Stable isotopic composition of atmospheric methane from the NOAA/ESRL Global Monitoring Division. *INSTAAR*, University of Colorado.
- Whiticar, M. J. (1999). Carbon and hydrogen isotope systematics of bacterial formation and oxidation of methane. *Chem. Geol.*, 161, 291–314.
- Zhao, Y., et al. (2020). On the role of trend and variability in the hydroxyl radical for the global methane budget. *ACP*, 20, 13011–13030.
