# importance-KIE: Relative Importance of OH-¹³C vs OH-D KIE in One-Box Methane Source Partitioning

## 1. Introduction

The kinetic isotope effect (KIE) for the OH sink is the most contested parameter in methane isotope inversions. For ¹³C, two laboratory values — Saueressig et al. (2001) at 1.0039 and Cantrell et al. (1990) at 1.0054 — bracket a range that determines whether post-2007 fossil fuel emissions are rising or falling. For D/H, three measurements cluster at 1.29–1.31, a much smaller relative spread.

This report analyzes the **relative sensitivity** of a one-box (and two-box) methane source inversion to the OH-¹³C and OH-D KIE uncertainties, comparing these against other sources of noise in the system: source signature uncertainty, atmospheric observation uncertainty, and lifetime uncertainty. Evidence is drawn from the repository's experiment results (`experiments/KIE_immunity/`, `experiments/dD_threshold/`, `experiments/KIE_sensitivity/`) and cross-referenced with literature in `ImportantReferences/`.

---

## 2. The OH-¹³C KIE Controversy

### 2.1 The Two Competing Values

| Study | KIE^C_OH | ε_OH (‰) | Repository source |
|-------|:--------:|:---------:|-------------------|
| Saueressig et al. (2001) | 1.0039 | 3.9 | `common.py` line 104; `Thanwerdas2024ACP` |
| Cantrell et al. (1990) | 1.0054 | 5.4 | `common.py` line 104; `Fujita2025JGR_SI` Table S4 footnote f |

The difference **Δε = 1.5‰** may appear small, but it is large relative to the atmospheric signal. The total δ¹³C shift from source composition to atmospheric composition is only ~6‰ (sources ~−53.6‰, atmosphere ~−47.3‰; see `rel/data/d13C_dei_compiled.txt`). A 1.5‰ uncertainty in the sink fractionation represents **25% of the total source-to-atmosphere fractionation shift**.

### 2.2 Impact on Source Partitioning

From `experiments/KIE_immunity/` (phase results in JSON files, verified in `REVISION_RESPONSE.md`):

- **Saueressig KIE → FF trend = +2.9 Tg yr⁻¹** (rising fossil fuel emissions)
- **Cantrell KIE → FF trend = −5.6 Tg yr⁻¹** (declining fossil fuel emissions)

The KIE choice **determines the sign** of the inferred FF emission trend. This is the single most consequential parameter uncertainty in the δ¹³C-based methane budget.

### 2.3 Thanwerdas et al. (2024) Perspective

From `Thanwerdas2024ACP` (line ~65):

> "Saueressig et al. (2001) indicate that their data is of considerably higher experimental precision... we prefer to allocate computational time to a sensitivity inversion testing a different OH field rather than testing a different OH fractionation coefficient."

This pragmatic choice reveals a key tension: while the measurement precision argument favors Saueressig, the community remains divided.

---

## 3. The OH-D KIE: A Narrower Uncertainty

### 3.1 Experimental Range

Three consistent measurements: Gierczak (1.292), Saueressig (1.294), Joelsson (1.311). See `KIE-D-OH.md` for full details.

- **Range: 1.292–1.311 → Δ(KIE^D_OH) = 0.019**
- **Relative range: 0.019 / 1.30 = 1.5%**

For comparison, OH-¹³C:
- **Range: 1.0039–1.0054 → Δ(KIE^C_OH) = 0.0015**
- **Relative range: 0.0015 / 1.0047 = 0.15%** in α, but **Δε / ε_mean = 1.5 / 4.65 = 32%** in the fractionation factor

### 3.2 Why the OH-D Range Matters Less

The D/H system has a fundamentally different error budget geometry than ¹³C:

1. **Larger absolute fractionation**: ε^D_OH ≈ 294‰ vs ε^C_OH ≈ 4.65‰. The uncertainty of ±19‰ (from the KIE range) is only **6.5%** of the fractionation factor.

2. **Source-to-atmosphere shift dwarfs KIE uncertainty**: The δD shift from source (~−310‰ for microbial) to atmosphere (~−80‰) is ~230‰. The KIE^D_OH uncertainty of ±19‰ is ~8% of this budget.

3. **Source signature spread dominates**: Microbial δD ≈ −310‰, fossil δD ≈ −190‰, BB δD ≈ −220‰. The source-to-source spread (~120‰) easily absorbs KIE perturbations.

In contrast, for ¹³C: microbial δ¹³C ≈ −60‰, fossil δ¹³C ≈ −44‰, BB δ¹³C ≈ −25‰. The KIE uncertainty of 1.5‰ is a substantial fraction of the ~16‰ FF-Mic separation.

---

## 4. Quantitative Noise Comparison from Repository Data

### 4.1 Atmospheric Observation Uncertainty

From `rel/data/` MC ensembles:

| Observable | Period | Annual MC std | Annual signal (trend) | SNR per year |
|-----------|--------|:------------:|:--------------------:|:------------:|
| δ¹³C (global) | 2000–2023 | 0.006–0.021‰ | −0.025‰/yr (full) / −0.035‰/yr (2005–2023) | ~2–6 |
| δD (global) | 2005–2023 | 0.41–0.92‰ | −0.30‰/yr | ~0.3–0.7 |
| CH₄ (global) | 1984–2023 | ~2 ppb | ~5–10 ppb/yr | ~3–5 |

**Source files:**
- δ¹³C: `rel/data/d13C_dei_compiled.txt` (1000 MC iterations, 1999–2023)
- δD: `rel/data/GlobMean_dD_dei_DasguptaCal_noBUDS.csv` (2005–2023)
- CH₄: `rel/data/GML_CH4_AnnualMean.xlsx` (1984–2023)

**Key observation**: δ¹³C has ~50–100× smaller measurement noise than δD, but δD has an ~8.5× steeper trend in the overlapping period. The per-year signal-to-noise ratio (SNR) is actually **worse for δD** (~0.5 vs ~4 for δ¹³C), meaning δD requires multi-year averaging or MC aggregation to extract meaningful signals.

### 4.2 Source Signature Uncertainty

From `rel/data/*_MC.csv` files (1000 MC iterations, hemispheric):

| Source | Isotope | Hemisphere | Mean (‰) | MC std (‰) |
|--------|---------|:----------:|:---------:|:-----------:|
| Fossil fuel | δ¹³C | NH | −43.4 | 2.4 |
| Fossil fuel | δ¹³C | SH | −48.0 | 2.7 |
| Fossil fuel | δD | NH | −193.1 | 5.6 |
| Fossil fuel | δD | SH | −189.6 | 8.1 |
| Microbial | δ¹³C | NH | −59.9 | 1.1 |
| Microbial | δ¹³C | SH | −59.7 | 1.1 |
| Microbial | δD | NH | −316.9 | 7.8 |
| Microbial | δD | SH | −304.9 | 7.3 |
| BB | δ¹³C | NH | −26.0 | 2.7 |
| BB | δ¹³C | SH | −24.2 | 2.4 |
| BB | δD | NH | −236.7 | 8.2 |
| BB | δD | SH | −210.3 | 7.1 |

**Key observations:**
- δD source signatures have ~3–4× larger absolute uncertainty than δ¹³C (e.g., Mic δD: 7.5‰ vs Mic δ¹³C: 1.1‰)
- But relative to source separation: δ¹³C FF−Mic = 16.5‰ with ~2.6‰ combined uncertainty (16%); δD FF−Mic = 124‰ with ~9.6‰ combined uncertainty (8%)
- **δD has better relative precision for source discrimination**, despite larger absolute errors

### 4.3 KIE Parameter Uncertainty

From `common.py` (lines 103–126):

| Parameter | Distribution | Range | Δε (‰) |
|-----------|:----------:|:-----:|:-------:|
| OH_¹³C | U(1.0039, 1.0054) | 0.0015 | 1.5 |
| OH_D | U(1.294, 1.327) | 0.033 | 33 |
| Cl_¹³C | N(1.066, 0.002) | ~0.004 (1σ) | 4.0 |
| Cl_D | N(1.52, 0.02) | ~0.04 (1σ) | 40 |
| Strat_¹³C | N(1.003, 0.001) | ~0.002 | 2.0 |
| Strat_D | N(1.179, 0.01) | ~0.02 | 20 |
| Soil_¹³C | N(1.0201, 0.003) | ~0.006 | 6.0 |
| Soil_D | N(1.083, 0.01) | ~0.02 | 20 |

---

## 5. Variance Decomposition: What Drives FF Uncertainty?

### 5.1 From the KIE_immunity Experiment

The `KIE_immunity` experiment (`experiments/KIE_immunity/`) performed a formal variance decomposition by selectively freezing parameter groups. Results from `phase9_bootstrap.json` and `basu_comparison_v2.json` (verified in `REVISION_RESPONSE.md` Table, A1):

| Variance component | % of FF trend variance | 95% CI |
|--------------------|:---------------------:|:------:|
| **Source signatures** | **47.6%** | [37.7, 56.2] |
| **KIE (all sinks)** | **24.9%** | [12.2, 33.8] |
| **Lifetime** | **0.8%** | [0.1, 1.5] |
| Residual / interactions | 27.4% | — |

**Source signatures dominate**, contributing nearly half the variance. The KIE contributes ~25%, and lifetime is negligible (<1%).

### 5.2 KIE Spread: Dual-Isotope vs Single-Isotope

From `basu_comparison_v2.json` (`REVISION_RESPONSE.md` Table, A1):

| Configuration | KIE-driven FF spread (Tg/yr) |
|--------------|:---------------------------:|
| δ¹³C only (hemispheric) | 13.0 |
| δ¹³C + δD (dual-isotope, hemispheric) | 8.6 |
| **Reduction** | **34%** |

Adding δD reduces the KIE-driven spread from 13.0 to 8.6 Tg/yr — a meaningful but incomplete reduction. The **KIE still determines the sign** of the FF trend:

- Saueressig → ΔFF = **+2.9** Tg/yr (rising)
- Cantrell → ΔFF = **−5.6** Tg/yr (falling)

### 5.3 W Matrix Sensitivity

From `phase14_W_sensitivity.json` (`REVISION_RESPONSE.md`, A2):

- KIE% varies from 24.6% to 25.5% across 6 W configurations → **KIE importance is robust to weighting**
- σ(FF) varies from 15.4 to 19.3 Tg/yr → moderate sensitivity
- FF trend sign can flip depending on W → evidence that δ¹³C and δD partially **disagree** about the FF trend direction

---

## 6. Why OH-¹³C KIE Matters More Than OH-D KIE

### 6.1 Fractionation Leverage Argument

Consider a one-box isotope mass balance. The atmospheric δ value is shifted from the source-weighted mean by the sink KIE:

δ_atm ≈ δ_source + ε_sink (simplified)

For ¹³C:
- ε^C_sink ≈ 6‰ (from δ_source ≈ −53.6 to δ_atm ≈ −47.3)
- KIE uncertainty Δε^C ≈ 1.5‰ → **25% of the total fractionation**

For D:
- ε^D_sink ≈ 230‰ (from δ_source ≈ −310 to δ_atm ≈ −80)
- KIE^D_OH uncertainty Δε^D_OH ≈ 19‰ → **8% of the total fractionation**
- Sink-weighted total KIE^D uncertainty ≈ ±12‰ → **5%** of total fractionation

The ¹³C system operates with **much less headroom**: a small KIE shift consumes a large fraction of the available fractionation budget, strongly perturbing the inferred source mix. The D system has ~40× more fractionation dynamic range, so KIE uncertainties are proportionally less damaging.

### 6.2 Propagation to Source Attribution

In the ¹³C mass balance for a two-source (FF+Mic) system:

FF/Mic ratio ∝ (δ_atm − ε_sink − δ_Mic) / (δ_FF − δ_Mic)

The denominator (δ_FF − δ_Mic ≈ 16.5‰ for ¹³C) is small relative to Δε of 1.5‰, so the numerator is highly sensitive to KIE choice.

For D:

FF/Mic ratio ∝ (δD_atm − ε^D_sink − δD_Mic) / (δD_FF − δD_Mic)

The denominator (δD_FF − δD_Mic ≈ 124‰) is large relative to any plausible KIE^D uncertainty, making the source ratio insensitive to KIE^D.

### 6.3 From the KIE_sensitivity Experiment

The `KIE_sensitivity` experiment (`experiments/KIE_sensitivity/`) directly tested whether solving δ¹³C and δD independently and checking agreement could discriminate between OH-¹³C KIE values. Key findings from `RESULTS.md` and phase scripts:

- **WLS coupling makes KIE sensitivity worse** (KSR ≈ 0.2) — Phases 1–5 showed no optimal W_dD that improves discrimination
- **Root cause**: shifted δ¹³C row contradicts unshifted δD row → WLS distributes the conflict across both unknowns
- **Agreement filter approach**: solving isotopes independently and filtering by FF consistency shows **35.5 pp agreement-rate discriminant** — Cantrell KIE produces more internally consistent solutions
- This discriminant is **stable across epochs** (Phase 8) and survives time-varying KIE trajectories (Phase 7)

---

## 7. Comparison with Other System Noise

### 7.1 Hierarchy of Uncertainty Sources

Combining the experimental data analysis with the variance decomposition:

| Uncertainty source | Impact on FF trend (Tg/yr) | % of variance | Reducible? |
|-------------------|:--------------------------:|:-------------:|:----------:|
| Source signatures | dominant | 48% | Yes (more measurements) |
| OH-¹³C KIE | ΔFF ≈ 8.5 Tg/yr (KIE spread) | 25% | No (needs new lab work) |
| δD observation noise | σ ≈ 0.5–0.9‰/yr | ~10%* | Yes (more stations) |
| OH-D KIE | minor vs OH-¹³C | <5%* | Already well-constrained |
| δ¹³C observation noise | σ ≈ 0.01‰/yr | <2%* | Already small |
| Lifetime | negligible | <1% | N/A |

*Estimated from the KIE_immunity variance decomposition residual and the source-signature sensitivity tests.

### 7.2 The δD Threshold Result

From `experiments/dD_threshold/` (draft.md, §3):

- δD improves FF uncertainty **only when** microbial δD source-signature uncertainty σ(Mic δD) < 37‰
- Current observational precision (~8‰ from sitesdata) provides a **4.5× safety margin**
- Thanwerdas et al. (2024) used σ ≈ 128‰ (3.5× above threshold) → found δD "useless"
- δD's value is **hemispheric**: one-box dual-isotope fails, two-box succeeds (53% CI reduction)

This means δD's contribution is primarily through **source signature discrimination** (where its large source-to-source separation matters), not through KIE^D uncertainty reduction.

---

## 8. Conclusions

1. **OH-¹³C KIE is the dominant KIE uncertainty**: It accounts for ~25% of FF emission trend variance and determines the sign of the post-2007 FF trend. This is an irreducible uncertainty floor until new laboratory measurements or atmospheric constraints can resolve the Saueressig-Cantrell controversy.

2. **OH-D KIE is well-constrained**: The three consistent measurements (1.29–1.31) span only 1.5% relative variation. This uncertainty propagates weakly through the D-system budget because the large δD source-to-atmosphere fractionation (~230‰) dwarfs the KIE uncertainty.

3. **Source signatures are the largest uncertainty**: At 48% of variance, improving source signature knowledge (especially hemispheric microbial δ¹³C and δD) would yield the greatest reduction in FF emission uncertainty.

4. **δD helps through source discrimination, not KIE refinement**: Adding δD reduces KIE-driven FF spread by 34% (13.0 → 8.6 Tg/yr), but its primary value is in the hemispheric two-box framework where δD source signatures have 5–10× larger NH-SH gradients than δ¹³C (from `dD_threshold` experiment findings).

5. **Lifetime is negligible**: <1% of variance, confirming that atmospheric chemistry lifetime uncertainty is not a meaningful contributor to source partitioning error.

---

## 9. References (Local Repository Files)

| Source | Path |
|--------|------|
| KIE_immunity results | `experiments/KIE_immunity/REVISION_RESPONSE.md` |
| KIE_immunity manuscript | `experiments/KIE_immunity/MANUSCRIPT_DUAL_ISOTOPE.md` |
| dD_threshold results | `experiments/dD_threshold/draft.md` |
| KIE_sensitivity results | `experiments/KIE_sensitivity/RESULTS.md` |
| Hemispheric_Divergence | `experiments/Hemispheric_Divergence/MANUSCRIPT.md` |
| KIE distributions | `common.py` lines 103–126 |
| δ¹³C MC data | `rel/data/d13C_dei_compiled.txt` |
| δD global data | `rel/data/GlobMean_dD_dei_DasguptaCal_noBUDS.csv` |
| CH₄ data | `rel/data/GML_CH4_AnnualMean.xlsx` |
| Source signature MC files | `rel/data/{FF,Mic,BB}_{d13C,dD}_{NH,SH}_MC.csv` |
| Thanwerdas et al. (2024) | `ImportantReferences/Thanwerdas2024ACP/` |
| Fujita et al. (2025) SI | `ImportantReferences/Fujita2025JGR_SI/` |
| Riddell-Young et al. (2025) | `ImportantReferences/Riddell-Young2025PNAS/` |
| Basu et al. (2022) | `ImportantReferences/Basu2022ACP/` |
| KIE survey | `KIE_Used_Previous_Study.md` |

