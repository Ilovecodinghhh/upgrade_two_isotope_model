# Master Inventory of Methane Isotope Data

## Compiled from 16 Papers + Model Repository
**Date:** 2026-05-11 | **Repository:** Ilovecodinghhh/upgrade_two_isotope_model

---

## 1. ISOTOPIC SOURCE SIGNATURES

### 1.1 δ¹³C Source Signatures (‰, vs. VPDB)

| Source Category | Subcategory | δ¹³C (‰) | Uncertainty | Reference | Geographic Provenance | Temporal Coverage | Method |
|---|---|---|---|---|---|---|---|
| **Fossil Fuel (Thermogenic)** | Global mean (EDGAR-weighted) | MC draws from CSV | ±3–5‰ (1σ from MC) | Sherwood et al. (2017, ESSD); EDGAR v8.0 weighting | Global, country-weighted | Time-invariant | Compilation of field measurements + inventory weighting |
| Fossil Fuel | Global mean (CT-CH₄-weighted) | MC draws from CSV | ±3–5‰ (1σ from MC) | Bruhwiler et al. (2014); CarbonTracker-CH₄ posterior | Global, posterior-weighted | Time-invariant | Inversion posterior + field measurements |
| Fossil Fuel | Overall range | −37 to −45 | — | Schaefer et al. (2016); Schwietzke et al. (2016) | Global average | — | Literature compilation |
| Fossil Fuel | Coal subsector | −41.2 | country-specific available | Sherwood et al. (2017) | Global mean; country-level in He 2026 JGR | Time-invariant | Field sampling |
| Fossil Fuel | Oil and gas subsector | −44.6 | country-specific available | Sherwood et al. (2017) | Global mean; country-level in He 2026 JGR | Time-invariant | Field sampling |
| Fossil Fuel | Extreme low end (biogenic-origin gas) | < −60 | — | Lu et al. (2021); Menoud et al. (2022) | Specific basins | — | Field measurements |
| Fossil Fuel | Schwietzke revised global | −44.0 | ±0.7 (1σ) | Schwietzke et al. (2016, Nature) | Global, inventory-weighted | Pre-industrial to present | 12,000+ measurements compiled |
| **Microbial (Biogenic)** | Overall range | −70 to −50 | — | Chang et al. (2019); Sherwood et al. (2017); Whiticar & Schaefer (2007) | Global | — | Literature range |
| Microbial | Wetlands (spatially resolved) | Mapped grid-cell values | Spatially variable | Ganesan et al. (2018) | Global gridded | Time-invariant | Process model + observations |
| Microbial | Wetlands (global mean) | −62 to −65 | ±5‰ | Oh et al. (2022); Ganesan et al. (2018) | Global | — | Spatially-varying maps |
| Microbial | Rice cultivation | −63 | country-specific | Sherwood et al. (2017); Whiticar & Schaefer (2007) | Global mean; country-level available | Time-invariant | Field/lab measurements |
| Microbial | Ruminants (C3-fed) | −67.9 | ±3‰ | Sherwood et al. (2017) | Global | Time-invariant | Field measurements |
| Microbial | Ruminants (C4-fed) | −54.5 | ±3‰ | Sherwood et al. (2017) | Global | Time-invariant | Field measurements |
| Microbial | Agriculture sector (global emission-weighted) | varies by grid cell | — | He et al. (2026, JGR): rice + ruminant + C3/C4 weighting | Global, spatially-resolved | Time-varying (inventory-driven) | Inventory-weighted compilation |
| **Pyrogenic (Biomass Burning)** | Overall range | −17 to −26 | — | Dlugokencky et al. (2011) | Global | — | Literature range |
| Pyrogenic | C3 vegetation | −25 | ±2‰ | Lassey et al. (2007) | Global | — | Field measurements |
| Pyrogenic | C4 vegetation | −12 | ±2‰ | Lassey et al. (2007) | Global | — | Field measurements |
| Pyrogenic | Spatially-resolved (C3/C4 mix) | Grid-cell values | — | He et al. (2026, JGR): C4 fraction from Still et al. (2009) | Global, spatially-resolved | Time-invariant | C3/C4 vegetation mapping |

### 1.2 δD Source Signatures (‰, vs. VSMOW)

| Source Category | Subcategory | δD (‰) | Uncertainty | Reference | Geographic Provenance | Temporal Coverage | Method |
|---|---|---|---|---|---|---|---|
| **Fossil Fuel** | Global (EDGAR-weighted) | MC draws from CSV | ±15–20‰ (1σ from MC) | Sherwood et al. (2017); EDGAR v8.0 weighting | Global, country-weighted | Time-invariant | Compilation + inventory weighting |
| Fossil Fuel | Global (CT-CH₄-weighted) | MC draws from CSV | ±15–20‰ (1σ from MC) | Bruhwiler et al. (2014); CarbonTracker-CH₄ | Global, posterior-weighted | Time-invariant | Inversion + field measurements |
| Fossil Fuel | Typical range | −150 to −200 | — | Whiticar (1999); Sherwood et al. (2017) | Global | — | Literature compilation |
| Fossil Fuel | Thanwerdas prior | −183 | ±37 (20% of value) | Thanwerdas et al. (2024, ACP) Table 2 | Global average | — | Literature prior |
| **Microbial** | Global MC draws | MC from CSV | 1σ = 7–8.2‰ | Bao et al. (this work); Riddell-Young (2025) Table 1 | Global | Time-invariant | MC sampling |
| Microbial | Typical range | −250 to −400 | — | Whiticar (1999) | Global | — | Literature |
| Microbial | Wetlands (Thanwerdas prior) | −320 | ±128 (40%) | Thanwerdas et al. (2024, ACP) Table 2 | Global | — | Literature prior |
| Microbial | Agriculture (Thanwerdas prior) | −310 | ±93 (30%) | Thanwerdas et al. (2024, ACP) Table 2 | Global | — | Literature prior |
| **Pyrogenic** | Annual MC draws | MC from CSV (BB_dD_annual.csv) | variable | Bao et al. (this work) | Global | Annual values | Compiled |
| Pyrogenic | Thanwerdas prior | −200 | ±70 (35%) | Thanwerdas et al. (2024, ACP) Table 2 | Global | — | Literature prior |

---

## 2. SINK FRACTIONATION FACTORS (KIE = α = k_light/k_heavy)

### 2.1 δ¹³C Kinetic Isotope Effects

| Sink | KIE (α) | Uncertainty | ε (‰) = (1/α − 1)×1000 | Reference | Method | Notes |
|---|---|---|---|---|---|---|
| **OH + CH₄** | 1.0039 | fixed | −3.9 | Saueressig et al. (2001, JGR) | Laboratory (298K) | Lower bound |
| **OH + CH₄** | 1.0054 | fixed | −5.4 | Cantrell et al. (1990, JGR) | Laboratory (298K) | Upper bound |
| **OH + CH₄** | U(1.0039, 1.0054) | full range | −3.9 to −5.4 | Literature range | — | **Default: sampled uniformly** |
| **OH + CH₄** (He 2026 JGR "medium") | ~1.0039 | — | ~−3.9 | He et al. (2026, JGR) EXP1 | 3D model | Used in optimization |
| **OH + CH₄** (He 2026 JGR "higher") | ~1.0054 | — | ~−5.4 | He et al. (2026, JGR) EXP3 | 3D model | Sensitivity test |
| **Cl + CH₄** | 1.066 | ±0.002 (1σ) | −61.9 | Saueressig et al. (1995, GRL) | Laboratory | Normal distribution |
| **Stratosphere** | 1.003 | ±0.001 (1σ) | −3.0 | Saueressig et al. (2001); Lassey et al. (2007, ACP) | Lab + atmospheric obs | Normal distribution |
| **Soil uptake** | 1.0201 | ±0.003 (1σ) | −19.7 | Snover & Quay (2000); Tyler et al. (1994); Reeburgh et al. (1997) | Field measurements | Average of 3 studies |

### 2.2 δD Kinetic Isotope Effects

| Sink | KIE (α) | Uncertainty | ε (‰) = (1/α − 1)×1000 | Reference | Method | Notes |
|---|---|---|---|---|---|---|
| **OH + CH₃D** | 1.294 | fixed | −227 | Saueressig et al. (2001, JGR) | Laboratory | Lower bound |
| **OH + CH₃D** | 1.327 | fixed | −246 | Whitehill-Joelson average | Laboratory | Upper bound |
| **OH + CH₃D** | U(1.294, 1.327) | full range | −227 to −246 | Literature range | — | **Default: sampled uniformly** |
| **Cl + CH₃D** | 1.52 | ±0.02 (1σ) | −342 | Saueressig et al. (2001, JGR) | Laboratory | Normal distribution |
| **Stratosphere (D)** | 1.179 | ±0.01 (1σ) | −152 | Dyonisius et al. (2020, Nature); Beck et al. (2018) | Ice core + atmospheric | Normal distribution |
| **Soil uptake (D)** | 1.083 | ±0.01 (1σ) | −77 | Snover & Quay (2000, GBC) | Field measurements | Normal distribution |

---

## 3. SINK FRACTIONAL CONTRIBUTIONS

| Configuration | OH | Cl | Stratosphere | Soil | Total | Reference | Notes |
|---|---|---|---|---|---|---|---|
| **Global default** | 0.835 | 0.035 | 0.070 | 0.060 | 1.000 | CarbonTracker documentation; composite | Used in 1-box models |
| **Thanwerdas** | 0.899 | 0.006 | 0.030 | 0.065 | 1.000 | Thanwerdas et al. (2024, ACP) | Very low Cl; sensitivity test |
| **NH default** | 0.825 | 0.040 | 0.070 | 0.065 | 1.000 | Estimated (higher Cl in NH tropics) | 2-box NH |
| **SH default** | 0.850 | 0.028 | 0.070 | 0.052 | 1.000 | Estimated (lower Cl, less land) | 2-box SH |

---

## 4. METHANE FLUX ESTIMATES (Tg CH₄ yr⁻¹)

### 4.1 Total Global Emissions

| Period | Total Emissions | Reference | Method | Notes |
|---|---|---|---|---|
| 2019 | 571 | He et al. (2026, Science) | TROPOMI inversion (GEOS-Chem) | IMI posterior |
| 2021 (peak) | 601 | He et al. (2026, Science) | TROPOMI inversion | Peak year |
| 2022–2024 | 571–575 | He et al. (2026, Science) | TROPOMI inversion | Returned to 2019 level |
| 2024 | 575 (568–578) | He et al. (2026, Science) | TROPOMI inversion ensemble | 27-member range |
| 1980–1989 | 521 ± 14 | He et al. (2026, JGR) EXP1 | GFDL-AM4.1 + obs constraint | |
| 1990–1998 | 555 ± 15 | He et al. (2026, JGR) EXP1 | GFDL-AM4.1 + obs constraint | |
| 1999–2006 | 568 ± 13 | He et al. (2026, JGR) EXP1 | GFDL-AM4.1 + obs constraint | |
| 2007–2017 | 607 ± 11 | He et al. (2026, JGR) EXP1 | GFDL-AM4.1 + obs constraint | |
| 2000–2017 | ~550–590 | Basu et al. (2022, ACP) | 3D inversion + δ¹³C | |

### 4.2 Sectoral Emissions (Multi-year means, Tg yr⁻¹)

| Source Sector | 1999–2006 | 2007–2017 | Change | Reference | Notes |
|---|---|---|---|---|---|
| Agriculture + Waste | 206 ± 5 | 225 ± 6 | +19 | He et al. (2026, JGR) EXP1 | CEDS-based |
| Energy (CEDS) | 120 ± 10 | 150 ± 6 | +30 | He et al. (2026, JGR) EXP1 | |
| Energy (EDGAR v5.0) | 88 ± 6 | 111 ± 8 | +23 | He et al. (2026, JGR) EXP2 | Lower than CEDS |
| Wetlands | 180 ± 10 | 177 ± 7 | −3 | He et al. (2026, JGR) EXP1 | Minor decrease |
| Biomass Burning | 12 ± 5 | 5 ± 5 | −7 | He et al. (2026, JGR) EXP1 | Strong decline |
| Wetlands (EXP3, higher KIE) | 171 ± 10 | 168 ± 7 | −3 | He et al. (2026, JGR) EXP3 | |
| Biomass Burning (EXP3) | 21 ± 5 | 15 ± 6 | −6 | He et al. (2026, JGR) EXP3 | Less decline than EXP1 |

### 4.3 He 2026 Science — Sectoral Breakdown (2024)

| Source Sector | Fraction of Total (%) | Emissions ~Tg/yr (from 575 total) | Reference |
|---|---|---|---|
| Wetlands | 28% | ~161 | He et al. (2026, Science) |
| Livestock | 25% | ~144 | He et al. (2026, Science) |
| Waste | 16% | ~92 | He et al. (2026, Science) |
| Oil and gas | 9% | ~52 | He et al. (2026, Science) |
| Other (coal, rice, etc.) | 22% | ~126 | He et al. (2026, Science) |

### 4.4 Regional Emission Trends (2019–2024)

| Region | Trend (Tg yr⁻²) | Range (ensemble) | Dominant Sector | Reference |
|---|---|---|---|---|
| East Africa | +1.5 | 1.4–1.7 | Wetlands + Livestock | He et al. (2026, Science) |
| South America | +1.3 | −0.4 to 2.6 | Oil/gas + waste | He et al. (2026, Science) |
| Europe | +1.1 | 0.9–1.5 | Livestock + waste | He et al. (2026, Science) |
| China | +0.4 | 0.1–0.6 | Coal (northward shift) | He et al. (2026, Science) |
| Africa total | +0.9 | 0.8–1.3 | Mixed | He et al. (2026, Science) |
| US + Canada | ~0 | — | Flat | He et al. (2026, Science) |
| Hudson Bay Lowlands | −0.2 | −0.1 to −0.2 | Wetlands (drought) | He et al. (2026, Science) |
| W. Siberian Lowlands | −0.3 | −0.3 to −0.5 | Wetlands (drought) | He et al. (2026, Science) |
| Amazon | −2.7 | −2.1 to −3.5 | Wetlands (2023 drought) | He et al. (2026, Science) |

### 4.5 Oil/Gas Emissions Trend

| Metric | Value | Reference |
|---|---|---|
| Global oil/gas 2019–2024 change | −9% | He et al. (2026, Science) |
| Rice 2019–2024 change | −17% | He et al. (2026, Science) |
| Livestock 2019–2024 change | +15% | He et al. (2026, Science) |
| Waste 2019–2024 change | +11% | He et al. (2026, Science) |

---

## 5. METHANE LIFETIME AND OH

| Parameter | Value | Uncertainty | Reference | Notes |
|---|---|---|---|---|
| **Total lifetime (IPCC)** | 9.0 yr | ±0.5 | IPCC AR5/AR6; Prather et al. (2012) | Fixed default |
| **Time-varying (He 2026 Science)** | τ(t) = 9.0 − 0.017×(t−2010) | — | He et al. (2026, Science): TROPOMI-derived OH trend | τ(1999)=9.19, τ(2022)=8.80 |
| **OH lifetime (He 2026 Science)** | 11.1 yr (mean 2019–2024) | — | He et al. (2026, Science) | Against OH only; 6% longer than prior |
| **OH interannual variability** | ±2% | — | He et al. (2026, Science); methylchloroform constraint | |
| **NH/SH OH ratio** | 1.01 (0.96–1.03) | ensemble | He et al. (2026, Science) | Prior was 1.07 |
| **Interhemispheric exchange τ_ex** | 1.0 yr | ±0.1 (1σ) | Patra et al. (2011, ACP); Naus et al. (2019) | Normal distribution |
| **Lifetime ratio NH/SH** | NH: 0.95, SH: 1.05 | — | Lawrence et al. (2001) | Shorter τ in NH (more OH) |

---

## 6. ATMOSPHERIC OBSERVATIONS (Model Inputs)

| Observable | File | Source | Temporal Coverage | Spatial |
|---|---|---|---|---|
| Global CH₄ (ppb) | GML_CH4_AnnualMean.xlsx | Lan et al. (2024), NOAA GML | Annual means | Global |
| NH/SH δ¹³C-CH₄ | ch4c13_nh_sh_mean.xlsx | White et al. (2023), INSTAAR/NOAA | Monthly | Hemispheric |
| δ¹³C MC iterations | d13C_dei_compiled.txt | Bao et al. (this work) | Annual | Global (1000 MC) |
| δD global (Umezawa cal.) | GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx | Umezawa et al. (2012); Rice et al. (2016) | Annual | Global (MC) |
| CarbonTracker CH₄ | CarbonTracker_CH4.xlsx | NOAA CarbonTracker-CH₄ | — | — |

---

## 7. HEMISPHERIC PARAMETERS (2-Box Model)

| Parameter | Value | Uncertainty | Reference | Notes |
|---|---|---|---|---|
| NH–SH CH₄ gradient | 80–100 ppb | — | NOAA GML; Dlugokencky et al. (2011) | Start: 80, End: 100 |
| δD interhemispheric offset | 6.0 ‰ | — | Riddell-Young et al. (2025, PNAS) | NH δD ~12‰ lower than SH |
| δD offset (old estimate) | 1.5 ‰ | — | Original v3.0 estimate | Superseded |
| BB NH/SH split | NH: 55%, SH: 45% | — | van der Werf et al. (2017), GFED4 | Fire distribution |

---

## 8. SOURCE-WEIGHTED SIGNATURES AND SINK FRACTIONATION (He 2026 JGR Table 2)

| Period | Source-weighted δ¹³C (‰) | Sink-weighted ε (‰) | Reference |
|---|---|---|---|
| 1980–1989 | −54.72 ± 0.27 | −6.21 ± 0.09 | He et al. (2026, JGR) EXP1 |
| 1990–1998 | −55.32 ± 0.37 | −6.43 ± 0.02 | He et al. (2026, JGR) EXP1 |
| 1999–2006 | −55.25 ± 0.20 | −6.32 ± 0.05 | He et al. (2026, JGR) EXP1 |
| 2007–2017 | −54.96 ± 0.24 | −6.20 ± 0.02 | He et al. (2026, JGR) EXP1 |
| 1999–2006 (high KIE, EXP3) | −54.64 ± 0.20 | −5.70 ± 0.05 | He et al. (2026, JGR) EXP3 |
| 2007–2017 (high KIE, EXP3) | −54.37 ± 0.24 | −5.57 ± 0.02 | He et al. (2026, JGR) EXP3 |

**Key insight:** Switching from "medium" to "higher" OH–¹³C fractionation shifts the sink-weighted ε by ~0.6‰ — demonstrating the KIE sensitivity your Title 2 paper would quantify.

---

## 9. BIOMASS BURNING EMISSIONS (Fixed input for 2×2 models)

| Mode | Description | Reference |
|---|---|---|
| CT_GFED4_mean | Fixed global mean from CarbonTracker/GFED4 | van der Werf et al. (2017) |
| CT_GFED4_annual | Time-varying annual from GFED4 | van der Werf et al. (2017) |
| Declining | 9% end fraction (strong decline) | Worden et al. (2017) sensitivity |

---

# DATA DENSITY SUMMARY

## Regions/Sources with MOST Evidence

| Category | Density | Key Sources | Confidence |
|---|---|---|---|
| **FF δ¹³C global** | ★★★★★ | 12,000+ measurements (Sherwood 2017); country-level; 3 inventory weightings (EDGAR, CT-CH₄, GlobUnc) | Very High |
| **OH–¹³C KIE** | ★★★★☆ | 2 independent lab measurements (Cantrell 1990, Saueressig 2001); heavily cited | High (but unresolved discrepancy) |
| **Global total emissions** | ★★★★★ | Multiple inversions (He, Basu, Saunois); TROPOMI + GOSAT + surface | Very High |
| **East Africa/Tropical trends** | ★★★★☆ | TROPOMI, GOSAT, multiple studies (Zhang 2021, He 2026, Chandra 2024) | High |
| **Cl–¹³C KIE** | ★★★★☆ | Single lab study (Saueressig 1995) but well-constrained (±0.002) | High |
| **Wetland δ¹³C (spatial)** | ★★★☆☆ | Ganesan (2018) map; Oh (2022) alternative; limited validation | Moderate |

## Regions/Sources with LEAST Evidence

| Category | Density | Gap Description | Impact on Model |
|---|---|---|---|
| **δD source signatures (all categories)** | ★★☆☆☆ | Far fewer measurements than δ¹³C; large uncertainties (±30–130‰ in Thanwerdas) | **CRITICAL** — directly limits Title 1 (δD threshold) |
| **Hemispheric δD observations** | ★☆☆☆☆ | No published NH/SH δD-CH₄ time series; relies on sparse station data | **CRITICAL** — limits 2-box δD validation |
| **Soil sink KIE (both isotopes)** | ★★☆☆☆ | Single study (Snover & Quay 2000); ecosystem variability poorly characterized | Moderate — soil is only 6% of total sink |
| **Stratospheric sink δD KIE** | ★★☆☆☆ | Ice core-based (Dyonisius 2020); Beck (2018); limited direct measurements | Moderate — strat is 7% of sink |
| **Tropical wetland δD** | ★☆☆☆☆ | Almost no direct δD measurements from tropical wetlands | High — tropics dominate microbial budget |
| **Time-varying source signatures** | ★★☆☆☆ | Most studies assume time-invariant signatures; Chang (2019) showed ruminant δ¹³C evolves with C3/C4 diet | Moderate — affects trend interpretation |
| **Southern Hemisphere observations (δ¹³C)** | ★★☆☆☆ | Far fewer stations; Cape Grim, SPO main sources | Moderate — limits 2-box SH constraint |
| **Cl sink fraction** | ★★☆☆☆ | Ranges from 0.6% (Thanwerdas) to 3.5% (default) — factor of 6 disagreement | High — affects δD budget via large Cl–D KIE |
| **Post-2020 emissions** | ★★★☆☆ | He 2026 Science covers 2019–2024; limited isotopic constraint after 2022 | Growing gap |

## Critical Uncertainties for Your Box Model

1. **OH–¹³C KIE (1.0039 vs. 1.0054):** Single largest structural uncertainty. No new lab measurement since 2001. Your Title 2 directly addresses this.

2. **Cl sink fraction (0.6% vs. 3.5%):** Factor of ~6 disagreement between Thanwerdas and your default. Cl has the largest δD KIE (α=1.52), so this massively affects the δD budget.

3. **Microbial δD signatures:** Uncertainty of ±93–128‰ (Thanwerdas 2024) vs. your MC draws with σ≈7–8‰. This discrepancy is central to Title 1 — which uncertainty specification is correct?

4. **NH/SH δD gradient:** Your model assumes 6‰ offset (Riddell-Young 2025) but this has never been directly measured from hemispheric station data. A key testable prediction from Title 3.

5. **Time-varying lifetime:** He 2026 suggests τ declining at 0.017 yr/yr since 2010, implying increasing OH. Most other studies assume constant. Major impact on trend attribution.
