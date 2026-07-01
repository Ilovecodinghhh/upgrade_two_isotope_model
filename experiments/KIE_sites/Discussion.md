# KIE_sites — Discussion: Improving the Seasonal Amplitude Ratio via Wetland Source Correction

## Goal

Extract a cleaner observational constraint on the OH ¹³C KIE from the seasonal amplitude ratio:

```
R = Δδ¹³C_seasonal / ΔδD_seasonal ≈ (α_13C − 1) / (α_D − 1)
```

OH is the dominant atmospheric CH₄ sink and peaks in summer. Both δ¹³C and δD are enriched by OH destruction, so when OH dominates, R directly reflects the KIE ratio. The current analysis (phases 1–5) shows R is heavily contaminated by wetland source seasonality at NH sites, inflating R far above the pure-OH prediction (0.013–0.017). The plan below describes how to subtract the wetland contribution and recover a corrected ratio.

---

## 1. The Source Contamination Problem

The observed seasonal cycle of each isotope is a sum of sink and source contributions:

```
Δδ_observed = Δδ_sink (OH enrichment) + Δδ_source (wetland dilution)
```

The contamination is **asymmetric**:
- δ¹³C: source–atm gap = −62 − (−47.3) = −14.7 ‰ (small gap, but ε_13C ≈ 5.8 ‰ is also small)
- δD:   source–atm gap = −310 − (−86) = −224 ‰ (huge gap, but ε_D ≈ 280 ‰ is also huge)

Relative contamination: source/sink ratio is **2.5 for δ¹³C** vs **0.8 for δD**. So source seasonality contaminates the δ¹³C amplitude 3.2× more than δD (relative to the OH signal). This inflates R above the true KIE ratio at any site with significant wetland seasonality.

---

## 2. Site-Specific δD Source Signatures — Completed

We built a database of wetland δD-CH₄ for each of the 12 KIE sites, using two independent methods:

| Method | Basis | Strengths | Weaknesses |
|--------|-------|-----------|------------|
| **OIPC + Douglas 2021 regression** | Predict δD-CH₄ from site precipitation δ²H | Site-specific; captures local climate | Extrapolates at extreme latitudes; assumes local wetland |
| **Douglas 2021 zonal mean** | Emission-weighted mean from 129-site compilation | Empirical; physically grounded | Coarse (3 latitude bands) |

**Key insight**: These stations measure background air, not local emissions. The relevant δD source is the emission-weighted average of upstream wetlands, not what a hypothetical local wetland would emit.

**Recommendations per site**:

| Site | Lat | OIPC (‰) | Zonal (‰) | Recommended | Reasoning |
|------|-----|-----------|-----------|-------------|-----------|
| ALT | +82.5 | −431 | −374 ± 10 | **Zonal** | No local wetlands; regression extrapolates |
| ZEP | +78.9 | −344 | −374 ± 10 | **Zonal** | No local wetlands; Gulf Stream anomaly |
| BRW | +71.3 | −384 | −374 ± 10 | **OIPC** | Near North Slope tundra; values agree |
| CBA | +55.2 | −326 | −324 ± 14 | **OIPC** | Near boreal wetlands; excellent agreement |
| MHD | +53.3 | −306 | −324 ± 14 | **Zonal** | Maritime; sources are continental |
| AZR | +38.8 | −284 | −324 ± 14 | **Zonal** | Excluded (non-MBL); no local sources |
| MLO | +19.5 | −320 | −301 ± 15 | **Zonal** | Excluded (non-MBL); altitude artefact |
| KUM | +19.6 | −287 | −301 ± 15 | **Zonal** | No local sources |
| ASC | −8.0 | −277 | −301 ± 15 | **Zonal** | Excluded; remote island |
| SMO | −14.2 | −299 | −301 ± 15 | **Zonal** | Excluded; methods agree |
| CGO | −40.7 | −294 | −301 ± 15 | **Zonal** | SH background; source via transport |
| SPO | −90.0 | −556 | −301 ± 15 | **Zonal** | OIPC meaningless for Antarctica |

Largest shifts from old global −310 ‰: NH high-lat sites move −64 to −74 ‰; SH sites move +9 ‰. The choice matters most at NH sites (heavy source contamination) and least at SH sites (minimal sources) — fortunately the SH sites dominate the KIE constraint.

Files: `data/dD_source_comparison.md`, `data/dD_source_database.json`, `data/oipc_precipitation_dD.csv`, `analysis/build_dD_source_db.py`.

---

## 3. δ¹³C Source Signatures — Low Priority

The δ¹³C of wetland CH₄ has a much smaller latitude dependence than δD:
- Tropical (C₄ influence): ≈ −56.7 ‰ (Ganesan et al., 2018)
- Boreal (C₃ dominated): ≈ −67 ‰
- Range: ~10 ‰ (vs 73 ‰ for δD)

This is primarily a C₃/C₄ vegetation effect, not a smooth latitude gradient. For our sites — all at high latitudes or remote marine — the relevant wetlands are overwhelmingly C₃-dominated. The current global constant (−62 ‰) is adequate. This is a second-order correction compared to δD and wetland flux seasonality.

---

## 4. Wetland Emission Seasonality — The Key Missing Piece

### What we need

An independent estimate of the seasonal amplitude of wetland CH₄ emissions (**ΔS_wetland/S**) in the source region of each site. This removes the need to solve for source seasonality algebraically from the isotope equations.

### Li et al. (2026) ESSD dataset

Li et al. (2026) provide monthly 1°×1° gridded natural wetland CH₄ emissions for 2000–2025, based on 35 Global Methane Budget model estimates (22 bottom-up + 13 top-down) emulated with XGBoost. Data available at: https://doi.org/10.5281/zenodo.18870108

Key features:
- Monthly resolution at 1°×1° — can extract seasonal cycle for any region
- Ensemble uncertainty from 35 model runs × 10 ERA5 ensemble members
- Covers our analysis period (2005–2010)
- Already provides latitude-band aggregates matching our site groupings

From this dataset, we extract for each site's source region:
- **Amplitude** of seasonal wetland emission cycle
- **Phase** (peak month) of wetland emissions

---

## 5. Vector Subtraction — Handling the Source–Sink Phase Shift

### The problem

Both OH and wetlands peak in summer, but not in the same month:
- **OH**: peaks ~July (driven by solar UV)
- **Boreal wetlands**: peak ~Aug–Sep (temperature + water table lag)

If the phase offset is ~1–2 months, simple amplitude subtraction is wrong. The observed amplitude is a **vector sum**, not a scalar sum:

```
A_observed² = A_sink² + A_source² + 2·A_sink·A_source·cos(Δφ)
```

### The solution: phasor decomposition

The harmonic fit gives B and C coefficients (or equivalently amplitude and phase) for each observed signal. Treat seasonal cycles as complex phasors:

```
For each isotope (δ¹³C or δD):

  Z_observed = B_obs + i·C_obs                     ← from harmonic fit (Phase 2)
  Z_source   = A_source × e^(i·φ_source)           ← from Li2026 + isotope signatures
  Z_sink     = Z_observed − Z_source                ← vector subtraction

  A_sink     = |Z_sink|
  φ_sink     = arg(Z_sink)
```

Then:

```
R_corrected = A_sink(δ¹³C) / A_sink(δD)  ≈  (α_13C − 1) / (α_D − 1)
```

### Built-in consistency check

After vector subtraction, the sink phases for δ¹³C and δD should agree — both should peak when OH peaks (~July). If they diverge, it signals residual contamination from non-wetland sources or non-OH sinks.

### How to compute the source phasor

For each isotope, the source phasor amplitude is:

```
A_source(δ¹³C) = |δ¹³C_wetland − δ¹³C_atm| × (ΔS_wetland / S_total)
A_source(δD)   = |δD_wetland − δD_atm|     × (ΔS_wetland / S_total)
```

Where:
- δ¹³C_wetland ≈ −62 ‰ (global), δ¹³C_atm = −47.3 ‰ → gap = 14.7 ‰
- δD_wetland = site-specific from our database
- ΔS_wetland / S_total = fractional seasonal amplitude of wetland emissions from Li2026
- Phase = peak month of wetland emissions from Li2026

---

## 6. Implementation Plan

### Step 1: Download Li2026 data
- Source: https://doi.org/10.5281/zenodo.18870108
- Monthly 1°×1° wetland CH₄ emissions, 2000–2025

### Step 2: Extract wetland seasonality per site
- Define source region for each site (latitude band or regional box)
- Compute monthly climatology (2005–2010 to match INSTAAR period)
- Fit annual harmonic → amplitude (ΔS_wetland/S) and phase (peak month)

### Step 3: Build source phasors
- Combine wetland flux seasonality with isotopic signatures (δD from our database, δ¹³C ≈ −62 ‰)
- Construct Z_source for each site and each isotope

### Step 4: Vector-subtract and compute corrected ratio
- Z_sink = Z_observed − Z_source for δ¹³C and δD
- R_corrected = |Z_sink(δ¹³C)| / |Z_sink(δD)|
- Propagate uncertainties via Monte Carlo (ensemble spread from Li2026 + isotope signature uncertainties + harmonic fit uncertainties)

### Step 5: Extract α_13C_OH
- R_corrected ≈ ε_13C_bulk / ε_D_bulk
- Correct for non-OH sinks (Cl, soil, stratosphere) as in current Phase 5
- Compare with Saueressig (1.0039) and Cantrell (1.0054)

### Data requirements summary

| Parameter | Source | Status |
|-----------|--------|--------|
| A_obs, φ_obs for δ¹³C and δD | Phase 2 harmonic fits | ✅ Done |
| δD_wetland per site | OIPC + Douglas 2021 | ✅ Done |
| δ¹³C_wetland | Global −62 ‰ | ✅ Adequate |
| ΔS_wetland/S and φ_wetland per site | Li et al. (2026) ESSD | ❌ To extract |
| Non-OH sink KIEs | Literature values | ✅ In Phase 5 |

---

## 7. Expected Impact

- **NH high-latitude sites** (ALT, ZEP, BRW): large source correction. Currently R ≈ 0.06–0.14 (3–10× above pure-OH). After correction, these sites should move toward R ≈ 0.013–0.020, becoming usable for KIE constraint.
- **SH sites** (CGO, SPO): small correction. Already near pure-OH prediction. Vector subtraction will tighten uncertainties.
- **Overall**: the corrected ratio from multiple sites should converge, providing a multi-site observational KIE constraint independent of lab measurements.
