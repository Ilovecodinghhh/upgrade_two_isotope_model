# Improving KIE Constraints at NH High-Latitude Sites — Discussion Notes

*2026-05-18 — working notes, not yet implemented*

---

## The Problem

At NH high-latitude sites (ALT, ZEP, BRW, etc.), the observed seasonal amplitude ratio R = A(δ¹³C)/A(δD) is 3–10× above the pure-OH prediction. Source seasonality (wetlands, rice) peaks in the same season as OH, inflating the δ¹³C cycle disproportionately.

Our current source-correction approach (Phase 5 Approach 2) uses a two-step method:

1. **Step 1:** Solve for ΔS/S and ΔL/L using CH₄ ppb + δD (independent of α_13C_OH)
2. **Step 2:** Subtract the source contribution from δ¹³C to isolate the sink KIE

This works in principle, but δD is **not** assumed free of source influence — the source effect on δD is explicitly modeled via `(δ_src_D − δ_atm_D) × ΔS/S`. The method trades the α_13C_OH unknown for an assumed wetland δD signature (−310 ± 30‰), which is the dominant uncertainty source and widens the 95% CI by ~2.4× compared to the SH direct approach.

## Why the δD Source Signature Matters So Much

| Parameter | Value | Uncertainty (1σ) | Impact |
|-----------|-------|-------------------|--------|
| δD_source (wetland) | −310‰ | ±30‰ | Dominant — changes source–atm gap by ~13% |
| δD_atm | −86‰ | small | Well-measured |
| Source–atm gap | 224‰ | ~30‰ | Directly scales ΔS/S estimate |

A 30‰ shift in δD_source → ~13% shift in ΔS/S → propagates into the residual sink amplitude → shifts the extracted α_13C_OH.

## Possible Directions to Improve NH Constraints

### 1. Better constrain source δD signature

- **Wintertime baseline:** In winter, OH is weak and wetland emissions near zero. The winter-to-summer δD change is more directly attributable to the summer pulse. Could extract site-specific source δD from seasonal anomalies.
- **Latitude-specific signatures:** Boreal wetlands (high-lat) may have δD ≈ −330 to −350‰ (more depleted precipitation) vs tropical wetlands ≈ −280 to −300‰. Using latitude-dependent source δD instead of a single global −310‰ could reduce bias.
- **Precipitation δD data:** Use GNIP/IAEA global precipitation isotope network to estimate local meteoric water δD → infer local wetland CH₄ δD.

### 2. Add a fourth observable

Currently: 3 observables (CH₄, δ¹³C, δD), 2 unknowns (ΔS/S, ΔL/L). A fourth constraint could allow solving for more unknowns (e.g., let δD_source be free, or separate wetland vs rice vs fossil source seasonality).

Candidates:
- **Ethane (C₂H₆):** Purely fossil/biomass burning, no microbial source. Seasonal cycle traces OH seasonality without source contamination. Could isolate OH component.
- **¹⁴CH₄ (radiocarbon methane):** Distinguishes fossil (¹⁴C-dead) from biogenic. Directly constrains fossil vs microbial source seasonality if co-located data exists.
- **CH₃Cl or CH₃Br:** Different OH-KIEs, could provide independent constraint on OH seasonality magnitude.

**Bottleneck:** Co-located data availability is limited for all of these.

### 3. Use transport model output to prescribe source seasonality

Instead of solving for ΔS/S from our box model, prescribe it from:
- Wetland process models (WetCHARTs, LPJ-wsl)
- Global CH₄ inversions (GOSAT-based, TROPOMI-based)

This makes the isotope equations over-determined, and both δ¹³C and δD can constrain the sink KIE.

**Risk:** Model-dependent. Wetland models disagree on high-latitude seasonality by factors of ~2.

### 4. Exploit phase information

Currently we only use **amplitudes** from harmonic fits. The **phases** carry additional information:
- OH peaks ~July (NH)
- Wetlands peak slightly later (~August–September at high latitudes due to soil thaw lag)
- A two-component harmonic model (OH at one phase + source at another phase) could separate them even when they partially overlap temporally
- Requires sufficient intra-annual sampling (~20–40 flask samples/year) — marginal but possible

### 5. Multi-site joint inversion (most promising?)

Currently each site is treated independently. A joint fit across all 12 sites with:
- **Shared parameters:** α_13C_OH, α_D_OH, OH seasonality latitude profile
- **Site-specific parameters:** local source seasonality (ΔS/S), local source signatures

The SH sites anchor the sink side, NH sites anchor the source side, and the latitude gradient provides additional leverage. This would have far more statistical power than site-by-site analysis.

### 6. Wait for more δD data

INSTAAR overlap is only 2005–2010 (~5 years). IMAU and other labs are now measuring δD at more sites from ~2022 onward. Extended records in a few years may improve harmonic fits and tighten constraints.

## Tentative Ranking (by promise vs effort)

| Rank | Direction | Promise | Effort | Notes |
|------|-----------|---------|--------|-------|
| 1 | Multi-site joint inversion | High | Medium | Best leverage from existing data |
| 2 | Latitude-specific source δD | Medium | Low | Low-hanging fruit, reduces dominant uncertainty |
| 3 | Phase-based separation | Medium | Medium | Clever but needs sufficient temporal resolution |
| 4 | Adding ethane as 4th observable | High | High | Powerful if co-located data exists |
| 5 | Transport model source constraint | Medium | High | Model-dependent, but could over-determine system |
| 6 | Wait for more data | Guaranteed | None | Time is the only cost |

## Status

Discussion only — no implementation yet. To be revisited after further consideration.
