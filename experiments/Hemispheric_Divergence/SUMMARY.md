# Title 3: Hemispheric Divergence Reconciliation

**Full Title:** "Hemispheric Asymmetry in Methane Source Trends Revealed by a Two-Box Dual-Isotope Model: Northern Fossil Stability vs. Southern Microbial Surge"

**Target Journal:** Nature Communications / Science Advances

---

## Scientific Question

Do Northern and Southern Hemisphere methane source trends diverge when constrained by hemispheric δ¹³C and δD observations simultaneously, and does this divergence reconcile the disagreement between 3D inversions (which find FF increases) and global box models (which find FF stable/declining)?

## Core Contradiction Addressed

Contradiction A — the biggest split in the field:
- **FF increasing:** Basu 2022, Thanwerdas 2024, Skeie 2023 (all 3D δ¹³C inversions)
- **FF stable/declining:** Riddell-Young 2025, He 2026, Fujita 2025 (all global box models with dual isotopes)

**These cannot both be correct** — unless the disagreement is a spatial aliasing artifact.

## The Reconciliation Hypothesis

A global box model averages over hemispheres where trends differ in sign:
- NH: FF may be increasing (consistent with 3D inversions seeing regional FF growth in China/Middle East)
- SH: Microbial sources surging (tropical wetlands + livestock)
- Global aggregate: net "microbial dominance" because SH growth > NH-FF growth

A 1-box model sees: "FF stable, Mic increasing" ✓ (Riddell-Young)
A 3D inversion sees: "FF increasing in NH" ✓ (Basu, Thanwerdas)
**Both could be correct — they're measuring different things.**

## Novelty

- **Only existing dual-isotope hemispheric model** (v3.0b/v3.1b)
- First study to explicitly test whether spatial aliasing explains the 1-box vs. 3D disagreement
- Can show whether NH-FF trend from a simple 2-box matches 3D inversion NH-FF posteriors

## Key Datasets and Their Roles

| Dataset | Role |
|---------|------|
| `ch4c13_nh_sh_mean.xlsx` | Hemispheric δ¹³C — core observational constraint for 2-box |
| `GML_CH4_AnnualMean.xlsx` | Global CH₄ concentration |
| NH/SH sink fractions in `inputs.py` | Asymmetric sinks drive different hemispheric KIEs |
| Basu 2022 station .nc files | NH/SH station coverage for validation |
| Dasgupta 2025 EGU SI | Tropical wetland vs. agriculture — validate SH-Mic |
| Zhang 2021 (GOSAT tropical) | Independent: tropical livestock + wetlands |
| Chandra 2024 (post-2019 surge) | Extend to check 2019–2022 SH-dominated growth |
| CarbonTracker_CH4.xlsx | Independent emission estimates for cross-checking |

## Falsifiable Prediction

The 2-box model shows:
1. NH-FF trend is *positive* (agreeing with 3D inversions)
2. SH-Mic trend is *strongly positive* (consistent with tropical wetland studies)
3. Global aggregate from 2-box shows FF stable/declining (agreeing with 1-box models)

If all three hold, the field's biggest disagreement is resolved as a spatial aliasing effect.

## Impact If Confirmed

- Reconciles the last decade's most contentious debate in atmospheric methane science
- Demonstrates that simple multi-box models can bridge the gap between 1-box and 3D approaches
- Has direct policy implications: NH fossil mitigation targets may need to be stronger than global-average analyses suggest
