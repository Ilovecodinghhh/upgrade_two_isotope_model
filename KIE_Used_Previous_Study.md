# OH-¹³C KIE Values Used in Previous Studies

A systematic survey of 16 publications from the `ImportantReferences/` directory, checking which OH-¹³C Kinetic Isotope Effect (KIE) value each study adopted in their model.

The two contested values:
- **Saueressig et al. (2001)**: α = 1.0039
- **Cantrell et al. (1990)**: α = 1.0054

---

## Studies That Use an Isotope Mass-Balance Model (KIE Matters Directly)

| Study | OH-¹³C KIE Used | Model Type | Notes |
|-------|-----------------|------------|-------|
| **Thanwerdas et al. 2024 (ACP)** | **1.0039 (Saueressig)** | 3D inversion (LMDz), ¹²CH₄ + ¹³CH₄ | Explicit in KIE table. Also includes Cl, O(¹D), soil KIEs all from Saueressig. |
| **Dasgupta et al. 2025 (EGU)** | **1.0054 (Cantrell)** | Two-box model, δ¹³C + δD | Explicit in Table 2. **The only study that defaults to Cantrell.** Conference abstract (not full peer-reviewed paper). |
| **Basu et al. 2022 (ACP)** | **1.0039 (Saueressig)** default; also tests **1.0054 (Cantrell)** | 1-box mass balance, δ¹³C | Most transparent about the controversy. States: "this earlier result has not been refuted in the literature, nor is there any independent evidence supporting one set of coefficients over another." |
| **Chandra et al. 2024 (Comm Earth Env)** | **Tests both**: L-KIE = 1.0039, H-KIE = 1.0054 | 3D ACTM | Found ~1.2‰ offset between the two, but "impact on the temporal δ¹³C-CH₄ trend was negligible." Does not pick a winner. Also mentions Whitehill et al. (2017) value of 1.0061. |
| **He et al. 2026 (JGR)** | **Tests two**: EXP1 (medium) vs EXP3 (higher fractionation) | 3D AM4.1, ¹²CH₄ + ¹³CH₄ | EXP1: sink-weighted ε ≈ −6.32‰ (Saueressig range). EXP3: ε ≈ −5.70‰ (Cantrell range or above). Found higher KIE requires negative BMB emissions → considered less plausible. |
| **Yu et al. 2026 (Nature Comm)** | **Tests three**: 1.0039, 1.0046, 1.0061 | 4D-Var (GEOS-Chem adjoint), CH₄ + δ¹³C + δD + C₂H₆ | Most thorough KIE sensitivity analysis. Notes literature range spans 1.0039–1.0061. Finds "sensitivity of model to KIE of OH has limited seasonal dependencies." |
| **Riddell-Young et al. 2025 (PNAS)** | **Net sink KIE = 1.0082** (fixed, zero uncertainty) | Two-box model, δ¹³C + δD | Net value consistent with **Saueressig 1.0039 for OH** weighted with Cl (1.066), soil (1.020), stratosphere (1.014). Does not vary the KIE at all. |
| **Schwietzke et al. 2016 (Nature)** | **~Saueressig range** (sampled in Monte Carlo) | Global isotope mass balance | Fractionation sampled as uncertain parameter in MC ensemble. Effectively Saueressig-anchored. |
| **Rice et al. 2016 (PNAS)** | **Saueressig** (base case) | Box model, δ¹³C | Sensitivity scenario S7 tests "decreased KIE of OH sink." Claims effect on decadal trends is "minimal." |
| **Worden et al. 2017 (Nature Comm)** | **ε = −6.8‰** (sink-weighted, fixed) | One-box δ¹³C model | Does not decompose into individual reaction KIEs. −6.8‰ is broadly consistent with **Saueressig**. |
| **Fujita et al. 2025 (JGR)** | **Net KIE_C = 1.0065 [1.005–1.008]** | Bayesian inversion, CH₄ + δ¹³C + δD + Δ¹⁴C | Treats net sink KIE as a free Bayesian parameter. Prior range spans from slightly above Saueressig to slightly above Cantrell territory. |

## Studies That Do Not Use a ¹³C KIE (Concentration-Only or OH-Focused)

| Study | Approach | Why KIE Is N/A |
|-------|----------|---------------|
| **He et al. 2026 (Science)** | 3D inversion of CH₄ mixing ratios + TROPOMI | No isotope model; constrains total emissions only. |
| **Maasakkers et al. 2019 (ACP)** | GOSAT XCH₄ inversion | Optimizes gridded emissions + OH. No isotope tracers. |
| **Naus et al. 2019 (ACP)** | MCF-based OH inversion (two-box) | Constrains OH lifetime via methyl chloroform. No ¹³C-CH₄. |
| **Zhang et al. 2021 (Nature Comm)** | GOSAT CH₄ inversion (GEOS-Chem) | Optimizes emissions, wetlands, OH. No isotopic tracers. |
| **Zhao et al. 2023 (ACP)** | Bottom-up vs top-down CH₄ sink reconciliation | Focuses on OH precursor observations to correct chemistry models. No isotope budget. |

---

## Summary

Of the **11 studies that use a ¹³C KIE** in their isotope model:

| Category | Count | Studies |
|----------|-------|---------|
| **Default to Saueressig (1.0039)** | 5 | Thanwerdas 2024, Basu 2022, Schwietzke 2016, Rice 2016, Worden 2017 |
| **Default to Cantrell (1.0054)** | 1 | Dasgupta 2025 (EGU abstract) |
| **Fixed net KIE consistent with Saueressig** | 1 | Riddell-Young 2025 |
| **Explicitly test multiple values** | 3 | Chandra 2024, He 2026 (JGR), Yu 2026 |
| **Treat as free Bayesian parameter** | 1 | Fujita 2025 |

**Saueressig (1.0039) is the dominant community default.** No study has used a systematic observational discriminant to choose between the two values — they either pick one, test both without resolving the dispute, or marginalize over the uncertainty.
