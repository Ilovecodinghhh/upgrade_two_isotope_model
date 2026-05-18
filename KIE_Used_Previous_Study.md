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
| **Yu et al. 2026 (Nature Comm)** | **Default 1.0039 (Saueressig)**; tests 1.0039, 1.0046, 1.0061 | 4D-Var (GEOS-Chem adjoint), CH₄ + δ¹³C + δD + C₂H₆ | The three test values span the full literature range: 1.0039 (Saueressig, lower bound), 1.0046 (midpoint), 1.0061 (upper bound, likely from Gierczak et al. 1997). **1.0061 is not from Cantrell (1.0054)** — it is the upper end of refs 76–84 range. After testing, they adopt 1.0039 because it gives best agreement with observed δ¹³C (Supplementary Fig. S16). |
| **Riddell-Young et al. 2025 (PNAS)** | **1.0054 (Cantrell)** → net sink KIE = 1.0082 | Two-box model, δ¹³C + δD | SI Table S3 explicitly lists OH KIE = 1.0054 (ref 1 = Cantrell). Weighted with Cl (1.066, 3.5%), stratosphere (1.003, 6%), soil (1.020, 7%), OH (83.5%) → net = 1.0082. Fixed, zero uncertainty. **Previously misidentified as Saueressig-consistent; corrected.** |
| **Schwietzke et al. 2016 (Nature)** | **Net ε = −6.3‰ [−7.1, −5.5]** (sampled distribution) | Global isotope mass balance, δ¹³C | Does NOT use an OH-specific KIE. Uses a single net fractionation factor from Miller (2005), sampled as a normal distribution in MC. Net ε = −6.3‰ → net KIE ≈ 1.0063. The range (−7.1‰ to −5.5‰, i.e. net KIE 1.0055–1.0071) spans values consistent with Saueressig OH and partially with Cantrell OH, but central value leans closer to Saueressig when back-calculated. **Not a direct Saueressig adopter — effectively sidesteps the choice by treating net fractionation as uncertain.** |
| **Rice et al. 2016 (PNAS)** | **Saueressig** (base case) | Box model, δ¹³C | Sensitivity scenario S7 tests "decreased KIE of OH sink." Claims effect on decadal trends is "minimal." |
| **Worden et al. 2017 (Nature Comm)** | **ε = −6.8‰** (sink-weighted, fixed) | One-box δ¹³C model | Does not decompose into individual reaction KIEs. −6.8‰ is broadly consistent with **Saueressig**. |
| **Fujita et al. 2025 (JGR)** | **Net KIE_C = 1.0065 [1.005–1.008]** | Bayesian inversion, CH₄ + δ¹³C + δD + Δ¹⁴C | Treats net sink KIE as a free Bayesian parameter. Prior range spans from slightly above Saueressig to slightly above Cantrell territory. |

## Theoretical Calculations

| Study | OH-¹³C KIE Calculated | Method | Notes |
|-------|----------------------|--------|-------|
| **Melissas & Truhlar 1993 (J. Chem. Phys.)** | **1.005 at 273–353 K** | Ab initio VTST (IVTST/SCT), MP-SAC2//MP2/adj-cc-pVTZ | Interpolated variational transition state theory with multidimensional tunneling. Calculated k₁₂/k₁₃ = 1.005 at all atmospheric temperatures (250–353 K), in "excellent agreement" with Cantrell et al. (1990). Temperature-independent. **Provides independent first-principles support for the Cantrell value.** Predates Saueressig (2001) by 8 years. Also showed that the ¹³C KIE's near-unity value arises from cancellation of inverse vibrational/rotational contributions against normal translational/tunneling contributions. See `ImportantReferences/Melissas1993ACP/`. |

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
| **Default to Saueressig (1.0039)** | 4 | Thanwerdas 2024, Basu 2022, Rice 2016, Worden 2017 |
| **Default to Cantrell (1.0054)** | 2 | Dasgupta 2025 (EGU abstract), Riddell-Young 2025 |
| **Net ε distribution (sidesteps the choice)** | 1 | Schwietzke 2016 |
| **Explicitly test multiple values** | 3 | Chandra 2024, He 2026 (JGR), Yu 2026 |
| **Treat as free Bayesian parameter** | 1 | Fujita 2025 |
| **Theoretical support for Cantrell** | 1 | Melissas & Truhlar 1993 (ab initio VTST → 1.005) |

### Key Observations

1. **Saueressig (1.0039) remains the most common default** (4 of 11 studies), but the gap is smaller than previously assessed.
2. **Riddell-Young 2025 uses Cantrell**, not Saueressig — their SI Table S3 explicitly lists OH KIE = 1.0054 (corrected from earlier analysis).
3. **Yu et al. 2026's third test value (1.0061)** is NOT Cantrell — it is the upper bound of the literature range from Gierczak et al. (1997). Despite testing three values, they adopt Saueressig (1.0039).
4. **Schwietzke 2016 does not pick an OH KIE** — they use a net fractionation distribution from Miller (2005), effectively marginalizing over the uncertainty.
5. **Melissas & Truhlar (1993) provide independent theoretical evidence favoring Cantrell.** Their ab initio variational transition state theory calculation yields k₁₂/k₁₃ = 1.005 at 273–353 K, in excellent agreement with Cantrell (1.0054) and notably higher than Saueressig (1.0039). This contradicts Basu et al.'s (2022) statement that "there is no independent evidence supporting one set of coefficients over another." The theoretical calculation predates Saueressig's measurement by 8 years.
6. **No study has used a systematic observational discriminant** to choose between the two values — they either pick one, test both without resolving the dispute, or marginalize over the uncertainty.
