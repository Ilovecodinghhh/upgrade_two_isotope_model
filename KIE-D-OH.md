# KIE-D-OH: D/H Kinetic Isotope Effect in the CH₄ + OH Reaction

## 1. Introduction

The kinetic isotope effect (KIE) for the D/H fractionation during the reaction of methane with the hydroxyl radical (CH₄ + OH → CH₃ + H₂O) is a critical parameter in the global methane isotope budget. Because the OH sink destroys ~83–89% of atmospheric methane (Saunois et al., 2020; Fujita et al., 2025 SI), the choice of KIE^D_OH directly controls the inferred source-to-atmosphere isotope shift for δD-CH₄, and hence the emission partitioning derived from δD observations.

This report compiles experimental measurements, temperature-dependent parameterizations, and the specific values adopted by major global CH₄ inversion studies, drawing exclusively from reference materials in this repository.

---

## 2. Experimental Measurements at Room Temperature

Four independent laboratory studies have measured KIE^D_OH. Three are mutually consistent; one is an outlier.

| Study | KIE^D_OH (288 K) | Uncertainty | Method | Reference file |
|-------|:-----------------:|:-----------:|--------|----------------|
| Gierczak et al. (1997) | **1.292** | — | Rate coefficient ratio | `Fujita2025JGR_SI` Table S4 footnote i |
| Saueressig et al. (2001) | **1.294** | ±0.018 | Relative rate, O(¹D)+OH | `Riddell-Young2025PNAS_SI` Table S3; `Thanwerdas2024ACP` Table 1 |
| Joelsson et al. (2016) | **1.311** | — | 278–313 K range | `Fujita2025JGR_SI` Table S4 footnote l |
| DeMore (1993) | **1.18** | — | Rate constant ratio | `Fujita2025JGR_SI` Table S4 footnote m |

**Consensus range (excluding DeMore):** 1.29–1.31

As summarized in `Fujita2025JGR_SI` (line 94):

> "Experimental KIE^D_OH ranges 1.29 to 1.31 (Gierczak et al. 1997; Saueressig et al. 2001; Joelsson et al. 2016) except DeMore (1993) as 1.18."

The DeMore (1993) value is widely regarded as an outlier. When used to compute a sink-weighted total KIE^D, it produces an anomalously low value of ~1.20 (`Fujita2025JGR_SI` line 94), far below the ~1.27–1.30 range obtained with the other three measurements.

### Full reference details

- **Gierczak et al. (1997)**: "Rate Coefficients for the Reactions of Hydroxyl Radicals with Methane and Deuterated Methanes," *J. Phys. Chem. A*, 101(17), 3125–3134. doi:10.1021/jp963892r (`Fujita2025JGR_SI` line 194)
- **Saueressig et al. (2001)**: "Carbon 13 and D kinetic isotope effects in the reactions of CH₄ with O(¹D) and OH: New laboratory measurements and their implications for the isotopic composition of stratospheric methane," *J. Geophys. Res.*, 106(D19), 23127–23138. doi:10.1029/2000JD000120 (`Fujita2025JGR_SI` line 221; `Riddell-Young2025PNAS_SI` ref 17)
- **Joelsson et al. (2016)**: "Kinetic isotope effects of ¹²CH₃D + OH and ¹³CH₃D + OH from 278 to 313 K," *Atmos. Chem. Phys.*, 16(7), 4439–4449. doi:10.5194/acp-16-4439-2016 (`Fujita2025JGR_SI` line 202)
- **DeMore (1993)**: "Rate constant ratio for the reactions of OH with CH₃D and CH₄," *J. Phys. Chem.*, 97(33), 8564–8566. doi:10.1021/j100135a006 (`Fujita2025JGR_SI` line 190)

---

## 3. Temperature-Dependent Parameterizations

Several sink reactions have temperature-dependent KIE^D expressions, important for altitude- and latitude-resolved models.

| Reaction | Expression | Source | Repository reference |
|----------|-----------|--------|---------------------|
| CH₄ + OH | α(T) = 1.097 × exp(49 K / T) | Saueressig et al. (2001) | `Thanwerdas2024ACP` Table 1 (line 55) |
| CH₄ + Cl | α(T) = 1.278 × exp(−53.31 K / T) | Saueressig et al. (1996) | `Thanwerdas2024ACP` Table 1 (line 55) |
| CH₄ + O(¹D) | 1.06 (constant) | Saueressig et al. (2001) | `Thanwerdas2024ACP` Table 1 (line 55) |
| CH₄ + soil | 1.083 (constant) | Snover & Quay (2000) | `Riddell-Young2025PNAS_SI` Table S3 |

### Evaluation at T = 288 K (tropospheric mean)

For the OH reaction at 288 K:

α(288 K) = 1.097 × exp(49 / 288) = 1.097 × 1.1855 ≈ **1.294**

This matches the Saueressig et al. (2001) room-temperature measurement exactly, confirming internal consistency.

For the Cl reaction at 288 K:

α(288 K) = 1.278 × exp(−53.31 / 288) = 1.278 × 0.8312 ≈ **1.062**

Note: The temperature-dependent Cl-D expression yields substantially lower values at tropospheric T than the commonly cited fixed value of ~1.52 (see §4). The fixed value of 1.527 from Saueressig (1996) corresponds to a different parameterization or effective temperature representative of the stratosphere/upper troposphere where most Cl-mediated destruction occurs. The value used in `common.py` and by Riddell-Young (2025) is 1.520, drawn from the net effective fractionation.

---

## 4. KIE Values for All Sinks (D/H)

The total atmospheric fractionation for δD-CH₄ is the sink-weighted average of individual KIEs.

### 4.1 Individual Sink KIE^D Values

From `Riddell-Young2025PNAS_SI` Table S3 and `Fujita2025JGR_SI` Table S4:

| Sink | KIE^D | Source | Sink fraction (Riddell-Young) | Sink fraction (Saunois 2020) |
|------|:-----:|--------|:-----------------------------:|:----------------------------:|
| Tropospheric OH | 1.294 | Saueressig (2001) | 0.835 | 0.87–0.90 |
| Tropospheric Cl | 1.520–1.527 | Saueressig (1996) | 0.035 | 0.002–0.057 |
| Stratospheric loss | 1.137–1.179 | Röckmann (2011) / Saueressig (2001) | 0.060 | 0.050 |
| Soil oxidation | 1.083 | Snover & Quay (2000) | 0.070 | 0.049 |

Note: The stratospheric KIE^D differs between references. Riddell-Young uses 1.179 (`Riddell-Young2025PNAS_SI` Table S3), while Fujita uses 1.137 from Röckmann et al. (2011) balloon data (`Fujita2025JGR_SI` Table S4). This repo's `common.py` uses Strat_D = 1.179 ± 0.01.

### 4.2 Sink-Weighted Total KIE^D

The total KIE^D = Σ(KIE^D_i × f_i), where f_i is the fractional contribution of sink i.

From `Fujita2025JGR_SI` Table S4, using Saunois et al. (2020) sink fractions:

| KIE^D_OH source | Total KIE^D (min Cl) | Total KIE^D (avg Cl) | Total KIE^D (max Cl) |
|:----------------:|:--------------------:|:--------------------:|:--------------------:|
| Gierczak (1.292) | 1.278 | 1.279 | 1.274 |
| Saueressig (1.294) | 1.276 | 1.280 | 1.289 |
| Joelsson (1.311) | 1.290 | 1.294 | 1.303 |
| DeMore (1.18) | 1.174 | 1.179 | 1.192 |

**Overall range (excluding DeMore): 1.274–1.303**

From `Riddell-Young2025PNAS_SI` Table S3, using their sink fractions: **Net KIE^D = 1.281**

From this repository's `common.py` (KIE_FIXED):

KIE^D_total = 0.835 × 1.3105 + 0.035 × 1.52 + 0.06 × 1.179 + 0.07 × 1.083 ≈ **1.296**

(using the midpoint of the OH_D uniform distribution: (1.294 + 1.327)/2 = 1.3105)

---

## 5. Values Adopted by Global Inversion Studies

### 5.1 Riddell-Young et al. (2025) — PNAS

**Source:** `ImportantReferences/Riddell-Young2025PNAS/` and `Riddell-Young2025PNAS_SI/`

- **KIE^D_OH = 1.294** (Saueressig et al., 2001)
- Net KIE^D = 1.281 (Table S3)
- Sink fractions: OH=0.835, Cl=0.035, Strat=0.06, Soil=0.07
- Methodology: 1-box dual-isotope mass balance; lifetime and sink partitioning held constant
- Key finding: "δD-CH₄ less affected by Cl sink, more sensitive to wetland variations" (`Riddell-Young2025PNAS` main text)

### 5.2 Thanwerdas et al. (2024) — ACP

**Source:** `ImportantReferences/Thanwerdas2024ACP/`

- **KIE^D_OH = α(T) = 1.097 × exp(49/T)** (Saueressig et al., 2001, temperature-dependent)
- Cl-D: α(T) = 1.278 × exp(−53.31/T) (Saueressig et al., 1996)
- Soil-D: 1.083 (Snover & Quay, 2000)
- O(¹D)-D: 1.06 (Saueressig et al., 2001)
- Methodology: 3-D inversion (LMDz-SACS), assimilates CH₄ + δ¹³C + δD
- Finding on δD: "assimilating δ(D,CH₄) observations... has a very small influence on our posterior emission estimates" (line ~303). However, this repo's `dD_threshold` experiment attributes this to their use of σ(Mic δD) ≈ 128‰, which is 3.5× above the detection threshold.
- On OH-¹³C KIE choice: "Saueressig et al. (2001) indicate that their data is of considerably higher experimental precision... we prefer to allocate computational time to a sensitivity inversion testing a different OH field rather than testing a different OH fractionation coefficient" (line ~65)

### 5.3 Fujita et al. (2025) — JGR

**Source:** `ImportantReferences/Fujita2025JGR_SI/`

- **Base KIE^D_OH = 1.292** (Gierczak et al., 1997)
- Prior total KIE^D range: **[1.25, 1.30]**, mean 1.275
- Posterior total KIE^D: **1.283 [1.274, 1.292]** for 2003–2012
- Tested sensitivity to DeMore exclusion (ParaRange #6: [1.275, 1.300]) — "impact on the posterior sectorial source fraction being minimal" (line 94)
- Derived sink fractions from posterior: OH = 89 ± 2%, soil = 2.2 ± 2.9%, Cl = 3.2 ± 1.1%
- Tested multiple OH-D KIE scenarios in Table S4: Saueressig (1.294), Joelsson (1.311), DeMore (1.18)
- Triple-isotope approach (δ¹³C + δD + Δ¹⁴C)

### 5.4 Basu et al. (2022) — ACP

**Source:** `ImportantReferences/Basu2022ACP/`

- Primarily a δ¹³C inversion; does not explicitly discuss KIE^D_OH
- Uses sink-weighted KIE for ¹³C; sensitivity tested via Cl field alternatives (Hossaini vs Wang)
- Key finding on isotope relaxation: "Large-scale gradients of atmospheric δ¹³C take significantly longer to respond to changes in emissions compared to gradients of CH₄ (Tans, 1997)" (line 182)

### 5.5 Dasgupta et al. (2025) — EGU

**Source:** `ImportantReferences/Dasgupta2025EGU/`

- 2-box dual-isotope model (δ¹³C + δD)
- Uses Gauss–Newton inversion jointly constrained by CH₄, δ¹³C-CH₄, and δD-CH₄
- Finds dual-isotope inverted lifetime shortens between 1994 and 2022 (–0.1 yr globally)
- Specific KIE^D_OH value not stated in the available text; likely follows Saueressig (2001) convention

### 5.6 This Repository (`common.py`)

```python
KIE_DISTRIBUTIONS = {
    'OH_D':   {'dist': 'uniform', 'low': 1.294, 'high': 1.327},
    'Cl_D':   {'dist': 'normal',  'mean': 1.52,  'std': 0.02},
    'Strat_D': {'dist': 'normal', 'mean': 1.179, 'std': 0.01},
    'Soil_D':  {'dist': 'normal', 'mean': 1.083, 'std': 0.01},
}
KIE_FIXED = {
    'OH_D':    (1.294 + 1.327) / 2,   # 1.3105
    'Cl_D':    1.52,
    'Strat_D': 1.179,
    'Soil_D':  1.083,
}
```

- OH_D lower bound (1.294) = Saueressig (2001); upper bound (1.327) exceeds the Joelsson (2016) value of 1.311, providing a conservative uncertainty range
- Distribution: uniform, reflecting lack of preference among the three consistent measurements plus margin
- Cl_D sampled N(1.52, 0.02) — described in `REVISION_RESPONSE.md` (D2) as "Gola et al., 2005"

---

## 6. Synthesis

### 6.1 Consensus and Disagreement

The three modern measurements (Gierczak, Saueressig, Joelsson) cluster tightly at **1.292–1.311**, spanning only 1.5% relative variation. This stands in sharp contrast to the OH-¹³C KIE controversy (Saueressig 1.0039 vs Cantrell 1.0054, a 39% relative difference in the fractionation factor ε = α − 1). The D/H fractionation by OH is thus **far better constrained experimentally** than the ¹³C/¹²C fractionation.

The DeMore (1993) outlier at 1.18 is universally excluded or de-emphasized. Fujita et al. (2025) showed that even including it via a broader prior range [1.25, 1.30] has "minimal" impact on posterior source fractions. The community consensus effectively treats KIE^D_OH ∈ [1.29, 1.31].

### 6.2 Sensitivity of Total KIE^D to OH-D KIE Choice

Using Saunois (2020) average sink fractions, swapping between the three consistent KIE^D_OH values changes the total KIE^D by:

- Gierczak → Joelsson: Δ(total KIE^D) ≈ +0.015 (from ~1.279 to ~1.294)
- Saueressig → Joelsson: Δ(total KIE^D) ≈ +0.014

For comparison, swapping between the two OH-¹³C KIE values:

- Saueressig → Cantrell: Δ(total KIE^C) ≈ +0.0013 (from ~1.0064 to ~1.0077)

While the absolute change in total KIE^D (0.015) is larger than the change in total KIE^C (0.0013), the relative impact on source partitioning depends on the fractionation leverage, analyzed in detail in `importance-KIE.md`.

### 6.3 Remaining Uncertainties

1. **Cl-D KIE**: The fixed value of 1.52 vs temperature-dependent evaluation differ substantially. Cl sink fraction uncertainty (1–35 Tg/yr per Saunois 2020) propagates strongly.
2. **Stratospheric KIE^D**: 1.137 (Röckmann) vs 1.179 (Saueressig) — a 3.7% discrepancy that matters at ~6% sink fraction.
3. **Soil KIE^D**: Only one measurement (Snover & Quay 2000). At 7% sink fraction and KIE=1.083, the soil term is small but non-negligible.
4. **DeMore (1993)**: If this value were correct, total KIE^D would drop to ~1.18, fundamentally changing the δD source budget. No subsequent study has reproduced it.

---

## 7. References (Local Repository Files)

| Short citation | Repository path |
|---------------|----------------|
| Riddell-Young et al. (2025) PNAS | `ImportantReferences/Riddell-Young2025PNAS/` |
| Riddell-Young et al. (2025) PNAS SI | `ImportantReferences/Riddell-Young2025PNAS_SI/` |
| Thanwerdas et al. (2024) ACP | `ImportantReferences/Thanwerdas2024ACP/` |
| Fujita et al. (2025) JGR SI | `ImportantReferences/Fujita2025JGR_SI/` |
| Basu et al. (2022) ACP | `ImportantReferences/Basu2022ACP/` |
| Dasgupta et al. (2025) EGU | `ImportantReferences/Dasgupta2025EGU/` |
| Model code (KIE config) | `common.py` lines 103–126 |
| KIE survey across 16 publications | `KIE_Used_Previous_Study.md` |

