# The D/H Kinetic Isotope Effect: The Dominant Uncertainty in Dual-Isotope CH₄ Source Partitioning

## Executive Summary

Cross-study comparison reveals that the **effective D/H kinetic isotope effect (KIE) of the atmospheric CH₄ sink is the single largest uncertainty** in using δD-CH₄ for source attribution. Three major studies — Rice et al. (2016), Thanwerdas et al. (2024), and Riddell-Young et al. (2025, "Ben") — all use the **same laboratory measurements** (Saueressig et al. 2001) but arrive at different effective KIEs due to how they apply these measurements to the atmosphere. The resulting ~57‰ shift in derived source δD **exceeds the isotopic separation between fossil fuel and biomass burning sources**, fundamentally limiting δD's power for FF vs BB partitioning.

**Practical recommendation**: use δD primarily for **microbial vs non-microbial separation** (where the isotopic gap is ~100-130‰), not for FF vs BB separation (where the gap is only ~10-50‰ and smaller than the KIE uncertainty).

---

## 1. The Three Studies

| Feature | Rice et al. 2016 (PNAS) | Thanwerdas et al. 2024 (ACP) | Riddell-Young 2025 (GRL) |
|---------|-------------------------|------------------------------|--------------------------|
| **Model type** | 3D GEOS-Chem Bayesian | 3D LMDz-SACS variational | 1-box Monte Carlo |
| **δD in inversion?** | **NO** (validation only) | YES (but found unhelpful) | YES (MC constraint) |
| **Effective KIE_D** | 1.264 (Fig S6 box model) | Emergent from 3D (not stated) | 1.281 (weighted sum) |
| **Period** | 1984–2009 | 1998–2018 | 2005–2018 |
| **δD data used** | Cape Meares (1 station) | NOAA + INSTAAR network | Ben's own global mean |

## 2. All Three Use the Same Lab Data

The OH + CH₃D reaction KIE comes from **Saueressig et al. (2001)**:

$$
\frac{k(\text{CH}_4)}{k(\text{CH}_3\text{D})} = 1.097 \times \exp\left(\frac{+49}{T}\right)
$$

> **Note**: The Thanwerdas 2024 Table 1 (from OCR) shows `exp(-49K/T)` but this is an extraction error. The correct sign is **positive**, confirmed by: at T = 296K (lab), the formula gives `1.097 × exp(49/296) = 1.294`, matching Ben's cited value of 1.294 from the same Saueressig paper.

| Reaction | 13C KIE | D/H KIE | Source |
|----------|---------|---------|--------|
| OH + CH₄ | 1.0039 | 1.097 × exp(49/T) | Saueressig 2001 |
| O(¹D) + CH₄ | 1.013 | 1.06 | Saueressig 2001 |
| Cl + CH₄ | 1.043 × exp(6.455/T) | 1.278 × exp(53.31/T) | Saueressig 1995/1996 |
| Soil uptake | 1.020 | 1.083 | Snover & Quay 2000 |

## 3. Where the Discrepancy Comes From

### 3.1 Temperature Dependence

The D/H KIE has **strong temperature dependence** (unlike ¹³C KIE which is nearly constant):

| T (K) | OH KIE_D | Context |
|--------|----------|---------|
| 230 | 1.358 | Stratosphere |
| 250 | 1.335 | Upper troposphere |
| 270 | 1.315 | Tropical troposphere |
| 277 | 1.309 | Flux-weighted mean? |
| 296 | 1.294 | **Lab temperature** |

**Ben evaluates at lab temperature (296K)**, giving OH_KIE = 1.294. The actual tropospheric CH₄ destruction occurs at ~270K, where OH_KIE = 1.315. This means Ben **underestimates** the OH fractionation.

### 3.2 Sink Partitioning

| Sink | Ben (fraction) | Rice (fraction) |
|------|---------------|----------------|
| OH | 0.835 | ~0.87 |
| Cl | 0.035 | (not explicit) |
| Stratosphere | 0.07 | ~0.077 |
| Soil | 0.06 | ~0.058 |

### 3.3 Stratospheric KIE

- **Ben**: 1.179 (Dyonisius 2020, Beck 2018 — empirical from firn air)
- **Rice**: 1.153 (earlier estimates)
- **Thanwerdas**: computed from OH + O(¹D) + Cl at stratospheric temperatures

### 3.4 Net Effect

| Study | Net KIE_D | α_D = 1/KIE | ε_D (‰) |
|-------|-----------|-------------|---------|
| Ben (lab T = 296K) | 1.281 | 0.781 | -219 |
| Ben (corrected T = 270K) | 1.300 | 0.769 | -231 |
| Rice (effective) | 1.264 | 0.791 | -209 |
| Thanwerdas (3D emergent) | ~1.26–1.30 | ~0.77–0.79 | ~-210 to -230 |

## 4. Impact on Source Attribution

Using the 1-box model to derive source δD from atmospheric observations:

| Parameter set | Derived source δD (2005–2017 mean) |
|--------------|------------------------------------|
| Ben (KIE = 1.281) | ~ -265‰ |
| Rice (KIE = 1.264) | ~ -210‰ |

**The 57‰ shift is larger than the source separations:**
- FF vs BB: 10–50‰ (depending on whose signatures)
- Microbial vs non-microbial: 100–130‰

## 5. Source δD Signature Disagreements

| Source | Rice 2016 | Ben 2025 | Thanwerdas 2024 |
|--------|-----------|----------|-----------------|
| FF (fossil fuel) | -175‰ | -186‰ | -183‰ |
| Microbial (weighted) | -316‰ | -299‰ | -310‰ (AGW) |
| BB (biomass burning) | -169‰ | -217‰ | -200‰ |
| Wetlands | -322‰ | (in Mic) | -320 to -360‰ |

The BB δD disagreement (169 vs 217 vs 200‰) is itself ~50‰ — comparable to the KIE uncertainty.

## 6. Recommendations for Our Model

### What δD CAN do:
- **Separate microbial from non-microbial** sources (δD gap ~ 100-130‰, robust across all studies)
- Provide an independent check on total microbial trends
- Detect gross inconsistencies in the δ¹³C-based partitioning

### What δD CANNOT reliably do (currently):
- Distinguish FF from BB (δD gap ~ 10-50‰, smaller than KIE uncertainty)
- Provide precise absolute source strengths for individual categories

### Model design implications:
1. **Use a 2×2 system**: δD constrains Microbial vs Non-Microbial total; δ¹³C handles FF vs BB within non-microbial
2. **Include KIE uncertainty**: sample effective T_OH from N(270, 15) K → propagates to OH_KIE
3. **Do NOT fix KIE at lab temperature**: always evaluate at tropospheric-mean T
4. **Use consistent T for Cl and OH**: both have strong T-dependence

---

## Appendix: Verification

- **Output plot**: `Output_dD_comparison/rice2016_replication.png`
- **Ben's code**: `Ben-BoxModel/Riddell-Young_2025_MassBalancePackage/.../dD_MassBalance_MC.py` lines 71–82
- **Thanwerdas Table 1**: `research-methane-isotope-db/methane_isotope_db/mineru_extractions_trimmed/Thanwerdas2024ACP.md`
- **Rice Fig S6 parameters**: τ = 9.7 yr, α_13C = 1.0056, α_D = 1.264 (from SI PDF)
