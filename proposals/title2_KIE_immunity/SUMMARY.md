# Title 2: Dual-Isotope KIE Immunity

**Full Title:** "Dual-Isotope Constraints Halve the KIE-Driven Ambiguity in Post-2007 Fossil vs. Microbial Methane Attribution"

**Target Journal:** Global Biogeochemical Cycles

---

## Scientific Question

Does simultaneously assimilating δ¹³C and δD observations reduce the sensitivity of methane source partitioning to the unresolved OH–¹³C kinetic isotope effect (KIE = 1.0039 vs. 1.0054), and if so, by how much?

## Core Contradiction Addressed

Contradictions A + E:
- The OH–¹³C KIE — a single lab parameter — shifts FF/Mic attribution by 20–40%
- Cantrell (1990): KIE = 1.0054 → more FF attribution
- Saueressig (2001): KIE = 1.0039 → more Mic attribution
- Basu 2022 demonstrated this sensitivity explicitly
- **No study has tested whether adding δD immunizes the budget against this choice**

## Novelty

- Your model already samples the full KIE range (unique in the field)
- The 2×2 framework solves δ¹³C and δD *independently* → can directly compare how much each isotope's answer depends on KIE
- The 3×3 framework solves *simultaneously* → can show the synergistic constraint
- No other study has structured this as a formal sensitivity experiment

## Key Datasets and Their Roles

| Dataset | Role |
|---------|------|
| KIE sampling framework (`inputs.py`) | Experimental design: full KIE space |
| `FF_d13C_GlobMC_EDGAR.csv` + `FF_d13C_GlobMC_CTCH4.csv` | Two FF signature sets — test inventory dependence |
| `d13C_dei_compiled.txt` | δ¹³C observational uncertainty |
| `GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx` | δD observations |
| Basu 2022 station .nc files | Validation against 3D inversion posteriors |
| Schwietzke 2016 (revised FF signatures) | Test signature sensitivity × KIE interaction |

## Falsifiable Prediction

The dual-isotope approach reduces the KIE-driven spread in FF emissions by ≥40% compared to δ¹³C-only, making the FF/Mic partition robust to the Saueressig-vs-Cantrell choice.

## Impact If Confirmed

The field can stop waiting for a new lab measurement of the OH–¹³C KIE and instead invest in δD monitoring networks. Policy-relevant methane source attribution no longer hinges on a 35-year-old unresolved lab discrepancy.
