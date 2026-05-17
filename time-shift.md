# time-shift: Temporal Offset Between δ¹³C-CH₄ and δD-CH₄ Atmospheric Records

## 1. Introduction

Atmospheric δ¹³C-CH₄ and δD-CH₄ both respond to changes in the methane source mix and sink fractionation, but they do not track each other in lockstep. Differences in their temporal evolution — phase lags, differential response rates, and divergent sub-period behavior — provide clues about which source categories are driving CH₄ growth and whether sink changes (particularly OH) play a role.

This report examines temporal offsets between δ¹³C and δD time series using (1) the observational data in `rel/data/`, (2) hemispheric records in `rel/data/HemMean_*.csv`, and (3) the theoretical and empirical basis from literature in `ImportantReferences/`.

---

## 2. Data Description and Coverage

### 2.1 Available Records

| Observable | File | Period | N_MC | Frequency |
|-----------|------|--------|:----:|:---------:|
| δ¹³C global | `rel/data/d13C_dei_compiled.txt` | 1999–2023 | 1000 | Annual |
| δD global | `rel/data/GlobMean_dD_dei_DasguptaCal_noBUDS.csv` | 2005–2023 | — | Annual (smoothed) |
| δD hemispheric | `rel/data/HemMean_dD_annual_DasguptaCal_noBUDS.csv` | 2005–2023 | — | Annual |
| δ¹³C hemispheric | `rel/data/ch4c13_nh_sh_mean.xlsx` | 1998.5–2023.7 | — | Monthly (deseasonalized) |
| CH₄ global | `rel/data/GML_CH4_AnnualMean.xlsx` | 1984–2023 | — | Annual |

**Critical limitation**: The δ¹³C record begins in 1999 and the δD record begins in 2005, giving only **19 years of overlap** (2005–2023). Sub-period analyses are necessarily short.

### 2.2 Global Mean Time Series (Overlapping Period)

| Year | δ¹³C (‰ VPDB) | δD (‰ VSMOW) |
|:----:|:-------------:|:------------:|
| 2005 | −47.093 | −76.76 |
| 2006 | −47.104 | −76.03 |
| 2007 | −47.091 | −76.43 |
| 2008 | −47.063 | −76.55 |
| 2009 | −47.137 | −76.31 |
| 2010 | −47.179 | −77.34 |
| 2011 | −47.223 | −77.16 |
| 2012 | −47.241 | −76.88 |
| 2013 | −47.263 | −76.11 |
| 2014 | −47.269 | −77.09 |
| 2015 | −47.303 | −77.41 |
| 2016 | −47.298 | −77.85 |
| 2017 | −47.356 | −78.17 |
| 2018 | −47.404 | −78.42 |
| 2019 | −47.487 | −79.00 |
| 2020 | −47.479 | −79.18 |
| 2021 | −47.579 | −80.88 |
| 2022 | −47.675 | −81.97 |
| 2023 | −47.781 | −82.52 |

---

## 3. Observed Temporal Patterns

### 3.1 Long-Term Trends

| Period | δ¹³C trend (‰/yr) | δD trend (‰/yr) | δD / δ¹³C ratio |
|--------|:------------------:|:----------------:|:----------------:|
| 2005–2012 | −0.023 | −0.11 | 4.6 |
| 2013–2023 | −0.050 | −0.60 | 11.8 |
| 2005–2023 (full) | −0.035 | −0.30 | 8.6 |

Both isotopes become more negative over time, consistent with increasing microbial emissions (microbial sources are isotopically light in both systems). However, the **δD decline accelerates dramatically after ~2013**, while δ¹³C accelerates more gradually. The δD-to-δ¹³C trend ratio increases from 4.6 (2005–2012) to 11.8 (2013–2023) — a 2.6× steepening of the δD response relative to δ¹³C.

### 3.2 Year-to-Year Correlation

Correlating the year-to-year changes (first differences):

**r = 0.41** (N = 18 years of differences)

This moderate positive correlation means the two isotopes move in the same direction about 56% of the time (10 out of 18 years show same-sign changes). However, the **early period (2005–2012) shows almost no coherence**: only 2 of 7 years have same-sign changes. The **late period (2015–2023)** is more coherent: 7 of 8 years agree.

This pattern — incoherent early, coherent late — is consistent with a transition from a period where different sources affected the two isotopes differently (e.g., fossil fuel changes affecting δ¹³C more than δD) to a period where a single dominant source change (microbial) drives both.

### 3.3 Normalized Trajectory Analysis

Normalizing both records to [0, 1] over the 2005–2023 span reveals the phasing:

- **2005–2013**: δ¹³C has **progressed further** through its total decline than δD. The normalized δ¹³C leads δD by up to 0.36 normalized units (in 2013). This means δ¹³C was already responding to something (likely the post-2007 microbial increase) that δD had not yet fully registered.

- **2014–2020**: The gap narrows. Both isotopes are progressing through their declines at roughly similar (normalized) rates.

- **2021–2023**: δD **catches up and briefly overtakes** δ¹³C. By 2022, δD_norm = 0.905 vs δ¹³C_norm = 0.847. The steep δD decline in 2021 (−1.7‰, the largest single-year change in the record) drives this convergence.

**Interpretation**: δ¹³C responded earlier to the microbial source increase (detectable from ~2007 as a more negative trend), while δD's response was delayed and then abruptly intensified ~2019–2021. This phasing is consistent with the hypothesis that the early microbial signal was dominated by **tropical wetland sources** with δ¹³C signatures distinct from fossil fuels but δD signatures that overlap more with other sources, making the δD signal harder to detect until the microbial contribution became large enough.

---

## 4. Physical Mechanisms for Temporal Offsets

### 4.1 Isotopic Relaxation Time (Tans, 1997)

The most fundamental reason for temporal offsets between isotopic records is the **longer relaxation time of isotope ratios** compared to total CH₄ concentration.

From `Thanwerdas2024ACP` (line 324):

> "the relaxation time for isotopic composition in the atmosphere in response to a perturbation is much larger (decades; Tans, 1997) than that for CH₄ itself"

From `Basu2022ACP` (line 182):

> "Large-scale gradients of atmospheric δ¹³C take significantly longer to respond to changes in emissions compared to gradients of CH₄ (Tans, 1997), requiring multi-decade spin-ups for models trying to simulate atmospheric δ¹³C"

The CH₄ lifetime is ~9 years, but the effective relaxation time for δ values is much longer because the isotope ratio depends on the **ratio of two species** (e.g., ¹³CH₄ / ¹²CH₄), each of which has slightly different lifetimes due to the KIE. A perturbation to the source mix shifts the isotope ratio, which then relaxes back on a timescale of τ / (α − 1), where α is the KIE.

For ¹³C: τ_relax ≈ 9 yr / 0.005 ≈ 1800 yr (for OH alone)
For D: τ_relax ≈ 9 yr / 0.30 ≈ 30 yr (for OH alone)

This enormous difference — **δ¹³C relaxes ~60× slower than δD** — means:

- δ¹³C responds sluggishly to source changes but integrates them over longer periods
- δD responds more rapidly to source changes but also reflects shorter-term variability

In practice, both isotopes respond to source changes on the CH₄ lifetime timescale (~9 yr), but the asymptotic adjustment to a new steady state takes much longer for ¹³C than for D. This means that **for sustained source changes, δD reaches its new equilibrium faster than δ¹³C**.

### 4.2 OH Decline and the ¹²CH₄–¹³CH₄ Time Lag

From `Thanwerdas2024ACP` (line 303):

> "If the OH sink is the only sink, a decline in OH concentrations has no effect on δ(¹³C,CH₄) in the long term (several decades) because the mean fractionation is not affected. However, in the short term (a decade), as OH concentrations decrease, ¹²CH₄ and ¹³CH₄ atmospheric lifetimes increase. Due to the fractionation effect, there is a time lag between increases in ¹²CH₄ and ¹³CH₄ amount fractions. ¹²CH₄ accumulates faster than ¹³CH₄, leading to a decrease in δ(¹³C,CH₄)."

This mechanism produces a **transient δ¹³C shift** during OH changes that mimics a source change. Critically, the same mechanism applies to δD but with ~60× larger KIE, meaning the transient D shift from an OH change would be proportionally much larger. However, the δD shift also relaxes faster, so the net observable effect depends on the timescale of the OH perturbation.

### 4.3 Differential Source Sensitivity

From `Riddell-Young2025PNAS` (main text):

> "δD-CH₄ less affected by Cl sink, more sensitive to wetland variations"

This differential sensitivity means:

1. **Cl sink changes** (e.g., from marine boundary layer chemistry) preferentially shift δ¹³C (Cl-¹³C KIE = 1.066 → ε = 66‰) but have a relatively smaller effect on δD (Cl-D KIE = 1.52, but Cl fraction is only ~3.5% of the total sink).

2. **Wetland emission changes** produce larger δD shifts because microbial δD sources (−310 to −320‰) are far from atmospheric δD (~−80‰), while microbial δ¹³C (−60‰) is closer to atmospheric δ¹³C (~−47‰). A given Tg increase in wetland emissions shifts δD by ~3× more than it shifts δ¹³C (in terms of fraction of the source-to-atmosphere offset).

3. **Fossil fuel changes** shift δ¹³C more efficiently because FF δ¹³C (−44‰) is far from microbial δ¹³C (−60‰), providing a 16‰ lever. For δD, FF (−190‰) is closer to atmospheric (−80‰), giving less leverage relative to microbial (−310‰).

### 4.4 Hemispheric Asymmetry

From `rel/data/HemMean_dD_annual_DasguptaCal_noBUDS.csv`:

| Hemisphere | δD trend (‰/yr, 2005–2023) | δD value 2005 (‰) | δD value 2023 (‰) |
|:----------:|:--------------------------:|:------------------:|:------------------:|
| NH | −0.46 | −85.0 | −90.0 |
| SH | −0.40 | −69.4 | −78.0 |

The NH declines ~15% faster than SH. The NH–SH gradient (NH more negative) averages ~−14.5‰ but is not constant — it widens from −11.4‰ (2009) to −17.1‰ (2020), then narrows sharply to −10.2‰ (2022) as the SH δD drops precipitously (possibly a data artifact or a real SH microbial surge).

From `experiments/Hemispheric_Divergence/MANUSCRIPT.md` findings:

> "Post-2006 CH₄ growth is driven by asymmetric hemispheric microbial trends: NH microbial emissions increase at +6.6 Tg/yr² (100% of MC positive), while SH microbial emissions are stable."

This hemispheric asymmetry naturally produces different δ¹³C and δD temporal patterns because the two isotopes have different NH–SH source signature gradients (δD has 5–10× larger hemispheric gradients than δ¹³C, per `dD_threshold` experiment findings).

---

## 5. The "Early-Late" Phase Transition (~2013)

The data reveal a distinct regime change around 2013:

### 5.1 Before 2013: δ¹³C Leads

- δ¹³C trend: −0.023‰/yr (moderate decline)
- δD trend: −0.11‰/yr (slight decline, noisy)
- Annual change correlation: poor (2/7 same-sign)
- δ¹³C-normalized has progressed further than δD-normalized

**Interpretation**: The early post-2007 microbial increase was detectable in δ¹³C because: (a) the ¹³C-based source attribution has a longer observational record with tighter measurement precision (σ ~ 0.01‰ vs ~0.5‰ for δD), and (b) the early microbial signal may have included a fossil fuel component (rising coal-mine emissions in China?) that shifts δ¹³C distinctly but has a δD signature similar to other fossil sources.

### 5.2 After 2013: δD Catches Up

- δ¹³C trend: −0.050‰/yr (accelerating)
- δD trend: −0.60‰/yr (strongly accelerating)
- Annual change correlation: strong (7/8 same-sign)
- δD-normalized convergence with δ¹³C-normalized

**Interpretation**: By 2013, the microbial source dominance was sufficiently large that both isotopes track the same dominant signal. The δD acceleration is disproportionately steep because (a) δD has more leverage for microbial sources (the δD source-to-atmosphere separation is larger), and (b) the potential contribution of wetland-specific δD depletion (very negative microbial δD) amplifies the signal.

### 5.3 2021 Anomaly

The year 2021 shows the largest single-year δD decline (−1.7‰) in the entire record, while δ¹³C also has a large decline (−0.10‰). Both are consistent with a massive microbial pulse — potentially related to increased tropical wetland emissions driven by La Niña conditions.

---

## 6. Literature Context

### 6.1 Inversion Spin-Up Implications

From `Thanwerdas2024ACP` (line 158):

> "posterior source signatures slowly move away from the prior value over time. After 2–3 years, the posterior value finally reaches a new and rather stable state. In other words, as the influence of the initial conditions on the isotopic composition decreases, the system prefers to optimize the source signatures, hence slowly reaching the posterior value."

This 2–3 year inversion spin-up reflects the same isotopic relaxation physics: the model needs time for isotopic gradients to adjust. The longer relaxation time for δ¹³C vs δD means that initial condition errors in δ¹³C persist longer in inversions, potentially affecting the early years of any time series comparison.

### 6.2 Thanwerdas et al. on δD Information Content

From `Thanwerdas2024ACP` (line ~303):

> "assimilating δ(D,CH₄) observations... has a very small influence on our posterior emission estimates"

However, the `dD_threshold` experiment in this repository demonstrated that this finding was an artifact of overly large δD uncertainty assumptions (σ ≈ 128‰). With proper uncertainty quantification (σ ≈ 8‰), δD provides substantial additional constraint, particularly in the two-box (hemispheric) framework.

### 6.3 Dasgupta et al. (2025) Dual-Isotope Lifetime

From `Dasgupta2025EGU` (line 197):

> "the dual-isotope inverted lifetime of CH₄ shortens between 1994 and 2022 (–0.1 years), with hemispheric adjustments diverging (NH: –0.2 years; SH: +0.3 years)"

A shortening lifetime (increasing OH) would produce a transient isotopic effect: δ¹³C and δD both shift toward heavier values (less negative) due to the increased fractionation, partially offsetting the source-driven trend toward more negative values. If this lifetime change is real, it could partially explain why the δD record appears "delayed" relative to δ¹³C in the early period — the OH-driven transient isotope effect is proportionally larger for D (larger KIE), temporarily counteracting the source-driven decline.

---

## 7. Quantitative Summary

### 7.1 Key Numbers

| Metric | δ¹³C | δD | Ratio |
|--------|:----:|:--:|:-----:|
| Overall trend (2005–2023) | −0.035‰/yr | −0.30‰/yr | 8.6 |
| Trend 2005–2012 | −0.023‰/yr | −0.11‰/yr | 4.6 |
| Trend 2013–2023 | −0.050‰/yr | −0.60‰/yr | 11.8 |
| Annual MC std | 0.006–0.021‰ | 0.41–0.92‰ | ~50× |
| Source-to-atm shift | ~6‰ | ~230‰ | ~40× |
| KIE relaxation time | ~1800 yr* | ~30 yr* | ~60× |
| Year-to-year correlation | r = 0.41 | — | — |

*Theoretical asymptotic; practical response occurs on the ~9 yr CH₄ lifetime.

### 7.2 Is There a Discrete Time Lag?

A simple cross-correlation of the detrended annual residuals does not reveal a statistically significant lag at annual resolution with only 19 data points. The data are more consistent with:

1. **A gradual phase transition** (not a fixed lag): δ¹³C responds first to a moderate microbial increase (2007–2013), then δD responds strongly once the microbial signal is large enough to overcome its higher noise floor (2013–2023).

2. **Different sensitivity profiles**, not a time delay: Both isotopes respond simultaneously to the same source change, but δ¹³C has better early detection (lower noise) while δD has larger eventual amplitude (more leverage).

3. **A possible OH contribution** that transiently suppresses the δD signal in the early period through the larger D/H fractionation effect (per Thanwerdas 2024, §4.2 above).

---

## 8. Implications for Dual-Isotope Inversions

1. **Sub-period sensitivity**: One-box inversions should not assume δ¹³C and δD always carry independent information at annual timescales. Before ~2013, the two isotopes show poor year-to-year coherence, suggesting they respond to different aspects of the source/sink system. After 2013, they become redundant for the dominant (microbial) signal but complementary for resolving secondary source changes.

2. **Two-box advantage**: The hemispheric framework captures the NH–SH asymmetry that drives much of the temporal offset. The NH leads the microbial trend, and the 5–10× larger δD hemispheric gradient means δD adds the most value in the two-box (not one-box) configuration — consistent with the `dD_threshold` experiment finding.

3. **Running means are critical**: Given δD's annual std of ~0.5‰ against a trend of ~0.3‰/yr, the 5-year running mean used in this repository's models (noted in `CLAUDE.md` conventions) is appropriate and necessary to extract the δD signal.

4. **The 2021 anomaly deserves caution**: The −1.7‰ single-year δD shift is 3.4σ from the annual noise level and could reflect either a real wetland pulse or a data processing artifact. Its influence on any trend analysis covering 2005–2023 is substantial.

---

## 9. References (Local Repository Files)

| Source | Path |
|--------|------|
| δ¹³C global MC data | `rel/data/d13C_dei_compiled.txt` |
| δD global data | `rel/data/GlobMean_dD_dei_DasguptaCal_noBUDS.csv` |
| δD hemispheric data | `rel/data/HemMean_dD_annual_DasguptaCal_noBUDS.csv` |
| δ¹³C hemispheric data | `rel/data/ch4c13_nh_sh_mean.xlsx` |
| CH₄ global data | `rel/data/GML_CH4_AnnualMean.xlsx` |
| Thanwerdas et al. (2024) | `ImportantReferences/Thanwerdas2024ACP/` |
| Basu et al. (2022) | `ImportantReferences/Basu2022ACP/` |
| Dasgupta et al. (2025) | `ImportantReferences/Dasgupta2025EGU/` |
| Riddell-Young et al. (2025) | `ImportantReferences/Riddell-Young2025PNAS/` |
| dD_threshold experiment | `experiments/dD_threshold/` |
| Hemispheric_Divergence experiment | `experiments/Hemispheric_Divergence/` |
| Model conventions | `CLAUDE.md` |

