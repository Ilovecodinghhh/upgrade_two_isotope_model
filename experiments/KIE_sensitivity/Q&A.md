# Q&A — KIE Sensitivity Experiment Figures

---

## Fig 1: KSR Summary (2×2 histograms)

**Q: How is the p-value calculated?**

A: It's a **two-sample Kolmogorov–Smirnov (KS) test** comparing the distribution of Δ Emission values from the Saueressig run vs the Cantrell run. The KS test measures the maximum distance between the two empirical CDFs and returns a p-value for whether they could have been drawn from the same underlying distribution. Computed using `scipy.stats.ks_2samp`.

**Q: What is "Δ Emission"?**

A: It's the **trend change** in emissions (Tg/yr) over the analysis period. Specifically: the smoothed (5-year running mean) emission for the final window (2020–2022) minus the initial window (2005–2007). So each iteration gives one Δ value representing "how much did FF (or Mic) emissions change from early to late period." The histograms show the distribution of these Δ values across 1000 Monte Carlo iterations.

---

## Fig 2: Uncertainty Timeseries

**Q: How is the uncertainty calculated?**

A: It's the **2σ (2 × standard deviation) across the 1000 MC iterations** for each year. For each year *j*, we compute `2 × np.nanstd(emissions[j, :])` where `emissions[j, :]` is the vector of 1000 iteration estimates for that year. The "δ¹³C-only (2σ)" line uses the sampled-KIE Run C from Phase 1; the "Dual isotope (2σ)" line uses the same Run C from Phase 2. The green fill = reduction between the two.

---

## Fig 3: Emission Timeseries (median + 2.5/97.5 CI)

**Q: Why is the Dual Isotope median for FF so low?**

A: Because the WLS dual-isotope system is poorly conditioned. When you add the δD equation as a hard constraint, the solver (bounded `lsq_linear`) is pulled by the large δD source-signature uncertainties (δD_FF ≈ −180 ± 50 ‰, δD_Mic ≈ −310 ± 30 ‰). The δD balance often "prefers" low FF solutions because the δD-inferred source is dominated by the OH-D KIE uncertainty (α_D = 1.294–1.327) which systematically shifts the inferred source-weighted δD toward more negative values, mimicking a microbial-dominated mixture. The WLS median for FF gets dragged down as a result — this is exactly the failure mode that Phases 1–5 diagnose (KSR < 1 = coupling makes things *worse*).

---

## Fig 4 (obsolete) / Fig 5: 2-Box Fixed Model

**Q: Why is the KSR for Mic-2box so differently high?**

A: In the 2-box (NH/SH) framework the Mic signal is amplified because:
1. **SH is microbial-dominated** — most tropical wetlands are in or near the SH tropics; splitting hemispheres isolates a hemisphere where Mic >> FF.
2. **Exchange isotopic terms** partially cancel the δ¹³C KIE effect for FF (because NH↔SH exchange fluxes carry isotopic information that stabilises the NH FF estimate) but do *not* cancel for Mic (because Mic dominates the SH budget and the exchange correction adds noise there).
3. So the Cantrell–Saueressig *spread* for Mic is damped more in the 2-box than in the 1-box (the exchange isotope flux acts as a weak additional constraint), giving a higher KSR_Mic. Meanwhile KSR_FF stays similar because FF is mainly an NH quantity where the exchange term is smaller relative to the source.

---

## Fig 5 (was fig4 in some versions): 2-Box Fixed — x-axis span

**Q: Why does the Δ-Emission span go from −300 to +500? Can you make it smaller?**

A: The wide span comes from the dual-isotope WLS in the 2-box model — the poorly-conditioned system generates some extreme trend estimates when all the KIE + source signature + exchange uncertainties compound. The large tails reflect the fundamental instability of the coupled WLS approach (which is the whole point of Phases 1–5: WLS coupling doesn't work). If you want to truncate the axis for visual clarity, you could set `ax.set_xlim(-200, 200)` in `phase4b_two_box_fixed.py` (in the bottom-row histogram panels) to focus on the bulk of the distribution while clipping outlier tails.

---

## Fig 6: Weight Sweep

**Q: What does "KIE spread" mean here?**

A: **KIE spread = |mean(Δ_emission)_Cantrell − mean(Δ_emission)_Saueressig|** (in Tg/yr). It measures how much the *emission trend* changes when you switch from the Saueressig OH-¹³C KIE to the Cantrell value. A small spread means the inversion is insensitive to the KIE choice (good). A large spread means the KIE matters a lot (bad). The x-axis is the δD weight `w_dD` in the WLS; the y-axis shows KIE spread. The figure shows that increasing `w_dD` from 0 (pure δ¹³C) monotonically *increases* KIE spread — i.e. no optimal weight exists.

---

## Fig 7: Cl × Weight Interaction

**Q: Where do `high_Cl` and `default` values come from? Why is KIE spread in Tg/yr?**

A: The Cl sink fraction configurations come from the literature:
- **`default` (Cl = 3.5%):** standard global budget (e.g. Saunois et al. 2020, Table S4: 3–5% Cl fraction)
- **`thanwerdas_low` (Cl = 0.6%):** Thanwerdas et al. (2024) lower-end estimate
- **`high_Cl` (Cl = 6.5%):** represents an upper-bound scenario (Allan et al. 2007-era estimates)

The remainder of the sink is reallocated to OH so the total sums to 1.

"KIE spread" is in **Tg/yr** because it's defined the same way as in Fig 6: the difference in *emission trends* (Tg/yr) when switching between Saueressig and Cantrell. It has physical units because it propagates a dimensionless KIE difference through the mass-balance into emission estimates.

---

## Fig 8: Agreement Framework

**Q: What is "spread = 1.98 Tg/yr"?**

A: That's the **KIE spread for FF** under the agreement-filtered ensemble = |mean(Δ_FF)_Cantrell − mean(Δ_FF)_Saueressig| after filtering. Compared to the δ¹³C-only baseline spread of ~2.0 Tg/yr (Phase 1), it means the agreement filter gives a similar-magnitude spread — i.e. KSR ≈ 1 at the Phase 6 default threshold of 100 Tg/yr. (Phase 6b and 8a later showed that at threshold = 50 Tg/yr, this drops to 0.62 Tg/yr → KSR = 3.21.)

---

## Fig 9: Threshold Sweep

**Q: What is "sampled"?**

A: "Sampled" means the OH-¹³C KIE is **drawn uniformly from [1.0039, 1.0054]** at each MC iteration — i.e. it's neither fixed to Saueressig nor to Cantrell, but treats the KIE itself as an uncertain parameter. This represents the "agnostic" case where you don't know which lab value is correct. Its agreement rate falls between the two fixed values, as expected.

---

## Fig 10: Agreement Timeseries + Lifetime Effect

**Q: Where is the "varying lifetime" and what is its source?**

A: The varying lifetime is defined in `common.py` as:

```python
τ(t) = 9.0 − 0.017 × (t − 2010)     # years
```

Source: **He et al. (2026), Science** — they estimated a declining CH₄ lifetime of approximately −0.017 yr/yr from combined δ¹³C + δD + ¹⁴CH₄ constraints. At year 2010, τ = 9.0 yr (standard IPCC value); by 2022 it's ~8.8 yr; back in 1999 it was ~9.19 yr. The "fixed" mode simply uses τ = 9.0 yr for all years.

Fig 10(b) shows bars comparing the agreement rate under varying τ vs fixed τ for both Saueressig and Cantrell. The varying lifetime slightly changes the absolute rates but **does not change the discriminant** (Cantrell still > Saueressig by ~24 pp).

---

## Fig 11: OSSE Recovery

**Q: In panel (a), where do "TRUE" and "TRUE KIE" come from?**

A: 
- **"TRUE" (black line)** = synthetic ground-truth emissions computed from the observed CH₄ mass balance, partitioned using He 2026 Science fractions: FF = 24%, Mic = 71%, BB = 5% (with BB decreasing from 6%→4% over 1999–2022).
- **"True KIE" (green)** = inversion using the same KIE that was used to generate the synthetic atmosphere (OH-¹³C = 1.0046, OH-D = 1.310). This is the "perfect knowledge" case — it shows how well the inversion could do if we knew the correct KIE.

**Q: Compared to (b), does filtering make FF emissions grow?**

A: No — filtering **removes iterations where δ¹³C and δD disagree**. What you see is that:
- In (a) (unfiltered), the Saueressig and Cantrell medians bracket the true line from different sides.
- In (b) (filtered), the surviving iterations are those where both isotopes gave consistent answers — which preferentially *keeps* iterations where the KIE is closer to being correct. For Cantrell (red), filtering doesn't change the level much because Cantrell was already closer to truth (1.0046 is closer to 1.0054 than to 1.0039). For Saueressig, filtering slightly shifts the median upward because it removes the most extreme under-estimates that failed the δD consistency check.

So it's not that filtering "grows" FF — it *selects* the more internally consistent subset of the Monte Carlo ensemble.

**Q: What is "recovery" in this context?**

A: Recovery = how well the inversion **recovers the known true emissions** from synthetic observations. We measure it with:
- **Bias** = mean(estimated − true) across all years
- **RMSE** = root-mean-square of (estimated − true) across all years

Better recovery = smaller bias and RMSE. The OSSE shows that agreement-filtering gives ~7% bias reduction and ~5% RMSE reduction vs unfiltered, but cannot fully eliminate the ~18 Tg/yr systematic bias from using the wrong KIE.

---

## Fig 12: Time-Varying OH KIE

**Q: What is "convergent" in this context?**

A: "Convergent" = a scenario where the Saueressig KIE **drifts from 1.0039 (1999) toward 1.0046 (2022)** — i.e. it's hypothetically "converging" toward the midpoint of the Saueressig–Cantrell range over the 23-year record. This is the aggressive test: if the *effective* KIE is becoming more Cantrell-like by 2022, does the discriminant that was built on *constant* KIE assumptions still work? Answer: yes, Δ = +18.7 pp, still significant.

The term "convergent" here does NOT mean the two trajectories meet each other — it means the low-endpoint (Saueressig) converges toward the midpoint. The Cantrell trajectory it's compared against remains constant.

---

## Fig 13: Fine Threshold Sweep

**Q: What is the "95% CI"?**

A: It's a **bootstrap 95% confidence interval on the agreement-rate difference Δ** (= rate_Cantrell − rate_Saueressig). Computed by:
1. Flatten the (years × iterations) agreement matrix into a 1D binary array (agree / not-agree)
2. Resample 2000 times with replacement
3. For each bootstrap sample, compute rate_Cantrell and rate_Saueressig independently
4. Take the 2.5th and 97.5th percentiles of the 2000 difference values

If the 95% CI is entirely above zero, the discriminant is statistically significant at that threshold.

**Q: How to choose between panel (b) discriminant power and panel (c) KSR when picking a threshold?**

A: They optimise for different things:
- **(b) Maximum discriminant power (τ = 90 Tg/yr):** maximises the *observational* separation between Cantrell and Saueressig. Use this if your goal is **diagnosing which KIE is correct** — i.e. the "KIE discriminant" use case. This is the primary scientific result.
- **(c) Maximum KSR (τ = 50 Tg/yr):** maximises how much the agreement filter *shrinks* the emission-trend sensitivity to the KIE choice. Use this if your goal is **reducing KIE uncertainty in the final budget** — i.e. reporting the tightest possible FF/Mic estimates.

In practice, 50–100 Tg/yr is the physically meaningful range. Lower thresholds give better KIE damping but retain fewer iterations (stricter filter). For publication, **report the discriminant at τ = 90 Tg/yr** (strongest signal) and **report the KSR at τ = 50 Tg/yr** (best sensitivity reduction). They serve complementary purposes.

---

## FigM6: Summary KSR across all methods

**Q: For the agreement filter, did you use the 1-box or 2-box model? Why?**

A: The agreement filter results in FigM6 (and in Phases 6, 6b, 6c, 7, 8) all use the **1-box (global mean) model**. Reasons:
1. **Current data loading:** Our `common.py` loads only the global-mean δD product (`GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx`). The 2-box scripts (`2x2_two.py`, `3x3_two.py`) construct hemispheric δD by applying a **fixed offset hack**: `dD_NH = dD_global − 6‰`, `dD_SH = dD_global + 6‰` (parameter `DD_IH_OFFSET = 6.0` in `common.py`). This is a rough approximation, not real hemispheric observations.
2. **Simplicity and interpretability:** The 1-box keeps the agreement test clean — one δ¹³C budget and one δD budget, no inter-hemispheric exchange corrections confounding the consistency check.
3. **Phase 4b showed the 2-box doesn't help for WLS:** KSR(FF) was 0.22 in the 2-box WLS, essentially the same as the 1-box WLS (0.20). Adding spatial resolution doesn't fix the fundamental problem, so the 1-box is sufficient to demonstrate the agreement filter works.

**However — hemispheric δD data DO exist** and can be constructed from Riddell-Young's own data package. The `dD_globmean.py` script in Riddell-Young's supplementary data (path: `ImportantReferences/Riddell-Young2025PNAS_DS/Riddell-Young_2025_dD_GlobMean/`) computes and saves:
- **Weekly NH and SH means** with 1000 MC iterations (`smoothedNH_matrix`, `smoothedSH_matrix`)
- **Semi-hemispheric annual means** for 4 latitude bands: PN (Polar North), TN (Tropical North), TS (Tropical South), PS (Polar South) — available in `data/dD_semihem_results.xlsx` (2005–2024, all data / no INSTAAR / no INSTAAR+BUDS variants)
- **Hemispheric CSV output** (`HemMean_dD_dei_DasguptaCal_noBUDS.csv`) with columns `Glob_smooth_mean`, `NH_smooth_mean`, `SH_smooth_mean`

The station-level data files (33 sites from ~82°N to ~90°S across NOAA/INSTAAR, MPI, Tohoku/NIPR, and IMAU networks) are all included. A proper 2-box agreement filter could use these real NH/SH δD products instead of the ±6‰ offset hack — this would be a meaningful extension.

**Note on the 2-box scripts (`2x2_two.py`, `3x3_two.py`):** These *do* run δD hemispheric inversions, but they use the offset hack (`dD_NH = global − 6`, `dD_SH = global + 6`), NOT real hemispheric observations. The ±6‰ offset was taken from the approximate NH–SH gradient reported in Riddell-Young et al. (2025), but it's a constant applied to a single global time series — it doesn't capture year-to-year hemispheric variability, which is the whole point of going to 2-box.

---

## Data Source for δD

**Q: What is the data source for your δD? From Riddell-Young or from `Old_files_before_organize/improved_dD_pipeline.py`?**

A: The current model uses **`rel/data/GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx`** loaded in `common.py` (line 412). This is:

- A **global-mean δD-CH₄** reconstruction with 1000 Monte Carlo iterations
- Calibrated on the **Umezawa (2012) scale**, excluding BUDS (Boulder Unified Data Set) measurements
- Constructed following the same methodology as the Riddell-Young (2025) δD reconstruction

The `Old_files_before_organize/improved_dD_pipeline.py` was an *earlier, improved version* of this reconstruction (adding more stations, updated calibration, NOAA 2024 updates). However, the current `common.py` loads the **original Riddell-Young-style Excel file**, NOT the output of `improved_dD_pipeline.py`. The improved pipeline generated `improved_dD_global_mean.csv` and `improved_dD_MC_iterations.csv`, but these are not yet integrated into the experiment.

**In short:** The δD source is effectively the **Riddell-Young (2025) product** in its original form (Umezawa-calibrated, no BUDS). The `improved_dD_pipeline.py` exists as an upgrade path but hasn't been plugged into the model yet.
