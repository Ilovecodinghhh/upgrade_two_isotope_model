# KIE_sites — Figure Explanations

Brief descriptions of each figure produced by the KIE_sites experiment.

---

## Phase 1: Data Extraction & Pairing

### `fig1_timeseries_grid.png`
Monthly-mean time series of **δ¹³C-CH₄** (left column, blue) and **δD-CH₄** (right column, orange) for all 12 co-located sites, sorted from north (ALT, 82°N) to south (SPO, 90°S). Each row is one site. This is the raw data overview — it shows the available record length, the long-term trend (δ¹³C becoming more negative over time), and hints of seasonal oscillation at each station.

### `fig1_data_coverage.png`
Heatmap of **paired data availability** by month and site. Green cells indicate months where both δ¹³C and δD observations exist simultaneously. Sites are sorted by latitude (N→S). Most paired coverage comes from the INSTAAR period (~2005–2010); some sites (ALT, ZEP, MLO) have extended records from other labs. Shows that SPO and high-latitude sites have the densest coverage, while tropical sites (ASC, SMO) are sparser.

---

## Phase 2: Seasonal Harmonic Fitting

### `fig2_seasonal_cycles.png`
**Folded seasonal cycles** for all 12 sites (rows) across three quantities (columns): δ¹³C (blue), δD (red), and CH₄ concentration in ppb (green). Monthly anomalies (trend removed) are plotted against month-of-year, with the best-fit annual harmonic overlaid as a smooth curve. This is the core data for extracting seasonal amplitudes and phases. NH sites show clear sinusoidal patterns; SH and tropical sites have smaller amplitudes.

### `fig2_harmonic_summary.png`
Two-panel summary of the harmonic fit results across all sites. **(a)** Seasonal amplitude ratio R = A(δ¹³C)/A(δD) vs latitude, with the pure-OH prediction band (Saueressig–Cantrell range, ~0.013–0.017) shown in green. **(b)** Phase difference (δ¹³C peak minus δD peak, in months) vs latitude, with a ±1 month "good alignment" band in grey. Key takeaway: most NH sites have ratios far above the pure-OH band, signalling source contamination. SMO is a phase-offset outlier.

---

## Phase 3: Cross-site Synthesis

### `fig3_ratio_vs_latitude.png`  *(same as fig2_harmonic_summary panel a, but standalone)*
Duplicate/companion of the harmonic summary left panel — amplitude ratio vs latitude with error bars and site labels. Included here for the Phase 3 standalone context.

### `fig3_site_classification.png`
Amplitude ratio vs latitude with sites classified as **clean** (blue circles, |Δφ| < 2 months and sufficient data) or **excluded** (grey triangles: SMO, AZR, MLO, ASC). The red horizontal line and pink band show the **weighted mean ± 1σ of clean sites** (R = 0.0397 ± 0.0066). The green band is the pure-OH prediction. The clean-site mean is 3–4σ above the pure-OH band, confirming that source seasonality inflates the observed ratio.

### `fig3_phase_diagnostic.png`
Phase difference (δ¹³C − δD peak timing, in months) vs latitude for all 12 sites. Blue dots are sites within the ±2-month threshold (green band); grey dots are flagged outliers. If OH alone drove both isotope cycles, the phases would align perfectly (Δφ ≈ 0). Most sites cluster near zero, but **SMO** (−5 months) is a clear outlier, and **MLO** and **ASC** also show offsets suggesting source interference.

---

## Phase 4: Source Deconvolution

### `fig4_decomposition.png`
Stacked bar charts decomposing the observed seasonal amplitude into **sink-driven** (blue, OH fractionation) and **source-driven** (orange, microbial seasonality) components for **(a)** δ¹³C and **(b)** δD at each site. NH high-latitude sites (ALT, ZEP, BRW) have the largest total amplitudes with a roughly 40–60% source contribution. SH sites (CGO, SPO) are almost entirely sink-driven. This explains why the raw amplitude ratio at NH sites exceeds the pure-OH prediction.

### `fig4_sink_ratio.png`
Two-panel figure. **(a)** The implied α\_13C\_OH if all seasonal variation were attributed to OH sink alone, plotted vs latitude. SH sites (CGO ≈ 1.005, SPO ≈ 1.003) fall near the lab values; NH sites imply absurdly high α values (1.01–1.04), confirming source contamination. Dashed lines mark the Saueressig (1.0039) and Cantrell (1.0054) lab values. **(b)** Bar chart comparing the **observed** amplitude ratio (blue) to the **sink-only** predicted ratio under Saueressig (green) and Cantrell (orange) assumptions. At every site, observed > predicted, consistent with an additive source component inflating δ¹³C more than δD.

---

## Phase 5: KIE Extraction

### `fig5_kie_constraint.png`
The final result figure, three panels. **(a)** Probability density of α\_13C\_OH from the **SH direct constraint** (blue, using CGO+SPO where sources are minimal): point estimate 1.0034, 95% CI [0.9975, 1.0092]. **(b)** PDF from the **source-corrected** approach using all clean sites after deconvolution (orange): point estimate 1.0162, 95% CI [1.0026, 1.0312]. **(c)** Summary comparison — both estimates shown as box-and-whisker alongside the Saueressig (green dashed, 1.0039) and Cantrell (orange dashed, 1.0054) lab values. The SH direct estimate favours Saueressig but cannot definitively exclude Cantrell; both lab values fall within the 95% CI.

---

## Phase 6: Phasor Source Correction

### `fig8_phasor_decomposition.png`
Polar "phasor clock" diagrams for BRW, CBA, CGO, and SPO. Angle is peak month and radius is seasonal amplitude, so the figure shows both timing and magnitude in the same visual grammar as Fig 11. Green arrows are observed phasors, red dashed arrows are wetland source phasors, and blue arrows are the corrected sink phasors. Rows separate d13C and dD. Each isotope row uses a common radius scale across sites to preserve the NH-vs-SH magnitude contrast; CGO and SPO include zoom insets because their vectors are much smaller than the NH corrections.

---

## Phase 7: Individual-Year Stability

### `fig12_yearly_stability.png`
Site-by-site yearly amplitude ratios from separate annual harmonic fits for 2005-2010, requiring at least 8 paired months per year. The dashed line is the full-period Phase 2 ratio and the grey band is its 95% CI. Point color gives paired-month coverage. Red panel titles flag sites where year-to-year scatter exceeds the full-period Phase 2 uncertainty, highlighting where the 5-year overlap may understate harmonic-fit uncertainty.
