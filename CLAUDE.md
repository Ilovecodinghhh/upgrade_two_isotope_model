# CLAUDE.md

## Project Overview

Methane dual-isotope (δ¹³C + δD) box model suite for source partitioning of global CH₄ emissions into fossil fuel (FF), microbial (Mic), and biomass burning (BB) components. The models solve isotope mass balance equations using atmospheric observations to infer annual source contributions from 1988–2024.

## Repository Structure

```
├── common.py                  # Shared utilities (isotope math, data I/O, KIE, lifetime)
├── models/
│   ├── core.py                # Shared model infrastructure (solver, plotting, quality monitor)
│   └── inputs.py              # Input choice catalog (KIE, signatures, lifetime options)
├── 2x2_one.py                 # BB-fixed, separate δ¹³C & δD, global 1-box
├── 2x2_two.py                 # BB-fixed, separate δ¹³C & δD, NH/SH 2-box
├── 2x2_three.py               # BB-fixed, separate δ¹³C & δD, NHext/Trop/SHext 3-box
├── 3x3_one.py                 # Simultaneous δ¹³C+δD, global 1-box (3×3 solve)
├── 3x3_two.py                 # Simultaneous δ¹³C+δD, NH/SH 2-box
├── 3x3_three.py               # Simultaneous δ¹³C+δD, NHext/Trop/SHext 3-box
├── rel/                       # Input data files (observations, source signatures, MC samples)
│   ├── data/                  # CSV/XLSX data: CH₄, δ¹³C, δD, CarbonTracker BB, source sigs
│   ├── build_*_sources.py     # Scripts to build regional source signature inputs
│   └── MASTER_DATA_INVENTORY.md
├── experiments/               # Focused analysis experiments
│   ├── dD_threshold/          # δD detection threshold analysis
│   ├── Hemispheric_Divergence/ # NH vs SH isotope divergence study
│   ├── KIE_immunity/          # KIE immunity analysis (dual-isotope manuscript + review)
│   └── KIE_sensitivity/       # KIE parameter sensitivity tests (multi-phase analysis pipeline)
├── sitesdata/                 # Site-level measurement data (δD, δ¹³C, CH₄ ppb)
│   ├── isotope_dD/            # 35 stations, Riddell-Young 2025
│   ├── isotope_d13C/          # 25 stations, NOAA/INSTAAR 2023
│   ├── methane_ppb/           # 91 stations, NOAA GML 2025
│   └── DATA_COLLECTION_REPORT.md
├── ImportantReferences/       # Reference papers and their supplement datasets
├── Old_files_before_organize/ # Legacy scripts and notes (reference only, paths may be stale)
├── dD_Gap_Fixed.md            # Documentation: hemispheric δD data gap fix (2020–2023)
├── KIE_Used_Previous_Study.md # Survey of OH-¹³C KIE values across 16 publications
├── Primary_Annual_Database.md # Data sources, processing, citations for the model inputs
└── .gitignore                 # Ignores __pycache__/, Output/, .DS_Store, Excel temp files
```

## Model Variant Naming Convention

**`{system}_{boxes}.py`** where:
- **System**: `2x2` = BB fixed from CarbonTracker, solve FF+Mic per isotope independently; `3x3` = all three sources free, solve simultaneously using both isotopes
- **Boxes**: `one` = global 1-box; `two` = NH/SH 2-box; `three` = NHext/Trop/SHext 3-box

## Key Scientific Concepts

- **Source partitioning**: Total CH₄ = FF + Mic + BB, each with distinct isotopic signatures
- **δ¹³C** (‰ VPDB): Separates thermogenic/fossil (~−44‰) from microbial (~−62‰) and BB (~−22‰)
- **δD** (‰ VSMOW): Separates fossil (~−180‰) from microbial (~−310‰) and BB (~−215‰)
- **KIE** (Kinetic Isotope Effect): Fractionation by OH, Cl, soil, stratospheric sinks
- **Isotope mass balance**: Source δ derived from atmospheric δ + sink KIE + lifetime

## Running Models

```bash
# Run a single model variant with default config
python3 2x2_one.py

# Run with custom config via command-line args (see --help)
python3 3x3_one.py --kie sampled --lifetime time_varying

# Output goes to Output/ (gitignored) as config.json, results.csv, results.png
```

## Dependencies

- Python 3.8+
- numpy, pandas, matplotlib, scipy (lsq_linear)
- openpyxl (for .xlsx reading)
- No environment file exists — install manually: `pip install numpy pandas matplotlib scipy openpyxl`

## Data

- **Atmospheric observations**: `rel/data/` — CH₄ global mean, δ¹³C NH/SH means, δD global MC iterations
- **Source signatures**: `rel/data/` — Monte Carlo CSV files for FF/Mic/BB × δ¹³C/δD × region
- **Station-level data**: `sitesdata/` — raw flask measurements from NOAA, INSTAAR, Riddell-Young
- **Reference datasets**: `ImportantReferences/` — published supplement data (Basu2022 ObsPack, Riddell-Young2025, etc.)

## Branches

- `master` — main development branch with model code
- `sites` — site-level data collection and organization
- `three-box` — three-box model development

## Experiments

Each experiment under `experiments/` is a self-contained analysis with its own `figures/`, `results/`, `analysis/` subdirectories and manuscript/results markdown files.

### 1. `dD_threshold/` — When Does δD Help? (Target: PNAS)

**Question:** At what microbial δD source-signature uncertainty does adding δD improve (or degrade) methane source attribution vs δ¹³C-only?

**Core result:** A sharp threshold exists at **σ(Mic δD) ≈ 37‰**. Below it, δD halves FF emission uncertainty (53% CI reduction); above it, δD actively degrades the solution. Current observational precision (~8‰) is well within the threshold — a 4.5× safety margin.

**Key findings:**
- Resolves the Riddell-Young (2025) vs Thanwerdas (2024) contradiction: Thanwerdas used σ≈128‰ (3.5× above threshold), inevitably finding δD useless
- δD's value is hemispheric: one-box dual-isotope *fails* (CI widens), two-box *succeeds* (53% improvement) because δD source signatures have 5–10× larger NH-SH gradients than δ¹³C
- Robust across 9 KIE×lifetime configs, 5 data versions (v1–v5), all year ranges
- Full draft manuscript in `draft.md`; shared model runner in `analysis/core.py` (all phases import from it)

**Status:** v5 complete (Luo 2024 C4 map). 6 analysis phases + comprehensive figures done. Code refactored into shared `core.py`.

### 2. `Hemispheric_Divergence/` — NH vs SH Microbial Trends (Target: ACP)

**Question:** Does hemispheric resolution reveal spatial structure in CH₄ growth that one-box models miss?

**Core result:** Post-2006 CH₄ growth is driven by **asymmetric hemispheric microbial trends**: NH microbial emissions increase at +6.6 Tg/yr² (100% of MC positive), while SH microbial emissions are stable (−1.1 Tg/yr²). This rules out globally symmetric mechanisms (e.g., uniform OH decline).

**Key findings:**
- NH drives >73% of the dual-isotope improvement; SH contributes a steady ~15%
- FF absolute levels corrected from ~50 → ~115 Tg/yr (now matches EDGAR) after Phase A fixes: observed IH gradient + uncertainty-based W matrix
- δD contributes 33% to cost function with proper weighting (was ~2% with old W)
- Global FF trend: −2.49 Tg/yr² (declining), driven by NH; contradicts some 3D inversions
- Validated against EDGAR: FF(2010) = 115 Tg/yr (EDGAR: 110), NH FF share = 72% (EDGAR: 72%)

**Status:** v4 results complete. Manuscript draft v3 in `MANUSCRIPT.md`. Critical assessment identifies remaining limitations (prescribed BB, noisy FF temporal CV, 3D inversion discrepancy).

### 3. `KIE_immunity/` — Dual-Isotope KIE Sensitivity & Variance Decomposition (Target: ACP)

**Question:** Does adding δD reduce sensitivity to the OH-¹³C KIE controversy (Saueressig 1.0039 vs Cantrell 1.0054)?

**Core result:** δD reduces FF trend uncertainty by **38%** (σ: 31→19.2 Tg/yr) and KIE-driven spread by **34%** (13.0→8.6 Tg/yr). But the KIE still determines the **sign** of the post-2007 FF trend (Saueressig: +2.9; Cantrell: −5.6 Tg/yr) — it is an irreducible uncertainty floor.

**Key findings:**
- Variance decomposition: Source signatures dominate (48%), KIE = 25%, lifetime < 1%, residual = 27%
- Earlier "KIE immunity" finding was an artifact of global-mean source signatures; hemispheric heterogeneity (v3→v4) restored KIE sensitivity
- Extensive structural tests: W matrix sensitivity, BB perturbation, MC convergence, solver diagnostics, EDGAR validation
- 18 analysis phases (phase5–phase18) with JSON result files as single source of truth
- Revision response addresses all 28 reviewer comments

**Status:** Manuscript revision 1 complete in `MANUSCRIPT_DUAL_ISOTOPE.md`. All numbers verified against JSON ground truth. `analysis/run_all.py` reproduces everything.

### 4. `KIE_sensitivity/` — Agreement Filter as KIE Discriminant (Multi-Phase Pipeline)

**Question:** Can the δ¹³C-δD agreement rate discriminate between the competing OH-¹³C KIE values?

**Core result:** WLS coupling of δ¹³C+δD makes KIE sensitivity **5× worse** (KSR ≈ 0.2). But solving the two isotopes *independently* and filtering by consistency reveals a **35.5 pp agreement-rate discriminant** — Cantrell's KIE produces more internally consistent δ¹³C/δD solutions than Saueressig's.

**Key findings:**
- Phases 1–5: WLS coupling always degrades KIE sensitivity (no optimal weight exists; even w_dD=0.01 amplifies spread 4×)
- Root cause: shifted δ¹³C row contradicts unshifted δD row → WLS distributes the conflict across both unknowns
- Phase 6–6c: Agreement filter framework — solve isotopes independently, keep iterations where FF estimates agree within threshold
- Phase 7: Discriminant survives time-varying KIE trajectories
- Phase 8: Stable across three independent 8-year epochs (1999–2006, 2007–2014, 2015–2022)
- Phase 9 correction: KSR=2.5–3.2 at N=1000 was inflated; at N=5000, KSR stabilizes at 1.12 — filter's value is as a discriminant, not sensitivity reducer
- `make_manuscript_figs.py` generates all publication figures

**Status:** Phases 1–9 complete. RESULTS.md has full per-phase tables. All phase scripts in top-level directory; results in `results/` subdirectories.

### Cross-Experiment Relationships

The four experiments share the same 2-box model infrastructure and data but answer distinct questions:
- **dD_threshold** asks: *under what uncertainty does δD help?* (answer: σ < 37‰)
- **Hemispheric_Divergence** asks: *what does 2-box reveal that 1-box misses?* (answer: NH microbial drives growth)
- **KIE_immunity** asks: *does δD reduce KIE sensitivity?* (answer: partially — 34% reduction, but KIE still determines FF trend sign)
- **KIE_sensitivity** asks: *can δ¹³C-δD agreement discriminate between KIE values?* (answer: yes — Cantrell is more internally consistent)

Common data: `rel/data/` MC CSVs, `common.py` / `models/core.py` model engine. Each experiment has its own `analysis/core.py` for experiment-specific shared code.

## Code Conventions

- All isotope values in permil (‰); δ¹³C vs VPDB, δD vs VSMOW
- CH₄ in ppb; conversion factor PT = 2.815 ppb→Tg
- Monte Carlo: 1000 iterations standard for uncertainty propagation
- Results smoothed with 5-year running mean where noted
- Constants: C13Std = 0.011113, DStd = 0.00015576

## Notes on Old_files_before_organize/

Contains legacy versioned scripts (v1.0 through v4.0) and past analysis notes. Paths inside may be stale. Consult only when you need historical context on earlier model iterations. The README.txt inside has guidance.
