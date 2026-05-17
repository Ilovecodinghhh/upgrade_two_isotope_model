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
├── Output_2x2_{one,two,three}/ # Model outputs (config.json, results.csv, results.png)
├── Output_3x3_{one,two,three}/
├── rel/                       # Input data files (observations, source signatures, MC samples)
│   ├── data/                  # CSV/XLSX data: CH₄, δ¹³C, δD, CarbonTracker BB, source sigs
│   ├── build_*_sources.py     # Scripts to build regional source signature inputs
│   └── MASTER_DATA_INVENTORY.md
├── experiments/               # Focused analysis experiments
│   ├── dD_threshold/          # δD detection threshold analysis
│   ├── Hemispheric_Divergence/ # NH vs SH isotope divergence study
│   ├── KIE_sensitivity/       # KIE parameter sensitivity tests
│   └── KIE_immunity/          # KIE immunity analysis
├── sitesdata/                 # Site-level measurement data (δD, δ¹³C, CH₄ ppb)
│   ├── isotope_dD/            # 35 stations, Riddell-Young 2025
│   ├── isotope_d13C/          # 25 stations, NOAA/INSTAAR 2023
│   └── methane_ppb/           # 91 stations, NOAA GML 2025
├── ImportantReferences/       # Reference papers and their supplement datasets
├── Old_files_before_organize/ # Legacy scripts (kept for reference)
└── md_files/                  # Analysis notes and summaries
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

# Output goes to Output_{variant}/ as config.json, results.csv, results.png
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

## Code Conventions

- All isotope values in permil (‰); δ¹³C vs VPDB, δD vs VSMOW
- CH₄ in ppb; conversion factor PT = 2.815 ppb→Tg
- Monte Carlo: 1000 iterations standard for uncertainty propagation
- Results smoothed with 5-year running mean where noted
- Constants: C13Std = 0.011113, DStd = 0.00015576
