# Hemispheric δD-CH₄ Data Gap Fix (2020–2023)

## Problem

The hemispheric (NH/SH) δD-CH₄ annual means from `dD_hemispheric_MC.py` were **NaN for 2020–2023**, even though multiple stations had smoothed δD data extending well past 2020.

## Root Cause

In `dD_hemispheric_MC.py`, the weekly date grid was built from the **longest station file**:

```python
max_len = max(s['data'].shape[0] for s in stations)
ref_station = [s for s in stations if s['data'].shape[0] == max_len][0]
date_full = ref_station['data'][:, 0]
```

The longest station was **gvnIMAU** (Neumayer Station, Antarctica), which spans 1988–2020.0 (1640 weekly points). Although gvnIMAU had the most rows, it **ended at exactly 2020.0**.

When other stations (gvnMPI, namMPI, syoNIPR, smoIMAU, cgoIMAU, sisMPI, etc.) were aligned to this grid, any data past 2020.0 was **silently clipped** because the reference grid had no indices for those dates:

```python
offset = np.searchsorted(date_full, s_dates[0])
end = min(offset + slen, n_weeks)  # clipped to n_weeks = 1640
```

For gvnMPI (2018–2024, 313 points), only 104 points fell within the grid (2018–2020); the remaining 209 points covering 2020–2024 were discarded.

## Station Coverage

### SH MBL stations with data past 2020

| Station | Network | Latitude | Time Range | Flag |
|---------|---------|----------|------------|------|
| gvnIMAU | IMAU | −70.67° | 1988.6–2020.0 | MBL |
| gvnMPI | MPI | −70.67° | 2018.0–2024.0 | MBL |
| namMPI | MPI | −23.56° | 2013.2–2022.4 | MBL |
| syoNIPR | NIPR | −69.01° | 1995.2–2023.0 | MBL |
| smoIMAU | IMAU | −14.25° | 2022.1–2024.2 | MBL |
| cgoIMAU | IMAU | −40.68° | 2023.1–2024.2 | MBL |

### SH stations contributing per year (with mean δD)

| Year | N stations | Stations |
|------|-----------|----------|
| 2018 | 4 | gvnIMAU(−72.9), gvnMPI(−72.1), namMPI(−71.0), syoNIPR(−71.5) |
| 2019 | 4 | gvnIMAU(−73.8), gvnMPI(−71.4), namMPI(−73.7), syoNIPR(−72.0) |
| 2020 | 3 | gvnMPI(−71.5), namMPI(−73.5), syoNIPR(−72.5) |
| 2021 | 3 | gvnMPI(−72.4), namMPI(−74.8), syoNIPR(−73.6) |
| 2022 | 4 | gvnMPI(−73.5), namMPI(−75.1), smoIMAU(−78.8), syoNIPR(−76.4) |
| 2023 | 3 | cgoIMAU(−78.6), gvnMPI(−75.0), smoIMAU(−80.0) |

Note: NOAA/INSTAAR δD-CH₄ flask data ended ~2009-2010. Post-2010 coverage relies entirely on MPI (Max Planck), IMAU (Utrecht), and NIPR (Japan) networks.

## Fix

Replaced the single-station reference grid with a **unified date grid** spanning all stations:

```python
all_t0 = min(s['data'][0, 0] for s in stations)
all_t1 = max(s['data'][-1, 0] for s in stations)
dt = 7.0 / 365.25  # weekly spacing
date_full = np.arange(all_t0, all_t1 + dt/2, dt)
```

Also updated station-to-grid alignment to use `searchsorted` with `side='left'` and `np.clip` for robust nearest-point mapping.

## Results

| Year | NH δD (‰) | SH δD (‰) | NH–SH (‰) |
|------|----------|----------|-----------|
| 2018 | −87.8 | −72.2 | −15.6 |
| 2019 | −88.0 | −72.5 | −15.5 |
| **2020** | **−89.2** | **−72.1** | **−17.1** |
| **2021** | **−90.4** | **−73.9** | **−16.5** |
| **2022** | **−90.3** | **−80.2** | **−10.2** |
| **2023** | **−90.0** | **−78.0** | **−12.0** |

## Caveats

1. **SH 2022–2023 sparse coverage:** Only 3–4 stations. The sharp SH shift (−73.9 → −80.2‰ in 2022) coincides with smoIMAU and cgoIMAU entering the network, which report substantially more negative values (−78 to −80‰) than the established polar stations (−73 to −75‰).

2. **Station composition change:** The 2022–2023 SH mean is qualitatively different from 2018–2021 because the station mix changed (lower-latitude stations smoIMAU at −14° and cgoIMAU at −41° replaced higher-latitude syoNIPR at −69° in 2023). This is real data but should be interpreted cautiously.

3. **Global δD was unaffected** — the `dD_globmean.py` script uses a different 4-block gap-filling approach that already covered 2020–2023.

## Commit

`cd64c7a` — "Fix hemispheric dD gap 2020-2023: extend date grid to cover all stations"
