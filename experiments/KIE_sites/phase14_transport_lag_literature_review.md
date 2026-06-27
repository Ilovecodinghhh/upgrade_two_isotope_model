# Phase 14: NH Wetland Transport Lag and Attenuation Literature Review

## Purpose

This note summarizes literature relevant to the Phase 14 sensitivity experiment for northern-hemisphere wetland influence on southern-hemisphere background sites, especially `CGO` and `SPO`.

The immediate question is how to justify the two transport parameters used for the diagnostic NH wetland residual:

```text
NH_high + NH_mid scale = 0.05 / 0.12 / 0.20
NH_high + NH_mid lag   = +2.8 months
```

The central interpretation is:

```text
scale = effective annual phasor attenuation
lag   = effective annual phasor phase lag
```

These should not be described as direct emission fractions, mass transport fractions, or literal travel times.

## Main Conclusion

The current Phase 14 nominal setting is defensible as a first-order annual-harmonic approximation to interhemispheric transport.

Literature based on SF6 and related long-lived tracers commonly gives interhemispheric exchange or NH-midlatitude-to-SH transit times of order `1.3-1.5 yr`. If this transport is approximated as a first-order linear mixing process, the annual-cycle response has:

```text
amplitude attenuation ~0.1-0.12
phase lag             ~2.8 months
```

This matches the current nominal Phase 14 choice:

```text
nominal scale = 0.12
lag           = +2.8 months
```

Thus the current settings are best described as a sensitivity experiment anchored to observed interhemispheric tracer transport, not as a station-specific footprint model.

## Relevant Literature

| Study | Relevant result | Use for Phase 14 |
|---|---:|---|
| [Geller et al. 1997](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/97GL00523), SF6 latitudinal distribution | Derived interhemispheric exchange time about `1.3 +/- 0.1 yr`. | Supports a nominal exchange timescale near `1.3 yr`. |
| [Levin and Hesshaimer 1996](https://d-nb.info/1205371907/34), SF6 and Kr-85 | Surface-based two-box estimates are about `1.5-1.7 yr`; model flux-based estimate can be closer to `1.1 yr`. | Supports a plausible range rather than one exact transport time; also warns that surface gradients can bias exchange-time estimates high. |
| [Patra et al. 2011](https://acp.copernicus.org/articles/11/12813/2011/acp-11-12813-2011.pdf), TransCom-CH4 | Multi-model average interhemispheric exchange time about `1.39 yr`, consistent with observational estimates near `1.3-1.5 yr`. | Strong support for using `1.3-1.5 yr` as a literature-anchored transport scale. |
| [Holzer and Waugh 2015](https://sites.krieger.jhu.edu/waugh/files/2015/07/Holzer_Waugh_2015.pdf), transit-time distributions | Mean transit time from NH midlatitude surface increases from about `1.1 yr` at Samoa to about `1.4 yr` at South Pole; Tasmania/Cape Grim and South Pole are similar; TTD width is about `1.3 yr`. | Most directly relevant to `CGO` and `SPO`. Supports similar parameters for the two sites, but also shows that real transport is a broad distribution, not a single delay. |
| [Yang et al. 2019](https://elib.dlr.de/126249/2/Yang_et_al-2019-Geophysical_Research_Letters.pdf), SF6 exchange time versus SF6 age | Models can match interhemispheric exchange time while still having biased NH-midlatitude-to-SH age; transport from northern extratropics into the tropics is especially important. | Justifies keeping Phase 14 as a sensitivity envelope rather than claiming a precise transport correction. |
| [Schuck et al. 2024](https://acp.copernicus.org/articles/24/689/2024/), upper-tropospheric SF6 gradient | At the most southern latitudes, lag times relative to northern midlatitude surface exceed `1 yr`; typical NH-to-SH mean transit or age estimates span roughly `1.1-2.6 yr`, depending on metric and location. | Supports a low-impact case for older or more dispersed air, and emphasizes uncertainty in transport metrics. |
| [CSIRO Kennaook / Cape Grim greenhouse gas data](https://www.csiro.au/en/research/natural-environment/atmosphere/latest-greenhouse-gas-data) | Cape Grim methane seasonality is influenced by wetland release, OH destruction, and global transport from source regions; baseline air is well mixed Southern Hemisphere background. | Supports including transport as a relevant process for Cape Grim, while keeping the experiment background-site focused. |
| [East et al. 2024](https://repository.library.noaa.gov/view/noaa/68686/noaa_68686_DS1.pdf), methane seasonality | Southern Hemisphere methane seasonality is comparatively smooth and largely OH-driven; NH wetland magnitude, distribution, and timing are critical for NH methane seasonality. | Supports treating NH wetland residual at SH sites as a diagnostic sensitivity term, not the default dominant SH seasonal driver. |
| [Dowd et al. 2023](https://acp.copernicus.org/articles/23/7363/2023/), CH4 seasonal-cycle amplitude | Changes in observed methane seasonal-cycle amplitude can reflect non-local emissions and transport. | Supports the broader idea that remote source regions can affect seasonal-cycle diagnostics. |

## Why a 1.3-1.5 yr Transport Scale Gives a 2.8 Month Annual Phasor Lag

The Phase 14 `lag` is not a literal air-parcel transit time. It is the phase shift of the annual harmonic after transport and mixing.

For a simple first-order mixing response:

```text
H(omega) = 1 / (1 + i * omega * tau)
omega    = 2*pi yr^-1
```

The annual-cycle attenuation and phase lag are:

```text
attenuation = 1 / sqrt(1 + (omega * tau)^2)
lag_years   = atan(omega * tau) / omega
lag_months  = 12 * lag_years
```

Representative values are:

| tau, yr | Annual phasor attenuation | Annual phasor lag, months |
|---:|---:|---:|
| 1.1 | 0.14 | 2.73 |
| 1.3 | 0.12 | 2.77 |
| 1.4 | 0.11 | 2.78 |
| 1.5 | 0.11 | 2.80 |
| 2.0 | 0.08 | 2.85 |

This explains why a literature-supported transport timescale near `1.3-1.5 yr` maps naturally to:

```text
scale ~0.1-0.12
lag   ~2.8 months
```

It also explains why the phase lag saturates near `3 months` for very slow first-order mixing. The lag is a harmonic phase response, not the same thing as the mean transit time modulo one year.

## Interpretation of Current Phase 14 Scenarios

The current NH wetland sensitivity settings can be interpreted as:

```text
scale = 0.05
```

Low-impact case. Represents stronger dilution, broader transit-time spreading, older effective air, or weaker station sensitivity to NH wetland annual phasors.

```text
scale = 0.12
```

Nominal case. Consistent with a first-order annual response for an interhemispheric exchange or transit scale of about `1.3-1.5 yr`.

```text
scale = 0.20
```

High-impact case. Represents faster pathways, less phase dispersion, or a stronger-than-nominal NH annual residual at `CGO` and `SPO`.

The common lag:

```text
lag = +2.8 months
```

is consistent with the annual-harmonic response of the nominal interhemispheric transport timescale.

## Methodological Cautions

1. The Li 2026 wetland product is a local gridded emission field. It does not include atmospheric transport, station footprints, or lagged station influence.

2. The Phase 14 NH contribution is an imposed sensitivity term:

```text
Z_wetland_SH_site =
  Z_SH_extra
+ scaled_lagged(Z_NH_high)
+ scaled_lagged(Z_NH_mid)
```

It is not a CTM-derived source-receptor relationship.

3. `scale` should not be described as a transported mass fraction. It is an effective annual phasor attenuation, combining dilution, transport dispersion, sampling geometry, and station representativeness.

4. A single lag is an approximation. Holzer and Waugh show that NH-to-SH transport has a broad transit-time distribution. A true treatment would convolve the wetland emission time series with a site-specific TTD or use a CTM footprint.

5. Applying the same lag and scale to `NH_high` and `NH_mid` is a simplifying assumption. Real transport from high-latitude and midlatitude wetland regions may differ.

6. Applying the same NH transport parameters to `CGO` and `SPO` is reasonable for a first sensitivity experiment because literature suggests weak differences between Tasmania/Cape Grim and South Pole transit times, but this should not be over-interpreted.

## Recommended Text for Manuscript or Notes

Suggested concise wording:

```text
For the SH wetland-source sensitivity test, delayed NH wetland influence was represented as an effective annual-harmonic residual rather than a footprint-resolved transport calculation. Based on SF6-constrained interhemispheric exchange and NH-midlatitude-to-SH transit times of order 1.3-1.5 yr, a first-order annual response implies an amplitude attenuation of about 0.1-0.12 and a phase lag of about 2.8 months. We therefore used +2.8 months as the nominal NH lag and explored effective phasor attenuation factors of 0.05, 0.12, and 0.20 as low, nominal, and high structural envelopes.
```

Suggested wording to avoid:

```text
12% of NH wetland emissions are transported to CGO/SPO.
```

Better wording:

```text
The nominal case allows a 0.12 effective annual phasor residual from NH wetland seasonality at CGO/SPO.
```

## Recommended Next Step

The current Phase 14 setup is acceptable as a diagnostic envelope. If we want to strengthen it without building a full CTM workflow, the next useful extension would be a small two-parameter sensitivity grid:

```text
scale = 0.05, 0.10, 0.12, 0.15, 0.20
lag   = 2.5, 2.8, 3.0 months
```

This would separate sensitivity to annual phasor amplitude from sensitivity to the assumed transport phase shift.

For a stronger physical treatment, the right next step would be to replace the single `scale + lag` approximation with a site-specific transit-time distribution or CTM-derived source-receptor footprint for `CGO` and `SPO`.
