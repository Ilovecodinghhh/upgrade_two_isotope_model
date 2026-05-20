# δD-CH₄ Source Signature Comparison: OIPC Regression vs Douglas Zonal Mean

## Purpose

For the KIE_sites seasonal-cycle analysis, each site needs a δ²H-CH₄ source
signature to correct for microbial source seasonality (phases 4–5). Two
independent approaches are compared here:

| Method | What it represents | Strengths | Weaknesses |
|--------|--------------------|-----------|------------|
| **OIPC + Douglas regression** | Predicted δD of wetland CH₄ **emitted at the station coordinates** | Site-specific; captures local climate (maritime, altitude, Gulf Stream) | Assumes the station itself is a wetland; extrapolates at extreme latitudes |
| **Douglas (2021) zonal mean** | Emission-weighted mean δD of wetlands **within a latitude band** | Empirical average over real wetland measurements; physically grounded | Coarse (3 bands); ignores within-band climate variability |

### The fundamental question

These are atmospheric monitoring stations, not wetland sites. The δD source
signature that matters for the seasonal cycle is **the δD of the microbial
CH₄ that actually reaches the station** — not what a hypothetical local
wetland would emit. At stations near extensive wetlands (BRW, CBA), the two
coincide. At remote ocean stations (ASC, SMO, SPO), the source CH₄ has
travelled thousands of kilometres and its δD reflects the **upstream emission
region**, not local precipitation.

---

## Regression details

**Primary regression** (Douglas 2021 Table S2, wetlands, growing-season δ²Hp):

    δ²H-CH₄ = 0.705 (±0.074) × δ²Hp_gs − 284.5 (±6.1)
    R² = 0.633,  RMSE = 22.4 ‰,  n = 55 wetland sites

**Douglas (2021) Table 1 zonal means** (wetlands only):

| Band | δ²H-CH₄ (‰) | ±1σ | Flux (Tg yr⁻¹) |
|------|-------------|-----|-----------------|
| <30° (tropics) | −301 | 15 | 115 |
| 30–60°N | −324 | 14 | 25 |
| >60°N | −374 | 10 | 9 |

---

## Site-by-site comparison and recommendation

### ALT — Alert, Nunavut (82.5°N)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −431 ± 28 ‰ | −374 ± 10 ‰ |
| **Basis** | Growing-season δ²Hp = −207 ‰ at 82.5°N | >60°N wetland mean |

**→ Recommend: Douglas zonal (−374 ± 10 ‰)**

Alert is the world's northernmost permanently inhabited settlement. There are
no significant wetlands at 82°N — the nearest major wetland sources are the
Hudson Bay Lowlands (~55°N) and Siberian/Scandinavian boreal zones (55–70°N).
The OIPC regression is extrapolating well beyond its training data (Douglas
2021 has almost no sites above 70°N). The −431 ‰ prediction is driven by
Alert's extremely depleted precipitation, but this precipitation does not feed
any local wetlands. The zonal mean (−374 ‰) represents the emission-weighted
average of real >60°N wetlands that produce the CH₄ reaching Alert.

### ZEP — Ny-Ålesund, Svalbard (78.9°N)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −344 ± 24 ‰ | −374 ± 10 ‰ |
| **Basis** | Growing-season δ²Hp = −84 ‰ at 78.9°N | >60°N wetland mean |

**→ Recommend: Douglas zonal (−374 ± 10 ‰)**

Svalbard precipitation is anomalously warm for 79°N due to the Gulf Stream /
North Atlantic Drift (δ²Hp only −84 ‰ in summer, comparable to 50°N
continental sites). This makes the OIPC prediction look "too enriched" for
a high-Arctic station. But the key point is that Svalbard has negligible local
wetland CH₄ emissions. The seasonal CH₄ signal at ZEP comes from Eurasian
boreal wetlands (55–70°N), whose δD is well represented by the >60°N zonal
mean. The OIPC value, while technically a valid regression output, describes
a hypothetical Svalbard wetland that does not exist at meaningful scale.

### BRW — Barrow/Utqiaġvik, Alaska (71.3°N)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −384 ± 25 ‰ | −374 ± 10 ‰ |
| **Basis** | Growing-season δ²Hp = −141 ‰ at 71.3°N | >60°N wetland mean |

**→ Recommend: OIPC regression (−384 ± 25 ‰), or either — they agree.**

BRW is the one high-latitude site where OIPC adds genuine value. It sits at
the edge of the North Slope tundra, one of the most productive Arctic wetland
regions. The local precipitation δD (−141 ‰ growing season) is a plausible
proxy for the water feeding nearby wetlands. The OIPC prediction (−384 ‰) is
consistent with the zonal mean (−374 ± 10) and with Chanton et al. (2006)
direct measurements at Alaskan wetlands (−308 to −394 ‰). Using the OIPC
value captures the specific North Slope climate. The larger uncertainty
(±25 ‰) appropriately reflects regression scatter.

### CBA — Cold Bay, Alaska (55.2°N)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −326 ± 24 ‰ | −324 ± 14 ‰ |
| **Basis** | Growing-season δ²Hp = −58 ‰ at 55.2°N | 30–60°N wetland mean |

**→ Recommend: OIPC regression (−326 ± 24 ‰), or either — excellent agreement.**

Cold Bay is on the Alaska Peninsula, downstream of extensive boreal wetlands.
The OIPC prediction matches the zonal mean almost exactly. The site-specific
value is marginally preferred because it captures the actual subarctic maritime
climate at 55°N.

### MHD — Mace Head, Ireland (53.3°N)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −306 ± 23 ‰ | −324 ± 14 ‰ |
| **Basis** | Growing-season δ²Hp = −30 ‰ at 53.3°N | 30–60°N wetland mean |

**→ Recommend: Douglas zonal (−324 ± 14 ‰)**

Mace Head is a clean Atlantic background station on Ireland's west coast.
Its local precipitation δD (−30 ‰ growing season) is very enriched for 53°N
due to the maritime climate. However, MHD has minimal local wetland CH₄
sources. The seasonal CH₄ signal arriving at MHD originates from boreal
wetlands across the broader 30–60°N band (Scandinavia, Russia, Canada). The
zonal mean better represents this spatially integrated source. The OIPC value
would be appropriate if there were productive wetlands in coastal western
Ireland feeding CH₄ directly to the station, but the dominant signal is
hemispheric/continental.

### AZR — Terceira, Azores (38.8°N)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −284 ± 23 ‰ | −324 ± 14 ‰ |
| **Basis** | Growing-season δ²Hp = +1 ‰ at 38.8°N | 30–60°N wetland mean |

**→ Recommend: Neither — site is excluded (non-MBL). If used: zonal (−324 ± 14 ‰).**

The Azores are subtropical mid-Atlantic with no significant local CH₄ sources.
The OIPC prediction (−284 ‰) reflects the warm, isotopically enriched
subtropical precipitation, but no wetlands exist here. Any source seasonality
at AZR reflects long-range transport from NH continental wetlands (30–60°N),
so the zonal mean is more physically appropriate. Note AZR is already excluded
from the KIE analysis (non-MBL station).

### MLO — Mauna Loa, Hawaii (19.5°N, 3397 m)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −320 ± 24 ‰ | −301 ± 15 ‰ |
| **Basis** | Growing-season δ²Hp = −50 ‰ at 3397 m | <30°N wetland mean |

**→ Recommend: Neither — site is excluded (non-MBL). If used: zonal (−301 ± 15 ‰).**

MLO is at 3397 m elevation, making its precipitation δD much more depleted
than sea-level tropical stations. The OIPC regression (−320 ‰) is an artefact
of the altitude effect on precipitation δD — there are no wetlands at 3400 m
on a volcanic summit. MLO samples the free troposphere and is already excluded
(non-MBL). The <30°N zonal mean (−301 ‰) is appropriate for the tropical
microbial sources that dominate the CH₄ reaching this latitude.

### KUM — Cape Kumukahi, Hawaii (19.6°N, 3 m)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −287 ± 23 ‰ | −301 ± 15 ‰ |
| **Basis** | Growing-season δ²Hp = −3 ‰ at 3 m | <30°N wetland mean |

**→ Recommend: Douglas zonal (−301 ± 15 ‰)**

KUM is at sea level on Hawaii's coast, with tropical oceanic precipitation
(δ²Hp ≈ −3 ‰). The OIPC regression (−287 ‰) describes what a hypothetical
Hawaiian wetland would emit — more enriched than the <30°N mean. But there are
no significant tropical wetlands on Hawaii. The CH₄ reaching KUM comes from
Asian rice paddies, SE Asian wetlands, and other tropical sources. The <30°N
zonal mean (−301 ‰) better represents this mixture. The 14 ‰ difference
(−287 vs −301) is within uncertainty of either estimate.

### ASC — Ascension Island (8.0°S)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −277 ± 23 ‰ | −301 ± 15 ‰ |
| **Basis** | Growing-season δ²Hp = +11 ‰ (year-round) | <30°N wetland mean |

**→ Recommend: Douglas zonal (−301 ± 15 ‰)**

Ascension is a tiny volcanic island in the equatorial Atlantic. Precipitation
δD is positive (marine tropical, enriched by warm evaporation). The OIPC
prediction (−277 ‰) describes a hypothetical equatorial marine-island wetland.
In reality, ASC has essentially zero local CH₄ sources. Its atmospheric CH₄
signal reflects long-range transport from African and South American tropical
wetlands, whose emission-weighted δD is well approximated by the <30°N zonal
mean (−301 ‰). ASC is already excluded from the clean analysis (low δ¹³C
amplitude).

### SMO — Tutuila, American Samoa (14.2°S)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −299 ± 23 ‰ | −301 ± 15 ‰ |
| **Basis** | Growing-season δ²Hp = −21 ‰ (year-round) | <30°N wetland mean |

**→ Recommend: Either — excellent agreement. Use zonal (−301 ± 15 ‰) for consistency.**

SMO shows near-perfect agreement between the two methods. Samoa's tropical
oceanic precipitation gives an OIPC prediction almost identical to the <30°N
zonal mean. For consistency with other tropical/SH stations, the zonal mean
is recommended. SMO is already excluded from the clean analysis (phase offset
+ low amplitude), but if included, either value works.

### CGO — Cape Grim, Tasmania (40.7°S)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −294 ± 23 ‰ | −301 ± 15 ‰ |
| **Basis** | Growing-season δ²Hp = −14 ‰ (Oct–Apr) | <30°N wetland mean |

**→ Recommend: Douglas zonal (−301 ± 15 ‰), with a caveat.**

CGO is one of the two critical SH stations for the KIE constraint. Its mild
maritime climate (Tasmania) gives an OIPC prediction (−294 ‰) close to the
tropical/global mean. However, CGO has minimal local wetland sources. The
small source component in its seasonal cycle likely originates from
interhemispheric transport of NH tropical/boreal wetland CH₄. The <30°N zonal
mean (−301 ‰) is a reasonable approximation — but note that the "correct"
source δD for CGO is really a transport-weighted average of global wetland
emissions, which is closer to the flux-weighted global mean (−310 ± 15 ‰
from Douglas 2021). Using −301 ‰ is slightly more enriched than the global
mean, but the difference (9 ‰) is within uncertainty. The key point: CGO's
seasonal cycle is sink-dominated, so the source correction is small regardless.

### SPO — South Pole (90.0°S)

| | OIPC regression | Douglas zonal |
|---|---|---|
| **Value** | −556 ± 37 ‰ | −301 ± 15 ‰ |
| **Basis** | Growing-season δ²Hp = −386 ‰ at South Pole | <30°N wetland mean |

**→ Recommend: Douglas zonal (−301 ± 15 ‰), or global mean (−310 ± 15 ‰).**

The OIPC regression gives a physically meaningless value (−556 ‰) because it
applies a wetland regression to Antarctic Plateau precipitation (δ²Hp ≈ −400 ‰).
There are no wetlands anywhere on Antarctica. SPO's CH₄ comes entirely from
the global atmosphere, delayed by ~1 year of interhemispheric mixing. The
relevant source δD is the **global flux-weighted mean** of all wetland
emissions, which Douglas (2021) estimated at −310 ± 15 ‰. Since the <30°N
band dominates global wetland flux (115 of 149 Tg yr⁻¹), using the <30°N
mean (−301 ‰) is also defensible. Like CGO, SPO's seasonal cycle is heavily
sink-dominated, so source correction sensitivity is low.

---

## Summary table

| Site | Lat | OIPC (‰) | ±1σ | Zonal (‰) | ±1σ | **Recommended** | **Value (‰)** | Reasoning |
|------|-----|-----------|-----|-----------|-----|-----------------|---------------|-----------|
| ALT | +82.5 | −431 | 28 | −374 | 10 | **Zonal** | **−374 ± 10** | No local wetlands; regression extrapolates |
| ZEP | +78.9 | −344 | 24 | −374 | 10 | **Zonal** | **−374 ± 10** | No local wetlands; Gulf Stream anomaly in OIPC |
| BRW | +71.3 | −384 | 25 | −374 | 10 | **OIPC** | **−384 ± 25** | Near North Slope tundra wetlands; matches zonal |
| CBA | +55.2 | −326 | 24 | −324 | 14 | **OIPC** | **−326 ± 24** | Near boreal wetlands; excellent agreement |
| MHD | +53.3 | −306 | 23 | −324 | 14 | **Zonal** | **−324 ± 14** | Maritime; sources are continental, not local |
| AZR | +38.8 | −284 | 23 | −324 | 14 | **Zonal** | **−324 ± 14** | Excluded (non-MBL); no local sources |
| MLO | +19.5 | −320 | 24 | −301 | 15 | **Zonal** | **−301 ± 15** | Excluded (non-MBL); altitude artefact |
| KUM | +19.6 | −287 | 23 | −301 | 15 | **Zonal** | **−301 ± 15** | No local sources; tropical zonal mean |
| ASC | −8.0 | −277 | 23 | −301 | 15 | **Zonal** | **−301 ± 15** | Excluded; remote island, no local sources |
| SMO | −14.2 | −299 | 23 | −301 | 15 | **Zonal** | **−301 ± 15** | Excluded; methods agree; zonal for consistency |
| CGO | −40.7 | −294 | 23 | −301 | 15 | **Zonal** | **−301 ± 15** | SH background; source via transport |
| SPO | −90.0 | −556 | 37 | −301 | 15 | **Zonal** | **−301 ± 15** | Antarctic; OIPC meaningless; global mean source |

### Old vs new comparison (clean sites only)

| Site | Old (global) | New (recommended) | Shift | Impact on source correction |
|------|-------------|-------------------|-------|-----------------------------|
| ALT | −310 | −374 | **−64** | Large: source–atm gap increases 29% |
| ZEP | −310 | −374 | **−64** | Large: same as ALT |
| BRW | −310 | −384 | **−74** | Large: source–atm gap increases 33% |
| CBA | −310 | −326 | −16 | Moderate |
| MHD | −310 | −324 | −14 | Moderate |
| KUM | −310 | −301 | +9 | Small (opposite sign) |
| CGO | −310 | −301 | +9 | Small (opposite sign) |
| SPO | −310 | −301 | +9 | Small (opposite sign) |

---

## Key takeaways

1. **OIPC regression works well at mid-latitudes (BRW, CBA, SMO)** where
   stations are near real wetlands and the regression is within its training
   domain.

2. **Douglas zonal means are more appropriate for remote sites** (ALT, ZEP,
   MHD, KUM, ASC, CGO, SPO) where the "source δD" should represent the
   emission-weighted upstream wetland region, not local precipitation.

3. **The largest corrections are at NH high latitudes** (ALT, ZEP, BRW):
   shifting from −310 to −374/−384 ‰ increases the source–atmosphere δD gap
   by 29–33%. This will substantially change the Phase 4 decomposition and
   Phase 5 source-corrected KIE at these sites.

4. **SH sites shift only slightly** (CGO, SPO: +9 ‰). Since these are the
   cleanest sites for KIE extraction (Approach 1), the impact on the final
   α_13C_OH is small — reassuring for our main result.

5. **The choice between methods matters most where it matters least** (NH sites
   with heavy source contamination) and **matters least where it matters most**
   (SH sites with minimal sources). This is fortunate for the KIE constraint.
