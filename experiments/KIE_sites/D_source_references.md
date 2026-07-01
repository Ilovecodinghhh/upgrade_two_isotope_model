# δD Source Signatures of Wetland/Microbial CH₄ — Reference List

References relevant to understanding the **latitude-dependent δD of microbial methane sources**, with emphasis on wetlands, rice paddies, and the relationship between source water δD and emitted CH₄ δD. Organized by topic.

**Status of each reference:**
- ✅ = Full paper + data available in this repo
- 📄 = Markdown summary in repo (no raw data)
- ❌ = Not in repo — to be obtained

---

## 1. Global Compilations of CH₄ Isotopic Source Signatures

### ✅ Douglas et al. (2021) — **KEY PAPER: Latitude-resolved freshwater δD-CH₄**
- **Citation:** Douglas, P. M. J., Stratigopoulos, E., Park, S., and Phan, D.: Geographic variability in freshwater methane hydrogen isotope ratios and its implications for global isotopic source signatures, *Biogeosciences*, 18, 3505–3527, https://doi.org/10.5194/bg-18-3505-2021, 2021.
- **In repo:** `Other References/Isotope/Douglas2021BG.md` (summary) + `ImportantReferences/Douglas2021EGU/` (full paper) + `ImportantReferences/Douglas2021EGU_DS.xlsx` (dataset) + `ImportantReferences/Douglas2021EGU_SI/` (supplement)
- **Key finding:** Compiled 129 freshwater sites globally. Provides latitude-resolved δ²H-CH₄ for three bands:

  | Source | Latitude | δ²H-CH₄ (‰) | ±1σ | Flux (Tg/yr) |
  |--------|----------|-------------|-----|-------------|
  | Wetlands | <30°N | −301 | 15 | 115 |
  | Wetlands | 30–60°N | −324 | 14 | 25 |
  | Wetlands | >60°N | **−374** | 10 | 9 |
  | Inland waters | <30°N | −301 | 12 | 80 |
  | Inland waters | 30–60°N | −308 | 18 | 64 |
  | Inland waters | >60°N | **−347** | 9 | 16 |
  | Rice paddies | <30°N | −324 | 8 | 19 |
  | Rice paddies | 30–60°N | −325 | 8 | 12 |

- **Regression:** δ²H-CH₄ vs δ²H-H₂O slope ≈ 0.52–0.69 (flatter than Waldron 1999's 0.675). δ²H-H₂O explains ~42% of the variance in δ²H-CH₄.
- **Flux-weighted global freshwater mean:** −310 ± 15‰ (natural wetlands alone: −310 ± 25‰).
- **Bottom-up global source δ²H-CH₄:** −278 ± 15‰ — overlaps with top-down estimates from Rice et al. (2016).
- **Critical for KIE_sites:** The ~73‰ spread from tropical (−301‰) to boreal (−374‰) wetlands is much larger than our assumed 30‰ uncertainty on a global mean.

### ✅ Sherwood et al. (2017) — Global inventory v2017
- **Citation:** Sherwood, O. A., Schwietzke, S., Arling, V. A., and Etiope, G.: A global inventory of gas geochemistry data from fossil fuel, microbial and burning sources, version 2017, *Earth Syst. Sci. Data*, 9, 639–656, https://doi.org/10.5194/essd-9-639-2017, 2017.
- **In repo:** `Other References/Isotope/Sherwood2017ESSD.md` + `ImportantReferences/Sherwood2017ESSD/` + `ImportantReferences/Sherwood2017ESSD_DS.xlsx`
- **Key content:** Most comprehensive database of δ¹³C and δ²H source signatures by category. Used by most global budget studies. Does **not** provide latitude-resolved signatures for microbial sources — this is the gap Douglas 2021 fills.

### 📄 Milkov & Etiope (2018) — Revised gas geochemistry classification
- **Citation:** Milkov, A. V. and Etiope, G.: Revised genetic diagrams for natural gases based on a global dataset of >20,000 samples, *Org. Geochem.*, 125, 109–120, 2018.
- **In repo:** `Other References/Isotope/Milkov2018OrgGeochem.md`
- **Key content:** Updated Whiticar-Schoell δ¹³C vs δD classification diagrams using >20,000 samples. Provides distributions for thermogenic, microbial, and abiotic CH₄.

### ✅ Etiope et al. (2019) — Gridded geological CH₄ isotopes
- **Citation:** Etiope, G., et al.: Gridded maps of geological methane emissions and their isotopic signature, *Earth Syst. Sci. Data*, 11, 1–22, 2019.
- **In repo:** `Other References/Isotope/Etiope2019ESSD.md` + `Etiope2019ESSD/` (full paper with images) + `Etiope2019ESSD_DS/` (gridded data) + `Etiope2019ESSD_SI/` (supplement)
- **Key content:** Geological CH₄ δD and δ¹³C signatures gridded globally. Useful for constraining non-microbial background.

---

## 2. The δD(CH₄) vs δD(H₂O) Relationship — Foundational Studies

### ✅ Waldron et al. (1999) — Original global δD(CH₄)–δD(H₂O) regression
- **Citation:** Waldron, S., Lansdown, J. M., Scott, E. M., Fallick, A. E., and Hall, A. J.: The global influence of the hydrogen isotope composition of water on that of bacteriogenic methane from shallow freshwater environments, *Geochim. Cosmochim. Acta*, 63, 2237–2245, 1999.
- **In repo:** `Other References/D_source/Waldron1999GCA/` (full paper + images)
- **Key finding:** First established the global regression from 46 sites:
  ```
  δD(CH₄) = 0.675(±0.1) × δD(H₂O) − 284(±6) ‰
  ```
  ~50% of variation in natural δD-CH₄ explained by δD-H₂O. Also showed laboratory incubations give a steeper relationship (slope 0.444, intercept −321‰), suggesting post-production processes (oxidation, diffusion) modify the in-situ signal. Argued that δD-CH₄ may be largely **independent of methanogenic pathway** due to enzyme-mediated H isotope exchange.

### ✅ Chanton et al. (2006) — Alaska wetlands, 60–70°N transect
- **Citation:** Chanton, J. P., Fields, D., and Hines, M. E.: Controls on the hydrogen isotopic composition of biogenic methane from high-latitude terrestrial wetlands, *J. Geophys. Res.*, 111, G04004, 2006.
- **In repo:** `Other References/D_source/Chanton2006JGR/` (full paper + images)
- **Key finding:** Measured δD-H₂O (−108 to −161‰) and δD-CH₄ (−308 to −394‰) along a N–S Alaska transect from 60°N to 70°N. Found:
  - δD-CH₄ and δD-H₂O are significantly correlated along the latitude gradient
  - Apparent αD (H₂O → CH₄) varies from 1.26 to 1.42, inversely correlated with αC (CO₂ → CH₄)
  - Combined Waldron+Chanton dataset gives flatter slope: **0.54 ± 0.05**
  - Demonstrates that **both** water δD and methanogenic pathway influence the emitted CH₄ δD at high latitudes

### ✅ Whiticar et al. (1986) — Foundational isotope classification
- **Citation:** Whiticar, M. J., Faber, E., and Schoell, M.: Biogenic methane formation in marine and freshwater environments: CO₂ reduction vs. acetate fermentation — Isotope evidence, *Geochim. Cosmochim. Acta*, 50, 693–709, 1986.
- **In repo:** `Other References/D_source/Whiticar1986GCA/` (full paper + images)
- **Key finding:** The foundational "Whiticar diagram" paper. Proposed that δ¹³C and δD together can separate:
  - **CO₂ reduction** (marine): δ¹³C = −110 to −60‰, δD = −250 to −170‰
  - **Acetate fermentation** (freshwater): δ¹³C = −65 to −50‰, δD = −400 to −250‰
  - H₂O → CH₄ hydrogen fractionation: −180 ± 20‰ for CO₂ reduction
  - Now known to oversimplify (Waldron 1999 and Milkov 2018 showed the boundaries are blurred)

### ✅ Whiticar (1999) — Comprehensive δD fractionation review
- **Citation:** Whiticar, M. J.: Carbon and hydrogen isotope systematics of bacterial formation and oxidation of methane, *Chem. Geol.*, 161, 291–314, 1999.
- **In repo:** `Other References/D_source/Whiticar1999CG/` (full paper + images)
- **Key finding:** Comprehensive review of fractionation factors for ¹³C and D in methanogenesis and methanotrophy:
  - Acetoclastic: δD-CH₄ can be as negative as −531‰
  - CO₂ reduction: δD-CH₄ ≈ −170 to −250‰
  - Aerobic oxidation εH = 95–285 (enriches residual CH₄ in D)
  - Anaerobic oxidation also fractionates H isotopes
  - Key reference values for KIE during oxidation

---

## 3. Methanogenic Pathway Controls on δD

### 📄 Okumura et al. (2016) — H₂ availability controls δD fractionation
- **Citation:** Okumura, T., et al.: Hydrogen and carbon isotope systematics in hydrogenotrophic methanogenesis under H₂-limited and H₂-enriched conditions, *Prog. Earth Planet. Sci.*, 3, 14, 2016.
- **In repo:** `Other References/Isotope/Okumura2016PEPS.md`
- **Key finding:** CH₄–H₂O hydrogen isotope fractionation during hydrogenotrophic methanogenesis: αH ≈ 0.67–0.69 under low-pH₂ (syntrophic) conditions. The "δD-H₂ effect" — substrate H₂ isotopic composition can also influence CH₄ δD. Contributes to the ~58% of δD-CH₄ variance NOT explained by δD-H₂O.

### 📄 Conrad & Claus (2009) — Pathway-specific fractionation in lake sediments
- **Citation:** Conrad, R. and Claus, P.: Characterization of stable isotope fractionation during methane production in the sediment of a eutrophic lake, *Limnol. Oceanogr.*, 54(2), 457–471, 2009.
- **In repo:** `Other References/Isotope/Conrad2009LO.md`
- **Key finding:** Used methyl fluoride inhibitor to specifically block acetoclastic methanogenesis. Found αCO₂/CH₄ = 1.03–1.09 for hydrogenotrophic pathway. Pathway mix varies with depth (35–60% hydrogenotrophic). Pathway controls both δ¹³C and δD.

### ❌ Valentine et al. (2004) — H₂ thermodynamics and isotope fractionation
- **Citation:** Valentine, D. L., et al.: Carbon and hydrogen isotope fractionation by moderately thermophilic methanogens, *Geochim. Cosmochim. Acta*, 68, 1571–1590, 2004.
- **Key finding:** "Differential reversibility" — carbon (and hydrogen) isotope fractionation depends on H₂ availability. Key mechanistic paper for understanding why fractionation varies between environments.

### ❌ Sugimoto & Wada (1995) — δD fractionation in rice paddies
- **Citation:** Sugimoto, A. and Wada, E.: Hydrogen isotopic composition of bacterial methane: CO₂/H₂ reduction and acetate fermentation, *Geochim. Cosmochim. Acta*, 59, 1329–1337, 1995.
- **Key finding:** Early pathway-dependent δD fractionation measurements in rice paddies.

---

## 4. Precipitation δD — The Underlying Latitude Gradient

### ❌ Bowen & Revenaugh (2003) — Global precipitation δD model
- **Citation:** Bowen, G. J. and Revenaugh, J.: Interpolating the isotopic composition of modern meteoric precipitation, *Water Resour. Res.*, 39, 1299, 2003.
- **Key finding:** The Online Isotopes in Precipitation Calculator (OIPC). Provides gridded δD/δ¹⁸O in precipitation globally. Tropical δD-precip ≈ −10 to −30‰; high-latitude ≈ −100 to −200‰. Input to Douglas 2021 for estimating δD-H₂O at sites without measurements.

### ❌ IAEA/WMO GNIP Database
- **Citation:** IAEA/WMO: Global Network of Isotopes in Precipitation. https://nucleus.iaea.org/wiser
- **Key resource:** Direct measurements of precipitation δD at >1000 stations globally. Empirical foundation for the latitude effect.

---

## 5. Transport and Oxidation Effects on δD

### 📄 Chanton (2005) — Gas transport effects on wetland CH₄ isotopes
- **Citation:** Chanton, J. P.: The effect of gas transport on the isotope signature of methane in wetlands, *Org. Geochem.*, 36, 753–768, 2005.
- **In repo:** `Other References/Isotope/Chanton2005OrgGeochem.md`
- **Key finding:** Plant-mediated transport (aerenchyma) fractionates both ¹³C and D. Diffusive transport preferentially removes lighter isotopologues, enriching residual CH₄. The transport mechanism affects the **emitted** signature, complicating the link between production δD and atmospheric δD.

### ❌ Snover & Quay (2000) — Soil oxidation δD fractionation
- **Citation:** Snover, A. K. and Quay, P. D.: Hydrogen and carbon kinetic isotope effects during soil uptake of atmospheric methane, *Global Biogeochem. Cycles*, 14, 25–39, 2000.
- **Key finding:** αD ≈ 1.066 for soil uptake of CH₄. The fraction oxidized varies 10–90% across ecosystems, creating large variability in emitted δD.

### ❌ Hornibrook (2009) — Peatland isotope systematics
- **Citation:** Hornibrook, E. R. C.: The stable carbon isotope composition of methane produced and emitted from northern peatlands, AGU Geophysical Monograph Series, 184, 187–203, 2009.
- **Key finding:** Ombrotrophic bogs vs minerotrophic fens have distinct δ¹³C (and δD) signatures.

---

## 6. Regional/Latitude-Specific Measurements & Atmospheric δD

### ✅ Rice et al. (2016) — Atmospheric δD-CH₄ box model 1977–2005
- **Citation:** Rice, A. L., et al.: Atmospheric methane isotopic record favors fossil sources flat in 1980s and 1990s with recent increase, *Proc. Natl. Acad. Sci.*, 113, 10791–10796, 2016.
- **In repo:** `ImportantReferences/Rice2016PNAS/` (full paper) + `ImportantReferences/Rice2016PNAS_SI/` (supplement)
- **Key finding:** Applied box model to atmospheric δ¹³C and δD measurements (1977–2005). Used a **constant** −322‰ for both low- and high-latitude wetlands — Douglas (2021) showed this is inaccurate for both bands. Top-down global source δ²H estimate: −258 to −289‰.

### ❌ Bock et al. (2010) — Ice core δD, paleo source signatures
- **Citation:** Bock, M., et al.: Hydrogen isotopes preclude marine hydrate CH₄ emissions at the onset of Dansgaard-Oeschger events, *Science*, 328, 1686–1689, 2010.
- **Key finding:** Used differentiated signatures: tropical wetlands −320‰, boreal −370‰. Their tropical estimate is lower than Douglas (2021)'s −301‰.

### ❌ Fisher et al. (2011) — Arctic atmospheric δ¹³C constraint
- **Citation:** Fisher, R. E., et al.: Arctic methane sources: Isotopic evidence for atmospheric inputs, *Geophys. Res. Lett.*, 38, L21803, 2011.
- **Key finding:** Atmospheric δ¹³C in Arctic suggests high-latitude wetland δ¹³C ≈ −68 ± 4‰ (more depleted than in situ measurements).

---

## 7. Spatially Resolved Signatures for Atmospheric Models

### 📄 Ganesan et al. (2018) — Spatially resolved wetland δ¹³C
- **Citation:** Ganesan, A. L., et al.: Spatially resolved isotopic source signatures of wetland methane emissions, *Geophys. Res. Lett.*, 45, 3737–3745, 2018.
- **In repo:** `Other References/Isotope/Ganesan2018GRL.md`
- **Key finding:** Used JULES land surface model to generate gridded wetland δ¹³C. Low-latitude δ¹³C ≈ −56.7‰ (C₄ influence), high-latitude ≈ −67‰. Methodology transferable to δD.

### 📄 Röckmann et al. (2016) — In situ δD at Cabauw tower
- **Citation:** Röckmann, T., et al.: In situ observations of the isotopic composition of methane at the Cabauw tall tower site, *Atmos. Chem. Phys.*, 16, 10469–10487, 2016.
- **In repo:** `Other References/Isotope/Rockmann2016ACP.md`

### 📄 Feinberg et al. (2018) — Spatially varying δ¹³C in global models
- **Citation:** Feinberg, A. I., et al.: Isotopic source signatures: Impact of regional variability on the δ¹³CH₄ trend, *J. Geophys. Res. Atmos.*, 123, 4841–4856, 2018.
- **In repo:** `Other References/Isotope/Feinberg2018JGR.md`
- **Key content:** Showed regional source signature variability significantly affects simulated atmospheric δ¹³C. Analogous work needed for δD.

### 📄 Chen et al. (2025) — Updated atmospheric δ¹³C and δD
- **In repo:** `Other References/Isotope/Chen2025ACP.md`

---

## 8. Summary: Latitude-Resolved δD-CH₄ Source Signatures

Best-estimate values from Douglas et al. (2021) Table 1:

| Source category | Latitude band | δD-CH₄ (‰ VSMOW) | ± (1σ) | Flux (Tg/yr) |
|----------------|---------------|---------------------|--------|---------------|
| **Wetlands** | **<30°N** | **−301** | 15 | 115 |
| **Wetlands** | **30–60°N** | **−324** | 14 | 25 |
| **Wetlands** | **>60°N** | **−374** | 10 | 9 |
| Inland waters | <30°N | −301 | 12 | 80 |
| Inland waters | 30–60°N | −308 | 18 | 64 |
| Inland waters | >60°N | −347 | 9 | 16 |
| Rice paddies | <30°N | −324 | 8 | 19 |
| Rice paddies | 30–60°N | −325 | 8 | 12 |
| Enteric ferm. | Global | −308 | 28 | 111 |
| Landfills | Global | −297 | 6 | 65 |
| Coal mining | Global | −232 | 5 | 42 |
| Oil & gas | Global | −189 | 2 | 79 |
| Biomass burning | Global | −211 | 15 | 17 |
| Termites | Global | −343 | 50 | 9 |
| Permafrost | Global | −374 | 15 | 1 |

**Key takeaway for KIE_sites experiment:** At our high-latitude sites (ALT 82°N, ZEP 79°N, BRW 71°N), the relevant wetland δD is **−374 ± 10‰**, not our assumed global mean of **−310 ± 30‰**. This means the source–atmosphere gap is 374 − 86 = **288‰** (not 224‰) — a 29% increase that would substantially change the Phase 5 source correction results.

---

## 9. Papers Still Needed (❌) — Priority Ranking

| Priority | Reference | Topic |
|----------|-----------|-------|
| ⭐⭐⭐ | Valentine et al. (2004) *GCA* | H₂ thermodynamics and δD fractionation |
| ⭐⭐ | Bock et al. (2010) *Science* | Ice core δD, tropical vs boreal |
| ⭐⭐ | Snover & Quay (2000) *GBC* | Soil oxidation δD KIE |
| ⭐⭐ | Bowen & Revenaugh (2003) *WRR* | Precipitation δD model (OIPC) |
| ⭐ | Sugimoto & Wada (1995) *GCA* | Rice paddy δD fractionation |
| ⭐ | Fisher et al. (2011) *GRL* | Arctic atmospheric δ¹³C |
| ⭐ | Hornibrook (2009) *AGU Mono.* | Peatland bog vs fen systematics |
