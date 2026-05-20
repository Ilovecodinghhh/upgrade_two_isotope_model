# δD Source Signatures of Wetland/Microbial CH₄ — Reference List

References relevant to understanding the **latitude-dependent δD of microbial methane sources**, with emphasis on wetlands, rice paddies, and the relationship between source water δD and emitted CH₄ δD. Organized by topic.

---

## 1. Global Compilations of CH₄ Isotopic Source Signatures

These papers provide the databases and summary statistics (including δD) used in global CH₄ budget models.

### Douglas et al. (2021) — **KEY PAPER: Latitude-resolved freshwater δD-CH₄**
- **Citation:** Douglas, P. M. J., Stratigopoulos, E., Park, S., and Phan, D.: Geographic variability in freshwater methane hydrogen isotope ratios and its implications for global isotopic source signatures, *Biogeosciences*, 18, 3505–3527, https://doi.org/10.5194/bg-18-3505-2021, 2021.
- **In repo:** `Other References/Isotope/Douglas2021BG.md`
- **Key finding:** Compiled 129 sites globally. Provides latitude-resolved δ²H-CH₄ for three bands:
  - Wetlands (<30°N): **−301 ± 15‰**
  - Wetlands (30–60°N): **−324 ± 14‰**
  - Wetlands (>60°N): **−374 ± 10‰**
  - Inland waters (<30°N): **−301 ± 12‰**
  - Inland waters (30–60°N): **−308 ± 18‰**
  - Inland waters (>60°N): **−347 ± 9‰**
  - Rice (<30°N): **−324 ± 8‰**
  - Rice (30–60°N): **−325 ± 8‰**
  - Flux-weighted global freshwater mean: **−310 ± 15‰**
- **Regression:** δ²H-CH₄ vs δ²H-H₂O slope ≈ 0.52–0.69 (flatter than Waldron 1999's 0.68 ± 0.1). δ²H-H₂O explains ~42% of the variance in δ²H-CH₄.
- **Critical for KIE_sites:** The ~73‰ spread from tropical to boreal wetlands is much larger than our assumed 30‰ uncertainty on a global mean.

### Sherwood et al. (2017) — Global inventory v2017
- **Citation:** Sherwood, O. A., Schwietzke, S., Arling, V. A., and Etiope, G.: A global inventory of gas geochemistry data from fossil fuel, microbial and burning sources, version 2017, *Earth Syst. Sci. Data*, 9, 639–656, https://doi.org/10.5194/essd-9-639-2017, 2017.
- **In repo:** `Other References/Isotope/Sherwood2017ESSD.md`
- **Key content:** The most comprehensive database of δ¹³C and δ²H source signatures by category. Used by most global budget studies. However, does **not** provide latitude-resolved signatures for microbial sources — this is the gap Douglas 2021 fills.

### Milkov & Etiope (2018) — Revised gas geochemistry classification
- **Citation:** Milkov, A. V. and Etiope, G.: Revised genetic diagrams for natural gases based on a global dataset of >20,000 samples, *Org. Geochem.*, 125, 109–120, https://doi.org/10.1016/j.orggeochem.2018.09.002, 2018.
- **In repo:** `Other References/Isotope/Milkov2018OrgGeochem.md`
- **Key content:** Updated Whiticar-Schoell classification diagrams (δ¹³C vs δD) using >20,000 samples. Provides distributions for thermogenic, microbial (acetoclastic vs hydrogenotrophic), and abiotic CH₄.

### Etiope & Sherwood (2019) — Etiope ESSD database update
- **Citation:** Etiope, G., Ciotoli, G., Schwietzke, S., and Schoell, M.: Gridded maps of geological methane emissions and their isotopic signature, *Earth Syst. Sci. Data*, 11, 1–22, https://doi.org/10.5194/essd-11-1-2019, 2019.
- **In repo:** `Other References/Isotope/Etiope2019ESSD.md`
- **Key content:** Geological CH₄ δD and δ¹³C signatures gridded globally. Useful for constraining non-microbial background.

---

## 2. The δD(CH₄) vs δD(H₂O) Relationship — Foundational Studies

The core physical relationship governing latitude dependence: source water δD controls microbial CH₄ δD.

### Waldron et al. (1999) — Original δD(CH₄) vs δD(H₂O) regression
- **Citation:** Waldron, S., Lansdown, J. M., Scott, E. M., Fallick, A. E., and Hall, A. J.: The global influence of the hydrogen isotope composition of water on that of bacteriogenic methane from shallow freshwater environments, *Geochim. Cosmochim. Acta*, 63, 2237–2245, https://doi.org/10.1016/S0016-7037(99)00192-1, 1999.
- **Not in repo**
- **Key finding:** First established the global regression: δD-CH₄ ≈ 0.68(±0.1) × δD-H₂O − 284(±13)‰. Showed that precipitation δD is the primary global predictor of wetland CH₄ δD. Douglas (2021) updated this with a flatter slope (0.52–0.69) and more high-latitude data.

### Chanton et al. (2006) — Alaska wetlands extend the regression
- **Citation:** Chanton, J. P., Fields, D., and Hines, M. E.: Controls on the hydrogen isotopic composition of biogenic methane from high-latitude terrestrial wetlands, *J. Geophys. Res.*, 111, G04004, https://doi.org/10.1029/2005JG000134, 2006.
- **Not in repo**
- **Key finding:** Added Alaskan wetland δD-CH₄ measurements to the Waldron dataset. Found slope ≈ 0.54 ± 0.05 (flatter than Waldron). Boreal/subarctic wetlands (δD-H₂O ≈ −100 to −150‰) produce CH₄ with δD ≈ −350 to −400‰. This is the paper that first demonstrated the steeper-than-expected depletion at high latitudes.

### Whiticar et al. (1986) — Isotope classification of microbial CH₄
- **Citation:** Whiticar, M. J., Faber, E., and Schoell, M.: Biogenic methane formation in marine and freshwater environments: CO₂ reduction vs. acetate fermentation — Isotope evidence, *Geochim. Cosmochim. Acta*, 50, 693–709, https://doi.org/10.1016/0016-7037(86)90346-7, 1986.
- **Not in repo**
- **Key finding:** The foundational isotope classification paper. Established the δ¹³C vs δD "Whiticar diagram" separating acetoclastic (δD ≈ −250 to −400‰) from hydrogenotrophic (δD ≈ −150 to −250‰) methanogenesis. Widely used but now known to oversimplify (see Milkov & Etiope 2018).

### Whiticar (1999) — Updated review of CH₄ isotope geochemistry
- **Citation:** Whiticar, M. J.: Carbon and hydrogen isotope systematics of bacterial formation and oxidation of methane, *Chem. Geol.*, 161, 291–314, https://doi.org/10.1016/S0009-2541(99)00092-3, 1999.
- **Not in repo**
- **Key finding:** Comprehensive review of fractionation factors for both ¹³C and D in methanogenesis and methane oxidation. Key reference values for acetoclastic (αD ≈ 1.22–1.44) and hydrogenotrophic (αD ≈ 1.14–1.51) H/D fractionation.

---

## 3. Methanogenic Pathway Controls on δD

The pathway of methanogenesis (acetoclastic vs hydrogenotrophic) modulates the δD(CH₄)–δD(H₂O) relationship, adding scatter beyond the simple latitude trend.

### Okumura et al. (2016) — H₂ availability controls δD fractionation
- **Citation:** Okumura, T., Kawagucci, S., Saito, Y., Matsui, Y., Takai, K., and Imachi, H.: Hydrogen and carbon isotope systematics in hydrogenotrophic methanogenesis under H₂-limited and H₂-enriched conditions, *Prog. Earth Planet. Sci.*, 3, 14, https://doi.org/10.1186/s40645-016-0088-3, 2016.
- **In repo:** `Other References/Isotope/Okumura2016PEPS.md`
- **Key finding:** δD fractionation between CH₄ and H₂O during hydrogenotrophic methanogenesis depends on H₂ availability. Low-pH₂ (syntrophic) cultures: αH ≈ 0.67–0.69. This variability contributes to the ~58% of δD-CH₄ variance NOT explained by δD-H₂O.

### Conrad (2009) — Pathway-specific fractionation in lake sediments
- **Citation:** Conrad, R. and Claus, P.: Characterization of stable isotope fractionation during methane production in the sediment of a eutrophic lake, *Limnol. Oceanogr.*, 54(2), 457–471, 2009.
- **In repo:** `Other References/Isotope/Conrad2009LO.md`
- **Key finding:** Determined pathway-specific ¹³C fractionation factors using inhibitors. Hydrogenotrophic methanogenesis increases with depth. The pathway mix affects both δ¹³C and δD of produced CH₄.

### Valentine et al. (2004) — H₂ thermodynamics and isotope fractionation
- **Citation:** Valentine, D. L., Chidthaisong, A., Rice, A., Reeburgh, W. S., and Tyler, S. C.: Carbon and hydrogen isotope fractionation by moderately thermophilic methanogens, *Geochim. Cosmochim. Acta*, 68, 1571–1590, https://doi.org/10.1016/j.gca.2003.10.012, 2004.
- **Not in repo**
- **Key finding:** Demonstrated "differential reversibility" — carbon isotope fractionation depends on H₂ availability. Also measured hydrogen isotope fractionation for thermophilic methanogens. Key mechanistic paper for understanding why δD fractionation varies between environments.

### Sugimoto & Wada (1995) — δD fractionation in rice paddies
- **Citation:** Sugimoto, A. and Wada, E.: Hydrogen isotopic composition of bacterial methane: CO₂/H₂ reduction and acetate fermentation, *Geochim. Cosmochim. Acta*, 59, 1329–1337, 1995.
- **Not in repo**
- **Key finding:** Early measurement of pathway-dependent δD fractionation. Rice paddy methanogenesis has a significant acetoclastic component, producing CH₄ with different δD than purely hydrogenotrophic environments.

---

## 4. Precipitation δD — The Underlying Latitude Gradient

The latitude dependence of source water δD drives the latitude dependence of wetland CH₄ δD.

### Bowen & Revenaugh (2003) — Global precipitation δD model
- **Citation:** Bowen, G. J. and Revenaugh, J.: Interpolating the isotopic composition of modern meteoric precipitation, *Water Resour. Res.*, 39, 1299, https://doi.org/10.1029/2003WR002086, 2003.
- **Not in repo**
- **Key finding:** The Online Isotopes in Precipitation Calculator (OIPC). Provides gridded estimates of δD and δ¹⁸O in precipitation globally. Input to Douglas 2021 for estimating δD-H₂O at sites without measurements. Tropical δD-precip ≈ −10 to −30‰; high-latitude ≈ −100 to −200‰.

### IAEA/WMO GNIP Database — Global Network of Isotopes in Precipitation
- **Citation:** IAEA/WMO: Global Network of Isotopes in Precipitation. The GNIP Database. https://nucleus.iaea.org/wiser
- **Not in repo** (database, not a paper)
- **Key resource:** Direct measurements of precipitation δD at >1000 stations globally. The empirical foundation for the latitude effect.

---

## 5. Transport and Oxidation Effects on δD

Soil oxidation and plant-mediated transport modify the δD of CH₄ between production and atmospheric emission.

### Chanton (2005) — Gas transport effects on wetland CH₄ isotopes
- **Citation:** Chanton, J. P.: The effect of gas transport on the isotope signature of methane in wetlands, *Org. Geochem.*, 36, 753–768, https://doi.org/10.1016/j.orggeochem.2004.10.007, 2005.
- **In repo:** `Other References/Isotope/Chanton2005OrgGeochem.md`
- **Key finding:** Plant-mediated transport (aerenchyma) fractionates both ¹³C and D. Diffusive transport preferentially removes lighter isotopologues, enriching the remaining CH₄. The transport mechanism affects the **emitted** signature, complicating the link between production δD and atmospheric δD.

### Snover & Quay (2000) — Soil oxidation δD fractionation
- **Citation:** Snover, A. K. and Quay, P. D.: Hydrogen and carbon kinetic isotope effects during soil uptake of atmospheric methane, *Global Biogeochem. Cycles*, 14, 25–39, https://doi.org/10.1029/1999GB900089, 2000.
- **Not in repo**
- **Key finding:** Measured the KIE for D/H during soil uptake of CH₄: αD ≈ 1.066. Soil oxidation before emission enriches the residual CH₄ in D (shifts δD positive). The fraction oxidized varies by 10–90% across ecosystems, creating large variability in emitted δD.

### Hornibrook (2009) — Peatland δ¹³C and δD systematics
- **Citation:** Hornibrook, E. R. C.: The stable carbon isotope composition of methane produced and emitted from northern peatlands, in: *Carbon Cycling in Northern Peatlands*, AGU Geophysical Monograph Series, 184, 187–203, https://doi.org/10.1029/2008GM000828, 2009.
- **Not in repo**
- **Key finding:** Reviews how methanogenic pathway, transport, and oxidation interact to determine the isotopic composition of peatland CH₄ emissions. Ombrotrophic bogs vs minerotrophic fens have distinct δ¹³C (and δD) signatures.

---

## 6. Regional/Latitude-Specific Wetland δD Measurements

Field measurements from specific latitude bands.

### Walter Anthony et al. (2012) — Arctic thermokarst lake ebullition
- **Citation:** Walter Anthony, K. M. and Anthony, P.: Constraining spatial variability of methane ebullition seeps in thermokarst lakes using point process models, *J. Geophys. Res. Biogeosci.*, 118, 1–19, https://doi.org/10.1002/jgrg.20087, 2013.
- **In repo:** `Other References/Isotope/WalterAnthony2008JGR.md`
- **Key content:** Characterizes CH₄ ebullition in panarctic thermokarst lakes. Provides flux and spatial variability data for Arctic sources. δD of Arctic lake ebullition CH₄ is typically very depleted (< −350‰).

### Fisher et al. (2011) — Arctic δ¹³C constraint on wetland emissions
- **Citation:** Fisher, R. E., Sriskantharajah, S., Lowry, D., Lanoisellé, M., Fowler, C. M. R., James, R. H., Hermansen, O., Lund Myhre, C., Stohl, A., Greinert, J., Nisbet-Jones, P. B. R., Mber, J., and Nisbet, E. G.: Arctic methane sources: Isotopic evidence for atmospheric inputs, *Geophys. Res. Lett.*, 38, L21803, https://doi.org/10.1029/2011GL049319, 2011.
- **Not in repo**
- **Key finding:** Atmospheric δ¹³C measurements in the Arctic suggest high-latitude wetland δ¹³C ≈ −68 ± 4‰ (more depleted than in situ measurements). Implies atmospheric sampling may capture a different mix of wetland types than field campaigns.

### Bock et al. (2010) — Ice core δD constraints on paleo wetland emissions
- **Citation:** Bock, M., Schmitt, J., Möller, L., Spahni, R., Blunier, T., and Fischer, H.: Hydrogen isotopes preclude marine hydrate CH₄ emissions at the onset of Dansgaard-Oeschger events, *Science*, 328, 1686–1689, https://doi.org/10.1126/science.1187651, 2010.
- **Not in repo**
- **Key finding:** Used ice core δD-CH₄ to constrain source changes. Applied differentiated source signatures: tropical wetlands −320‰, boreal wetlands −370‰. Their tropical estimate is significantly lower than Douglas (2021)'s −301‰, illustrating the ongoing uncertainty in these values.

### Rice et al. (2016) — Box model analysis of atmospheric δD-CH₄
- **Citation:** Rice, A. L., Butenhoff, C. L., Teama, D. G., Röger, F. H., Khalil, M. A. K., and Rasmussen, R. A.: Atmospheric methane isotopic record favors fossil sources flat in 1980s and 1990s with recent increase, *Proc. Natl. Acad. Sci.*, 113, 10791–10796, https://doi.org/10.1073/pnas.1522923113, 2016.
- **Not in repo**
- **Key finding:** Applied a box model with atmospheric δD-CH₄ measurements (1977–2005) to constrain source changes. Used a **constant** −322‰ for both low- and high-latitude wetlands — Douglas (2021) showed this is inaccurate for both bands.

---

## 7. Spatially Resolved Isotopic Source Signatures for Atmospheric Models

Papers that implement or argue for latitude-dependent source signatures in global models.

### Ganesan et al. (2018) — Spatially resolved wetland δ¹³C signatures
- **Citation:** Ganesan, A. L., Stell, A. C., Gedney, N., Comyn-Platt, E., Hayman, G., Rigby, M., Poulter, B., and Hornibrook, E. R. C.: Spatially resolved isotopic source signatures of wetland methane emissions, *Geophys. Res. Lett.*, 45, 3737–3745, https://doi.org/10.1002/2018GL077536, 2018.
- **In repo:** `Other References/Isotope/Ganesan2018GRL.md`
- **Key finding:** Used JULES land surface model + wetland isotope parameterization to generate gridded wetland δ¹³C-CH₄. Found large regional differences: low-latitude δ¹³C ≈ −56.7‰ (C₄ plant influence), high-latitude ≈ −67‰. Argues spatially resolved signatures are critical for atmospheric inversions. (Focuses on δ¹³C, not δD, but methodology is transferable.)

### Röckmann et al. (2016) — In situ δD measurements at Cabauw
- **Citation:** Röckmann, T., et al.: In situ observations of the isotopic composition of methane at the Cabauw tall tower site, *Atmos. Chem. Phys.*, 16, 10469–10487, https://doi.org/10.5194/acp-16-10469-2016, 2016.
- **In repo:** `Other References/Isotope/Rockmann2016ACP.md`
- **Key content:** High-frequency δ¹³C and δD measurements at a European tower site. Demonstrates the potential for continuous δD monitoring to constrain local/regional sources with distinct δD signatures.

### Feinberg et al. (2018) — Improved δ¹³C source signatures
- **Citation:** Feinberg, A. I., Coulon, A., Stenke, A., Schwietzke, S., and Peter, T.: Isotopic source signatures: Impact of regional variability on the δ¹³CH₄ trend and spatial distribution, *J. Geophys. Res. Atmos.*, 123, 4841–4856, https://doi.org/10.1002/2017JD027784, 2018.
- **In repo:** `Other References/Isotope/Feinberg2018JGR.md`
- **Key content:** Implemented spatially varying δ¹³C source signatures in a global model. Showed regional variability in source signatures significantly affects simulated atmospheric δ¹³C fields. Analogous work is needed for δD.

### Chen et al. (2025) — Updated atmospheric δ¹³C and δD analysis
- **Citation:** Chen, Y., et al. (2025), *Atmos. Chem. Phys.*
- **In repo:** `Other References/Isotope/Chen2025ACP.md`
- **Key content:** Recent analysis using both δ¹³C and δD atmospheric measurements. Check for updated source signature assumptions.

---

## 8. Ruminant and Other Microbial δD — Latitude Dependence?

Ruminants are another major microbial source whose δD may be latitude-dependent (via feed water δD).

### Sherwood et al. (2017) — reports ruminant δD ≈ −308 ± 28‰
- See entry in Section 1 above.
- **Note:** Douglas (2021) Table 1 flags that enteric fermentation/manure δD has the **largest flux-weighted uncertainty** globally. No latitude-resolved data yet exist. Given that livestock drink local water (δD varies by latitude), there is a plausible but unstudied latitude gradient here too.

---

## 9. Atmospheric δD-CH₄ Measurements and Trends

Papers measuring atmospheric δD-CH₄ that could serve as top-down constraints.

### Sperlich et al. (2015) — δD-CH₄ in firn air
- **Citation:** Sperlich, P., et al.: Carbon isotope ratios suggest no additional methane from boreal wetlands during the rapid Greenland Interstadial 21.2, *Atmos. Chem. Phys.*, 15, 7247–7257, 2015.
- **In repo:** `Other References/Isotope/Sperlich2015ACP.md`

### Ferretti et al. (2006) — Ice core δD-CH₄
- **Citation:** Ferretti, D. F., et al.: Kinetic isotope effects and their use in directly constraining the methane budget, *Atmos. Chem. Phys. Discuss.*, 6, 9823–9857, 2006.
- **In repo:** `Other References/Isotope/Ferretti2006ACPD.md`

### Kelly et al. (2022) — δ¹³C and δD trends
- **Citation:** Kelly, B. F. J., et al. (2022), *Atmos. Chem. Phys.*
- **In repo:** `Other References/Isotope/Kelly2022ACP.md`

---

## 10. Summary Table: Latitude-Resolved δD-CH₄ Source Signatures

Best-estimate values from Douglas et al. (2021) Table 1:

| Source category | Latitude band | δD-CH₄ (‰ VSMOW) | ± (1σ) | Flux (Tg/yr) | Primary reference |
|----------------|---------------|---------------------|--------|---------------|-------------------|
| Wetlands | <30°N | −301 | 15 | 115 | Douglas 2021 |
| Wetlands | 30–60°N | −324 | 14 | 25 | Douglas 2021 |
| Wetlands | >60°N | **−374** | 10 | 9 | Douglas 2021 |
| Inland waters | <30°N | −301 | 12 | 80 | Douglas 2021 |
| Inland waters | 30–60°N | −308 | 18 | 64 | Douglas 2021 |
| Inland waters | >60°N | **−347** | 9 | 16 | Douglas 2021 |
| Rice paddies | <30°N | −324 | 8 | 19 | Douglas 2021 |
| Rice paddies | 30–60°N | −325 | 8 | 12 | Douglas 2021 |
| Enteric ferm. | Global | −308 | 28 | 111 | Sherwood 2017 |
| Landfills | Global | −297 | 6 | 65 | Sherwood 2017 |
| Coal | Global | −232 | 5 | 42 | Sherwood 2017 |
| Oil & gas | Global | −189 | 2 | 79 | Sherwood 2017 |
| Biomass burning | Global | −211 | 15 | 17 | Sherwood 2017 |
| Termites | Global | −343 | 50 | 9 | Sherwood 2017 |
| Permafrost | Global | −374 | 15 | 1 | Douglas 2021 |

**Key takeaway for KIE_sites:** At our high-latitude sites (ALT 82°N, ZEP 79°N, BRW 71°N), the relevant wetland δD is **−374 ± 10‰**, not our assumed global mean of **−310 ± 30‰**. This means the source–atmosphere gap is 374 − 86 = **288‰**, not 224‰ — a 29% increase that would substantially change our source decomposition results.

---

## Papers Not in Repo — Priority for Future Reading

| Priority | Reference | Topic | Why needed |
|----------|-----------|-------|-----------|
| ⭐⭐⭐ | Waldron et al. (1999) *GCA* | Original δD-CH₄ vs δD-H₂O regression | Foundational — first global relationship |
| ⭐⭐⭐ | Chanton et al. (2006) *JGR* | Alaska wetlands, slope = 0.54 | High-latitude data that flattened the regression |
| ⭐⭐ | Whiticar et al. (1986) *GCA* | Isotope classification diagram | Foundational for δ¹³C–δD systematics |
| ⭐⭐ | Whiticar (1999) *Chem. Geol.* | Comprehensive δD fractionation review | Key fractionation factors |
| ⭐⭐ | Bock et al. (2010) *Science* | Ice core δD, paleo source signatures | Tropical vs boreal δD in paleo context |
| ⭐⭐ | Rice et al. (2016) *PNAS* | Atmospheric δD box model 1977–2005 | Used constant −322‰, now known to be wrong |
| ⭐⭐ | Snover & Quay (2000) *GBC* | Soil oxidation δD fractionation | αD ≈ 1.066 for soil uptake |
| ⭐⭐ | Valentine et al. (2004) *GCA* | H₂ thermodynamics and δD | Mechanistic pathway fractionation |
| ⭐ | Bowen & Revenaugh (2003) *WRR* | OIPC precipitation δD model | Underlying water δD latitude gradient |
| ⭐ | Fisher et al. (2011) *GRL* | Arctic atmospheric δ¹³C constraint | High-lat atmospheric perspective |
| ⭐ | Hornibrook (2009) *AGU Mono.* | Peatland isotope systematics | Bog vs fen δD differences |
| ⭐ | Sugimoto & Wada (1995) *GCA* | Rice paddy δD fractionation | Pathway-dependent δD in rice |
