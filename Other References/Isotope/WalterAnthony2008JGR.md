
--- Page 1 ---
Constraining spatial variability of methane ebullition seeps
in thermokarst lakes using point process models
Katey M. Walter Anthony1 and Peter Anthony 1
Received 28 August 2012; revised 14 June 2013; accepted 16 June 2013; published 19 July 2013.
[1]
Ebullition is an important but highly heterogeneous mode of methane emission in lakes.
Variability in both spatial distribution and temporal ﬂux creates difﬁculty in constraining
uncertainties in whole lake emission estimates. Analysis of short- and long-term ﬂux
measurements on 162 ebullition seeps in 24 panarctic lakes conﬁrmed that seep classes,
identiﬁed a priori according to bubble patterns in winter lake ice, have distinct associated
ﬂuxes irrespective of lake or region. To understand the drivers of ebullition’s spatial
variability and uncover ways to better quantify ebullition in ﬁeld work, we combined
point-process modeling with ﬁeld measurements of 2679 GPS-marked and classiﬁed
ebullition seeps in three Alaskan thermokarst (thaw) lakes that varied by region, permafrost
type, and seep distribution. Spatial analysis of ﬁeld data revealed that seeps cluster above
thawed permafrost soil mounds in lake bottoms. Seep density and clustering, determined
from ﬁeld observations, were used as parameters in a Poisson cluster process model to
simulate seeps across entire lake surfaces. Sampling results indicated that (1) applying
seep-class mean ﬂux values to unmeasured seeps counted on ice-bubble surveys does not
compromise accuracy of whole lake ﬂux estimates; (2) three distributed 50 m2 ice-bubble
survey transects more accurately estimate mean lake ebullition than 17 dispersed 0.2 m2
bubble traps; and (3) the uncertainty associated with whole lake mean ebullition estimated
by lake-ice survey transects is inversely related to seep density. Findings suggest that
transect ﬁeld data collected on a large number of widely distributed lakes can be combined
to provide a well-constrained, bottom-up estimate of regional lake ebullition.
Citation: Walter Anthony, K. M., and P. Anthony (2013), Constraining spatial variability of methane ebullition seeps in
thermokarst lakes using point process models, J. Geophys. Res. Biogeosci., 118, 1015–1034, doi:10.1002/jgrg.20087.
1.
Introduction
[2] Methane cycling in aquatic ecosystems involves
microbial methanogenesis in anaerobic sediments, diffusion
of aqueous methane through sediments and water, and
aerobic and anaerobic methane oxidation [Valentine and
Reeburgh, 2000; Valentine et al., 2001, Liikanen et al.,
2002]. Other important components of methane cycling and
transport include bubble formation and migration in sedi-
ments, bubble release at the sediment-water interface, partial
or full dissolution of bubbles during their rise through the
water column, and potential bubble release at the water
surface to the atmosphere [Martens and Klump, 1980;
Chanton et al., 1989; Joyce and Jewell, 2003; Boudreau
et al., 2001, 2005; McGinnis et al., 2006; Scandella et al.,
2011a, 2011b].
[3] In many lakes, ebullition (bubbling) is the dominant
mode of methane emissions [Crill et al., 1988; Keller and
Stallard, 1994; Casper et al., 2000; Bastviken et al., 2004,
2011; Walter et al., 2006; Del Sontro et al., 2010, 2011;
Schubert et al., 2012]. However, spatial and temporal
heterogeneity in the ebullition process gives rise to large
uncertainties in estimating ecosystem methane emissions
[Del Sontro et al., 2010; Wik et al., 2011a]. For instance,
ebullition is episodic and not representatively captured by
typical short-term measurements [Bastviken et al., 2004;
Varadharajan, 2009], a problem that is frequently encoun-
tered in sampling of rare ecological events [Thompson, 2012].
[4] In many northern lakes, the spatial clustering and
temporal variability in the ﬂux of discrete bubbling point
sources, herein called “ebullition seeps,” contribute to
heterogeneity in lake ebullition. While ebullition seeps
release bubbles episodically, with frequencies of several
minutes to weeks depending on the seep type, atmospheric
pressure dynamic, and season of year, the location of many
seeps in lakes with dense sediments remains stable over time
scales of at least seasons to years. Walter Anthony et al.
[2010] marked the location of 17 discrete ebullition seeps
in Siberian and Alaskan lakes and found that the majority
of seeps did not change in location over observation periods
of 2 and 8 years. Temporal consistency in the location of
ebullition seeps is likely the result of formation of fractures
Additional supporting information may be found in the online version of
this article.
1Water and Environmental Research Center, University of Alaska
Fairbanks, Fairbanks, Alaska, USA.
Corresponding author: K. M. Walter Anthony, Water and Environmental
Research Center, University of Alaska Fairbanks, 306 Tanana Loop, 525
Duckering Bldg., Fairbanks, AK 99775, USA. (kmwalteranthony@alaska.
edu)
©2013. American Geophysical Union. All Rights Reserved.
2169-8953/13/10.1002/jgrg.20087
1015
JOURNAL OF GEOPHYSICAL RESEARCH: BIOGEOSCIENCES, VOL. 118, 1015–1034, doi:10.1002/jgrg.20087, 2013


--- Page 2 ---
or “bubble tubes” in the ﬁne-grain sediments, which have
been described previously in other ﬁeld studies [Martens
et al., 1980.; Martens and Klump, 1980; Boudreau et al.,
2005;
Scandella
et al., 2011a] and lab experiments
[Scandella et al., 2011b]. When a sufﬁcient volume of gas
is produced or when the hydrostatic pressure drops enough
to dislocate large bubbles in sediments, the bubbles break
out creating preferential ﬂow channels in the sediments
[Varadharajan, 2009]. Once these fractures are formed,
escape of bubbles continues through the bubble tubes,
resulting in consistency in the location of bubble seepage
from the lake bottom.
[5] During ice-free conditions, methane-rich bubbles
released from lake sediments constitute direct emissions to
the atmosphere in lakes where water depth and bubble size
allow for bubbles to survive their ascent to the lake surface
[McGinnis et al., 2006], which was the case for all of the
>100 panarctic lakes where ebullition bubbles have been
recently investigated [Wik et al. 2011a, 2011b; Brosius
et al., 2012; Walter Anthony et al., 2010, 2012] and 20 lakes
in East Antarctica [Sasaki et al., 2009]. When lake ice forms,
rising bubbles are trapped by the bottom of the lake-ice sheet.
Throughout most of the winter, methane-rich bubbles freeze
into the ice, becoming sealed in the ice when the downward
growing lake ice thickens around them [Gow and Langston,
1977]. Walter et al. [2006] used distinctions in the patterns
that seep bubbles form in winter lake ice to classify discrete
ebullition seeps into four types: (A) kotenok, (B) koshka,
(C) kotara, and Hotspots. These relatively large ebullition
bubbles (1cm to >100 cm), rich in methane (up to 95%)
are distinct from tiny (millimeter scale), tubular ice-bubbles
with lower methane concentrations (usually <1%) that result
from the outfreezing of dissolved gases [Boereboom et al.,
2012]. A-type ebullition seeps are relatively small clusters
of bubbles in which individual bubbles stack on top of each
other in ice without merging laterally (Figure 1). Due to
relatively higher ebullition rates, individual bubbles released
from sediments in B-type seeps laterally merge into larger
bubbles under the ice prior to freezing in ice. Types A and
B seeps produce low gas-volume clusters of bubbles in lake
B
C
A
Hotspot
Seep 
class
Description
Stacks of small individual 
bubbles, >50% unmerged; 
usually <40 cm diameter
Bubbles >50% merged, 
clustered in multiple 
layers of ice; usually 
>25 cm diameter
Single large pockets of 
near 100% merged bubbles 
stacked in ice; occasional, 
isolated bubbles around 
<25% of perimeter; usually 
>40 cm diameter
Ice-free hole in lake ice 
due to frequent bubbling. 
Thin snow or ice crust, 
easily broken with fist or 
foot, often covers the hole.
Mean ± s.e.m., 
n seeps
22  ± 4, n=30
211  ± 39, n=42
1,726  ± 685, n=35
7,801  ± 764, n=55
Ebullition (ml gas d-1)
Min-max
2-96
(0-47)
32-1,081
(0-1,030)
525-2,649
(0-24,640)
3,110-16,918
(1,313-25,412)
Figure 1.
Ebullition seep classes (A, B, C, and Hotspot) distinguished according to visible patterns of
bubbles trapped in lake ice and the associated bubbling rates that generate these patterns. Error term repre-
sents the standard error of the mean (s.e.m.). The minimum and maximum daily ebullition rates are shown
for long-term measured seeps (no parentheses) and for short-term measured seeps (inside parentheses). The
top left photograph shows a 1 m wide bubble survey transect on lake ice cleared of snow. The diameters of
the example A, B, C-type bubble clusters in the photographs, as shown by the white lines, are 20 cm, 25 cm,
and 45 cm respectively. The open-hole Hotspot has a 30 cm diameter. This classiﬁcation scheme, presented
previously by Walter et al. [2006] and Walter Anthony et al. [2010], is shown here with additional new ﬂux
data from 20 lakes and greater detail in the distinguishing seep characteristics.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1016
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 3 ---
ice (Figure 1), whereas larger C seeps result in big pockets of
gas in ice separated vertically by ice layers containing few or
no bubbles. The solid ice layers in between the large gas
pockets of C-type seeps represent periods of relative quies-
cence in between large ebullition events. The frequency of
ebullition release from Hotspot seeps and the associated
convection in the water column created by rising bubble
plumes can be strong enough to maintain mostly ice-free
holes in winter lake ice (Figure 1). Not all lakes have
ebullition Hotspots. In thermokarst (thaw) lakes, Hotspot
seeps are typically concentrated along margins of most active
permafrost thaw [Walter et al., 2006], where labile organic
matter is made abundantly available to methanogenesis by
deep permafrost thaw in taliks [Kessler et al., 2012]. The
mean annual methane ebullition efﬂuxes from sediments
associated with each of the seep types were recently reported
based on year-round, long-term continuous ﬂux measure-
ments of ebullition seeps in four Siberian and Alaskan
thermokarst lakes [Walter Anthony et al., 2010].
[6] While
some
progress
has
been
made
toward
constraining heterogeneity in ebullition seep classiﬁcation
and temporal efﬂuxes from sediments, the pattern of spatial
variability in ebullition is unresolved. Taking advantage of
trapped bubbles in lake ice, researchers follow a relatively
new method outlined by the Pan-Arctic Lake Ice Methane
Monitoring Network (PALIMMN, http://ine.uaf.edu/werc/
palimmn/) and described in detail by Walter Anthony et al.
[2010] to examine the spatial distribution of bubbling events
in lakes. In ice-bubble surveys, transects are placed randomly
in different zones of the lake (stratiﬁed randomization) to ac-
count for different shore-type processes, such as thermokarst
erosion, and other sources of geographic variability that
could inﬂuence ebullition within lakes. How representative
these transect-derived ﬁeld survey data are of the true
whole-lake ebullition ﬂux from sediments is unknown
because a very limited fraction of the lake is surveyed.
Researchers who used this method found that the spatial
variability of seeps within individual lakes was as high as
102–103% in Alaskan, Siberian, and Swedish lakes [Walter
Anthony et al., 2010; Wik et al., 2011a]. Confounding this
problem, stratiﬁed randomization in the placement of bubble
traps in summer resulted in 4–5 times lower emission esti-
mates than those produced by ice-bubble surveys [Walter
et al., 2006; Zimov et al., 1997; Wik et al. 2011a, 2011b].
Thus, the accuracy of each method for estimating ﬂux and
the causes of discrepancy between methods need to
be resolved.
[7] Due to the spatial and temporal heterogeneity in
ebullition efﬂux from sediments and the tendency for total
ebullition ﬂux to be dominated by events that are rare in
space and time, accurate estimation of ebullition faces
problems similar to those faced when sampling other rare
ecological
events.
These
problems
include
the
near-
impossibility of systematic ground sampling to observe
every ebullition seep and the tendency of classical random
sampling methods (e.g., stratiﬁed randomization in the
placement of bubble traps) to produce sampling distributions
that underestimate the true mean [Thompson, 2012].
[8] Since we cannot know the entire real population of
ebullition seeps on a lake, ﬁeld estimates of ebullition ﬂux
cannot be compared to a known mean lake ﬂux. In this
situation, stochastic simulation provides the best method for
evaluating the effectiveness of various ﬁeld methods for
estimating
ebullition
[Thompson, 2012]. Spatial
point
process models allow characterization of the spatial patterns
of ebullition. In ice-bubble ebullition surveys, seep locations
are identiﬁed in lakes using global positioning systems (GPS)
and classiﬁed as A, B, C, or Hotspot based on ice bubble
morphology (Figure 1). Each seep has a location and a
classiﬁcation, and so may be considered as a marked point
process [Cressie, 1993; Bivand et al., 2008]. Coupling point
process models of seep location and class with ﬁeld measure-
ments of seep ﬂux allows simulation of realistic populations
of ebullition seeps and simulation-based evaluation of
ebullition sampling strategies.
[9] The objectives of this study were to use spatially
explicit, individual seep ﬁeld data combined with point-
process models to describe and interpret spatial patterns of
ebullition in three intensively studied thermokarst lakes
located in different permafrost types and regions of Alaska,
and to evaluate ﬁeld methods for quantifying ebullition.
First we test the hypothesis that speciﬁc seep classes (A, B,
C, and Hotspot) identiﬁed according to ice-bubble patterns
differ by ebullition rate irrespective of lake type or region.
Then, we investigate the nature and causes of spatial hetero-
geneity in methane ebullition in thermokarst lakes, particu-
larly as they relate to permafrost thaw geomorphology. We
combine our understanding of the spatial pattern of ebullition
with other ﬁeld observations of bathymetry and extent of
thermokarst erosion along lake margins to simulate ebullition
across the whole lake surfaces with point process modeling.
Using iterated and randomized sampling of these simulated
ebullition data sets, we compare the accuracies of the ice-
bubble survey transects and the classical random sampling
by bubble traps methods for estimating whole-lake ebullition
release from sediments. By quantifying the errors associated
with transect-based lake ebullition surveys, we provide a
solution for reducing uncertainty in regional-scale lake
ebullition estimates based on limited ﬁeld data.
2.
Methods
2.1.
Site Description
[10] This study focused on three thermokarst lakes located
in different regions of Alaska (Figure 2) that were
unglaciated during the late Pleistocene [Kaufman and
Hopkins, 1986; Manley and Kaufman, 2002]. The lakes
differed in size, depth, and their topographic position on the
landscape. Lake depths varied according to the distribution
of excess ground ice in permafrost in which they formed.
[11] Lake Claudi (informal name; 66.55°N, 164.45°W;
162,750 m2), located in the continuous permafrost zone on
the northern Seward Peninsula, is a 10 m deep lake formed
in yedoma-type permafrost. Yedoma was recognized ﬁrst as
an extensive permafrost type in Russia [Birkengof, 1933;
Tomirdiaro, 1980; Kaplina, 1981], characterized by ice-su-
persaturated Pleistocene-aged loess and loess-like deposits.
Thick yedoma deposits in Siberia were reported to store large
quantities of organic carbon (≈450 Gt) [Zimov et al., 2006],
which can contribute to production and emission of methane
and carbon dioxide upon permafrost thaw [Zimov et al.,
1997]. Yedoma-type permafrost also occurs extensively in
Alaska [Kanevskiy et al., 2011] and has been described on
the northern Seward Peninsula as late-Pleistocene-aged
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1017
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 4 ---
aeolian silt with deep syngenetic ice in the form of ice-
wedges, lenses and pore ice that comprise 60–80% of the
substrate by volume [Hopkins et al., 1955, Hopkins and
Kidd, 1988; Hoeﬂe et al., 2000]. Permafrost depth on
the northern Seward Peninsula generally exceeds 90 m
[Hopkins et al., 1955]. Thermokarst expansion was observed
along all margins of Lake Claudi during recent decades;
however, the most rapid expansion occurred along the NE
and SW shores [Jones et al., 2009, 2011], as shown in
Figure 3a. The bathymetry, morphology, and numerical
modeling of Lake Claudi’s expansion, parameterized by
site-speciﬁc ﬁeld data, were presented recently by West and
Plug [2008] and Kessler et al. [2012].
[12] Goldstream Lake (informal name; 64.92°N, 147.85°
W; 10,030 m2; 2.9 m deep), located in the discontinuous
permafrost zone of interior Alaska near Fairbanks, formed
in retransported late-Quaternary loess that is common on
many hill slopes and valley bottoms of interior Alaska
[Muhs and Budahn, 2006]. Unlike the majority of the
yedoma lowlands in Siberia [Zimov et al., 2006] and the
northern Seward Peninsula [Hopkins et al., 1955, Hopkins
and Kidd, 1988], where permafrost formed syngenetically
(concurrent with deposition) during the late Pleistocene,
colluvial forces and frost action in interior Alaska gradually
eroded loess down slope during the late Pleistocene and early
Holocene, forming icy, organic-rich deposits described as
“muck,” frequently several tens of meters to 100 m thick in
valley bottoms [Péwé, 1975; Hamilton et al., 1988; Muhs
and Budahn, 2006; Reyes et al., 2010]. Massive ice wedges,
2 to 4 m wide at their tops are common [Hamilton et al.,
1988]. Based on their high ice content, large syngenetic ice
wedges, cryostructure, presence of mammoth fauna remains,
often high organic carbon contents (0.38–6.8% C) [Hamilton
et al., 1988], and occurrence of distinct thermokarst
landforms, Kanevskiy et al. [2011] classiﬁed these interior
Alaska deposits as “yedoma-type.” Remote-sensing based
observation of lake extent shows that a partial drainage event
occurred
in
Goldstream
Lake
sometime
after
1949
(Figure 3b); however, we observed that thermokarst expan-
sion continues presently, predominately along the eastern
margin of Goldstream Lake.
[13] Ikroavik Lake (71.23°N, 156.63°W; 5,292,850 m2),
the largest of the study lakes, has a maximum depth of 2.4
m and is located on the Barrow Peninsula on the Arctic
Coastal Plain of northern Alaska. The region is underlain
by continuous permafrost with maximum thickness in excess
of 400 m [Hinkel et al., 2003]. Like most of the lakes on the
Barrow Peninsula, Ikroavik Lake is elliptical (Figure 3c),
with the major axis oriented a few degrees west of due north
and nearly perpendicular to the prevailing summer wind
direction [Carson and Hussey, 1962; Sellmann et al.,
1975]. Barrow Peninsula lakes developed in ice-rich perma-
frost [Livingstone et al., 1958; Black, 1969]; however, their
shallow depth is in part the consequence of numerous
iterations of the thermokarst-lake cycle in the region
[Hinkel et al., 2003, 2007; Eisner et al., 2005] that limits
the volume and distribution of excess ground ice in which
thermokarst lakes form today [Jorgenson and Shur, 2007].
Very few remnants of thicker ice-rich units remain [Eisner
et al., 2005; Frohn et al., 2005]. The shores types of
Ikroavik Lake
Lake Claudi
Goldstream Lake
Barrow
Kotzebue
Fairbanks
N
0
250 km
Figure 2.
Map of Alaska showing the location of study lakes (yellow squares) on the northern Seward
Peninsula (Lake Claudi), interior Alaska (Goldstream Lake), and the coastal plain of the Alaska North
Slope (Ikroavik Lake). The locations of nearby communities at Barrow (71.29°N, 156.79°W), Kotzebue
(66.90°N, 162.60°W), and Fairbanks (64.88°N, 147.82°W) are shown as white dots. The Alaska map is
the National Elevation Data Set 30 m hillshade raster.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1018
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 5 ---
(f)
(e)
(b)
0
0
250 m
30 m
(c)
0
750 m
T1
T3
SQ4
T2
T6
T7
T5
T11
T7
T3
T8
T2
T1 T10
T5
T4
T9
T6
24
25
1
23
2
N
(a)
(d)
SQ12
SQ13
26
8
28
4
29
30
3
7
5
SQ27
6
10
9
11
12
13
17
14
31
15
16
18
19
20
22
21
32
Figure 3.
The locations of 2007–2011 ﬁeld survey transects and simulated ebullition seeps on Lake Claudi,
Goldstream Lake, and Ikroavik Lake. Field survey transects are drawn to scale overlying recent satellite
images on (a) Lake Claudi, 2008 pan-sharpened, multispectral IKONOS© image with a resolution of 1 m,
and on (b) Goldstream Lake, 27 August 2008 Worldview image with a resolution of 0.5 m. The approximate
locations of transects for which precise GPS positions were unavailable are shown as dashed lines (Claudi T5,
T7). Field survey transects are not to scale on (c) Ikroavik Lake, CIR photo, 2.5 m resolution; actual survey
areas are provided in Table 2. The previous shoreline positions, shown as red lines, were from 1955, 1949,
and 1951 black and white aerial photographs with ≈1 m resolution for Lake Claudi, Goldstream Lake, and
Ikroavik Lake, respectively, based on Jones et al. [2009]. Baydjarakhs (conical thermokarst mounds), now
covered in shrub vegetation are seen across the bottom of a drained thermokarst lake basin adjacent to
Lake Claudi in the bottom left corner of Figure 3a. Baydjarakh spacing of approximately 10 m demonstrates
the morphology of the Pleistocene yedoma ice-wedge network. Baydjarakhs also cover the bottom of Lake
Claudi (Figure 7). (d–f) Lake Claudisim, Goldstreamsim Lake, and Ikroaviksim Lake simulated by point
process modeling using the ratios, distributions, and clustering of ebullition seep types A (green dots), B
(yellow dots), C (orange dots), and Hotspot (red dots) from ﬁeld observations are also shown. The randomly
placed transects and polygons on simulated lakes are drawn to scale.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1019
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 6 ---
Ikroavik Lake were described by Wohlschlag [1953] as
consisting of gravel that grades into sand at the northeastern
beach, a 1.5 m tundra bluff along the western margin that is
formed by ice movement in spring, and swampy shoreline
around the remainder of the lake. Expansion of the lake,
caused by wind and wave action and thermal erosion
occurred predominately along the eastern and southeastern
margins (Figure 3c).
2.2.
Field Sampling and Analysis
[14] The three lakes described above were surveyed
intensively for spatial bubble patterns. In addition to tran-
sects, we surveyed large polygons, referred to as “squares”
in these lakes, which allowed quantiﬁcation of statistical
parameters for spatial modeling the lakes. Extensive ice-
bubble surveys were conducted on lakes in the northern
Brooks Range region, near Toolik Field Station, Alaska.
Bubble-trap ﬂux measurements on 162 seeps from 24
panarctic lakes were used to scale the bubble surveys to
estimates of ﬂux. We also collected and radiocarbon dated
bubble gas to help infer its source.
2.2.1.
Ice-Bubble Surveys to Estimate Whole-Lake
Seep Ebullition
[15] We followed a stratiﬁed randomization sampling
design for surveying ebullition seeps in winter lake ice at
Lake Claudi, Goldstream Lake, and Ikroavik Lake. In
November 2008, we cleared snow from lake ice in a 639
m2 polygonal plot, called “Claudi Square,” offshore at the
northern margin of Lake Claudi. Each of the 1076 seeps in
the plot was categorized as A, B, C, or Hotspot using the seep
classiﬁcation method shown in Figure 1 and was mapped
using a centimeter-accuracy RTK differential GPS (Leica
Geosystems AG). As covariate data, we mapped the bathym-
etry of the lake bed with an intensity of ≈0.75 measurements
per m2 using sonar (Vexilar LPS-1 Hand-Held Depth Finder)
point measurements through ice. Approximately square
Table 1. Seep Density and Estimated Ebullition Flux from GPS Field Surveys During 2007–2011 on Lake Claudi, Goldstream Lake and
Ikroavik Lake
Seep Density (seeps m2)
Percent
of Lake
Surveyed
Ebullition
on
Transects
Whole-
Lake
Ebullition
Field Date
Transect
Area (m2)
A
B
C
HS
All Seeps
(mL gas m2 d1)
Lake Claudi
0.57%
Whole lake
1.43
0.28
0.04
0
1.74
152 ± 38
3 Nov 2008
T1
62
2.72
0.47
0.10
0
3.29
326
3 Nov 2008
T2
33
0.88
0.40
0
0
1.28
103
6 Nov 2008
T3
37
1.02
0.43
0
0
1.45
112
4 Nov 2008
SQ4
639
1.68
0.31
0.04
0
2.03
175
19 Apr 2009
T5
50
0
0
0
0
0
0
21 Apr 2009
T6
50
0.02
0.04
0
0
0.06
9
21 Apr 2009
T7
50
0.02
0
0
0
0.02
0
Goldstream Lake
10.74%
Thermokarst zone
0.88
0.42
0.14
0.09
1.53
1,064
Nonthermokarst zone
0.25
0.04
0.02
0.00
0.31
75
Whole lake, zone-weighted
0.36
0.11
0.04
0.02
0.53
0
253 ± 82
Whole lake, transects only
0.54
0.18
0.04
0.02
0.78
0
254 ± 62
12 Oct 2007
T1
71
0.40
0.43
0.06
0.07
0.95
750
12 Oct 2007
T2
87
0.01
0.05
0
0.01
0.07
99
12 Oct 2007
T3
61
0.13
0.05
0.02
0.02
0.21
169
20 Oct 2009
T4
50
0.94
0.24
0.08
0
1.25
208
20 Oct 2009
T5
58
0.70
0.17
0.03
0.02
0.92
244
20 Oct 2009
T6
68
1.36
0.22
0.01
0
1.60
102
26 Oct 2009
T7
30
0.43
0.10
0.03
0
0.57
88
26 Oct 2009
T8
50
0.67
0.14
0.04
0
0.85
112
19 Oct 2010
T9
48
0.84
0.23
0
0.04
1.11
393
19 Oct 2010
T10
38
0.72
0.37
0.11
0.03
1.22
485
20 Oct 2010
T11
54
0.04
0.06
0.07
0
0.17
140
29 Oct 2011
SQ12
225
1.03
0.41
0.17
0.10
1.71
1,162
30 Oct 2012
SQ13
236
0.30
0.01
0.00
0
0.31
16
Ikroavik Lake
0.03%
Whole lake
0.07
0.01
0
0
0.07
0
2.8 ± 1.3
28–29 Oct 2009
T1 and T2
37
0
0
0
0
0
0
27 Oct 2011
T3 - T22
703
0
0
0
0
0
0
29 Oct 2009
T23
72
0.01
0
0
0
0.01
0
28 Oct 2009
T24
34
0.09
0.06
0
0
0.15
14
28 Oct 2009
T25
33
0.12
0.12
0
0
0.24
28
25 Oct 2011
T26
136
0.01
0
0
0
0.01
0
27 Oct 2011
SQ27
276
0.03
0
0
0
0.03
1
26 Oct 2011
T28
161
0.03
0
0
0
0.03
1
26 Oct 2011
T29
361
0.27
0.02
0
0
0.29
9
27 Oct 2011
T30
26
0.04
0
0
0
0.04
1
27 Oct 2011
T31
19
0.05
0
0
0
0.05
1
27 Oct 2011
T32
17
0.06
0
0
0
0.06
1
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1020
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 7 ---
polygons were selected so that seep observations would be
present across a range of x and y axes values; polygon bound-
aries were also identiﬁed with RTK GPS for use in spatial
analysis. Ebullition seeps within six additional 1 m wide
survey transects that varied in length from 33 to 62 m
(Table 1) were mapped in Lake Claudi in either November
2008 or April 2009, the lake ice sheet of a single winter
season. Since lake ice freezes from the top down, surveying
ebullition seep bubbles in the top 30 to 40 cm of lake ice in
April reveals the ebullition patterns that formed in the fall
of the year prior. Therefore, the November 2008 and April
2009 data sets can be directly compared as spatial data sets
of contemporary ebullition activity. Transects in Lake
Claudi were placed so as to originate from different shore
types (virgin yedoma, T1; nonvirgin yedoma, T2), off-
shore (T3, SQ4), and in the central portions of the basin
(T5-7) (Figure 3a).
[16] In Goldstream Lake, we mapped the location of
ebullition seeps along eleven 1 m wide transects at a
frequency of three to ﬁve transects per year in October
2007, October 2009, and October 2010 (Figure 3b) using
Garmin 76CSx GPS with Wide Area Augmentation System
differential correction. Transects varied by length and seep
density (Table 1). Given minimal spatial overlap among
transects and observations of interannual seep-position
stability, spatial information from all transects in all years
was pooled into a single data set. In October 2011 we
mapped 303 ebullition seeps within two large polygons using
a
centimeter-accuracy
RTK
differential
GPS
(Leica
Geosystems AG). Goldstream Thermokarst Square (SQ12,
225 m2) was located an average of 7 m off shore along the
eastern
thermokarst
shore
of
the
lake
(Figure
3b).
Goldstream Center Square (236 m2) was located outside the
thermokarst zone in the approximate center of Goldstream
Lake, where the density of methane ebullition seeps was
signiﬁcantly lower. As covariate data, we mapped the
bathymetry of the lake bed with a 0.22 m2 point measurement
density on both Goldstream squares.
[17] We classiﬁed and mapped ebullition seeps with a
Garmin 76CSx GPS along ﬁve transects in Ikroavik Lake on
28–29 October 2009. In October 2011 we surveyed ebulli-
tion seeps and mapped bathymetry in 27 additional survey
plots, including Ikroavik Square, a large 276 m2 polygon.
Survey plots were randomly located across the lake, but
constrained to areas where wet snow had not frozen to the
lake-ice surface. We used centimeter-accuracy RTK differen-
tial GPS (Leica Geosystems AG) for 2011 surveys.
Locations of the 32 ﬁeld survey plots on Ikroavik Lake are
shown in Figure 3c.
[18] To relate ﬁeld observations of ebullition-seep cluster-
ing in lake ice to patterned permafrost ground, we mapped
the distribution of ice wedge polygons and baydjarakhs [con-
ical thermokarst mounds] along the eastern margin of
Goldstream Lake and along the margins of Lake Claudi in
summer. Field work at Ikroavik Lake took place only in
winter, at a time when snow cover inhibited accurate
mapping of ice wedge polygons in the surrounding tundra.
[19] Finally, we surveyed ebullition along 64 transects in
13 other lakes in the northern foothills of the Brooks Range
near Toolik Field Station (68.4°N, 149.4°W). The number
of randomly placed, dispersed transects per lake ranged from
two to nine, with median and mean values of four and ﬁve
transects per lake, respectively. These transects will be used
as an example of how to apply the results from the spatial
surveys and error estimation from the modeling.
2.2.2.
Seep Gas Collection and Radiocarbon Dating
of Methane
[20] We collected gas samples of ebullition bubbles from
49 seeps at the lake surface in Goldstream Lake during
2008–2011and Lake Claudi during a 2 week expedition in
Table 2. Geographic Information for 24 Lakes on Which Ebullition Flux Measurements Were Conducteda
Region
Lake
Permafrost
Ecology
Latitude
Longitude
Depth (m)
Number of Flux-Measured
Seeps
Alaska, Arctic Coastal Plain
Cake Eater*
continuous
tundra
71.279
156.637
1.8
8
Alaska, Arctic Coastal Plain
Ikroavik*
continuous
tundra
71.247
156.641
2.4
5
Alaska, Seward Peninsula
Rhonda
continuous
tundra
66.566
164.466
1.5
12
Alaska, Seward Peninsula
Cocker Gap
continuous
tundra
66.562
164.451
5.4
4
Alaska, Seward Peninsula
Fox Den
continuous
tundra
66.559
164.456
2.4
2
Alaska, Seward Peninsula
Claudi
continuous
tundra
66.552
164.450
11.0
7
Alaska, Seward Peninsula
Kim
continuous
tundra
66.516
164.251
3.3
10
Alaska, Seward Peninsula
Jaeger
continuous
tundra
66.501
164.260
13.3
19
Alaska, Seward Peninsula
Three Loon
continuous
tundra
66.497
164.254
9.0
5
Alaska, Interior
Vault
discontinuous
boreal forest
65.029
147.699
5.0
2
Alaska, Interior
Goldstream
discontinuous
boreal forest
64.916
147.848
2.9
28
Alaska, Interior
Doughnut
discontinuous
boreal forest
64.899
147.910
2.5
3
Alaska, Interior
Killarney
discontinuous
boreal forest
64.870
147.902
2.0
8
Alaska, Interior
Smith*
discontinuous
boreal forest
64.865
147.866
4.0
1
Alaska, Interior
Deuce*
discontinuous
boreal forest
64.864
147.940
6.0
3
Alaska, Interior
Stevens Pond
discontinuous
boreal forest
64.863
147.871
>0.8
2
Alaska, Interior
Ace*
discontinuous
boreal forest
64.862
147.936
9.0
8
Alaska, Interior
Rosie Creek beaver pond
discontinuous
boreal forest
64.770
148.079
3.9
2
Alaska, Southcentral
Tern
discontinuous
tundra
63.398
148.670
12.5
4
Canada, Yukon Territories
Dredge Pond
discontinuous
boreal forest
64.035
139.281
2.0
1
West Greenland
G11
continuous
tundra
67.056
50.441
>1
1
Russia, Kolyma Lowlands
Grass
continuous
boreal forest
68.752
161.377
12.0
5
Russia, Kolyma Lowlands
Shuchi*
continuous
boreal forest
68.746
161.393
11.0
20
Russia, Kolyma Lowlands
Tube Dispenser*
continuous
boreal forest
68.764
161.403
17.0
2
aLake names are informal, except those formal names indicated by asterisk.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1021
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 8 ---
October 2008 using submerged bubble traps. We followed
methods described previously for gas collection and labora-
tory analyses [Walter Anthony et al., 2012; Stuiver and
Polach, 1977], including bubble methane concentrations
and radiocarbon dating methane collected from two ice-
bubble pockets above baydjarakhs on Claudi Square.
2.2.3.
Ebullition Flux Measurements Using Submerged
Bubble Traps
[21] We measured ebullition ﬂux on short-term or long-term
time scales on 162 individual seeps classiﬁed a priori as A, B,
C, or Hotspot based on ice-bubble patterns (Figure 1) on 24
lakes located in Siberia, Alaska, northwest Canada, and west
Greenland (Table 2). Twenty-ﬁve percent of the 162 panarctic
seeps were located in our intensive study lakes, Lake Claudi
(seven seeps), Goldstream Lake (28 seeps), and Ikroavik
Lake (ﬁve seeps). This effort to quantify seep ebullition in
panarctic lakes is a component of the Pan Arctic Lake Ice
Methane Monitoring Network (PALIMMN). Flux data from
four of the 24 lakes reported on here were included in previous
publications: Three lakes in Siberia and one lake in Alaska
[Walter et al., 2006; Walter Anthony et al., 2010]. This study
reports new long-term and short-term ebullition ﬂux data from
20 additional lakes in Alaska, Canada, and Greenland.
Altogether this data set consists of ≈213,600 individual
ebullition ﬂux measurements made using submerged bubble
traps on 162 seeps in 24 panarctic lakes. These data represent
rates of bubble release from sediments, measured at the lake
surface during multiple seasons of the year. Since winter ice
impedes the majority of ebullition bubbles from reaching the
atmosphere during the ice cover season, the rates presented
here do not represent seasonal ebullition ﬂuxes to the
atmosphere. The latter requires consideration of gas exchange
between lake water and bubbles trapped under ice and lake-
ice phenology.
[22] Prior to measuring ﬂux, seeps were identiﬁed as
bubbles in ice and classiﬁed based on bubble morphology
(Figure 1). Then the ice was chipped away in order to place
a bubble trap directly over the identiﬁed seep. Short-term
measurements consisted of measuring the volume of gas
accumulated in submerged bubble traps (corrected for hydro-
static pressure) placed over seeps for periods of ≈20 min to
up to 6.1 days. The mean and median periods of short-term
measurements were 1.6 and 1.0 days per seep, respectively.
In most cases, short-term ﬂux measurements were made only
once per seep, but several seeps were measured twice over
the short term.
[23] Long-term ﬂux measurements consisted of continu-
ously
measuring
the
volume
of gas
accumulated
in
submerged traps per unit of time (minutes to weeks); how-
ever, these traps were left in place year-round and secured
over discrete seeps for periods of up to 359 days in
Siberian lakes and up to 700 days in Alaskan lakes. In the
Alaskan lakes, semi-automated bubble traps equipped with
wet-cup gas meters and event data loggers, similar to those
used in other ﬁelds of study [Massé et al., 1997], were placed
at various depths below the lake surface. We corrected bub-
ble volume data for hydrostatic pressure and lake water
temperature. These bubble traps allowed for the volumetric
measurement of all gas released from individual seeps over
the long period of measurement.
[24] We used unbalanced, type III two-way analysis of
variance (ANOVA) to test differences in measured ebullition
among the A, B, C, and Hotspot by lake, region, and
measurement period (short versus long) using the statistical
software SAS Institute. Degrees of freedom were adjusted
to account for any empty cells in an analysis.
2.3.
Spatial Statistical Modeling of Ebullition in Lakes
[25] We used data collected from the large survey squares
to test three hypotheses for spatial patterns of seeps. Lake
bed morphology was tested as a possible explanatory vari-
able for determining seep density. We used point-process
models to model seep locations regardless of explanatory
variables. Results of point process models were used to
create simulations of bubble occurrence in three lakes of
the same dimensions as the three intensively studied lakes,
Claudi, Goldstream and Ikroavik. Details and results of the
spatial statistical modeling are available in Figures S1, S2,
and S3 in the supporting information.
2.3.1.
Describing Spatial Patterns
[26] We used marked point process models to parameterize
the spatial distribution of ebullition seeps in lakes. Marked
point processes are widely used in ecological research to
analyze data in which the researcher is interested both in
the location of an object and also some attribute, or mark
[Cressie, 1993]. Here we are interested in the location of
seeps and the ebullition class, or mark, of each seep. Three
possible formulations exist for a marked point process
[Baddeley, 2008]: (1) seeps are equally likely to be located
anywhere in a region, and neither seeps of the same class
nor seeps of different classes exhibit the interactions of
attraction or repulsion (Figure 4a); (2) seeps of different clas-
ses are located through independent spatial processes
(Figure 4b); or (3) seeps of all classes are located through a
single spatial process, and seeps classes are conditionally
independent and identically distributed among locations
(Figure 4c). Determination of the correct formulation is
essential for accurate simulation of point processes [Illian
et al., 2008]. To determine which formulation was appropri-
ate for evaluation of ebullition seeps in thermokarst lakes, we
tested three null hypotheses [Baddeley, 2008]:
[27] H01: Seeps do not occur in clusters across the lake, and
seep class is random (Figure 4a).
[28] H02: The seep classes A, B, C, and Hotspot all follow
independent patterns of spatial distribution, and if clustering
occurs, the cluster patterns are different among seep classes
(Figure 4b).
[29] H03:
Ebullition seeps occur in clusters, but within
clusters seep class is random (Figure 4c).
[30] To determine which type of point process (H01, H02, or
H03) thermokarst-lake ebullition seeps follow, we analyzed
spatial seep patterns in the ﬁeld data sets of Claudi Square,
Goldstream Thermokarst Square, Goldstream Center Square,
and Ikroavik Square. We determined clustering (an estimate
of the spatial dependence between points at a range of scales)
by generating K-functions derived from the DGPS locations
of seeps [Baddeley, 2008; Van Lieshout and Baddeley,
1999]. All analyses were conducted using the spatstat
package [Baddeley and Turner, 2005] in R version 2.12.1
[R Development Core Team, 2010].
2.3.2.
Inhomogeneous Poisson Models
[31] Seeps may be considered to occur at higher densities
in some zones of a lake due to peculiar properties of that
zone. When information about these properties is known,
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1022
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 9 ---
inhomogeneous Poisson models may be used to describe
seep intensity as a function of one or more of these covari-
ates. To test the hypothesis that seep location is related to lake
bed morphology, we ﬁt the inhomogeneous Poisson model
with seep intensity modeled as a loglinear function of the
lake bed morphology covariate. Models were ﬁt to large
square data, where we had identiﬁed seep locations and had
measured lake bed morphology. We tested lakebed slope
and depth as potential morphology covariates. We ﬁt the
inhomogeneous Poisson process model using the Berman-
Turner algorithm [Berman and Turner, 1992], implemented
in the ppm procedure from the spatstat package in R. The ﬁt
of our model was determined using the likelihood ratio test
and the Akaike Information Criterion (AIC) [Baddeley, 2008].
2.3.3.
Point Process Models
[32] Seep locations were also modeled using a Poisson
cluster process model. Poisson cluster processes model the
clustering and intensity of points as an intrinsic part of their
behavior, without explanatory variables. This model allowed
us to describe and simulate seep clustering across whole
lakes where precise information on lake bed morphology
and other potential covariates is unknown. We ﬁt a Thomas
process variant of the generalized Poisson cluster process to
the data from each of Claudi Square, the two Goldstream
Squares, and Ikroavik Square using the kppm procedure from
the spatstat package in R. The ﬁt of the model was tested
using simulation to produce critical envelope values with a
95% signiﬁcance level for the model [Baddeley, 2008].
2.3.4.
Lake Simulation
[33] To test sampling strategies for estimating ebullition in
ﬁeld-based approaches, we created simulated lakes in which
the location and efﬂux from sediments of every ebullition
seep was exactly known. To cover a range of thermokarst
lake types in Alaska, we simulated seep ebullition across
the surfaces of three lakes, producing Lake Claudisim,
Goldstreamsim Lake, and Ikroaviksim Lake (Figures 3d–3f)
with the same areas and shapes as the lakes studied in ﬁeld
work. Each simulated lake consisted
of the known,
georeferenced lake surface populated with simulated seeps.
We used lake-speciﬁc ﬁeld survey data from squares and
transects to determine seep density and survey data from
squares to determine spatial parameters for Thomas process
models of seep clustering on each lake. We then used these
lake-speciﬁc models to simulate seep ebullition across the
entire surface of each lake.
[34] We created zones of differing seep intensities on lakes
according to generalization of ﬁeld observations. For
Goldstreamsim Lake, we deﬁned two zones of differing seep
intensity (Figure 3e): a high-intensity thermokarst zone
within ≈20 m wide belt along the eastern margin and a lower
intensity nonthermokarst zone, which was the rest of the
lake, using a seep intensity transition line observed in ice-
bubble ﬁeld surveys and early winter aerial photos of
Goldstream Lake. In aerial photos, the density of open-hole
Hotspots drops off sharply at this boundary, which was
approximately twice the distance from the 1949 shoreline
to the modern shoreline. If expansion rates in the past were
similar to that observed from 1949 to the present, these
observations suggest that enhanced methane production and
emission occur in response to thermokarst in Goldstream
Lake on the order of 100–120 years following the
thermokarst conversion of land to lake. Zone-speciﬁc ﬁeld
survey data of seep class densities and ratios collected on
squares and transects were used to create Thomas process
models for simulating seeps separately by zone.
[35] At Lake Claudi, thermokarst expansion is observed in
all directions (Figure 3a). For Lake Claudisim, we deﬁned
seep intensity as a continuous variable proportional to the
distance from the center of the lake (Figure 3b), consistent
with our ﬁeld observations of seeps [this study] and modeling
observations of methanogenesis [Kessler et al., 2012].
[36] At Ikroavik Lake, ebullition seep density was low
across the whole lake surface. The lake surface was consid-
ered one uniform zone for seep intensity on Ikroaviksim Lake.
[37] After simulating seep locations, we assigned a
simulated ebullition ﬂux value (SEF, mL seep1 d1) to each
simulated seep in each lake. From the PALIMMN ﬁeld data
set of ﬂux from 162 lake ebullition seeps, we determined that
seep ﬂux follows a lognormal distribution (Figure 5a). To
assign SEF values, we created random lognormal distribu-
tions for each class of simulated seeps on each simulated
lake, pooled these distributions by lake while maintaining a
class assignment, and randomly resampled them to produce
lake-speciﬁc ﬂux distributions for each simulated lake
(a)
(b)
(c)
Figure 4.
Examples of three hypothetical models for marked point patterns. (a) A point pattern that
follows complete spatial randomness (CSR), with marks assigned randomly to points. (b) A point pattern
in which two types of points (open and ﬁlled) are generated according to separate, independent spatial point
processes. (c) A point pattern in which all points are ﬁrst generated according to a single spatial point
process, and then marks are assigned to points through a random mechanism.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1023
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 10 ---
(Figure 5b–5d). Simulated ﬂuxes were randomly attached to
simulated seep locations.
[38] In addition to having a unique SEF value, each seep
was also assigned a PALIMMN ﬂux value, which was one
of the four mean ﬂuxes for each seep class (A, B, C, and
Hotspot) in Figure 1. PALIMMN values were assigned based
on which seep class the SEF value fell into. This means that
on simulated lakes, each seep had two ebullition values in
units of mL seep1 d1 assigned to it: (1) the PALIMMN
class average (Figure 1) and (2) a unique and random SEF
value from the lognormal distribution of seep ﬂuxes that
occur in lakes. Seep class assignment was retained in order
to compare sampling methodologies using mean ﬂux by
bubble class (PALIMMN) with methodologies measuring
ebullition ﬂux on each seep (SEF).
[39] We also created a separate set of 26 simulated lakes to
test for the effects of variation in the percentage of lake
surveyed, seep density, and clustering patterns on the accu-
racy and precision of lake-ice ebullition surveys. Eight of
these simulated lakes varied in area with constant seep
density and clustering pattern; nine varied in seep density
with constant area and clustering pattern; and nine varied in
clustering pattern with constant area and seep density
(Table 3). The same bubble class ratios and ebullition ﬂux
distributions were assigned to seeps in all of these 26
simulated lakes using the methodology described above.
[40] We
performed
all
lake
simulations
in
R
[R
Development
Core
Team,
2010],
using
the
rThomas
procedure from the R spatstat package [Baddeley and
Turner, 2005] to create seep locations and using the rlnorm
procedure in R to create the seep ﬂuxes.
2.4.
Sampling Evaluation
[41] Knowledge of the exact location, class, and associated
ﬂux of every seep on a simulated lake allowed us to test sam-
pling strategies commonly used in ﬁeld work. Speciﬁcally, we
compared the simulated lake ﬂuxes with estimated ﬂuxes from
Table 3. Lake Areas, Seep Densities, and Clustering Parameters
for 26 Simulated Lakes
Parameter
Lake
Radius
Area
Seep Density
(seeps m2)
Clustering
Parameters
(m)
(m2)
λ
κ
σ
μ
Lake area
1
30
2,815
0.30
0.041 0.758 7.46
2
50
7,833
0.30
0.041 0.758 7.46
3
70
15,364
0.30
0.041 0.758 7.46
4
100
31,375
0.30
0.041 0.758 7.46
5
200
125,581
0.30
0.041 0.758 7.46
6
400
502,488
0.30
0.041 0.758 7.46
7
800
2,010,285
0.30
0.041 0.758 7.46
8
1600 8,041,807
0.30
0.041 0.758 7.46
Seep density
1
200
125,581
0.02
0.041 0.758 0.50
2
200
125,581
0.04
0.041 0.758 1.00
3
200
125,581
0.08
0.041 0.758 2.00
4
200
125,581
0.12
0.041 0.758 3.00
5
200
125,581
0.19
0.041 0.758 5.00
6
200
125,581
0.33
0.041 0.758 8.00
7
200
125,581
0.52
0.041 0.758 13.00
8
200
125,581
0.80
0.041 0.758 20.00
9
200
125,581
2.06
0.041 0.758 50.00
Clustering pattern
1
200
125,581
0.50
0.010 0.758 50.00
2
200
125,581
0.50
0.020 0.758 25.00
3
200
125,581
0.50
0.040 0.758 12.50
4
200
125,581
0.50
0.060 0.758 8.33
5
200
125,581
0.50
0.080 0.758 6.25
6
200
125,581
0.50
0.100 0.758 5.00
7
200
125,581
0.50
0.150 0.758 3.33
8
200
125,581
0.50
0.200 0.758 2.50
9
200
125,581
0.50
0.300 0.758 1.67
Figure 5.
The distribution of ebullition ﬂux of 162 ﬁeld-
measured seeps weighted by the ratios of seep classes
observed among 16,364 seeps individually classiﬁed and
mapped on (a) 75 Alaskan lakes [Walter Anthony et al.,
2012] and simulated seeps on (b) Lake Claudisim, (c)
Goldstreamsim Lake, and (d) Ikroaviksim Lake. Colored dots
represent seep types: A, green; B, yellow; C, orange; and
Hotspot, red. Ebullition ﬂux (y axis) followed a lognormal
distribution for measured seeps and seeps on simulated lakes,
but the number of seeps and ratios of seep classes on each
simulated lake differed in accordance with ﬁeld observations.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1024
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 11 ---
both transect and ﬂoating bubble-trap methods to determine
the error associated with these two sampling approaches.
2.4.1.
Errors Associated With Estimating Ebullition
From Transects on Simulated Lakes
[42] Simulated transects were placed on simulated lakes
following
the
ﬁeld
survey
methodology
of
Walter
Anthony et al. [2010], whereby transects are evenly dis-
tributed and randomly placed along different shorelines
and different zones, including the approximate centers
of lakes. In the set of 26 simulated lakes, transect origins
were randomly placed on the shoreline or in the lake
center. Each shoreline transect was oriented toward the
center of the lake, and all transects were deﬁned as 50
m in length and 1 m in width.
[43] Evaluation of transects was conducted in ArcGIS 9.3
[ESRI]. Simulated seeps were converted to vector ﬁle formats
using the maptools package [Lewin-Koh and Bivand, 2011]
in R. Seeps were intersected with transects in ArcGIS using
the Intersect tool. For each iteration, we recorded the SEF
ﬂux (mL m2 d1), based on the sum of the precise simulated
ebullition ﬂux (SEF) of simulated seeps on transects divided
by transect area. To empirically determine errors associated
with the transect survey method, we calculated the mean
ebullition ﬂux for the simulated transects, calculated the
transect error as [E = P(predicted ﬂux from transects)  A(actual whole
lake mean ﬂux)], and determined the variance of the errors based
on 250 iterations of random transect placement for each lake.
We estimated the standard deviation of the errors as the
square root of the variance. We completed this process for
combinations of n = (1, 2, 3,…, 9, and 30) transects from each
of the three transect zones on each simulated lake. We
conducted analysis of variance to determine signiﬁcance of
density, clustering, area of lake surveyed and the number of
distributed transects surveyed on the relative error of
whole-lake ebullition estimates.
2.4.2.
Comparison of Floating Bubble-Trap
and Transect Methodologies
[44] We evaluated and compared tethered, ﬂoating bubble
traps with lake-ice transects as methodologies for estimating
ebullition. For evaluation of the ﬂoating bubble-trap method,
we conducted 1000 iterations of random placement of seven-
teen 0.2 m2 circular, simulated ﬂoating bubble traps on
Claudisim, Goldstreamsim, and Ikroaviksim lakes. For each it-
eration, we determined the total seep ebullition ﬂux captured
by the 17 ﬂoating bubble traps. For evaluation of the transect
method, we conducted 1000 iterations of random placement
of three 50 by 1 m transects on each of the simulated lakes.
The Goldstreamsim Lake shoreline was divided into one
thermokarst and two nonthermokarst segments, and one
transect origin was randomly placed in each segment. Each
simulated transect ran perpendicular to the shoreline toward
the lake center. For Ikroaviksim Lake, where simulated seep
density was uniform across the lake, the shoreline was
divided into three segments, and one transect origin was
randomly placed in each segment. For Claudisim Lake, where
simulated seep density was increased from the lake center
to the margins, one transect origin was randomly placed
within each of two shoreline segments, and one transect
origin was placed in the center of the lake, oriented with
a random heading. For each iteration, we determined both
the SEF and PALIMMN ebullition captured by the three
transects as follows:
[45] 1. SEF ﬂux (mL m2 d1) was the sum of the precise
simulated ebullition ﬂux (SEF) of simulated seeps on
transects divided by transect area.
[46] 2. PALIMMN ﬂux (mL m2 d1) was the sum of the
class-speciﬁc mean ebullition rates (PALIMMN, Figure 1) of
simulated seeps on transects divided by transect area.
[47] For Goldstreamsim Lake only, we repeated this
process using 5 by 10 m polygons instead of 1 by 50 m
transects (Figure 3e) to compare the spatial errors associ-
ated with the alternative survey method of 50 m2 polygons
instead of 50 m2 transects. Results from each iteration
were recorded and used to construct sampling distributions
for each method [Thompson, 2012]. Simulation of method-
ologies was conducted using Python scripting in ArcGIS
9.3 [ESRI]. We compared resulting sampling distributions
using the two-sample Kolomogorov-Smirnov test in R [R
Development Core Team, 2010].
3.
Results and Discussion
3.1.
Ice Bubble Patterns Represent Statistically
Distinct Fluxes
[48] Field-measured ebullition rates for 162 seeps in 24
panarctic lakes identiﬁed a priori according to ice-bubble
pattern (A, B, C, and Hotspot) varied by seep class
(ANOVA, F = 4.7246,115, p < 0.0001; Figure 1). Interaction
effects of seep class by lake or seep class by region were
not signiﬁcant. This indicates that the four different ice-
bubble patterns represent statistically distinct mean ebullition
ﬂuxes irrespective of the lake or region in which they are
measured. The range and magnitude of emissions covered
by the four seep classes in Figure 1 are also consistent with
observations in other lakes and reservoirs [Del Sontro
et al., 2010, 2011; Ramos et al., 2006; Schubert et al., 2012].
[49] While measurement period (long versus short) also
did not affect the mean ebullition rates by seep class at
α = 0.05 level, the range of ebullition values for individual
seeps was generally larger among the short-term measured
seeps (Figure 6). Short-term ebullition is strongly controlled
by hydrostatic pressure dynamics [Martens and Klump,
1980;
Mattson
and
Likens,
1990;
Fechner-Levy
and
Hemond, 1996; Scandella et al., 2011a], so measured ebulli-
tion rates vary by orders of magnitude over time scales of
hours to days [Glaser et al., 2004; Walter Anthony et al.,
2010; Goodrich, 2010; Goodrich et al., 2011]. In contrast,
ebullition values reported for the long-term monitored seeps
are the average of a large number of shorter-term ﬂux
measurements on individual seeps (up to ≈61,500 measure-
ments per seep), so the high and low ﬂuxes from individual
seeps measured over the long term converge on a more
representative ebullition rate for individual seeps. While a
single short-term measurement of a seep is highly unlikely
to represent the long-term mean ebullition rate from that
seep, the mean ebullition of a large number of seeps
measured over the short term at different times and in
different locations was not different from the mean ebullition
of different seeps measured over the long term (Figure 6).
Two caveats should be highlighted.
First, short-term
measurements of A and B seeps were almost entirely
conducted during the early winter ice-cover period following
removal of overlaying ice since these low-ﬂux seeps are
easiest to locate in lakes when seasonal ice cover reveals their
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1025
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 12 ---
locations as bubbles trapped in lake ice. Since low-ﬂux A and
B seeps originate from relatively shallow lake sediments
[Walter et al., 2008], cold winter temperatures in shallow
sediments slow rates of methanogenesis [Zimov et al.,
2001; Liikanen et al., 2002; Metje and Frenzel, 2005; Duc
et al., 2010], effectively lowering the seasonal efﬂux from
sediments from these seeps. The higher mean ebullition rate
reported for long-term monitored A and B seeps takes into
account summertime ebullition from these seep classes,
which was higher due to warmer temperatures in surface lake
sediments in summer [Walter et al., 2006]. Similar responses
to seasonal changes in solar bed warming and cooling have
been described for other temperature-dependent physical
and microbial processes, such as gas solubility and anaerobic
oxidation of methane coupled to sulfate reduction in shallow
marine sediments [Dale et al., 2008]. Second, while the
relative error (standard deviation divided by the mean) was
generally under one for A, B, and Hotspot seep classes,
regardless of long versus short measurement periods, mea-
surement period had a strong inﬂuence on relative error in
the C-class ebullition seeps. Long-term relative error of C-
class seeps was 0.4, while short-term relative error was 2.4.
This large discrepancy is most likely explained by the high-
amplitude ebullition dynamic characteristic of C-type seeps.
Long-term monitoring of C-type seeps revealed that C-seeps
emit bubbles at high rates typical of Hotspots over short
periods of time, but they also have long (days to weeks)
periods of quiescence. The lack of bubbling over these long
periods allows ice to thicken over C-seeps in winter.
Temporary interruptions of high-ﬂux ebullition in C-seeps
results in large, frequently isolated pockets of gas that stack
in between the ice layers. As with all seep classes, long-term
monitoring of C-seeps provides representative long-term
ebullition values for these seeps; however, since short-term
C-seep ebullition can vary from zero to 104 mL gas d1,
the standard deviation of short-term C-seep measurements
was high.
[50] We considered it appropriate to use the pooled ﬂux
data set consisting of short-term and long-term monitored
seeps to (1) populate seeps on simulated lakes with
PALIMMN ebullition values, and (2) during ﬁeld work,
apply the PALIMMN ebullition values to unmeasured seeps
identiﬁed by ice bubble patterns in survey transects for the
following reasons: (a) The mean ﬂuxes of A, B, C, and
Hotspots seeps were not statistically different between the
short-term and long-term data sets; (b) the pooled data set
represents a much larger sample size of seeps within each
seep class; and (c) seep classes had consistent ebullition rates
among a large number of widely distributed arctic lakes.
3.2.
Seep Clustering on Thermokarst Lakes
[51] Our tests of H01 and H02 rejected these two null hy-
potheses at α = 0.05 for Claudi Square, the two Goldstream
squares and Ikroavik Square (See Figures S1 and S2 in the
supporting information). We found (H1): ebullition seeps
exhibit signiﬁcant clustering for all seep classes; and (H2):
the clustering patterns were the same for classes A, B, C
and Hotspot. Our test of H03 failed to reject this null hypoth-
esis at α = 0.05 at all three lakes (See Figure S3 in the
supporting information), so we conclude that given the
locations of seeps, seep classes were conditionally indepen-
dent and identically distributed. In other words, seep density
is dependent on location, but seep type is not. The clustering
of
seeps
across
lake
surfaces
indicates
that
within
thermokarst lakes there are regular patterns of methane
production and transport of bubbles to the sediment-water
Seep ebullition (ml gas d-1)
Long
Short
Long
Short
Long
Short
Long
Short
(A)
(B)
(C)
Hotspot
Measurement period by seep class
1
10
100
1,000
10,000
100,000
Figure 6.
Ebullition ﬂuxes measured on 162 ebullition seeps (open symbols) pre-classiﬁed according to
ice-bubble patterns as A, B, C, and Hotspot. Fluxes were measured over long (up to 700 days) and short (up
to 6 days) periods for different seeps. Box plots show median, maximum, minimum, ﬁrst quartile, and third
quartile values.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1026
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 13 ---
interface and then to the atmosphere. This conclusion
allowed us to model seep location irrespective of ﬂux and
class; on simulated lakes, it allowed us to ﬁrst simulate seep
locations and then to randomly assign ﬂux to seeps
irrespective of location.
3.3.
The Spatial Relationship of Ebullition Seeps
and Permafrost Thaw
[52] Clustering patterns of seeps in lakes matched the
pattern of polygonal ground observed in the permafrost-
dominated regions. In ﬁeld work in the tundra and boreal
yedoma permafrost study regions, we observed that dense
clusters of methane seeps were distributed ≈10 m apart
across lake ice surfaces. The ice wedge polygons exhibited
identical ≈10 to 12 m spacing. As sedimentation and soil
formation progress in the polygon centers over millennia,
surrounding ice wedges syngenetically grow both vertically
with the sedimentation and in width by repeated frost crack-
ing and reﬁlling with vein ice. If a disturbance to the ground
surface occurs that causes permafrost to thaw, ice wedges
melt and the soils of the polygon center remain as slowly
eroding conical thermokarst mounds known as baydjarakhs
[Czudek and Demek, 1970; Kanevskiy et al., 2011].
Baydjarakhs form a regularly rough topography, similar in
morphology to the bottom of an egg carton, but with the
spacing and distribution of the former ice wedge polygons.
We observed 10 to 12 m spacing of subaerial baydjarakhs
near lake margins. Bathymetric maps created by sonar depth
measurements beneath Claudi Square (Figure 7), Goldstream
Thermokarst Square (not shown), and other linear survey
transects in these lakes revealed that this ≈10 m spacing of
baydjarakhs continued beneath the lake, on the yedoma-lake
bottoms.
[53] Were methane bubbles seeping out of the tops of these
baydjarakh mounds, or from the troughs in between mounds?
In most lakes, sediments and organic detritus are focused
toward the deep centers of basins [Lehman, 1975]. In
thermokarst lakes, sediment focusing occurs at multiple
0
50
100
150
(a)
(b)
Water depth (m)
Percent Slope (ΔY/ΔX * 100)
0
-8.5
0
40
0
25
20
Figure 7.
Methane ebullition seep locations from ﬁeld surveys (a) overlaying slope ratio and (b) bathymetry
of Claudi Square lake bed. The location of seeps in lake ice mapped with DGPS is shown as open circles in
Figure 7a and as blue dots in Figure 7b. In Figure 7b, lake-bed depth grades from white (2.5 m)
to brown, yellow, green, and turquoise (8.5 m). Methane seeps were disproportionately associated
with the slopes of baydjarakhs, the cone-shaped mounds of thawed permafrost soil distributed among
the low-elevation troughs formed by melting of massive ice wedges beneath the thermokarst lake.
Baydjarakhs contain the organic matter in thawed permafrost that fuels the majority methane production
in yedoma-type thermokarst lakes.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1027
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 14 ---
spatial scales, both toward the deep central basin of the lake
[Hopkins and Kidd, 1988; Murton, 1996], but also toward
the topographical low points where ice wedges melted
forming troughs in between baydjarakhs [personal observa-
tion from lake coring and paleolake exposure sampling,
unpublished data; Farquharson, 2012]. We used spatial
statistics and radiocarbon dating to determine whether seeps
were associated with baydjarakhs (consisting of organic-rich
thawed permafrost soil) or troughs (topographically low
points for focusing organic-rich surface sediments), both
potential sources of organic matter to fuel methanogenesis.
[54] The inhomogeneous Poisson model indicated that
lake bed morphology was highly associated with seep loca-
tion. Seep intensity in the ﬁeld data sets was signiﬁcantly
related to baydjarakh slope (χ2 test, p < 0.0001). Adding
depth to the model did not improve the model ﬁt. Seep inten-
sity was highest near ﬂatter tops of baydjarakhs and lowest
on the steep baydjarakh sideslopes (Figure 7). Baydjarakh
tops are very broad (6–8 m), while the trough bottoms are
quite narrow (0.5 m), so that >90% of ebullition occurs from
the tops of baydjarakhs. Radiocarbon ages of methane
collected from bubbles trapped in ice in the dense bubble
clusters above two separate baydjarakhs closest to the
thermokarst margin in Claudi Square were 40,400 and
41,300 years.
The
association
of
seep
clusters
with
baydjarakhs and the radiocarbon ages of methane that
matched those of the Pleistocene-aged organic matter in
yedoma permafrost are both pieces of evidence that support
the hypothesis that permafrost thaw fuels methane produc-
tion in thermokarst lakes, particularly lakes formed in icy
yedoma permafrost.
[55] Sonar surveys of Goldstream Lake and Lake Claudi
revealed consistent, lake-wide baydjarakh patterns in the lake
bed morphology. The signiﬁcant relationship between lake
bed morphology and seep density allowed us to apply the
seep clustering observed in polygon samples to the entire
lake surfaces.
3.4.
Modeling Seeps as a Point Process to Quantify
Uncertainty in Lake Emission Estimates
[56] Because we rejected H01 and H02, but failed to reject
H03, we modeled ebullition seep location as a single point
process and randomly assigned seep classes based on ﬁeld
observations of class ﬂux distributions (Figure 5) once seep
locations were determined. The Thomas process model also
provided excellent ﬁt of the seep location data over seep
separation distances of 0 to 8 m, indicating that the clustering
of seep locations could be modeled without explanatory data
from lake bed covariates. The resulting distribution of seeps
in the simulated lakes, Lake Claudisim, Goldstreamsim Lake,
and Ikroaviksim Lake (Figures 3d–3f), were consistent with
our ﬁeld observations of ebullition in the actual lakes;
however, unlike in the actual lakes, the location and ﬂux of
all ebullition seeps on simulated lakes was known exactly.
[57] Knowing the precise simulated ebullition ﬂux (SEF)
value for each individual seep on the 26 simulated lakes
(Table 3), we could distribute different numbers of transects
and polygons across simulated lakes to determine the error
associated with estimating whole-lake methane ebullition
based on limited spatial sampling. Increasing the number of
transects from one transect per zone to nine transects per zone
reduced the relative error associated with estimating whole
lake ebullition at the 95% conﬁdence level within all lakes
(Figure 8). Among lakes, contrary to our initial expectation,
differences in the uncertainty of ebullition estimates using
the survey transect method were not related to the fraction
of the lake surveyed or to the seep clustering parameter
(Figure 8 and Table 4). Rather, we found that uncertainty in
upscaling transect data to whole-lake ebullition ﬂux values
Relative error (%)
Clustering parameter (K)
Fraction of lake surveyed (%)
100
10
1
0.1
0.01
0.001
0
20
40
60
80
100
120
0
20
40
60
80
100
0.3
0.25
0.2
0.15
0.1
0.05
0
Seep density (# seeps m-2)
2
1.5
1
0.5
0
0
100
200
300
400
500
a
b
c
Figure 8.
The relative error (%) of the mean whole-lake
ebullition estimate based on one (open symbols) and nine
(closed symbols) transects per sector randomly distributed
on 26 simulated lakes that varied according to (a) seep
density, (b) seep clustering, and lake size, which translates
to fraction of the lake surveyed based on 50 m2 (1 by 50 m)
transects. Simulated ebullition ﬂuxes (SEF) were unique
ebullition values assigned to each seep randomly on the 26
simulated lakes following the distribution of ﬂuxes measured
on 162 seeps weighted by the ratios of these classes as they
occurred among 16,364 seeps individually classiﬁed and
mapped on 75 Alaskan lakes (Figure 5a) [Walter Anthony
et al., 2012].
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1028
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 15 ---
was inversely related to seep density (Figure 8a, R2 = 0.99).
Equation (1) deﬁnes the relative error RE at the 95%
conﬁdence level associated with the mean lake ebullition
estimate as a function of seep density D and the number of
distributed survey transects T:
RE ¼ 0:6866D0:481T0:51
(1)
[58] On lakes with low seep density, the probability was
low that an individual transect would intersect a representa-
tive number of seeps. In contrast, on lakes with high seep
density, the probability was much higher that an individual
transect would intersect a representative number of seeps.
Increasing the number of transects surveyed reduced uncer-
tainty in whole-lake ebullition estimates by increasing the
likelihood that transects represent the spatial heterogeneity
in the density of ebullition seeps across lakes (Table 5).
This concept is highlighted by the results shown in
Figure 8c: multiple dispersed transects result in lower errors
than a single transect, even when the fraction of the lake
surveyed is the same. These results are consistent with the
ﬁndings of Wik et al. [2011a], who also found that larger
numbers of widely distributed ebullition survey transects
reduce uncertainty in seep ebullition ﬂux estimates on two
lakes in northern Sweden.
3.5.
Comparison of Ice-Bubble Survey Transects
Versus Polygons
[59] In the comparison of an alternative survey method of
equal area, using 5 by 10 m polygons as opposed to 1 by
50 m linear transects, we found that surveying polygons on
Goldstreamsim nearly doubled the relative error associated
with the mean ﬂux estimate of lake ebullition. The high
variance associated with polygon surveys resulted from
similarity in the scale among polygons, seep clustering, and
baydjarakhs. While a 50 m long transect was likely to
intersect several baydjarakhs and associated seep clusters,
polygons were apt to fall on or between baydjarakhs and seep
clusters, giving rise to large variances. These results imply
that for ﬁeld surveys, the distributed 50 m transect method
is preferred to more nearly isometric polygons of equal area.
Long, continuous transects are more likely to cover decame-
ter-scale heterogeneity in seep distribution within and among
zones in the lake. Transects are also logistically easier to
survey in the ﬁeld because walking along a 1 m wide strip
is more efﬁcient and reduces the likelihood of missing or
double counting seeps, a challenge that is present in survey-
ing polygons.
3.6.
Comparison of Ice-Bubble Survey Transects
and Bubble-Trap Methods
[60] In the comparison of alternative survey methods for
estimating seep ebullition, using 0.2 m2 ﬂoating bubble traps
as opposed to 1 by 50 m linear transects, we found that 0.2 m2
ﬂoating bubble traps resulted in substantial underestimation
of the mean ﬂux. The magnitude of this underestimation
ranged from 1.5-fold to >100-fold, depending on seep
density. For Goldstreamsim Lake, where seep density was
near the middle of the observed range, ﬂoating bubble traps
produced a median estimate of whole-lake ebullition that
was 11.5 times less than the actual mean whole-lake ebulli-
tion and 10.9 times less than the estimate produced by the
linear transects. The classic bubble-trap method also
underestimated mean lake ﬂux in 63–90% of occasions
(Figure 9), also depending on seep density. The sampling
distribution from the ﬂoating bubble-trap method exhibited
extreme deviation from normality, indicating substantial
median-bias and associated difﬁculties with statistical infer-
ence from this method if it were used in the ﬁeld. In contrast,
the linear transect method underestimated mean lake
ebullition by a much smaller amount (5–34%) and produced
sampling distributions that approximated normality. We
attribute the poor performance of the ﬂoating bubble-trap
method to well-described problems encountered in sampling
of rare events [Thompson, 2012]. Given that the mean seep
density on the simulated lakes was 0.5 seeps m2 and that
each ﬂoating bubble trap covered 0.2 m2, there existed a
0.9 probability that individual bubble traps would capture
no ebullition seeps. As a result, the sampling distribution of
the ﬂoating bubble traps was highly skewed toward underes-
timation of ebullition ﬂux. This magnitude of underestima-
tion for ﬂoating bubble traps is similar to that observed by
Wik et al. [2011a, 2011b] and Walter et al. [2006] in
comparisons of these two methods in the ﬁeld.
Table 5. The Relative Error (RE, Percent of the Mean) for the 95%
Conﬁdence Level of Mean Lake Ebullition Based on Equation (1)a
Seep Density
Number of Transects
(# seeps m2)
3
6
9
15
21
27
90
0.01
343.8
234.8
187.9
141.9
117.9
102.7
53.0
0.03
202.7
138.4
110.8
83.6
69.5
60.5
31.2
0.05
158.5
108.3
86.6
65.4
54.4
47.3
24.4
0.07
134.8
92.1
73.7
55.6
46.2
40.3
20.8
0.09
119.5
81.6
65.3
49.3
41.0
35.7
18.4
0.2
81.4
55.6
44.5
33.6
27.9
24.3
12.5
0.4
58.3
39.8
31.9
24.1
20.0
17.4
9.0
0.6
48.0
32.8
26.2
19.8
16.5
14.3
7.4
0.8
41.8
28.5
22.8
17.2
14.3
12.5
6.4
1
37.5
25.6
20.5
15.5
12.9
11.2
5.8
1.2
34.4
23.5
18.8
14.2
11.8
10.3
5.3
1.4
31.9
21.8
17.4
13.2
10.9
9.5
4.9
1.6
29.9
20.4
16.4
12.4
10.3
8.9
4.6
1.8
28.3
19.3
15.5
11.7
9.7
8.4
4.4
2
26.9
18.4
14.7
11.1
9.2
8.0
4.1
aWhere RE is a function of seep density and number of distributed
transects per lake surveyed. The transect dimensions used to statistically
derive these relationships were 1 m × 50 m (50 m2), parameters that follow
the PALIMMN protocol [Walter Anthony et al., 2010] and capture the within-
zone variability in ebullition controlled by lake bed morphology (10 m scale;
section 3.3).
Table 4. Results of ANOVA for Prediction of Relative Error in
Whole-Lake Ebullition Estimation
Variable and Source of Variation
d.f.
MS
F
P
Seep density
1
11.3
9594.27
<0.001
Seep clustering
1
0.0
2.39
0.124
Percentage of lake surveyed
1
0.0
0.21
0.645
Number of transects
1
5.9
5048.06
<0.001
Seep density × number of transects
1
1.8
1568.10
<0.001
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1029
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 16 ---
Figure 9.
Sampling distributions for estimates of lake ebullition obtained through 1,000 iterations of
three sampling methodologies using (a–c) 17 tethered 0.2 m2 bubble traps and seep-speciﬁc ﬂux [Wik
et al., 2011a, 2011b]; (d–f) three distributed 50 m transects and seep-speciﬁc ﬂux, SEF; and (g–i) three dis-
tributed 50 m transects and PALIMMN-estimated ﬂux by class (Figure 1), applied to simulated ebullition
seeps on three lakes of different seep density (Claudisim, 1.23 seeps m2; Goldstreamsim, 0.53 seeps m2;
and Ikroaviksim, 0.03 seeps m2). Gray bars are a histogram of the sampling distribution of estimated ﬂux.
Solid black line indicates the cumulative sampling curve. Vertical dashed black line indicates the actual
mean ﬂux for the simulated lake; horizontal dashed black line indicates the intersection of actual mean ﬂux
with the sampling distribution. Vertical dashed red line indicates the median estimated ﬂux for the sampling
distribution. The bubble-trap method produced median ﬂux estimates that were 1.5, 11.5, and >100 times
lower than actual mean lake ﬂuxes on Claudisim, Goldstreamsim, and Ikroaviksim, respectively, and
underestimated actual mean ﬂux in 64, 80, and 90% of occasions for these lakes. The sampling distributions
of the bubble trap method exhibited extreme departure from normality and high median-bias toward under-
estimation of ﬂux. Ice-bubble survey transect methods based on SEF and PALIMMN seep values produced
median ﬂux estimates that were 1.09, 1.06, and 1.6 times lower than actual ﬂux for lakes Claudisim,
Goldstreamsim, and Ikroaviksim, respectively, and underestimated actual ﬂux in 60, 55, and 66% of
occasions for these lakes. The sampling distributions of ice-bubble survey transect methods were much
closer to normal distributions than those of the bubble trap method. All methods produced more narrowly
distributed sampling estimates with increasing seep density. The ice-bubble survey transect methods
produced results that were substantially more accurate than those of the bubble-trap method. Differences
between results based on SEF and PALIMMN transect methods were not signiﬁcant.
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1030
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 17 ---
3.7.
Comparison of Ice-Bubble Transect Estimates
Based on Class-Mean Fluxes With Exact Fluxes
[61] There
were
no
signiﬁcant
differences
between
sampling distributions based on applying simulated seep ﬂux
values (SEF) along transects versus applying the mean
PALIMMN ﬂux for seep classes on lakes Claudisim and
Goldstreamsim (Figure 9). On the very low seep-density
lake, Ikroaviksim, three transects intersected very few seeps
(mode = 4), resulting in a bimodal sampling distribution
when mean PALIMMN ﬂux values were applied. There
was no indication from values of the mean, median, or
variance that this resulted in decreased accuracy. When
sampling intensity on Ikroaviksim was increased to nine
transects, we observed no signiﬁcant differences between
SEF and PALIMMN sampling distributions (data not
shown). Because there is considerable variability in seep
ﬂux within class (Figures 5 and 6), applying the individual,
seep-speciﬁc SEF values induced as much error as it
resolved. This is a highly important result given that
determining the true ﬂux of individual seeps along transects
is logistically impractical in extensive ﬁeld work. The
implication of this ﬁnding is that PALIMMN ﬂux values
can be applied to seeps by class during ﬁeld surveys without
affecting the accuracy of the whole-lake mean seep ebulli-
tion estimate.
3.8.
Quantifying Uncertainty in Field-Based Lake
Ebullition Estimates
[62] In the past, upscaling spatially limited transect-level
ﬁeld data to whole-lake ebullition estimates was inhibited
by a lack of knowledge of the relative errors associated with
mean ﬂux estimates from transects. In this study, spatial
analysis of randomly placed transects across a variety of
simulated lakes, for which whole-lake ebullition was
precisely known, revealed the critical relationship between
the number of distributed survey transects and the density
of seeps on a lake to constrain the relative errors associated
with mean ebullition estimates from transects (section 3.4
and Figure 8a). Based on this relationship (equation (1)),
the relative error values predicted for different combinations
of seep densities and distributed 50 m2 survey transects are
provided in Table 5.
[63] If there is no prior knowledge about the areal distribu-
tion of processes that control ebullition in a lake, then survey
transects should be randomly distributed in lakes, and
ebullition estimated as the cumulative seepage from all
surveyed
areas
divided
by
the
total
area
surveyed.
Following this approach at Ikroavik Lake, the total amount
of ebullition observed among the 32 survey transects was
5227 mL d1. Dividing this ﬂux by the total area surveyed
(1873 m2) yielded a mean ﬂux of 2.8 mL m2 d1. To deter-
mine the uncertainty of this mean ﬂux estimate, we used the
mean density of seeps on Ikroavik (0.07 seeps m2) and the
number of transects surveyed in the lake (32 transects) in
equation
(1)
to
calculate
the
relative
error
(37%).
Propagating this spatial sampling error together with the
ﬂux-based relative error (26%) associated with seep-class
mean ebullition rates from Figure 1, yielded a mean whole-
lake seep ebullition estimate for Ikroavik Lake, 2.8 ± 1.3
mL m2 d1 at the 95% conﬁdence level. Following the
same approach for Lake Claudi, where the mean seep density
was 1.74 seeps m2 along 7 transects, the relative error on the
mean estimate of whole-lake ebullition was 25%, yielding a
whole-lake mean ebullition estimate of 152 ± 38 mL m2
d1 at the 95% conﬁdence level.
[64] In Goldstream Lake, since 11 transects were represen-
tatively distributed across two, well-constrained zones
(thermokarst and nonthermokarst), we determined the rela-
tive error (24%) based on the propagation of errors associated
with the mean seep density, 0.78 seeps m2 using equation
(1) and the ﬂux-based errors (Figure 1). This approach,
which excluded the two large surveyed polygons in the
thermokarst and nonthermokarst zones, yielded a whole-lake
ebullition estimate of 254 ± 62 mL m2 d1. Separately, we
used a weighted mean for only those transects that were fully
in the thermokarst or nonthermokarst zones together with the
large polygons that were fully in these zones, considering
that these zones occupied 18% and 82% of the lake area,
respectively. This weighted approach yielded a mean
whole-lake ebullition estimate, 253 ± 82 mL m2 d1, similar
to that of the transect-only method.
[65] Using equation (1) and the ﬂux-based errors of each
seep class to estimate the relative error in the mean estimate
of whole winter ﬂux for two lakes in northern Sweden stud-
ied by Wik et al. [2011a] (mean ﬂux 87 mg CH4 m2 d1,
mean seep density 1.63 seeps m2, surveyed area equivalent
to 5.5 1 m × 50 m transects), we calculate an estimate of
87 ± 25 mg CH4 m2 d1 at the 95% conﬁdence level.
Propagating the variances presented by Wik et al. [2011a]
to calculate a relative error in the estimate of the mean [Ku,
1966], we arrive at 87 ± 33 mg CH4 m2 d1 at the 95% con-
ﬁdence level. The close agreement between our simulation-
based error estimates and the estimates derived from Wik
et al. [2011a] provides corroborating evidence that ice-bubble
transect surveys yield well-constrained estimates of mean lake
ebullition on the Swedish study lakes.
3.9.
Mass-Based Methane Ebullition From Seeps
[66] To convert volumetric ebullition to a massed based
estimate of methane emissions to the atmosphere, the
methane concentration in ebullition bubbles at the surface
of the lake must be known. On Goldstream Lake, methane
concentrations (percent by volume) based on 246 measure-
ments of 42 ebullition seeps (type A = 82 ± 3%, n = 6; type
B = 83 ± 7%,
n = 3;
type
C = 85 ± 1%,
n = 14;
type
Hotspot = 89 ± 1%, n = 19; reported as mean ± standard error,
n number of seeps) were higher than those reported as the
mean from a much larger set of arctic and subarctic
PALIMMN lakes (A, 73%; B, 75%; C, 76%; and Hotspot,
78%) [Walter Anthony et al., 2010]. Applying these site-
speciﬁc methane concentrations measured in bubbles at the
surface of Goldstream Lake to the volumetric ebullition
survey data yielded a massed-based ﬂux estimate for the
whole lake, 157 ± 51 mg CH4 m2 d1 at the 95% conﬁdence
level. At Lake Claudi, where the bubble gas methane concen-
tration [B = 72 ± 4%, n = 4; C = 73 ± 4%, n = 3] was similar to
that of the PALIMMN values [Walter Anthony et al., 2010],
and in Ikroavik Lake, where we did not measure bubble
methane
concentration,
we
applied
PALIMMN
mean
concentrations. This resulted in mass-based ﬂux estimates
of 81 ± 20 mg CH4 m2 d1and 1.5 ± 0.7 mg CH4 m2 d1,
respectively. It should be noted that these estimates represent
the magnitude of mean annual methane ebullition that
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1031
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 18 ---
reaches the lake surface or bottom of the lake-ice depending
on season; so they do not represent actual emissions to the
atmosphere. The latter requires consideration of seasonal
variability in seep ebullition, ice phenology, and other
physical and biological processes that inﬂuence methane
cycling in lakes.
[67] One contributor to the higher methane concentration
values in Goldstream Lake is the relatively shallow water
column at this site (2.9 m), which limits the magnitude of
degassing from ascending bubbles [McGinnis et al., 2006].
Differences in bubble size could also inﬂuence the magnitude
of methane dissolution during ascent [McGinnis et al., 2006];
however, bubble sizes did not appear to be different when we
made ﬂux measurements or when individual bubbles were
trapped by ice. We did not collect quantitative bubble size
data at all lakes. Since methane solubility is temperature
dependent, sediment and water column temperature will also
inﬂuence bubble methane content seasonally [Del Sontro
et al., 2010]. Similarly, bubbles resting under winter ice prior
to entrapment when ice thickens around them are subject to
methane dissolution followed by potential microbial oxida-
tion in the water column [Walter et al., 2008; Boereboom
et al., 2012]. This process effectively lowers methane ebulli-
tion emissions when ice melts in spring since not all methane
released from sediments as bubbles reaches the atmosphere.
Other biogeochemical cycling process must also inﬂuence
seep-bubble methane concentrations, because bubbles of
similar size released from lakes with similar depths and tem-
perature regimes had highly variable methane concentrations
in the PALIMMN data set [Walter Anthony et al., 2010]. To
improve accuracy in mass-based ebullition ﬂux estimates,
ebullition seep gases should be sampled seasonally, analyzed
for methane concentration, and accounted for in whole-lake
ﬂux calculations.
3.10.
Implications for Estimating Regional Lake
Methane Ebullition Using the Example of Lakes Near
Toolik Field Station, Alaska
[68] If lakes in a region have commonalities in the physical
and biogeochemical processes that control methane produc-
tion and ebullition emission such as substrate availability,
water depth, permafrost characteristics, periglacial history,
and lake developmental cycles, then data from multiple lakes
can be pooled as if they were multiple transects belonging to
a single lake to determine a regional lake ﬂux estimate. For
example, the mean and standard error of seep density on
lakes near Toolik Field Station was 0.66 ± 0.11 seeps m2.
The total ebullition ﬂux from all 64 transects was 96.8 L
d1 surveyed over 2,996 m2, yielding a mean ebullition ﬂux
value of 32 mL m2 d1. Pooling these transect data as if
they were collected from a single large regional lake, the
relationship between seep density and the number of distrib-
uted transects surveyed (equation (1)) suggests that the
spatial relative error associated with this ebullition estimate
is 8.5%. Propagating this error with the ﬂux-based errors of
each of the seep classes in Figure 1 (24%), we conclude that
at the regional scale, there is 95% conﬁdence that the mean
seep ebullition in lakes in the northern foothills of the
Brooks Range near Toolik Field Station is 32 ± 8.4 mL m2
d1. This exercise demonstrates that pooling results of a
low number of transects per lake on a larger number of
widely distributed lakes reduces the error associated with
the mean, allowing for better-constrained regional lake
ebullition ﬂux estimates.
4.
Conclusions
[69] Constraining the magnitude of methane ebullition in
thermokarst lakes is important for understanding the role of
permafrost thaw in the arctic carbon cycle and feedbacks to
global climate change. Until now, the uncertainty in
estimates of mean lake ebullition associated with spatial
heterogeneity of ebullition seeps was unknown and hindered
the upscaling of transect-level seep surveys to whole-lake
ebullition estimates. We used short- and long-term ebullition
ﬂux measurements and point process modeling to quantify
the errors associated with ﬁeld survey methods and to reveal
spatial patterns of ebullition in a variety of thermokarst lakes
in Alaska. Our ﬁndings indicate that the classic randomized
bubble-trap method for estimating mean lake ebullition is
highly median-biased toward underestimation of ﬂux. The
ice-bubble survey transect method yielded more accurate
estimates of methane ebullition in lakes than the classic
bubble-trap method and were only slightly median-biased
toward underestimation of ﬂux. These results explain previ-
ously observed discrepancies between ebullition estimates
based on randomly placed bubble traps and estimates based
on ice-bubble transects [Wik et al., 2011a, 2011b; Walter
et al., 2006]. Class-mean ebullition ﬂuxes measured in the
ﬁeld by the Pan-Arctic Lake Ice Methane Monitoring
Network (PALIMMN) can be assigned to seeps identiﬁed
as A, B, C and Hotspot types according to their ice-bubble
patterns in other northern lakes without compromising the
accuracy or precision of ice-bubble survey transect estimates
of whole-lake mean ebullition emissions.
[70] Field surveys revealed that ebullition was higher and
better constrained in the yedoma-type thermokarst lakes that
had high seep densities. In these lakes, methane seepage was
greatest adjacent to thermokarst margins. Seeps were
clustered around conical thermokarst mounds (baydjarakhs)
protruding from the lake bottom with a spacing of ≈10 m.
Soil mounds are the ice wedge polygon centers and thus part
of a regular micro scale landscape pattern across permafrost
regions, which makes their distribution (including on lake
bottoms) also predictable. Radiocarbon ages of ebullition
methane, 40.4 kyr and 41.3 kyr, collected as bubbles emitted
from the mounds, revealed that Pleistocene-aged organic
matter contained in baydjarakhs fuels methanogenesis in
the lake bottom.
[71] Statistical results showed that increasing the number
of survey transects distributed widely across lake surfaces
signiﬁcantly reduces the uncertainty in whole-lake ebullition
estimates. Neither the fraction of the lake surveyed nor seep
clustering had a signiﬁcant effect on the accuracy of whole-
lake mean ebullition estimates. However, precision in
estimating ebullition seep ﬂuxes at the whole-lake scale
using the ice-bubble transect method was inversely related
to the density of seeps on lakes. A calculated uncertainty
can be applied to estimating ebullition from lakes by propa-
gating errors associated with seep-class mean ﬂuxes, the
number of distributed transects surveyed and seep densities.
Regional-scale lake ebullition estimates can also be better
constrained by pooling results of a low number of transects
per lake on a larger number of widely distributed lakes, as
WALTER ANTHONY AND ANTHONY: THERMOKARST LAKE METHANE EBULLITION
1032
 21698961, 2013, 3, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20087, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
