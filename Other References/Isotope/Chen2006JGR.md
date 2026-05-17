
--- Page 1 ---
Estimation of atmospheric methane emissions between 1996
and 2001 using a three-dimensional global chemical
transport model
Yu-Han Chen1,2 and Ronald G. Prinn1
Received 8 April 2005; revised 20 September 2005; accepted 20 January 2006; published 31 May 2006.
[1]
Using an atmospheric inversion approach, we estimate methane surface emissions for
different methane regional sources between 1996 and 2001. Data from 13 high-frequency
and 79 low-frequency CH4 observing sites have been averaged into monthly mean
values with associated errors arising from instrumental precision, mismatch error, and
sampling frequency. Simulated methane mole fractions are generated using the 3-D global
chemical transport model (MATCH), driven by NCEP analyzed observed meteorology
(T62 resolution), which accounts for the impact of synoptic and interannually varying
transport on methane observations. We adapted the Kalman filter to optimally estimate
methane flux magnitudes and uncertainties from seven seasonally varying (monthly
varying flux) and two aseasonal sources (constant flux). We further tested the sensitivity of
the inversion to different observing sites, filtered versus unfiltered observations,
different model sampling strategies, and alternative emitting regions. Over the 1996–2001
period the inversion reduces energy emissions and increases rice and biomass burning
emissions relative to the a priori emissions. The global seasonal emission peak is
shifted from August to July because of increased rice and wetland emissions from
southeast Asia. The inversion also attributes the large 1998 increase in atmospheric CH4 to
global wetland emissions. The current CH4 observational network can significantly
constrain northern emitting regions but not tropical emitting regions. Better estimates of
global OH fluctuations are also necessary to fully describe interannual methane
observations. This is evident in the inability of the optimized emissions to fully reproduce
the observations at Samoa.
Citation:
Chen, Y.-H., and R. G. Prinn (2006), Estimation of atmospheric methane emissions between 1996 and 2001 using a three-
dimensional global chemical transport model, J. Geophys. Res., 111, D10307, doi:10.1029/2005JD006058.
1.
Introduction
[2] Methane is the second most radiatively important
greenhouse gas attributable to human activity. It contributes
about 15% of the total 2.5 W m2 increase in radiative
forcing caused by the anthropogenic release of greenhouse
gases in the industrial age [Hansen and Sato, 2001]. Direct,
worldwide observations over the past two decades show the
atmospheric burden of methane increasing at approximately
0.5% per year. Sources of methane arise from both natural
(e.g., wetlands, termites) and anthropogenic (e.g., rice,
domesticated animals, natural gas) processes, but large
uncertainties exist in their individual magnitudes and tem-
poral variabilities on the global and regional scale. A better
understanding of the current methane budget is necessary
to predict possible changes and feedbacks due to climate
change, and to sensibly formulate emission reduction
strategies.
[3] Estimation of surface fluxes of methane and many
other trace gases on the global scale have relied on three
broad approaches. The first approach is extrapolation of few
direct flux (or proxy flux) measurements to larger regions
(e.g., Matthews and Fung [1987] for wetlands and Lerner et
al. [1988] for animals). The second approach uses process
models that represent the actual physical and biological
processes of methane production (e.g., Walter et al. [2001a]
and Cao et al. [1996b] for wetlands). These first and second
methods are considered to be ‘‘bottom-up’’ approaches, and
suffer from extrapolation and process model errors, respec-
tively. This study uses a third method known as the ‘‘top-
down’’ or atmospheric inverse modeling approach. Here,
worldwide methane observations and simulations from
atmospheric chemical transport models are combined to
estimate unknown CH4 fluxes and their uncertainties. This
approach has errors arising from insufficient observations
and model errors.
JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 111, D10307, doi:10.1029/2005JD006058, 2006
1Center for Global Change Science, Department of Earth, Atmospheric,
and Planetary Science, Massachusetts Institute of Technology, Cambridge,
Massachusetts, USA.
2Now at AAAS Science and Technology Policy Fellowship Program,
Washington, DC, USA.
Copyright 2006 by the American Geophysical Union.
0148-0227/06/2005JD006058
D10307
1 of 25


--- Page 2 ---
[4] Most previous methane inversion studies have used
weekly flask measurements, a single year of meteorolog-
ical transport, and solved for multiyear or seasonally
averaged methane fluxes or flux trends [e.g., Fung et
al., 1991; Hein et al., 1997; Houweling et al., 1999;
Dentener et al., 2003; Cunnold et al., 2002]. Over the
past two decades, however, high-frequency in situ mon-
itoring stations that measure atmospheric methane have
become active (e.g., Advanced Global Atmospheric Gases
Experiment (AGAGE)). These high-frequency measure-
ments capture the full variability of CH4 observations,
including the occurrence of subweekly synoptic transport
events [Chen and Prinn, 2005]. They offer, in principle,
much greater information about methane sources, sinks,
and transport compared to weekly measurements. In
addition to high data frequency, realistic driving meteo-
rology (as opposed to climatological or a single year of
winds) is required to accurately determine of interannual
and intra-annual fluxes. This is because interannual trans-
port can strongly affect the tracer observations. Chen and
Prinn [2005] showed the strong year-to-year influence of
El Nin˜o and the North Atlantic Oscillation (NAO) on
methane mole fractions at Samoa and Mace Head, re-
spectively, from transport alone. Warwick et al. [2002]
also examined the global influence of transport IAV on
atmospheric CH4 growth rates.
[5] In this study, methane emissions from different
regional sources and their uncertainties are optimally
estimated between 1996 and 2001 at monthly time
resolution. The inversion incorporates high- and low-
frequency observations and realistic meteorology. We
use the Model for Atmospheric Transport and Chemistry
(MATCH) driven by National Center for Environmental
Prediction (NCEP) analyzed observed winds at a resolu-
tion (1.8  1.8) higher than previous global studies.
We also adapt the Kalman filter to optimize methane
emissions at monthly time resolution. In addition to
computational convenience, the Kalman filter quantifies
the usefulness of each additional measurement in an
observational time series. Optimized interannual and av-
eraged monthly emissions and uncertainties are compared
to previous bottom-up and top-down estimates.
[6] We also test the sensitivity of the inversion results to
different observational choices and modeling scenarios. In
addition to varying the observational networks, we compare
the use of unfiltered (i.e., pollution events retained) versus
filtered (i.e., pollution events removed) observations. We
also compare the use of modeled monthly means using all
model time steps within a month versus time steps
corresponding to exact observational times. Finally, the
potential impact of model errors is discussed. In general,
the timescale of the results (e.g., multiyear average vs.
monthly anomalies) will determine which model error
dominates.
[7] This paper is organized as follows: Methane obser-
vations, the MATCH model (including OH sink), and
methane sources are described in sections 2, 3, and 4,
respectively. The Kalman filter inverse methodology is
described in section 5. Section 6 focuses on the optimized
seasonal cycle and annually averaged results, and includes
observational and model sensitivity tests. Section 7
describes the interannual methane fluxes for seven seasonal
sources at monthly time resolution. The conclusions are
contained in section 8.
2.
Observations
[8] An initial step in this study was the accumulation,
intercalibration, and organization of available methane time
series. Table 1 lists the sites names, locations, and labora-
tory affiliations, as well as error information described later
in section 5.1. The site locations are also plotted in Figure 1.
Although methane observations exist at other sites, we
include only those 92 sites that were active during our
1996 to 2001 period of study, are of reliable quality, and can
be adequately represented using a global CTM. As much in
situ high-frequency data as possible was included. The
13 high-frequency stations (sites 1–13), listed in Table 1
and shown in Figure 1, measure methane mole fractions in
situ between 24 and 36 times per day [e.g., Prinn et al.,
2000]. The 74 low-frequency sites (see Table 1 and Figure 1)
are locations where flask samples of air are collected, with
measurement at a later date [e.g., Climate Monitoring and
Diagnostics Laboratory (CMDL), 2001b]. Most sites are
located in the Northern Hemisphere; tropical land regions
which have important methane sources in particular are
undersampled (Figure 1). We have further divided the flask
measurements into sites that have more (41 sites) and less
(38 sites) than 42 out of 60 (i.e., 70%) of the monthly
averaged observations between July 1996 and June 2001,
which is the 5-year time period of the inversion. For brevity,
Table 1 includes only those 41 flask sites that were 70%
active (sites 14–54); the other flask sites are included in the
auxiliary material.1
[9] With 5 active sites, the Advanced Global Atmospheric
Gases Experiment (AGAGE) operates the greatest number of
high-frequency stations over the longest time period. Chen
and Prinn [2005] conducted a modeling analysis of the these
high-frequency observations using the identical version of
MATCH as used here. Peak fluctuations, including pollution
events that typically last a few days or less, that are well
characterized by high-frequency sampling can be reproduced
by MATCH. The high frequency of in situ measurements is a
significant advantage over flask sampling, which typically
has a temporal resolution of one observation per week.
However, flask sampling allows greater spatial coverage
because samples need only be collected, rather than mea-
sured, at a particular site. The Climate Monitoring and
Diagnostics Laboratory (CMDL, NOAA) operates the great-
est number of flask sites, in conjunction with several other
laboratories. In this study, duplicate and triplicate flask
measurements have been averaged to produce single CH4
mole fractions at each time.
[10] For both in situ and flask observations, we have
discarded samples flagged as having obvious contamina-
tion, but retained those samples labeled as polluted or
nonbaseline (e.g., AGAGE [Cunnold et al., 2002] and
CMDL [Conway et al., 1994]). The ‘‘Filt’’ column in
Table 1 lists whether or not observations at a particular
site include identified pollution events over 1996–2002.
A data set with pollution events removed is considered to
1Auxiliary material is available at ftp://ftp.agu.org/apend/jd/
2005jd006058.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
2 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 3 ---
Table 1. Methane Measuring Site Information for High-Frequency and Flask Observationsa
Site
Name
Location
Latitude
Longitude
Altitude,
m
Network
Affiliationb Filt.
Flask
Error
Sample
Frequency
Error
Mismatch
Error
Total
Error
13 High-Frequency Sites
1
alt
Alert, Greenland
82
62
210
4
N
0.0
0.8
4.9
6.5
2
brw
Barrow, Alaska
71
156
11
2
Y
0.0
0.8
21.5
21.9
3
mhd
Mace Head, Ireland
53
9
25
1
Y
0.0
1.2
14.8
15.4
4
frd
Fraserdale, Canada
49
81
250
4
N
0.0
2.0
13.9
14.6
5
coi
Cape Ochi-Ishi, Japan
43
145
100
5
N
0.0
1.9
6.8
8.1
6
thd
Trinidad Head, CA
41
124
140
1
Y
0.0
1.4
8.2
9.3
7
iza
Tenerife, Canary Islands
28
16
2360
2
N
0.0
2.1
6.0
7.6
8
mnm
Minamitorishima, Japan
24
153
8
6
N
0.0
2.3
5.4
7.2
9
hat
Hateruma
24
123
47
5
N
0.0
3.5
14.5
15.5
10
mlo
Mauna Loa, Hawaii
19
155
3397
2
Y
0.0
0.6
5.0
6.5
11
rpb
Barbados
13
59
42
1
Y
0.0
0.4
3.6
5.5
12
smo
Samoa
14
170
42
1
N
0.0
0.3
2.5
4.9
13
cgo
Cape Grim, Australia
41
145
94
1
Y
0.0
0.4
3.1
5.2
54 Flask Sites Operational for at Least 70% of the Months Used in the Inversion
14
zep
Zeppelin St., Norway
78
11
474
2
Y
1.8
11.7
3.0
13.3
15
stm
Atlantic Ocean, Norway
66
2
7
2
Y
1.4
11.5
4.2
13.2
16
ice
Storhofdi, Iceland
63
20
100
2
Y
1.6
14.6
4.8
16.2
17
Sis
Shetland Is., Scotland
60
1
30
3
Y
5.0
26.7
5.9
29.4
18
bal
Baltic Sea, Poland
55
16
7
2
Y
1.5
22.3
19.7
30.1
19
cba
Cold Bay, Alaska
55
162
25
2
Y
1.5
6.3
3.4
8.7
20
shm
Shemya Island, Alaska
52
174
40
2
Y
1.6
8.8
2.7
10.6
21
epc
Estevan Pt, BC, Canada
49
126
39
3
Y
1.8
18.2
9.2
21.1
22
hun
Hegyhatsal, Hungary
46
16
344
2
Y
1.7
38.2
93.6
101.2
23
lef
Park Falls, Wisconsin
45
90
868
2
Y
1.8
16.5
13.3
21.9
24
uum
Ulaan Uul, Mongolia
44
111
914
2
Y
1.6
15.0
6.8
17.3
25
kzd
Sary Taukum, Kazahkstan
44
77
412
2
Y
1.1
17.6
18.8
26.1
26
bsc
Black Sea, Romania
44
28
3
2
Y
2.0
24.4
39.0
46.4
27
kzm
Plateau Assy, Kazahkstan
43
77
2519
2
Y
1.4
14.7
16.2
22.4
28
nwr
Niwot Ridge, Colorado
40
105
3475
2
Y
1.4
11.7
10.3
16.3
29
uta
Wendover, Utah
39
113
1320
2
Y
1.7
11.6
24.6
27.7
30
azr
Azores
38
27
40
2
Y
1.8
12.2
3.5
13.8
31
tap
Tae-ahn Pen., Korea
36
126
20
2
Y
2.0
32.3
44.0
54.9
32
wlg
Mt. Wanliguan, China
36
100
3810
2
Y
1.6
14.9
16.0
22.5
33
lmp
Lampedusa, Italy
35
12
85
2
N
5.0
12.8
10.7
19.8
34
bme
Bermuda East
32
64
30
2
Y
1.7
15.0
7.4
17.5
35
bmw
Bermuda West
32
64
30
2
Y
1.7
15.5
7.9
18.2
36
wis
Sede Boker, Israel
31
34
400
2
Y
1.3
13.6
10.4
17.8
37
mid
Midway
28
177
4
2
N
1.4
10.6
4.9
12.7
38
key
Key Biscayne, FL
25
80
3
2
Y
1.7
25.4
39.0
46.8
39
ask
Assekrem, Algeria
23
5
2728
2
N
1.8
8.5
5.0
11.2
40
kum
Kumukahi, Hawaii
19
154
3
2
Y
2.1
8.2
3.4
10.6
41
gmi
Guam
13
144
2
2
Y
1.5
6.5
6.0
10.2
42
sey
Mahe Island, Seychelles
4
55
3
2
N
1.5
8.7
9.9
14.1
43
asc
Ascension Island
7
14
54
2
Y
1.6
1.9
1.7
5.8
44
cfa
Cape Ferguson, Aust.
19
147
2
3
Y
2.3
8.6
3.6
11.2
45
eic
Easter Island
27
109
50
2
Y
1.7
1.9
0.5
5.6
46
crz
Crozet, Indian Ocean
46
51
120
2
Y
1.8
2.4
0.3
5.9
47
mqa
Macquarie Island
54
158
12
3
Y
2.2
3.8
0.5
7.1
48
tdf
Teirra Del Fuego
54
68
20
2
N
1.6
4.7
2.4
7.5
49
psa
Palmer Station, Antarct
64
64
10
2
N
1.5
2.1
0.5
5.5
50
maa
Mawson St., Antarctica
67
62
32
3
Y
2.5
1.8
0.2
6.7
51
syo
Syowa Station, Antarct.
69
39
11
2
Y
1.8
1.7
0.3
5.7
52
hba
Halley Bay, Antarctica
75
26
10
2
N
1.3
1.4
0.3
5.1
53
spoc
South Pole
89
24
2810
3
Y
1.8
2.4
0.9
6.0
54
spo
South Pole
89
24
2810
2
Y
2.7
3.5
0.9
7.6
Laboratory Network
Definition
Intercalibration Factor
1
AGAGE Advanced Global
Atmospheric Gases
Experiment
1
2
CMDL
Climate Monitoring and
Diagnostics Laboratory
1.0119
3
CSIRO
Commonwealth Scientific
and Industrial Research
Organization (Australia)
1.0119
4
MSC
Meteorological Service of
Canada
1b
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
3 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 4 ---
be more easily modeled by global models, at the expense
of loss of near-station flux information. We test the use of
both unfiltered (polluted) and filtered data sets in the
inversion. Our study uses only actual observations, and
not smoothed/interpolated flask data such as from CMDL
[2001b], since we use MATCH’s capability to simulate
observed values at specific times [Chen and Prinn, 2005].
[11] The standards used for the absolute calibration of
methane mole fractions differ among most laboratories,
although most have been calibrated to either the AGAGE
or CMDL scale. We use the AGAGE standard, which is
based on the Tohoku gravimetric technique standard as
described by Cunnold et al. [2002]. The most recent
intercalibration factor of 1.0119 [Cunnold et al., 2002] is
used to convert other CH4 measurements based on the
CMDL scale, as done by Chen and Prinn [2005]. The
reported precision due to random instrumental error of most
methane measurements is between 0.07 and 0.2% [Cunnold
et al., 2002; CMDL, 2001a], equivalent to a 1–3 ppb error.
These absolute calibration and instrumental precision errors
are usually small compared to other types of observational
error, as described and quantified in section 5.1.
3.
Methane Sources
[12] Along with the compilation of the global methane
observations, a major preparatory effort in this study was
the accumulation of the most up-to-date, readily available
surface emission data sets for use in global atmospheric
models. A first step was the choice of emission source, since
often more than one emission distribution and magnitude is
available for a single process. Once chosen, individual
emission patterns were then (1) combined to create a
‘‘reference’’ CH4 flux, representing an initial guess of all
methane emissions, and (2) used to define the emission
patterns for individual sources whose magnitudes are opti-
mally estimated in the inversion (section 5). For certain
processes, this included further spatial aggregation or dis-
aggregation of the original data set.
[13] Surface flux fields of methane were obtained or
adapted from data sets described by Fung et al. [1991] for
wetland and rice, EDGAR 3.0 (Emission Database for
Global Atmospheric Research, available at http://
www.mnp.nl/edgar/) for aseasonal anthropogenic emis-
sions, and Hao and Liu [1994] for biomass burning.
Figure 2 shows the individual emission patterns of these
data sets at MATCH T62 resolution. The reference
emission magnitudes are listed in the first row of Table 2.
The original emission magnitudes for the different processes
have been scaled to more realistically reflect current bottom-
up estimates, with additional uniform scaling of approx-
imately 10% for all emissions to match the CH4 growth
rate, largely determined by the global OH concentration
(section 4). Note that these reference emission magnitudes
fall within the literature ranges (rows 1 and 2 in Table 2).
The relatively high total emission of 589 Tg yr1 is a
result of the prescribed global OH level. A plot repre-
sentative of the annual average pattern of the total
reference emission is given by Chen and Prinn [2005].
[14] The individual sources can also be divided into
seasonal (Wetlands, Rice, and Biomass Burning) and
aseasonal emissions (AnMSW, and Energy), as shown
in Table 2. Seasonal processes, such as wetlands, contain
significant variations in magnitude and location in suc-
cessive months. Aseasonal components, such as emissions
from fossil fuel methane leakage, are considered to vary
much more slowly, and are assigned a constant spatial
distribution over time. We have further divided the total
wetland distribution into northeastern (WetNE), north-
western (WetNW), and southern and tropical (WetS)
regions. This division allows each of these regions to
be solved separately in the inversion, rather than as a
whole. Several observational stations surround WetNE
and WetNW, allowing the measurement of distinct emission
signatures from these two regions. Further disaggregation
would be difficult due to the relative sparseness of methane
observations within these regions. Note that we have com-
bined wetlands in southeast Asia with rice emission to create a
unified emission pattern (RiceWet). These two processes
share substantially overlapping spatial and temporal charac-
teristics, especially in southeast Asia, which makes their
individual estimation difficult. As shown in Figure 2, greater
than 80% of the rice flux originates from China, India, and
southeast Asia. In the inversion, the rice region thus repre-
sents a large emitting region localized in southeast Asia, with
more diffuse emissions spread across other tropical and
temperate areas (Figure 2). The partitioning between rice
and wetland emissions in RiceWet is 21% wetland (located in
Southeast Asia) and 79% rice, based on reference values.
[15] Seasonal ITCZ shifts, which influence tropical pre-
cipitation, result in a strong seasonal variation in the
distribution of biomass burning in Asia, South America,
and Africa [see Hao and Ward., 1993; Hao and Liu, 1994;
Duncan et al., 2003b]. This results in a bimodal seasonal
behavior of tropical biomass burning emissions. In the
inversion, animal and waste (AnMSW) emissions are solved
as a combined emission since they have similar spatial
patterns. Figure 2 shows their combined distribution, with
the greatest emissions near highly populated centers. The
Energy emission distribution is a combination of gas and
Table 1. (continued)
Laboratory
Network
Definition
Intercalibration Factor
5
NIES
National Institute for
Environmental Studies
(Japan)
1
6
JMA
Japan Meteorological
Agency
1
aActive for at least 70% of the months between 1996 and 2001. See Figure 1 for locations and section 5.1 for description of error terms. All errors in ppb
CH4.
bMSC CH4 standard is within 0.3% of AGAGE.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
4 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 5 ---
coal emissions, which are dominated by Northern Hemi-
spheric sources.
[16] An additional 60 Tg yr1 source of methane included
in the reference run, but not solved for in the inversion, is
denoted as ‘‘Other’’ emissions in Table 2. This source
includes 36 Tg yr1 from industrial and biofuel combustion
sectors (EDGAR 3.0) and 23 Tg yr1 from termites [Fung
et al., 1991]. We do not solve for these ‘‘other’’ processes
because they relatively diffuse, i.e., distributed over large
space scales and therefore measurements are not therefore
very sensitive to them. We do not include geologic sources
of methane [Etiope and Klusman, 2002; Milkov et al., 2003]
and methane hydrates in this study, since their timing,
spatial distribution, and emission magnitude are highly
uncertain. For the similar reasons, we also ignore methane
surface sinks due to bacterial consumption [Ridgwell et al.,
1999].
4.
Chemical Transport Model and Atmospheric
Sink
[17] The MATCH model was developed to realistically
simulate atmospheric constituents using analyzed observed
meteorology [Rasch et al., 1997; Mahowald, 1996].
Throughout this work, MATCH is driven by NCEP reanal-
ysis meteorology at T62 spectral resolution, which corre-
sponds to approximately 1.8  1.8 horizontal resolution.
In the vertical, NCEP-driven MATCH has 28 sigma levels
between 1000 and 2.9 mb. The surface (bottom) layer
varies between 50 and 100 m in height. An identical version
of MATCH is described by Chen and Prinn [2005],
including the atmospheric OH sink, for a forward modeling
study of atmospheric methane. The spatial and temporal
pattern of the OH field corresponds to the output of a
T62 run of MATCH for 1997 using comprehensive
atmospheric chemistry [Lawrence et al., 1999; Jockel,
2000; von Kuhlmann et al., 2003]. As described by Chen
and Prinn [2005], we further adjusted the total magnitude
of this OH field to best fit high-frequency methlychloro-
form (MCF) observations between 1978 and 2001 at the
5 ALE/GAGE/AGAGE sites [Prinn et al., 2000]. The
simulation of MCF observations is a good diagnostic of
an OH field because the major sources of MCF are
relatively well known and its sink is dominated by
reaction with OH. The MCF mole fractions at all
AGAGE stations are fairly well reproduced, suggesting that
the spatiotemporal characteristics of the MATCH OH field are
broadly correct. The OH field has a tropospheric annual
average concentration of 1.1  106 molecules cm3
weighted by mass, which is within 5% of recent opti-
mized OH values between 1996 and 2001 using AGAGE
MCF data (J. Huang, personal communication, 2004).
5.
Inverse Methodology
[18] This section describes the inversion methodology
used to estimate methane emissions and their uncertainties.
The methodology is based on the Kalman filter, which has
been used in the study of different atmospheric trace gases
on the global scale [e.g., Prinn et al., 2001; Huang, 2000;
Figure 1.
Location of methane measuring sites. The large red letters denote high-frequency in situ
stations. Blue and green letters represent flask sites, with available data greater and less than 70% of the
60 months between July 1996 amd June 2001, respectively.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
5 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 6 ---
Figure 2.
Methane reference source spatial distributions. Annually averaged patterns are shown; note
that the seasonal sources (WetNE, WetNW, WetS, RiceWet, BBAF, BBAS, and BBAM) have significant
pattern changes from month to month. The aseasonal sources (AnMSW and Energy) have constant
emission patterns.
Table 2. Aggregated Methane Sourcesa
General Source
Number
of Sites
MATCH
Interpol.
Filt.
Obs.
Wetlandsb
Ricec
Biomass
Burningd
AnMSWe
Energyf
Otherg
Total
Reference total
151 ± 54
92 ± 38
29 ± 19
169 ± 169
89 ± 89
60
589
Literature Rangeh
92–237
25–126
23–55
120–180
53–159
500–600
Inversion cases
Control
54
Yes
No
145 ± 25
112 ± 29
43 ± 18
189 ± 3
48 ± 3
60
596
HF Only
13
Yes
No
148 ± 33
96 ± 29
37 ± 19
211 ± 4
46 ± 4
60
597
HF MBLi
14
Yes
No
146 ± 30
110 ± 31
43 ± 18
184 ± 3
54 ± 3
60
596
Control FullAv
92
Yes
No
146 ± 26
111 ± 30
46 ± 18
180 ± 3
54 ± 3
60
596
Control Filt
54
No
No
145 ± 25
114 ± 28
44 ± 18
187 ± 3
48 ± 3
60
597
Control FullAv Filt
92
No
No
145 ± 26
115 ± 29
48 ± 18
175 ± 3
55 ± 3
60
597
All Obs
92
Yes
Yes
143 ± 25
113 ± 29
48 ± 17
185 ± 3
48 ± 3
60
596
All Obs FullAv
54
Yes
Yes
144 ± 25
115 ± 28
47 ± 18
182 ± 3
49 ± 3
60
596
All Obs Filt
54
No
Yes
145 ± 26
115 ± 28
48 ± 18
178 ± 3
52 ± 3
60
597
aValues are in Tg CH4 yr1. The inversion results listed here are a summation of individual regional sources listed in Table 3. Section 5.5 describes the
nine inversion cases.
bWetlands = WetNE + WetNW + WetS + 0.21  RiceWet (see Table 3).
cRice = RiceWet  0.79 (see Table 3).
dBiomass Burning = BBAF + BBAM + BBAS.
eAnimal and waste emissions spatial distribution are adapted from EDGAR 3.0.
fEnergy = Natural Gas + Coal.
gIncludes 36 Tg yr1 industrial/biofuel combustion (EDGAR 3.0) and 23 Tg yr1 termites [Fung et al., 1991].
hRange from IPCC [2001] except for Wetlands [Walter et al., 2001a] and Energy [Walter et al., 2001a].
iIncludes 21 flask sites that sample marine boundary layer (MBL) air: asc, azr, bme, bmw, cba, chr, crz, eic, hba, ice, kum, mid, psa, shm, spo, syo, tdf,
cfa, cgo, maa, mqa.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
6 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 7 ---
Mahowald et al., 1997], but not for methane. On the
regional scale, the KF has been used to constrain emissions
of chlorine compounds [Kleiman and Prinn, 2000] and
methane [Janssen et al., 1999]. Here we adapt the Kalman
filter to estimate methane emissions from the regional
sources described in the previous section at monthly time
resolution for seasonally varying emissions and as time-
invariant fluxes for more steady emissions. Much of the KF
methodology is described by Prinn [2000] and Enting
[2002].
[19] Section 5.1 describes the formulation of the obser-
vational errors, which determine the relative importance
(i.e., weighting) of each observation in the Kalman filter,
and consequently, its influence on the final optimized
emissions. Sections 5.2–5.4 describe the several equations
which make up the full Kalman filter as used in this study.
For brevity, we have moved some of the methodological
description, including matrix information, to the auxiliary
material.
5.1.
Observational Errors
[20] The observational error (ek) associated with the
observed (yk
o) monthly mean at time step k is not known;
however, its probability distribution function can be
expressed by its standard deviation, sk. The overall error
in the observed mole fraction, in units of ppb, can be
estimated as the net effect of the errors associated with:
(1) the instrumental and sampling errors, (2) the sampling
frequency used to define the monthly mean, and (3) the
mismatch error between observations and model as shown
in equation (1). This assumes (reasonably) that these errors
are uncorrelated.
sk ¼
ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
s2
measurement þ s2
mismatch þ s2
sampling frequency
q
ð1Þ
The measurement error, smeasurement, arises from imperfec-
tions in instrumentation, sampling, and intercalibration.
smeasurement ¼
ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
s2
instrument precision þ s2
flask errors þ s2
intercalibration
q
ð2Þ
[21] The instrumental precision of CH4 measurements is
approximately 1–4 ppb CH4 at nearly all sites [Cunnold et
al., 2002; CMDL, 2001a]. A uniform and deliberately
generous instrumental precision error of 4 ppb is chosen
for all sites. After adjusting all measurements to the
AGAGE standard (section 2) and assuming perfect accu-
racy, the calibration error, sintercalibration, should be zero.
However, Cunnold et al. [2002] report a 1 ppb difference
at a number of colocated AGAGE and CMDL sites even
after intercalibration, which we take as sintercalibration. For
both in situ and flask measurements, we discard obser-
vations that are flagged for obvious contamination or
instrumental difficulties. In addition, we retain measure-
ments that are flagged as representing nonbackground (or
‘‘polluted’’) air, since we generally expect that MATCH
can simulate these measurements [see Chen and Prinn,
2005]. However, the inversion sensitivity to filtered (non-
polluted) observations is also examined.
[22] The sampling frequency error, ssamplingfrequency, quan-
tifies how well a monthly mean quantity is defined given a
finite number of (m) measurements. Assuming temporally
uncorrelated data, the error on a mean quantity due to
limited sampling frequency is best represented as the
standard error, following Wunsch [1996]:
ssampling frequency ¼
ﬃﬃﬃﬃﬃﬃﬃﬃﬃ
s2
mon
m
r
ð3Þ
where smon
2
is the monthly mean variance, and m is the
number of observations taken in that month. Since smon
2
is
not known at locations other the high-frequency stations, it
is approximated at most sites using MATCH high-frequency
CH4 output from the reference run. This validity of this
assumption is supported by the forward comparisons of
Chen and Prinn [2005], which shows similar variability
between MATCH modeled and observed high-frequency
values. As listed in Table 1, high-frequency observations
have a very low monthly standard error given the high
number of measurements each month (m  1000). The
standard error is much higher for the weekly flask
measurements (m  4). At the same site, the high-frequency
standard error would be about 15 to 16 times smaller than
the flask error. This quantifies the greater weight of high-
frequency measurements in the inversion, especially at
locations with a high monthly variance.
[23] The mismatch (or ‘‘representation’’) error, smismatch,
describes the difference between an observation made at a
single point in space and a model-simulated observation
representative of a large volume of air [e.g., Prinn, 2000].
The degree to which a point measurement fails to
represent this volume dictates the size of the mismatch
error. Among other factors, this difference depends on the
resolution of the model, and the observational method
and location. Most methane observing sites are situated to
sample large, well-mixed volumes of air, which can be
more accurately modeled. The mismatch error increases
significantly over continental sites near emitting regions,
as MATCH does not have the resolution to cope with
local influences. The mismatch error is difficult to quan-
tify and may include bias error, where the model sys-
tematically overestimates or underestimates the observed
mole fractions. Previous studies have used the temporal
variability at a particular site as a proxy for this error
(e.g., Prinn et al. [2000] for halocarbons at high-
frequency stations and Gurney et al. [2002] for CO2).
We have chosen to estimate the mismatch error at each
site using the standard deviation of the CH4 mole
fraction yik (mean yk) at the nine model grid cells i
which contain and surround each observing site,
smismatch =
ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
1
9
P
9
i¼1
yik  yk
ð
Þ2
s
. This assumes that the spatial
variability within a single grid cell is related to the
variability among the neighboring grid cells. The mis-
match error so defined is much larger over strongly
emitting continental sites compared to remote ocean
locations. The mismatch error at each site also varies
by month, consistent with seasonal changes in emissions
and transport.
[24] Table 1 lists individual and total monthly mean errors
calculated at each site, averaged between July 1996 and
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
7 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 8 ---
June 2001. The total error is the aggregate sum of the
individual errors, following equation (1). For high-frequency
observations, the measurement and mismatch errors have the
largest contributions, while the sampling frequency errors are
small. In addition to these errors, the flask observations
nearby to strongly emitting regions have very large sampling
frequency errors. At a given site, the errors on monthly mean
observations can vary due to monthly variations in the
mismatch and sampling frequency errors. Very large error
bars usually correspond to a monthly mean observation
defined by only 1 or 2 flask measurements.
5.2.
Measurement Equation
[25] The measurement equation relates modeled and
observed mole fractions through the following relationship:
yo
k ¼ ~yk þ Hkxk þ yadj
kT1 þ ek
ð4Þ
The time index k refers to a specific month of observation
and flux. The observation vector, yk
o, contains the monthly
mean observations at time k for all sites. ~yk contains the
corresponding modeled monthly means from the MATCH
reference run. The primary objective of the inversion is to
obtain estimates of the state vector, xk, whose individual
elements represent adjustments to the reference emissions
(section 5.4). These adjustments improve the modeled fit to
observed CH4 monthly means, weighted by observational
errors. The elements of xk are related to the simulated CH4
mole fractions through the time-dependent sensitivity
matrix Hk, which is generated through multiple MATCH
runs (described below). Although xk is synchronized to the
monthly observational vector, yk
o, it also includes elements
corresponding to emissions from T months previous to time
k, i.e., k  1, k  2, . . . k  T. This is because an
observation depends not only on the emission at a particular
month, but also on all previous emissions. The contribution
of optimized emissions prior to month k  T are contained
in ykT1
adj
.
[26] The sensitivity matrix, H, relates the CH4 mole
fraction at each site to the CH4 emissions in the state vector.
Its elements are generated from multiple MATCH runs of
the individual methane source of interest (e.g., WetNE,
WetNW, etc.) described in section 3. Denoting the reference
values with a tilde, the elements of H can be written as:
hijkk0 ¼ yik  ~yik
xjk0  ~xjk0 ﬃ@yik
@xjk0
ð5Þ
In the above, i represents site number and j methane
regional source. The sensitivity elements for the aseasonal
sources contain the cumulative emission influence of all
previous months (month 1 to month k). The seasonal
sensitivities only contain emission information for a single
specific month (k0) prior to the month of interest (k), i.e.,
k0 
 k. The seasonal and aseasonal sensitivities are thus
generated by two different types of MATCH simulations.
For
the
aseasonal
fluxes
(e.g.,
animal
emissions),
emissions from a single source are perturbed above their
reference levels by 20% in all months and the model then
run over the entire data period. The sensitivities are then
expressed in units of ppb (Tg yr)1. For the seasonal
sensitivities (e.g., wetland flux from a specific month), a
1-month pulse of methane is emitted and its subsequent
dispersion is calculated for all succeeding months.
[27] The Kalman filter, as applied here, requires a linear
relationship between sources and mixing ratios. Chen
[2003] confirmed the linearity of simulated methane pulses
using the same version of MATCH. The use of a prescribed
(but optimized) OH field also ignores possible CH4-OH
feedbacks. We expect this effect to be small, however, as the
optimized methane adjustments are relatively small and
should not impact global OH concentrations substantially.
This effect is discussed in the auxiliary material.
5.3.
State-Space Equation
[28] This section describes how the state vector (xk)
evolves from one time step to another, before the use of
the next monthly observational vector. The state-space
equation, equation (6), uses the state transition matrix, Mk
to propagate the state vector from one time step to another.
In equation (7), this matrix also propagates the state error
covariance matrix (Pk) associated with the state vector from
one time step to another.
xf
k ¼ Mk1  xa
k1 þ hk1 ﬃMk1  xa
k1
ð6Þ
Pf
k ¼ Mk1  Pa
k1  MT
k1 þ Qk1
ð7Þ
These equations are also known as the ‘‘forecast’’ and
‘‘forecast error’’ equations in other applications. The
superscripts f (forecast) and a (analysis) denote a quantity
before and after, respectively, the use of one month’s
observations to update the state. These equations allow the
emissions contained in the state vector to be solved
recursively as new observations are added (see section 5.4).
A full matrix description of Mk is included in the
auxiliary material, but basically shifts the elements of the
state vector down. The elements which exit the state
vector are considered to be the optimized (adjustment)
emissions.
[29] If the state vector included the total seasonal emis-
sions, rather than adjustments to the reference emissions,
then Mk would describe the periodic evolution of seasonal
emissions from one month to another. An approximation to
this seasonal variation is implicitly contained in the refer-
ence run, ~yk. The random forcing term, hk1, may include
multiple sources of error in Mk1. For example, the annually
repeating seasonality, implicit in the reference run, cannot
describe year-to-year changes in seasonal fluxes at monthly
time resolution. This error is not known exactly for each
month, but its statistics are contained in Qk1, the random
forcing covariance matrix. This error corresponds to an
initial error of the adjustment flux before the use of
observations.
5.4.
Kalman Filter for Time-Dependent Inversions
[30] The Kalman filter produces an optimal estimate of the
state vector with each new monthly observational data set.
The previously described state-space equation and measure-
ment equations are combined with the standard Kalman filter
gain matrix, state vector update, and error covariance update
expressions to yield the full KF equations:
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
8 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 9 ---
[31] State space extrapolation equations
xf
k ¼ Mk1  xa
k1
ð8Þ
Pf
k ¼ Mk1  Pa
k1MT
k1 þ Qk
ð9Þ
Model measurement equation
yk ¼ ~yk þ Hkxf
k þ yadj
kT1
ð10Þ
Kalman gain matrix
Kk ¼
Pf
kHT
k
HkPf
kHT
k þ Rk
ð11Þ
State update
xa
k ¼ xf
k þ Kk  yo
k  yk


ð12Þ
State error covariance update
Pa
k ¼ I  KkHk
ð
ÞPf
k
ð13Þ
These equations are repeated over all monthly observations
from k = 1 to k = N. The measurement equation has been
modified to include only modeled CH4 mole fractions,
shown as yk. The Kalman gain matrix combines the prior state
error Pk
f, the observational error Rk (derived in section 5.1),
and the sensitivity matrix Hk to compute a weighting matrix.
This matrix determines the degree to which observations will
modify the prior state through the state update equation,
equation (12). A large gain matrix K corresponds to a strong
sensitivity of the state vector (emissions) to observations,
while a small K corresponds to a very low sensitivity to
observations. The optimized state, xk
a, will almost always
undergo some change from the prior state, xk
f, based on the
size of the gain matrix, Kk, and the modeled-observational
difference, yk
o  yk. Equation (13) updates the previous
estimated flux uncertainties contained in Pk
f, to the new
estimated uncertainty in Pk
a. Note that the diagonal elements
of Pk
a are always less than or equal to those contained in Pk
f.
This error, or uncertainty, reduction is an indicator of how
effectively the observations (given particular observational
errors) constrain the fluxes.
[32] Figure 3 shows how the estimates of the WetNE
(northwest wetlands) subvector of x corresponding to May
1999 emissions change as observations from subsequent
months are used. The starting value for WetNE is zero since
it represents an initial adjustments to the reference run
emissions. The large initial errors correspond to the diago-
nal elements of the a priori error covariance for the new
emission elements. For most seasonal fluxes the a priori
uncertainties (contained in Q in Equation 7) are taken as
±100% of their reference monthly emission. A smaller
initial uncertainty of ±30% is used for the much larger
emissions from WetS and RiceWet. With each new monthly
observation the initial value changes and the error
decreases, following equations (8)–(13). Note that the
adjustment and error reduction for this source (and all the
other seasonal sources) is greatest in the first 3–4 months,
by which time most flux values have stabilized. This
stabilization arises from the decreased sensitivity of obser-
vations to a given monthly emission due to atmospheric
mixing after the first few months. The optimized emission
values for May 1999 are obtained after using 11 months of
data (i.e., after the March 2000 observation). At each month
k, a new subvector of optimized emissions is thus produced.
The optimized emissions that exit the state vector are used
to update ykT1
adj
, which contains the continued influence of
these emissions on the global background methane mole
fraction.
[33] The optimization of the aseasonal components are
more straightforward. Since constant emissions are as-
sumed, solution toward a single value occurs over all time
steps, as shown in Figure 4 for the animal and waste
(AnMSW) case. For the aseasonal errors, a large a priori
uncertainty equal to ±100% the reference magnitudes was
chosen. The rapid decrease in the uncertainty indicates that
the inversion is relatively insensitive to this starting value.
The optimized flux adjustment is reached with the final
observation in the inversion. The optimized total aseasonal
emission is the sum of the flux adjustment and the reference
value. Because the aseasonal emissions are only fully
optimized at the final step of the time series, earlier seasonal
emissions (whose estimates also depend on the aseasonal
values) at the beginning of the filter will not be fully
optimized. In order to fully optimize all seasonal monthly
values, the entire Kalman filter is repeated with aseasonal
values fixed to their optimized values from the first KF run.
The optimized seasonal fluxes show only small changes in
this second run.
[34] At each time step, an optimally estimated set of
seasonal emissions is produced and collected to create an
optimally estimated emission time series over the entire
inversion period. Figure 5 (left) shows the monthly emis-
sion time series (red) using methane observations between
July 1996 and June 2001 for a representative inversion case
described in section 5.5. Figure 5 also contains the annually
repeating reference emissions (blue). The inversion results
contain the average, seasonal, and interannual behavior of
methane emissions, which are discussed in sections 6 and 7.
The right side of Figure 5 superimposes the optimized (red)
uncertainties on top of the reference (blue) uncertainties.
Note that the inversion always acts to reduce the initial
uncertainty by amounts depending on the sensitivity of the
observations to particular source regions, as well as the
observational error. The error reduction can be large (e.g.,
WetNW) or small (e.g., biomass burning regions). For most
sources, the final optimized value lies inside the range of the
initial uncertainty. Poorly constrained regional fluxes, such
as BBAM, may for some months lie slightly outside the
initial error bar. The initial and final uncertainties overlap in
nearly all cases, however.
[35] A test of the inversion is that the model-observational
comparison improves when the using the optimized emis-
sions, compared to the reference. The optimized monthly
mean CH4 mole fractions (red) using the Control case
inversion (described in section 5.5) are shown in Figure 6,
which also includes the observations (black) and reference
run (blue). The errors on the monthly mean observations are
shown as vertical bars, as derived in this section 5.1. The
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
9 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 10 ---
Figure 3.
Evolution of estimated WetNE source for a single monthly emission (May 1997). The vertical
axis corresponds to emission adjustments from the reference value (Tg CH4/month), and the horizontal
axis corresponds to the addition of new observational data. The initial vertical line shows the assumed a
priori error for the May, 1997 emission. The emission adjustment and error reduction occurs with the use
of the first few months of data. The final optimized emission estimate is taken after 11 months of
observations have been used, by which time the emission estimate has stabilized.
Figure 4.
Kalman state vector evolution for the AnMSW methane source. The vertical axis corresponds
to emission adjustments from the reference value (Tg CH4/month), and the horizontal axis corresponds to
the addition of new observational data. The initial vertical error bar corresponds to the a priori error
associated with AnMSW. Unlike the seasonal sources the aseasonal sources are solved as constant fluxes
over the entire 5-year time period. The optimized flux adjustment is obtained after all observations have
been used.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
10 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 11 ---
improvements in model-observational comparison are appar-
ent at nearly all sites, and used for interpretation of the
inversion results in the following sections. As a further check
of the inversion, we also verified that the monthly means in
the optimized forward run and that produced by Equation 10
are in agreement.
5.5.
Inversion Cases
[36] We investigated several different inversion cases that
are based on different sets of actual and simulated obser-
vations as described in Table 2. The Control inversion uses
13 high-frequency (HF) and 41 flask sites as listed in
Table 1, using monthly mean observations with MATCH
output interpolated to exact observation times before
averaging. The control case incorporates those 41 flasks
sites active for at least 70% of the 60 months between
July 1996 and June 2001 (Table 1). The next case, HF
only, uses only observations from the 13 high-frequency
sites. This case is followed by HF MBL, which adds 21 flask
sites located in the clean marine boundary layer (see Table 2).
The next inversion, Control FullAv, uses the same observa-
tional sites as Control, but uses the full MATCH output to
determine the modeled monthly means. As mentioned earlier,
the model monthly means can be computed using modeled
mole fractions interpolated to the exact 4–5 flask sampling
times each month (as in Control), or using all model time steps
(FullAv cases) for that month.
[37] The Control Filt case also uses the same observa-
tional sites as Control, but uses filtered observations in
which the pollution events have been removed for applica-
ble sites (see Table 1). The Control FullAv Filt case
combines the Control FullAv and Control Filt cases by
using full MATCH monthly means and filtered observa-
tions. The next three cases, All Obs, All Obs Full, and All
Obs Filt are analogous to Control, Control FullAv, and
Control Filt, except that all 92 (rather than 54) available
observation stations within the time period are used. The
additional 38 observation sites have monthly data for less
than 70% of the inversion months (Table 1). The distinction
between sites with greater and less than 70% of the data is
made because sites with significant amounts of missing data
may affect the inversion results independently of their
observational values. Most previous atmospheric inversions
(e.g., Gurney et al. [2002] for CO2) use a constant obser-
vational network, avoiding possible effects from the abrupt
activity/inactivity of intermittent sites. These impacts are
sometimes difficult to separate from emission changes due
to the interannual variability of actual observations. The
Figure 5.
(left) Monthly optimized emissions in the control case and (right) corresponding uncertainties
for the seven seasonally optimized sources. Note the different scales for the different sources. The
reference emissions and uncertainties are shown in blue, and Control inversion case results are in red.
Figure 5 (right) shows the error reduction from the reference (a priori) uncertainties.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
11 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 12 ---
Figure 6.
Monthly mean observations for all sites following the listing in Table 1. The observations
(black) are shown with corresponding error bars described in section 5.1. The reference (blue) and
optimized (red) curves represent MATCH monthly means (determined from model output interpolated to
exact measurement time). The optimized values correspond to MATCH using the results of the All Obs
inversion (see section 5.4).
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
12 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 13 ---
Figure 6.
(continued)
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
13 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 14 ---
70% cutoff is taken as an approximate value above which
these impacts are minimized, but we also perform the
inversion with all observations. Many of the inversion
results for the 9 different cases are similar. We restrict the
comparison across all cases to the annual average results
(see section 6.2).
6.
Average Inversion Results
[38] Although atmospheric CH4 observations have been
available since the early 1980s, only by the mid-1990s were
sufficient stations operational to conduct an inversion using
only high-frequency data at the regional scale. The inver-
sion starts in July, 1996, one month before the 5th AGAGE
station at Cape Matatula, Samoa began its high-frequency
CH4 measurements. The Kalman filter is formally run using
nearly 6 years of observational data, until May, 2002. For
the averaged inversion results described below, we consider
the 5 year span between July 1996 and June 2001. This is
because the last 11 months of the inversion includes
monthly emissions which are suboptimized, as they do
not incorporate a full 11 months of observations.
[39] Because of the large amount of methane emission
information produced by the inversion, we divide the
discussion of the methane inversion results into three
temporal sets: average seasonal (i.e., annually repeating
monthly fluxes), average annual (i.e., average fluxes over
the entire inversion period), and interannually varying (i.e.,
monthly fluxes). This division also facilitates comparison to
bottom-up estimates, which are often made at different
temporal scales. The average seasonal cycles are derived
by averaging the interannual monthly fluxes into a single
annually repeating cycle for each source shown in Figure 5.
To obtain the average annual results for the seven seasonal
sources, we further average these monthly fluxes into
annual mean fluxes applicable over the 5-year inversion
period. The aseasonal sources are solved for as constant
fluxes over the entire inversion period, and hence do not
require averaging.
6.1.
Average Seasonal Results
[40] The multiyear monthly averaged optimized emis-
sions and emission uncertainties are shown in Figure 7.
The flux values (Figure 7, left) are generated by averaging
the monthly values between July 1996 and June 2001 of the
interannually varying fluxes shown in Figure 5. Three
inversion cases are shown for each source: Control (red),
All Obs (orange), and HF (green). The original reference
case is shown in blue. The relationship s ¼
ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
P
N
t
s2
t =N
s
,
where N = 5 years and st is uncertainty for a particular
month (e.g., January), is used to compute the average
uncertainties for each month (Figure 7, right).
[41] For the Control inversion, both WetNE and WetNW
show decreases compared to the reference, particularly
during the fall and spring. The overall decrease in the
northern wetland regions is consistent with the reference
run overestimate of the observed northern hemispheric mole
fractions (Figure 6). These effects are particularly evident at
the nearby stations of alt[1], brw[2], and frd[4] (Figure 6)
where the reference overestimate has been corrected by the
inversion. WetS emissions show a dip during the spring
months centered around May compared to the relatively
smooth reference emissions, although emissions are still
globally active throughout the year. The inversion strongly
increases RiceWet emissions and shifts the maximum from
August to July. The optimized emissions are decreased
below the reference after September, partially offsetting
the large increase in earlier months. At least part of the
large increase can be attributed to the effects of measure-
ment sites downwind of Asia, where the reference run
exaggerates the seasonal mole fraction trough during June
and July. This is most clearly seen in the reference mole
fraction underestimates during the summer at coi[5] and
hat[9] (Figure 6), which are two high-frequency stations
most sensitive to RiceWet. The inversion compensates by
shifting the RiceWet emissions to an earlier and more
intense maximum, followed by a sharper decline. Note that
an overestimate of tropical summertime OH in MATCH
values could also lead to an overestimate of the reference
seasonal mole fraction trough. However, the good represen-
tation of the methane seasonal cycle at Samoa, located in
the remote tropics far from emissions, does not suggest a
large error in the OH seasonality. The biomass burning
(BBAM, BBAF, BBAS) regions have stronger emissions
with greatly enhanced peak values relative to the reference.
This result is consistent with the need for tropical emissions
to exceed the reference values to better replicate the
observations. The total emission for all seasonal sources
(Figure 7) shows a significant increase above the reference.
The shift of the seasonal source total to a July peak is
mostly due to RiceWet, with smaller contributions from
WetS and WetNE. Overall, the All Obs case is very similar
to the Control, indicating that the 38 additional flask
measurements listed in Table 1 do not significantly change
the average seasonal results. For the more well-constrained
regions such as WetNE and WetNW, the difference between
Control and All Obs is indiscernible.
[42] Compared to the Control (and All Obs) inversion, the
HF only case includes nearly offsetting increases and
decreases in WetNE and WetNW, respectively. The presence
of flask sites sensitive to Eurasian wetland emitting regions
accounts for part of this difference. At flask sites that are
sensitive to WetNE such as zep[14], stm[15], and bal[18],
model simulations using the reference emissions overesti-
mate observed CH4. The decrease in WetNE emissions is
smaller upon removal of these sites in case HF. The overall
northern hemispheric decease is shifted to WetNW, whose
emissions are still sensitive to the high-frequency stations at
alt[1], brw[2], and frd[4]. The WetS emission change is
much smaller than the Control case, partially due to the
fewer numbers of high-frequency stations to constrain this
large and diffuse region. The strong increase in RiceWet is
also produced, although to a slightly less extent than for the
Control. The changes in the three biomass burning regions
in HF also generally follow Control, although to a lesser
degree. The HF seasonal total change is nearly identical to
the Control, indicating that the global seasonal solution is
largely independent of the addition of flask measurements to
the high-frequency network.
[43] Figure 7 (right) shows the average monthly uncer-
tainty reduction from the inversion. WetNE and WetNW
have the greatest uncertainty reduction; these sources have
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
14 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 15 ---
several observational sites within and downwind of these
regions. Moderate uncertainty reduction is seen for Rice-
Wet, largely due to sites downwind of Asia. The least well
constrained sources are WetS and the three biomass burning
regions, which are all located in the tropics. The additional
uncertainty reduction due to the flask network can be seen
by the differences between the Control and HF error bars.
The high-frequency sites account for most of the uncertainty
reduction for WetNE, WetNW, and RiceWet, since there are
sufficient nearby high-frequency sites to constrain these
sources. The flask measurements provide additional uncer-
tainty reduction for WetS, BBAF, BBAM, and BBAS.
BBAS is the least well constrained even with flask measure-
ments. Note that the additional 38 flask sites for All Obs
lead to only small additional uncertainty reductions.
6.2.
Average Annual Results
[44] The annually averaged optimized emissions (red)
and emission uncertainties from the control inversions for
aggregated processes are compared to reference emission
(blue) in Figure 8. The numerical results for individual
process emissions for all inversion cases were given in
Table 3. The aggregated processes, including Wetland
(WetNE, WetNW, WetS, and part of RiceWet), and Biomass
Burning (BBAM, BBAF, and BBAS) were listed in Table 2.
The use of these aggregated emissions facilitates compari-
son to bottom-up estimates, as well as comparison between
the nine inversion cases. The reference error bars in
Figure 8 are the a priori uncertainties of each aggregated
set of processes. The optimized case contains two sets of
error bars. The right error bar corresponds to the spread
in inversion values for all nine inversion cases in Table 2.
The left error bar represents the uncertainty from the
Control inversion. This uncertainty is always less than the
reference uncertainty due to the influence of the
measurements in the Kalman filter. The relation-
ship s ¼
ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
P
N
t
s2
t =N
s
, where N = 60 months and st is
monthly error, is used to compute the annual average
uncertainty for each seasonal process. The aseasonal
uncertainties are taken from the last step of the Kalman
filter as shown in Figure 4. That the uncertainties from
the aseasonal components are much smaller than the
seasonal components arises for two reasons. The first is
that the aseasonal sources have very broad geographical
Figure 7.
(left) Average seasonal cycle for optimized emissions and (right) corresponding uncertainties
for the seven seasonally varying sources. Note the different scales for the different sources. The reference
emissions and uncertainties are shown in blue, and inversion cases Control, HF, and All Obs are shown in
red, green, and orange, respectively. Figure 7 (right) shows the error reduction from the reference (a
priori) uncertainties.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
15 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 16 ---
distributions and, consequently, strong sensitivities at
many observing sites. The seasonal sources, in contrast, are
more localized regionally and/or have relatively low sensi-
tivities to the observing network (e.g., biomass burning). The
second reason is that the inversion solves the aseasonal
sources as constants over the entire time period, with flux
uncertainty reduction for each new observation. The seasonal
sources are solved as monthly fluxes which adds greater
uncertainty to their 5-year averages.
[45] For the aseasonal processes, the inversion spread
(right error bars) is larger than the inversion uncertainty (left
error bars). For the seasonal processes, in contrast, the
inversion uncertainties are much larger than the spread in
inversion cases. This results from the seasonal processes
being less well constrained than the aseasonal processes for
the two reasons given above. The rightmost (yellow hanging)
error bars in Figure 8 correspond to the range of estimates
found in the literature for each process. The optimized values
fall within this literature spread, although certain processes lie
at the far end of that range. The total emissions are at the high
end of the IPPC range of 500–600 Tg yr1, which is
dependent on the prescribed OH sink, as discussed earlier.
The higher total CH4 emission contributes to the higher
emission estimates for most processes when compared to
the literature range.
[46] Compared to the reference case, the inversions
increase AnMSW and decrease Energy emissions, respec-
tively. As discussed previously, the inversion acts to
decrease the reference run overestimate of the methane
interhemispheric gradient in the reference run. An in-
crease in AnMSW (globally distributed) and a decrease in
Energy (largely northern hemispheric distribution) results
in a net shift of emissions from northern to tropical and
southern regions. Wetlands show some reduction, mostly
from decreases in northern boreal emissions. RiceWet
emissions, which are dominated by sources between 0–
30N, increase by approximately 20% overall. Individual
optimized tropical biomass burning sources are increased
to between 30 and 80% above the reference. Table 2
indicates that the emission totals are similar among the
nine optimized cases, since the total is determined largely
by the prescribed global OH field. The addition of flask
data generally decreases emissions from AnMSW and
Wetlands, and increases Energy, Biomass Burning, and
Rice emissions. Addition of flask data also leads to
greater uncertainty reduction for each process. Note that
compared to the Control case, the addition of 38 flask
sites in All Obs leads to only slightly smaller uncertainty
reductions, partially because these additional observations
are active for less than 30% of the inversion period.
6.3.
Comparison to Bottom-up Estimates
6.3.1.
Wetlands
[47] Cao et al. [1996b] used a process model that incor-
porated substrate availability and climatic variables to
estimate global wetland emissions of 92 Tg yr1. Their
total is divided into northern, temperate, and tropical con-
tributions of 23.3, 17.2, and 51.4 Tg yr1, respectively.
Using another wetland process model, Walter et al. [2001b]
estimated total wetland emissions of 260 Tg yr1, with a
25% contribution from wetlands north of 30N. In both
studies, most of the wetland emissions are located in the
tropics. Our optimally estimated global wetland flux of
143–148 Tg yr1 (including a 21% contribution from
RiceWet) lies between these two estimates but is closer to
the Cao et al. [1996b] estimate. The dominance of tropical
and southern emissions (70%) over northern emissions
(30%) from the inversion is consistent with both process
model studies, as is our estimated peak wetland emissions
during July [Cao et al., 1996b].
Figure 8.
Annual average methane emissions. Shown are reference (blue) and optimized (red)
emissions using the Control inversion. The error on the reference is the assumed a priori inversion
uncertainty. The optimized values include two error bars: left bars corresponds to the inversion
uncertainty for the Control case, and right bars to the spread of inversion results from the different
inversion cases. The rightmost errors (yellow) represent the range of emission values found in the
literature. See Table 2 for numerical values.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
16 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 17 ---
6.3.2.
Rice
[48] Over the past few decades, there has been a down-
ward trend in the estimates of rice emissions using bottom-
up methods [Mosier et al., 1998]. This trend is largely due
to the results from the increasing numbers of in situ flux
measurements used for emission extrapolation, which col-
lectively indicate lower overall rice emissions. Sass [1994]
estimated global rice emissions between 25 and 54 Tg yr1
by combining rice cultivation area and flux estimates. Using
a process based model of CH4 emissions from rice, Cao et
al. [1996a] estimated a global flux of 53 Tg yr1. The
inverted emissions for the RiceWet emitting regions of 96–
115 Tg yr1 (79% of RiceWet) are about twice current
bottom-up rice estimates, but are closer to other top-down
studies [Lelieveld et al., 1998; Hein et al., 1997]. This
discrepancy may indicate a bias in the bottom-up extrap-
olations, but it may also partially arise from the presence of
emissions within rice emitting regions that are included in
our RiceWet emission region, but not accounted for in
bottom-up studies. Bottom-up estimates focus on rice cul-
tivation areas and may ignore emissions from nearby
inundated regions formed by either natural or anthropogenic
processes. Emissions during the nongrowing periods are
also excluded in some estimates [e.g., Yan and Cai, 2003].
Incorporation of these other emissions would likely lead to
larger bottom-up estimates from areas that include rice
cultivation.
[49] In our inversion, the increase in rice emissions above
the reference are most pronounced in July and August
(Figure 7). The global seasonality of rice paddy emissions
are difficult to estimate from bottom-up approaches because
different regions have different rice planting seasons. For
flooded rice fields, peak methane emissions usually occur
several weeks after the initial flooding. The timing and
magnitude of the fluxes depend strongly on soil character-
istics, plant type, and fertilizer composition, in addition to
climatological factors. Cao et al. [1996a] used a process
model to estimate that over half of rice emissions occur
between July and October globally, with peak emissions
occurring in August. In a study of methane emissions from
China, Yao et al. [1996] estimated peak methane emissions
in June using regional classification and emission rates from
six sites.
[50] Among the high-frequency sites, hat[9] and coi[5]
are two sites most sensitive to rice emissions, and we note
that these observations are reproduced well by the opti-
mized MATCH run. The measurements at the most sensitive
flask site, wlg[32] in north central China, are sometimes
overestimated by the optimized MATCH run during July
and August. However, the seven nearby South China Sea
observations (scsn) are not overestimated by MATCH
(auxiliary material), although these stations are less sensi-
tive to RiceWet emissions during the summer due to
transport. This discrepancy among stations suggests that
additional continental based measurements in this region
should be made. The large observational error bars of
wlg[32] (±23 ppb, Table 1) due to the sampling frequency
error make future high-frequency measurements very desir-
able in this region.
6.3.3.
Biomass Burning
[51] The optimization nearly doubles the reference bio-
mass burning fluxes for all tropical regions, although the
estimated emissions are still within the current literature
range. Biomass burning emissions are poorly constrained in
the inversion due to observational undersampling in the
tropics. Bottom-up studies also have large uncertainties,
which arise from the challenge of estimating total amounts
of biomass burning. In addition, the associated CH4 emis-
sion factor depends on the fire type and can be highly
variable [Hao and Ward, 1993]. The only bottom-up value
for global CH4 emissions corresponds to our reference case
from Hao and Ward [1993]. Our inversion confirms the
bimodal behavior of BB Africa (Figure 7), but doubles its
amplitude from the reference. This bimodal behavior is
caused by seasonal shifts in the ITCZ, which lead to
alternating dry and wet periods within a single region
[Duncan et al., 2003b]. Table 3 shows that total methane
emission decrease from BB Africa, BB America, to BB
Asia. This result is qualitatively consistent with the Hao and
Ward [1993] estimate applicable for the 1980s. Duncan et
al. [2003b] find greater CO emissions from Asia compared
to the Americas during the 1990s, although the applicability
of this result to CH4 is unclear. More exact determination of
the tropical biomass burning source for CH4 will require
more sensitive observations.
6.3.4.
Aseasonal Processes
[52] The Intergovernmental Panel on Climate Change
(IPCC) [2001] includes top-down and bottom-up estimates
which range from 120 to 180 Tg yr1 for animal and waste
(landfills, MSW, etc.) emissions [Fung et al., 1991; Hein et
al., 1997; Houweling et al., 1999; Lelieveld et al., 1998;
Olivier et al., 1999]. Our optimally estimated total of 178–
211 Tg yr1 for AnMSW is at the higher end of this IPCC
estimate. The increases in AnMSW in our inversions are
mostly offset by decreased Energy emission estimates
(although this compensation is not a constraint in the
Table 3. Reference and Nine Inversion Resultsa
Inversion Region
WetNE
WetNW
WetS
RiceWet
BBAF
BBAM
BBAS
AnMSW
Energy
Reference
29 ± 43
14 ± 22
83 ± 25
117 ± 38
11 ± 12
11 ± 13
7 ± 7
169 ± 169
89 ± 89
Control
22 ± 11
12 ± 6
81 ± 22
142 ± 29
20 ± 12
14 ± 12
9 ± 7
189 ± 3
48 ± 3
HF Only
29 ± 20
8 ± 9
86 ± 24
121 ± 31
16 ± 12
14 ± 13
7 ± 7
211 ± 4
46 ± 4
HF MBL
26 ± 17
10 ± 8
81 ± 23
139 ± 30
20 ± 12
14 ± 12
9 ± 7
184 ± 3
54 ± 3
Control FullAv
22 ± 12
12 ± 6
82 ± 22
141 ± 28
22 ± 12
14 ± 12
10 ± 7
180 ± 3
54 ± 3
Control Filt
22 ± 11
12 ± 6
81 ± 22
144 ± 29
21 ± 12
14 ± 12
9 ± 7
187 ± 3
48 ± 3
Control FullAv Filt
21 ± 12
12 ± 6
82 ± 22
145 ± 29
23 ± 12
15 ± 12
10 ± 7
175 ± 3
55 ± 3
All Obs
21 ± 11
12 ± 6
80 ± 22
143 ± 28
21 ± 11
17 ± 11
10 ± 7
185 ± 3
48 ± 3
All Obs FullAv
21 ± 11
12 ± 6
80 ± 22
146 ± 28
21 ± 11
17 ± 12
9 ± 7
182 ± 3
49 ± 3
All Obs Filt
21 ± 12
12 ± 6
82 ± 22
145 ± 28
21 ± 11
17 ± 12
10 ± 7
178 ± 3
52 ± 3
aValues are in Tg CH4 yr1. Section 5.5 describes the nine inversion cases.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
17 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 18 ---
inversion). The IPCC total energy emissions range is 75–
100 Tg yr1, based mostly on bottom-up studies using
economic data. Quay et al. [1999] used global measure-
ments of 14CH4 to estimate that fossil (energy) methane
emissions are 9–27% of the total methane emissions. This
percentage range has been applied to the reference emission
total of 590 Tg yr1 to compute the absolute range in
literature values shown in Figure 8. The optimally estimated
Energy emissions, which represent coal and natural gas
emissions, total 46–55 Tg yr1. To compare to the total
bottom-up energy estimates, an additional 18 Tg yr1 from
various industrial sources (included in the reference run as
constants but not solved by the inversion) are added, which
leads to an adjusted energy range of 64–73 Tg yr1. Note
that this optimized total is still at the lower end of the range
of bottom-up energy literature estimates.
[53] Most of the literature estimates of energy emissions
correspond to time periods before the 1996–2001 time
period of our inversion. The EDGAR 3.0 data set shows
decreases in the estimated gas and coal CH4 emissions
between 1990 and 1995. At least part of the decrease was
associated with the collapse of the centrally planned econ-
omies of Eastern Europe and the former Soviet Union
(FSU). Total gas production, and to a lesser extent coal,
are reported to have decreased from these countries starting
in the early 1990s [Energy Information Administration
(EIA), 2003]. Dlugokencky et al. [1994] also attribute
observed atmospheric methane reductions to decreased gas
leakage in FSU in the early 1990s. Coal production also
decreased in China in the late 1990s due to the centrally
planned closure of many small-scale mines and a switch to
imported oil [EIA, 2003]. Coal emissions in the United
States are also considered to have decreased by over 30%
between 1990 and 2000 partially due to increased methane
recovery efforts [EIA, 2002]. These reported energy trends
may explain some of the difference between the older
bottom-up literature values and the more current inversion
results reported here.
6.4.
Errors
[54] Although the observational mismatch error (see
section 5.1) accounts for some aspects of model transport
error, the inversion otherwise assumes a perfect atmospheric
model. Systematic errors in transport, OH sink, and emis-
sion geographic patterns can contribute to biased inversion
estimates. For example, deficiencies in the modeled inter-
hemispheric transport and hence mole fraction gradient
(IHG) can influence the annual average values. As men-
tioned in section 4, however, MATCH successfully repro-
duces the observed IHG of purely anthropogenic
compounds such as SF6 and CFC-11. The ability of the
model to reproduce observations at many high-frequency
CH4 sites adds further confidence to its synoptic-scale
transport. Errors in the spatial and temporal pattern of the
OH field may also contribute to a bias error, but large errors
are not detected from the methyl chloroform simulations
over the inversion time period [Chen and Prinn, 2005]. The
seasonal methane cycle at Samoa, dominated by tropical
OH values and interhemispheric transport, and relatively
insensitive to emission sources, is also well reproduced
(Figure 6).
[55] Another possible source of error lies in the assumed
emission spatial patterns used to solve for the emission
magnitudes. This error is potentially greater for the large
aseasonal sources that have many nearby and sensitive
observing sites. A different assumed spatial emission pat-
tern would alter the modeled sensitivity at these sites to
these sources and, consequently, alter the inverted emis-
sions. The sensitivity to the emission pattern is smaller for
localized regions that do not contain observing sites, such as
the tropical biomass burning regions. For these regions, a
somewhat different assumed emission pattern would still
likely result in a similar signal at a downwind site. The
proper test of this error would be to use alternative emission
patterns, especially for the aseasonal components. The
scarcity of alternative spatial distributions, as well as the
large computational burden associated with using even a
single set of assumed flux spatial patterns, makes these tests
difficult.
7.
Interannual Results
[56] This section describes the monthly, interannual meth-
ane emissions of the seven seasonal regional sources. The
interannual variability is illustrated by computing the
monthly anomalies (deviations) of each year’s seasonal
cycles in Figure 5 from the optimized average seasonal
cycles in Figure 7. The monthly anomalies for each emis-
sion source are shown in Figure 9 for the Control (blue) and
HF (red) cases. For each source, the summation of
anomalies over all months is zero. The global variability,
computed as the summation of the individual sources, is
also plotted. RiceWet and WetNE are the largest contrib-
utors to the interannual variability, followed by WetS and
WetNW. The absolute contributions from the three bio-
mass burning sources are small, although their variability
relative to their emission magnitudes are comparable to
the other processes. In addition to the actual amplitude of
emission variability from a particular source, the sensi-
tivity of the emission estimates to the observations also
contributes to the amplitude of the inverted emission
anomalies. For example, BBAS is poorly sampled by
the observing network; this may contribute to its relatively
low emission variability compared to WetNW despite having
a similar total flux magnitude.
[57] Figure 9 shows that the Control and HF cases are
generally similar, but have significant differences for certain
months. For example, HF shows a much larger increase for
WetNE in 1998 compared to the Control case. The anoma-
lies for WetS and all three biomass burning regions are
larger for the Control case, which incorporates many more
flask observations sensitive to these regions. The global
total variabilities, which depend on global changes in CH4
observations, are similar for the two cases. Note that Figure
9 does not suggest a strong linear trend in the total seasonal
methane emissions, but does indicate a significant emission
increase in 1998, as next discussed in section 7.1. Figure 10
compares the interannual differences between the Control
and FullAv cases. The month-to-month differences between
the two cases are usually small, although differences for
certain months are larger than the very small differences
between their annual average values (Table 3). For WetNE
and WetNW, the flux anomalies between the Control and
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
18 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 19 ---
Control FullAv differ the most during the strongly emitting
summer months. The Control anomaly is more negative
than Control FullAv by about 10 Tg yr1 for WetNE in May
1998. The Control FullAv peak anomaly is also centered in
June, rather than July, 1998. Figure 10 shows the conflicting
results that arise when using different sampling strategies
are used to compute flask monthly means in the model and
observations, and illustrates the additional uncertainty in-
troduced when using flask observations in the inversion.
[58] We further compared the Control case and the All
Obs case (Figure 11), showing the effect of adding the 38
additional flask stations which were active for less than the
70% of the months of the inversion (Table 1). Most of these
additional observations affect the inverted values for the
poorly constrained biomass burning emissions, rather than
the already well-constrained emissions such as WetNE and
WetNW. The shipboard measurements in the Pacific and the
South China sea, which account for most of the additional
observations in this 38-site database (Figure 6), are sensitive
to BBAM and BBAS emissions, respectively. The new flask
site in Namibia (nmb[63]) further constrains BBAF. These
additional flask measurements increase the interannual
variability (IAV) for biomass burning emissions. The peak
heights for BBAF, for example, are greater in the All Obs
case than the Control case (especially for BBAF). Note that
the additional stations generally enhance the magnitudes of
the monthly anomalies rather than alter their timing. This
suggests that although the smaller network can capture the
general characteristics of the IAV, more sensitive observing
sites are needed to accurately determine the full magnitude
of emission anomalies. This argues for more long-term
measuring sites directly downwind of these tropical regions.
The total emissions are nearly the same between the Control
and All Obs cases. We have also compared the Obs All and
Obs All FullAv cases (not shown), which compares the
effects of the different MATCH flask sampling strategies.
The differences between the two cases are similar in
magnitude to those between Control and FullAv (Figure 11.)
[59] Table 4 contain the flux anomalies averaged for each
year for the Control, HF, and All Obs inversion cases. Total
year-to-year emission changes fluctuate between +33 and
16 Tg yr1. This compares reasonably well with estimates
by Cunnold et al. [2002], who used a semi-inverse method
to compute global annual changes of up to ±37 Tg yr1.
The combined wetland values fluctuate between +19 and
15 Tg yr1. Walter et al. [2001b] used a wetland
process model based on changes in precipitation and
temperature to estimate annual changes of approximately
±20 Tg yr1. Using the difference between Alert and
Fraserdale mole fractions between 1990 and 1998, Worthy
et al. [2000] estimated a wetland emission variability in
the Hudson Bay Lowland (HBL), which is located within
our WetNW region, between ±0.23 to ±0.5 Tg yr1.
Assuming the HBL represents about 10% of northern
Figure 9.
Monthly mean anomalies (from 5-year mean value) for Control (blue) and HF (red) inversion
cases. Note the different vertical scales for different sources. Anomalies are repeated on the right side but
with identical vertical scales to emphasize contributions of individual regional sources to the total change.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
19 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 20 ---
hemispheric wetlands [Worthy et al., 2000], this corre-
sponds to a range of ±2.3 to ±5.0 Tg yr1 for Northern
wetlands, which is somewhat smaller than the ±7 Tg yr1
range for WetNE + WetNW emissions determined here.
There have been fewer studies of the interannual vari-
ability in rice emissions. However, rice emissions can be
expected to vary similarly to wetlands emissions, as both
are influenced by climatic conditions (although the rice
flooding stage is managed). Sass et al. [2002] measured
emissions at a single site in the U.S over nine years and
observed a year-to-year flux variability of approximately
±50% of the annual mean over the entire period. Biomass
burning also fluctuates significantly from year to year
[e.g., Duncan et al., 2003b], although much less is
known about associated methane emissions.
7.1.
Global Methane Increase in 1998
[60] The ENSO event that occurred in late 1997 and 1998
influenced climate on a global scale [e.g., Bell et al., 1999].
Global annual mean temperatures are considered to be
highest for 1998 since the advent of reliable, direct measure-
ments. In addition, both positive and negative precipitation
anomalies led to regional flooding (e.g., central China in
June–July 1998, Indonesia in the second half of 1998) and
drought (e.g., Indonesia in the first half of 1998), respec-
tively. These climate changes likely contributed to the
dramatic increase of global CH4 mole fractions during this
time [Dlugokencky et al., 2001]. This large increase is
clearly captured by the high-frequency observations at
Samoa, shown in Figure 12. Figure 12 also shows the
improvement in the model simulation when using Control
optimized emissions compared to the reference. The inter-
annual total methane flux anomaly in Figure 9 (right) shows
a global CH4 emission increase in 1998 followed by a steep
decline during early 1999. Table 4 lists the optimized CH4
anomalies between 1996 and 2001 for the Control inver-
sion, which shows the unusually large total anomaly of
+33 Tg yr1 during 1998. This total anomaly falls
between previously reported methane emission anomalies
of +24 Tg yr1 and +37 Tg yr1 as estimated by Dlugokencky
et al. [2001] and Cunnold et al. [2002], respectively, using
hemispherically averaged mole fractions.
[61] Average surface temperature anomalies of approxi-
mately +1C occurred over the northern and tropical land-
masses [Bell et al., 1999]. In general, increased wetland soil
Figure 10.
Monthly mean anomalies for Control and
Control FullAv inversions. This plot compares the effects of
using different flask sampling strategies.
Figure 11.
Monthly mean anomalies for Control and All
Obs inversions. This plot shows the effect of using the
extended flask stations.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
20 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 21 ---
temperatures enhance methanogenic activity, leading to
increased CH4 production and emission. Anoxic environ-
ments are also necessary for methanogens to survive, and
increased precipitation (expected during this El Nin˜o year)
generally increases anoxia in soils. Our inversely estimated
emission anomalies are consistent with increased methane
flux from wetlands in 1998. Rice emissions also show an
increase in 1998 in our inversions; this result is not
surprising because both natural and rice processes should
have similar sensitivities to temperature and, to a lesser
extent, precipitation. Table 4 contains the annual variability
for the HF and All Obs inversion cases. Across all three
inversion cases, the aggregated Wetlands and Rice 1998
anomaly flux range between 14 and 19 Tg yr1 and 8 and
18 Tg yr1, respectively.
[62] Dlugokencky et al. [2001] used the wetland CH4
emission model of Walter et al. [2001a] to estimate 1998
flux anomalies. The model was driven by NCEP soil
temperature and precipitation data from 1980 to 1999.
The spatial distribution of the wetland model was further
subdivided into northern and tropical wetlands based on the
Matthews et al. [1987] distribution. This distribution is
therefore spatially similar to our WetNE, WetNW, and WetS
distributions, as well as the wetland contribution to the
RiceWet distribution. In the bottom-up model, methane
fluxes increased by approximately 20% for a 1C temper-
ature increase, and 8% for a 20% precipitation increase. The
simulated 1998 CH4 emission anomalies for northern and
tropical wetlands were 12 and 13 Tg yr1, respectively,
compared to the 1980–1999 average. These values are
slightly larger than the northern and tropical wetland
anomalies of 5–10 and 8–10 Tg yr1 for the three
inversion cases in Table 4 (recall that 21% of RiceWet are
attributed to wetlands). The inversion also produces a rice
emission anomaly of 8–18 Tg yr1, which may include
further overlapping wetland emissions, as discussed in
section 3. Mikaloff Fletcher et al. [2004a], using an inverse
approach which incorporated 13CH4 observations, also
suggest that wetlands dominated the 1998 observed meth-
ane increase. Unfortunately, no bottom-up estimates of rice
emission anomalies for 1998 are available for comparison.
[63] Biomass burning is known to have increased for
certain regions during 1998. With the exception of BBAF
for the All Obs case, all biomass burning CH4 emission
anomalies are positive in 1998, although the magnitudes are
smaller than for Wetland and Rice. Southeast Asia experi-
enced ENSO related biomass burning between August 1997
and April 1998. Levine [1999] estimated an anomalous
Indonesian biomass burning emission of 1.9 Tg during the
late fall of 1997. The 1998 Control emissions in Table 4
indicates an anomaly of 0.8 Tg yr1 from BB Asia, with
most of the increase occurring in early 1998 as shown in
Figure 9. The Obs All case shows a similar flux behavior
(Figure 11), but with stronger peak and trough emissions.
Unfortunately, many of these stations in Obs All, such as
the ship tracks in the South China Sea and in the Pacific
were not active during the 1998 El Nin˜o event. The closest
high-frequency site, Hateruma (hat[9]), does not indicate an
obvious methane enhancement during the Southeast Asian/
Indonesian fires of August–December 1997, probably due
to its position upwind of these fire during these months.
Central America also experienced large amounts of burning
between April and June 1998 due to an ENSO related
drought [Duncan et al., 2003b]. The exact timing of this
event is only weakly reproduced by the inversion results
(Figure 9), which is likely due to the lack of nearby
observing sites. For the three inversion cases described in
Table 4, the overall range of 1998 from BBAM is between
+0.5 to +1.9 Tg yr1. Using satellite data, Duncan et al.
[2003b] and van der Werf et al. [2004] found smaller
positive anomalies of biomass burning in Africa compared
to Asia and America in 1998. The three inversion cases
show a 1998 methane BBAF emission anomaly ranging
between 1.1 to +2.6 Tg yr1. A small region in Eastern
Siberia experienced biomass burning between July and
September 1998 due to unusually warm and dry conditions.
Using AVHRR data, Kasischke and Bruhwiler [2002]
estimated methane releases of 2.9–4.7 Tg. Since these
particular wildfire areas are not explicitly represented, the
inversion would likely attribute this possible contribution to
enhanced emissions of WetNE.
[64] Figure 12 shows that the optimized MATCH forward
run compensates for about half the difference between the
reference run and observations. Given the remote location
of Samoa to strongly emitting regions, it is unlikely that any
of the emissions (including aseasonal sources) can further
compensate for this discrepancy without first adversely
influencing stations more sensitive to these emissions. The
mismatch between observed and optimized mole fractions
at Samoa suggests that decreased OH concentrations may
also have contributed to the increased 1998 methane growth
rate. Using CH3CCl3 measurements, Prinn et al. [2001]
Table 4. Inversion Annual Anomalies for Control, HF Only, and All Obs Casesa
Region
1996
1997
1998
1999
2000
2001
Control
HF
Only
All
Obs
Control
HF
Only
All
Obs
Control
HF
Only
All
Obs
Control
HF
Only
All
Obs
Control
HF
Only
All
Obs
Control
HF
Only
All
Obs
WetNE
2.7
2.6
2.8
1.9
0.3 1.3
1.5
8.9
1.6
2.8
6.7 3.5
2.0
0.3
2.3
2.5
0.5
2.8
WetNW
1.9
0.3 1.9
1.5
1.2
0.7
2.2
0.3
2.3
0.2
2.0
0.5
1.6
2.3
1.4
0.4
1.5
0.3
WetS
2.5
2.4
2.6
3.2
5.2
3.0
5.3
8.3
7.4
3.2
3.5 4.7
10.7
9.0 11.9
2.8
3.4
3.7
RiceWet
7.5
10.3 5.6
4.8
4.5 4.1
22.3
10.3
20.6
1.4
4.5 0.6
1.9
3.0
1.4
10.4
2.9 11.7
BBAF
0.8
1.0
1.4
1.1
0.8
1.5
0.5
2.6 1.1
0.5
0.0 0.2
0.6
1.6
1.5
2.3
2.8
3.1
BBAM
1.6
2.7
2.6
1.4
1.6 0.7
0.5
1.9
1.3
0.5
1.6
2.9
3.7
3.6
4.4
0.4
1.0
1.7
BBAS
0.1
0.1
0.0
1.2
0.8 0.5
0.8
0.5
0.3
0.3
0.0 0.4
0.9
0.5
0.8
0.2
0.2
0.2
Wet Total
1.7
2.1
2.3
1.8
3.3
1.5
13.7
19.0
15.5
6.1
7.1 7.9
13.8
10.9 15.3
2.7
2.1
3.7
Rice Total
5.9
8.2 4.4
3.8
3.6 3.2
17.6
8.1
16.3
1.1
3.5 0.5
1.5
2.4
1.1
8.2
2.3
9.2
BB Total
2.3
3.6
4.0
1.3
1.7
0.3
1.8
4.9
0.5
0.8
1.6
2.3
3.4
4.7
2.1
2.8
3.9
4.9
Total
1.9
6.7
1.9
0.7
1.4 1.3
33.1
32.0
32.3
6.5
5.2 6.1
15.7
13.2 16.3
8.3
8.3 10.4
aValues are in Tg yr1.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
21 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 22 ---
deduced anomalously low OH levels in 1997 and 1998. The
driver for reduced OH concentration may be increased
cloudiness (consistent with increased precipitation) which
would reduce the amount of sunlight necessary for OH
production. This effect would be pronounced for tropical
regions where OH concentrations dominate CH4 loss, such
as Samoa. Novelli et al. [2003] linked the strong 1997–
1998 increases in wildfires to large globally observed CO
mole fractions; they further calculated significant decreases
in global OH due to the large CO increase [see also Duncan
et al., 2003a].
7.2.
Model Errors
[65] The model errors that affect the interannual results
differ from those that affect the annual average results.
Errors in the interhemispheric exchange rate and the model
OH concentrations will bias the annual average modeled
atmospheric CH4 distribution, and consequently, the
inverted annually averaged emissions. In contrast, the
interannual emissions results depend less on average CH4
distributions and more on observed and simulated monthly
mole fraction anomalies. The deduced interannual emis-
sions are thus less affected by large-scale transport biases
which affect average CH4 distributions, and more affected
by errors in the short-term changes in the spatial flux
patterns. For example, the inversion may not have been
able to identify the boreal fires in Siberia during July–
September 1998 because this pattern was not explicitly
solved for. Instead, the inversion would have attributed this
possible increase to WetNE emissions. The individual
seasonal patterns are also assumed to be representative for
particular months, which becomes less true for very short
time periods. A solution is to use even smaller process
based emitting regions, ultimately at the grid resolution of
the model, but the sparsity of observations make such an
approach currently intractable at the global scale. A possible
solution is to incorporate mesoscale models for regions
sampled by nearby high-frequency observations [e.g.,
Kleiman and Prinn, 2000].
[66] The assumed constancy of aseasonal emissions is
another simplification in this inversion. Although bottom-up
studies indicate that the climate-driven seasonal processes
considered here dominate the interannual variability, ani-
mals, waste, and energy emissions can also change year to
year. Large fluctuations in these processes for the 1996–
2001 time period have not been published. Possible changes
in these aseasonal sources should, however, certainly be
modeled for inversions over longer time periods. It has been
hypothesized that reductions in gas and coal emissions led
to the global decline of the CH4 growth rate in the early
1990s [Dlugokencky et al., 1994; Law and Nisbet, 1996].
As a sensitivity test, the inversion was also modified to
solve for the linear aseasonal emission trends in addition to
their average values. The inversion could not adequately
constrain these small trends, because the trend sensitivities
to the global measuring network were small and not
sufficiently distinct. A trend of +1 Tg yr1 added to each
of the aseasonal source over the inversion period also had a
negligible impact.
Figure 12.
Observed (black) and simulated (red and blue) high-frequency mole fractions at Samoa. Red
and blue correspond to MATCH simulations using reference and Control optimized emissions,
respectively. Note the large increase in mole fraction observed in 1998. The optimized simulation (red)
captures most, but not all, of this increase.
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
22 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License


--- Page 23 ---
[67] Another potential model error involves using an
annually repeating OH field, as discussed for the 1998
observed CH4 mole fraction increase. Using methyl chlo-
roform (MCF) data, Prinn et al. [2001] determined an
increase followed by a decrease in the global average OH
field between the late 1970s and the late 1990s. Incorpora-
tion of an optimized interannually varying OH field, derived
from MCF studies for example, into the MATCH reference
and sensitivity runs could then account for interannual
changes in the sink.
8.
Summary and Conclusions
[68] We have carried out an inverse modeling study of
methane fluxes incorporating three new elements not
previously combined in methane inversions: (1) high-
frequency CH4 observations, (2) interannual transport in
the atmospheric transport model, and (3) the Kalman
filter solution of interannually varying monthly fluxes.
The first two elements are critical to determining methane
emissions at higher space and time resolution. The
development of new inversion techniques, such as done
here using the Kalman filter, will also play a role in
combining more realistic model and observational data.
[69] The annual average results of this study indicate a
reduction in northern hemispheric emissions and an increase
in tropical and southern hemispheric emissions compared to
the reference. The deduced average seasonal maximum in
emissions is shifted to July from the August reference, and
is more intense. This result is dominated by phase changes
and increases from rice emitting regions relative to the
reference. The inverted emissions for wetlands compare
well with other recent estimates. Inverted energy emissions
are on the low side of prior estimates, but are consistent
with expected recent decreases in methane emissions. The
inversion produces rice emissions nearly twice the current
bottom-up estimates. Part of the difference may be that the
inversion solves for emissions from the entire rice emitting
region, and includes those wetland emissions not formally
associated with rice cultivation. Global wetlands and rice
emissions dominate the observed CH4 increase in 1998,
consistent with a bottom-up study using a wetland processes
model [Dlugokencky et al., 2001]. The inverted methane
flux increases for 1998 cannot fully reproduce the CH4
growth rate at Samoa, suggesting that a decreased OH sink
may also have contributed to the observed increase.
[70] On the basis of different sensitivity tests, the impact
on the inversion results is influenced in order of decreasing
influence by the following: observational network, modeled
monthly means using exact observation time versus all
model time steps, and filtered versus unfiltered observa-
tions. Several sources of model error are not easily input
into the inversion framework. The use of alternative or
additional emission patterns (such as for boreal biomass
burning and soil uptake) can potentially affect the inversion
results. As bottom-up studies continue to better define the
spatial distributions of these fluxes, their inclusion into the
top-down approach should result in more accurate flux
estimates. This study also assumed an annually repeating
OH sink. Future studies should incorporate interannually
varying OH fields determined by other techniques, such as
inverse studies of methyl chloroform. This study and Chen
and Prinn [2005] have shown the importance of larger-scale
transport IAV. Another important aspect of model physics is
the planetary boundary layer (PBL) transport for sites in
strongly emitting regions. As more of these sensitive sites
become active, tests of the accuracy of the PBL transport
will become more crucial, perhaps using tracers that have
better-known emission magnitudes compared to methane.
[71] The inversion leads to significant uncertainty reduc-
tion for northern wetlands, moderate reduction for rice, and
small reductions for the tropical biomass burning and
swamp sources. This clearly suggests that more methane
observations are needed in the tropics to constrain biomass
burning and wetland emissions from this region. The
determination of the optimal placement of future observing
sites (i.e., network design) have been conducted for CO2
[Patra and Maksyutov, 2002; Glooret al., 2000]. The location
of optimal observing sites depends not only the species of
interest, but also the model to be used for the inversion. The
most sensitive sites for emission estimates are often the most
difficult to model, due to subgrid-scale effects. The network
design algorithm should strike a balance between emission
sensitivity, and model errors and resolution when determining
optimal sites. An inversion algorithm similar to that used here
could be used in principle to estimate the emission
uncertainty reductions, and hence usefulness, of new
observing sites.
[72] Methane isotopes (13CH4, 14CH4, CH3D) may also
aid in resolving isotopically distinct processes that have
otherwise similar spatial and temporal patterns [Mikaloff
Fletcher et al., 2004a, 2004b]. The use of multiple tracers
could aid in the attribution of methane sources, such as done
for emissions off the Asian continent with other anthropo-
genic gases [Bartlett et al., 2003]. Satellite data from
sources such as MOPPITT and TES can provide global
coverage of methane mole fractions, although information is
often coarse in the vertical, and is less precise than ground
based measurements.
[73] The use of a global, top-down approach to solve for
specific interannual CH4 flux regions/processes has reached
the stage of allowing approximately monthly time resolu-
tion. Further increasing the time resolution of the optimi-
zation may require the optimization of spatially smaller
emission regions, of which the smallest is ultimately the
model grid resolution. The appropriate tool in this case
would be an adjoint of the MATCH, which could in
principle efficiently optimize the flux from each model grid
cell. This technique has been applied to CH4 [Houweling et
al., 1999] at coarse resolution. However, this technique
would still require use of prior information in order to
successfully constrain specific sources to avoid severe ill
conditioning. Regional-scale modeling offers a complement
to the global approach; for example, Janssen et al. [1999]
and Wang and Bentley [2002] have estimated regional
methane emissions using mesoscale models and high-
frequency observations. Finally, future inverse studies
may seek to estimate model parameters that link methane
fluxes to climatological variables, as contained in CH4
process models. This will lead to a true coupling between
top-down and bottom-up approaches.
[74]
Acknowledgments.
We thank Phil Rasch and Brian Eaton for
help in using the MATCH model and Mark Lawrence for the OH field. We
D10307
CHEN AND PRINN: ATMOSPHERIC METHANE INVERSION
23 of 25
D10307
 21562202d, 2006, D10, Downloaded from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JD006058, Wiley Online Library on [03/05/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
