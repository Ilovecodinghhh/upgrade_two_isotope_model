#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026

@author: yufan bao
"""

# Code for Monte Carlo analysis of onebox delta13C-CH4 mass balance
from pathlib import Path
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# I/O paths
# Inputs are expected under ./rel (as provided by the user).
# Outputs are written to ./Output (created if missing).
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
REL_DIR = BASE_DIR / "rel"
DATA_DIR = REL_DIR / "data"
SRC_DIR = REL_DIR / "output"
OUT_DIR = BASE_DIR / "Output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# FOR SENSITIVITY TESTS: 
    # Make sure to smooth the output of each sensitivity test and save with appropriate filenames. 
    # The sensitivity tests are as follows
        # 1. Reduced Biomass burning (compiled as _BBDrop)
        # 2. Reduced Cl sink (compiled as _RedCl)
        # 3. Increase OH sink (compiled as _IncOH)
        # Tests with stable microbial and fossil fuel source signatures can be conducted manually by 
        # using mean of time series and saving output. The outputs of these tests have been added to the results folder.
        # The same is true for tests using the global mean FF source signature using CT-CH4 posterior flux-weighting, and saueressig OH fractionation.


#%% First, load data

# For dD
# Load updated annual mean DEI
glob_ann_dD_path = DATA_DIR / "GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx"
glob_ann_dD_df = pd.read_excel(glob_ann_dD_path)
glob_ann_dD_num = glob_ann_dD_df.apply(pd.to_numeric, errors="coerce")
glob_ann_dD_years = glob_ann_dD_num.iloc[:, 0].to_numpy(dtype=np.float64)
glob_ann_dD = glob_ann_dD_num.iloc[:, 1].to_numpy(dtype=np.float64) - 0.5
glob_ann_dD_unc = glob_ann_dD_num.iloc[:, 2].to_numpy(dtype=np.float64)
glob_ann_dD_uncR = glob_ann_dD_num.iloc[:, 4].to_numpy(dtype=np.float64)

# # Load all dD DEI iterations
# Expect metadata in the first 5 columns (year/mean/unc/..); MC iterations follow.
dD_AnnAvg_matrix = glob_ann_dD_num.iloc[:, 5:].to_numpy(dtype=np.float64)

# Load GML global annual means for CH4
CH4data = pd.read_excel(DATA_DIR / "GML_CH4_AnnualMean.xlsx").to_numpy()
CH4data = CH4data[21:,:]
# Note: This first load is not used later - saving for reference, but will be overwritten


# Load Carbon tracker methane
data2 = pd.read_excel(DATA_DIR / "CarbonTracker_CH4.xlsx")
years = data2.iloc[:,0].values
bbCT = data2.iloc[:,9].values #Prior emissions (GFED4)

# Load annual mean source signatures
BB_dD_data = pd.read_csv(SRC_DIR / "BB_dD_annual.csv", delimiter=',', header = None)
Mic_dD_data = pd.read_csv(SRC_DIR / "Mic_dD_AnnGlob.csv", delimiter=',', header = None)
Mic_dD_MC_trends = pd.read_csv(SRC_DIR / "Mic_dD_MC.csv", delimiter=',', header = None) 
Mic_dD_MC = Mic_dD_MC_trends.iloc[6:,1:]

# Load FF source sigantures
FF_dD_data = pd.read_csv(SRC_DIR / "FF_dD_GlobUnc.csv", delimiter=',')
FF_dD_MC_CTCH4_data = pd.read_csv(SRC_DIR / "FF_dD_GlobMC_CTCH4.csv", delimiter=',', header = None)
FF_dD_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_dD_GlobMC_EDGAR.csv", delimiter=',')
FF_dD_MC_EDGAR = FF_dD_MC_EDGAR_data.iloc[34:,1:]
FF_dD_MC_CTCH4 = FF_dD_MC_CTCH4_data.iloc[7:,1:]
FF_dD_MC_CTCH4 = pd.concat([FF_dD_MC_CTCH4, FF_dD_MC_CTCH4.iloc[[-1]]], ignore_index=True)

# Pad FF_dD_MC_EDGAR and FF_dD_MC_CTCH4 to 24 rows for d13C alignment
if FF_dD_MC_EDGAR.shape[0] < 24:
    pad_count = 24 - FF_dD_MC_EDGAR.shape[0]
    first_row = FF_dD_MC_EDGAR.iloc[0:1, :]
    pad_rows = pd.concat([first_row] * pad_count, ignore_index=True)
    FF_dD_MC_EDGAR = pd.concat([pad_rows, FF_dD_MC_EDGAR], ignore_index=True)
if FF_dD_MC_CTCH4.shape[0] < 24:
    pad_count = 24 - FF_dD_MC_CTCH4.shape[0]
    first_row = FF_dD_MC_CTCH4.iloc[0:1, :]
    pad_rows = pd.concat([first_row] * pad_count, ignore_index=True)
    FF_dD_MC_CTCH4 = pd.concat([pad_rows, FF_dD_MC_CTCH4], ignore_index=True)
# Calculate CTCH4 statistics
FF_dD_mean_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
FF_dD_mean_EDGAR = FF_dD_MC_EDGAR.iloc[:, 1:].mean(axis=1).to_numpy()
FF_dD_std_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()
R_FF_dD_mean_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
R_FF_dD_std_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy() 

# For d13C
# Load INSTAAR DEI d13C-CH4
C13data = pd.read_excel(DATA_DIR / "ch4c13_nh_sh_mean.xlsx").to_numpy()
# Calculate annual average
glob_dates = C13data[:,0]
glob_mean = C13data[:,1]
df = pd.DataFrame({'Date': glob_dates, 'Value': glob_mean})
years = np.floor(df['Date']).astype(int)
annual_avg = []
for year in np.unique(years):
    year_values = df.loc[years == year, 'Value']
    annual_avg.append({'Year': year, 'Mean': year_values.mean()})
Glob_annual_avg = pd.DataFrame(annual_avg)
d13C_glob = Glob_annual_avg.iloc[1:, 1].values # starting in 1999
years = Glob_annual_avg.iloc[1:, 0].values

# Load GML global annual means for CH4
CH4data = pd.read_excel(DATA_DIR / "GML_CH4_AnnualMean.xlsx").to_numpy()
CH4data = CH4data[11:,:]  # Skip early years
CH4 = CH4data[4:28,1] # Starting in 1999 (24 years: 1999-2022)
CH4year = CH4data[4:28,0]

# Load global annual mean DEI iterationslines
d13C_glob_iterations_data = np.loadtxt(str(DATA_DIR / "d13C_dei_compiled.txt"))
d13C_glob_iterations = d13C_glob_iterations_data[1:,1:]  # Start from row 1 (1999) to match d13C_glob, and skip first column (years)

# Load Carbon tracker methane
data2 = pd.read_excel(DATA_DIR / "CarbonTracker_CH4.xlsx")
yearss = data2.iloc[:,0].values
micCT = data2.iloc[:,7].values
ffCT = data2.iloc[:,3].values
bbCT = data2.iloc[:,9].values #Prior emissions (GFED4)

# Load annual mean Mic and BB source signatures(13C)
BB_d13C_data = pd.read_csv(SRC_DIR / "BB_d13C_annual.csv", delimiter=',', header = None) # Calculated from Luo C3 C4 map and CTCH4 BB emissions
Mic_d13C_data = pd.read_csv(SRC_DIR / "Mic_d13C_annual.csv", delimiter=',', header = None)
Mic_d13C_MC_trends = pd.read_csv(SRC_DIR / "Mic_d13C_MC.csv", delimiter=',', header = None) 
Mic_d13C_MC = Mic_d13C_MC_trends.iloc[:,1:]

# Load FF source sigantures(13C)
FF_d13C_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobUnc.csv", delimiter=',')
FF_d13C_MC_CTCH4_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobMC_CTCH4.csv", delimiter=',', header = None)
FF_d13C_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobMC_EDGAR.csv", delimiter=',')
FF_d13C_MC_EDGAR = FF_d13C_MC_EDGAR_data.iloc[28:,1:] # starting in 1999
FF_d13C_MC_CTCH4 = FF_d13C_MC_CTCH4_data.iloc[1:,1:]
FF_d13C_MC_CTCH4 = pd.concat([FF_d13C_MC_CTCH4, FF_d13C_MC_CTCH4.iloc[[-1]]], ignore_index=True)
# Calculate CTCH4 statistics
FF_d13C_mean_EDGAR = FF_d13C_MC_EDGAR.iloc[:, 1:].mean(axis=1).to_numpy()
FF_d13C_mean_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
FF_d13C_std_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()
R_FF_d13C_mean_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
R_FF_d13C_std_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()


#%% Load model components

# Sink KIEs for carbon (13C)
OH_KIE_13C = 1.0054    # Saueressig et al., 2001 is 1.0039; Cantrell is 1.0054
Cl_KIE_13C = 1.066     # Saueressig
Strat_KIE_13C = 1.003  # Saueressig; Lassey et al., 2007
Soil_KIE_13C = 1.0201  # Average of Snover & Quay; Tyler; Reeburgh

# Sink KIEs for hydrogen (D)
OH_KIE_D = 1.294       # Saueressig et al., 2001 (Whitehill-Joelson Avg ~ 1.327)
Cl_KIE_D = 1.52        # Saueressig
Strat_KIE_D = 1.179    # Dyonisius et al., 2020; Beck et al., 2018
Soil_KIE_D = 1.083     # Snover and Quay

# Sink fractional contributions (sum to ~1; CarbonTracker documentation / Thanwerdas)
OH_Sink = .835
OH_Sink_Than = .899
Cl_Sink = .035
Cl_Sink_Than = .006
Strat_Sink = .07
Strat_Sink_Than = .03
Soil_Sink = .06
Soil_Sink_Than = .065

# Effective (bulk) KIE for the total sink (per isotope system)
Sink_13C = OH_KIE_13C * OH_Sink + Cl_KIE_13C * Cl_Sink + Strat_KIE_13C * Strat_Sink + Soil_KIE_13C * Soil_Sink
Sink_D = OH_KIE_D * OH_Sink + Cl_KIE_D * Cl_Sink + Strat_KIE_D * Strat_Sink + Soil_KIE_D * Soil_Sink

# Calculate change in KIE if Cl sink decreases
Cl_Sink_red = np.concatenate((np.full(6, Cl_Sink), np.linspace(Cl_Sink, .011, 18)))
Sink_Rem = Cl_Sink - Cl_Sink_red
OH_Sink_RedCl = OH_Sink + Sink_Rem * OH_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Strat_Sink_RedCl = Strat_Sink + Sink_Rem * Strat_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Soil_Sink_RedCl = Soil_Sink + Sink_Rem * Soil_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Check = OH_Sink_RedCl + Strat_Sink_RedCl + Soil_Sink_RedCl + Cl_Sink_red
Sink_13C_RedCl = (
    OH_KIE_13C * OH_Sink_RedCl
    + Cl_KIE_13C * Cl_Sink_red
    + Strat_KIE_13C * Strat_Sink_RedCl
    + Soil_KIE_13C * Soil_Sink_RedCl
)
Sink_D_RedCl = (
    OH_KIE_D * OH_Sink_RedCl
    + Cl_KIE_D * Cl_Sink_red
    + Strat_KIE_D * Strat_Sink_RedCl
    + Soil_KIE_D * Soil_Sink_RedCl
)

# Calculate change in KIE if Cl sink increases (Van Herpen et al., 2023)
Cl_Sink_inc_frac = np.array([1, 1, 1, 1, 1, 1, 1, 0.97, 0.94, 0.58, 0.70, 1.13, 0.79, 1.09, 0.80, 1.48, 1.86, 2.08, 1.62, 1.85, 1.85, 1.85, 1.85, 1.85])
Cl_Sink_inc = Cl_Sink_inc_frac * Cl_Sink
Sink_13C_Cl_inc = (OH_KIE_13C * OH_Sink + Cl_KIE_13C * Cl_Sink_inc + Strat_KIE_13C * Strat_Sink + Soil_KIE_13C * Soil_Sink) / (
    OH_Sink + Cl_Sink_inc + Soil_Sink + Strat_Sink
)
Sink_D_Cl_inc = (OH_KIE_D * OH_Sink + Cl_KIE_D * Cl_Sink_inc + Strat_KIE_D * Strat_Sink + Soil_KIE_D * Soil_Sink) / (
    OH_Sink + Cl_Sink_inc + Soil_Sink + Strat_Sink
)

# Create increasing OH following Olaf et al., trends
OH_ann_change = OH_Sink*.003
OH_Sink_inc = np.concatenate((np.full(6, OH_Sink), np.linspace(OH_Sink, OH_Sink + OH_ann_change * (18 - 1), 18)))
# Calculate net effective KIE under increasing OH (renormalize by total sink)
Sink_13C_OH_inc = (OH_KIE_13C * OH_Sink_inc + Cl_KIE_13C * Cl_Sink + Strat_KIE_13C * Strat_Sink + Soil_KIE_13C * Soil_Sink) / (
    OH_Sink_inc + Cl_Sink + Strat_Sink + Soil_Sink
)
Sink_D_OH_inc = (OH_KIE_D * OH_Sink_inc + Cl_KIE_D * Cl_Sink + Strat_KIE_D * Strat_Sink + Soil_KIE_D * Soil_Sink) / (
    OH_Sink_inc + Cl_Sink + Strat_Sink + Soil_Sink
)


# Create sink fractionation scenario based on varying OH (no chemistry)
# Constants for model
C13Std = 0.011113                   # standard for 13C-12C isotopic ratio (IUPAC 2024). Old: .0112020968 
DStd = 0.00015576      # standard for D-H isotopic ratio
Watm = 28.96           # Molecular weight of atmosphere (g/mole)
Matm = 5.15 * 10**21   # Mass of atmosphere (g)
Lifetime = 9
Lifetime_OH_inc = Lifetime / (OH_Sink_inc + Cl_Sink + Strat_Sink + Soil_Sink)
Lifetime_Cl_inc = Lifetime / (OH_Sink + Cl_Sink_inc + Strat_Sink + Soil_Sink)
PT = 2.815 # Conversion factor for ppb to Tg using the molar mass of the atmosphere. 
# ---------------------------------------------------------------------------
# Isotope conventions (unified)
# - R is the isotope ratio: 13C/12C for carbon, D/H for hydrogen
# - delta (‰): δ = (R/Rstd - 1) * 1000
# - f is the heavy-isotope fraction: f = R / (1 + R)
# Mass-balance is linear in heavy-isotope *amount* (or fraction), not in δ.
# ---------------------------------------------------------------------------
def delta_to_R_d13C(delta_permil):
    return (delta_permil / 1000.0 + 1.0) * C13Std

def delta_to_R_dD(delta_permil):
    return (delta_permil / 1000.0 + 1.0) * DStd

def R_to_fraction(R):
    return R / (1.0 + R)

def fraction_to_R(f):
    return f / (1.0 - f)

def delta_to_fraction_d13C(delta_permil):
    return R_to_fraction(delta_to_R_d13C(delta_permil))

def delta_to_fraction_dD(delta_permil):
    return R_to_fraction(delta_to_R_dD(delta_permil))

# Conversion from R ratio back to delta notation
def R_to_delta_d13C(R):
    """Convert 13C/12C ratio to delta-13C notation"""
    return ((R - C13Std) / C13Std) * 1000

def R_to_delta_dD(R):
    """Convert D/H ratio to delta-D notation"""
    return ((R - DStd) / DStd) * 1000

def fraction_to_delta_d13C(f):
    return R_to_delta_d13C(fraction_to_R(f))

def fraction_to_delta_dD(f):
    return R_to_delta_dD(fraction_to_R(f))


#%% Run mixing ratio model forward to determine global source strength (no sink uncertainty)

# 1 box model (calculates 1999 source and on)
SumSource = np.zeros(len(CH4) - 1)
SumSource_OH_inc = np.zeros(len(CH4) - 1)
SumSource_Cl_inc = np.zeros(len(CH4) - 1)
SumCH4change = np.zeros(len(CH4) - 1)
SumSink = np.zeros(len(CH4) - 1)

for i in range(len(CH4) - 1):
    # Source stable sink
    Source = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / Lifetime
    # Source increased sink due to increased OH Sink 
    Source_OH_inc = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / Lifetime_OH_inc[i]
    # Source increased sink due to increased Cl sink
    Source_Cl_inc = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / Lifetime_Cl_inc[i]
    CH4change = CH4[i + 1] - CH4[i] 
    SumSource[i] = Source
    SumSource_OH_inc[i] = Source_OH_inc
    SumSource_Cl_inc[i] = Source_Cl_inc
    SumCH4change[i] = CH4change 
    SumSink[i] = CH4[i] * PT / Lifetime
    
# 1 box model to calculate change in sink assuming stable source
SumSink_StableSource = np.zeros(len(CH4) - 1)
Source = SumSource[0]
for i in range(len(CH4) - 1):
    Sink_StableSource = -(CH4[i + 1] * PT - CH4[i] * PT - Source)
    SumSink_StableSource[i] = Sink_StableSource
    
#%% Calculate sensitivity test change in sink
    
# Sink strengths in tg
Cl_Tg = 19
Cl_Tg_double = np.linspace(19, 38, 23)  # 23 elements to match SumSink length (len(CH4) - 1)
Strat_Tg = 39
Soil_Tg = 33
OH_Tg = SumSink - Cl_Tg - Strat_Tg - Soil_Tg
OH_Tg_doubleCl = SumSink - Cl_Tg_double - Strat_Tg - Soil_Tg
OH_Tg_StableSource = SumSink_StableSource - Cl_Tg - Strat_Tg - Soil_Tg


# Calculate KIEs
# C13 first
C13KIE_OnlyOH = (OH_KIE_13C*OH_Tg + Cl_KIE_13C*Cl_Tg + Strat_KIE_13C*Strat_Tg + Soil_KIE_13C*Soil_Tg) / SumSink
C13KIE_DoubleCl = (OH_KIE_13C*OH_Tg_doubleCl + Cl_KIE_13C*Cl_Tg_double + Strat_KIE_13C*Strat_Tg + Soil_KIE_13C*Soil_Tg) / SumSink
C13KIE_StableSource = (OH_KIE_13C*OH_Tg_StableSource + Cl_KIE_13C*Cl_Tg + Strat_KIE_13C*Strat_Tg + Soil_KIE_13C*Soil_Tg) / SumSink_StableSource





#%% Begin Monte Carlo Box model analysis

# Define constants and uncertainties
# Microbial
mic_d13C =  Mic_d13C_data.iloc[:, 1].mean() # Currently, no trend in microbial d13C-CH4. Could add this to test sensitivity
mic_d13C_U =  Mic_d13C_data.iloc[:, 2].mean() 
# Fossil
ff_d13C = FF_d13C_data.iloc[28:,1] #Starting in 1999
ff_d13C_U = FF_d13C_data.iloc[28:,2] 
# Biomass burning
bb_d13C = BB_d13C_data.iloc[1:,1]  # Skip first row to match bb_d13C_U
bb_d13C_U = BB_d13C_data.iloc[1:,2]
# Extract mean and standard deviation. The final row is the mean/std of all years and MC iterations (24000 annual estimates)
mean_bb = bb_d13C.iloc[-1]
mean_bb_U = bb_d13C_U.iloc[-1]
# Extend length of BB data to match full range of years
bb_d13C = np.concatenate((bb_d13C, (np.full(1, mean_bb))))
bb_d13C_U = np.concatenate((bb_d13C_U, (np.full(1, mean_bb_U))))
BB = np.mean(bbCT) # Biomass burning emissions. Does not change over time
BB_Thanwerdas = np.array([31] * 10 + [28] * 15) # Poor mans thanwerdas biomass burning
BBneg = np.concatenate((np.full(6, BB), np.linspace(BB,BB*.09, 18)))# Scenario with reduced biomass burning emissions starting in 2005

# Define constants and uncertainties
# Microbial
mic_dd = Mic_dD_data.iloc[:, 1].mean()
mic_dd_U = 7 # Mic_dD_data.iloc[:, 2].mean() 
# Fossil
ff_dd = FF_dD_data.iloc[34:,1] #Starting in 2005
ff_dd_U = FF_dD_data.iloc[34:,2] 
# Biomass burning
bb_dd = BB_dD_data.iloc[:,1]
bb_dd_U = BB_dD_data.iloc[:,2]
# Extract mean and standard deviation. The final row is the mean/std of all years and MC iterations (24000 annual estimates)
mean_bb = bb_dd.iloc[-1]
mean_bb_U = bb_dd_U.iloc[-1]
# Extend length of BB data to match full range of years
bb_dd = np.concatenate((np.full(3, mean_bb), bb_dd))
bb_dd = np.concatenate((bb_dd, (np.full(1, mean_bb))))
bb_dd_U = np.concatenate((np.full(3, mean_bb_U), bb_dd_U))
bb_dd_U = np.concatenate((bb_dd_U, (np.full(1, mean_bb_U))))
BB = np.mean(bbCT) # Biomass burning emissions. Does not change over time
#Create negative trend in BB for sensitivity test
BBneg = np.linspace(BB, BB*.09, 18)

# Data alignment and preparation for MC loop
print("\nPreparing data for 3-source Monte Carlo model...")
ff_d13C = np.array(ff_d13C).flatten()
ff_d13C_U = np.array(ff_d13C_U).flatten()
bb_d13C = np.array(bb_d13C).flatten()
bb_d13C_U = np.array(bb_d13C_U).flatten()
ff_dd = np.array(ff_dd).flatten()
ff_dd_U = np.array(ff_dd_U).flatten()
bb_dd = np.array(bb_dd).flatten()
bb_dd_U = np.array(bb_dd_U).flatten()

# Align dD data to match d13C timeframe (1999-2022)
# The MC loop uses SumSource which has length len(CH4)-1, so all source arrays should match that length
# Extend dD arrays with first values repeated to match d13C length (1999-2004 uses 2005 value)
target_length = len(SumSource)  # This is len(CH4) - 1

# Check and pad ff_dd to match target_length
while len(ff_dd) < target_length:
    ff_dd = np.concatenate([np.array([ff_dd[0]]), ff_dd])
    ff_dd_U = np.concatenate([np.array([ff_dd_U[0]]), ff_dd_U])
    bb_dd = np.concatenate([np.array([bb_dd[0]]), bb_dd])
    bb_dd_U = np.concatenate([np.array([bb_dd_U[0]]), bb_dd_U])

# Trim to exact target length if any are too long
ff_dd = ff_dd[:target_length]
ff_dd_U = ff_dd_U[:target_length]
bb_dd = bb_dd[:target_length]
bb_dd_U = bb_dd_U[:target_length]

# Pad dD_AnnAvg_matrix to match d13C length
pad_length_dD = len(ff_d13C) - dD_AnnAvg_matrix.shape[0]
if pad_length_dD > 0:
    first_row = dD_AnnAvg_matrix[0:1, :]  # Get first row as 2D array
    pad_rows = np.repeat(first_row, pad_length_dD, axis=0)
    dD_AnnAvg_matrix = np.vstack([pad_rows, dD_AnnAvg_matrix])

# Pad Mic_dD_MC to match MC loop length (len(SumSource) = 23)
# The original Mic_dD_MC has 18 rows and needs to be 23 for the loop
mc_target_length = len(SumSource)
if Mic_dD_MC.shape[0] < mc_target_length:
    pad_count = mc_target_length - Mic_dD_MC.shape[0]
    first_row_MC = Mic_dD_MC.iloc[0:1, :]
    pad_rows_MC = pd.concat([first_row_MC] * pad_count, ignore_index=True)
    Mic_dD_MC = pd.concat([pad_rows_MC, Mic_dD_MC], ignore_index=True)
elif Mic_dD_MC.shape[0] > mc_target_length:
    Mic_dD_MC = Mic_dD_MC.iloc[:mc_target_length, :]

# Verify data shapes
print(f"SumSource length: {len(SumSource)}")
print(f"ff_d13C length: {len(ff_d13C)}")
print(f"ff_dd length: {len(ff_dd)}")
print(f"Mic_d13C_MC shape: {Mic_d13C_MC.shape}")
print(f"Mic_dD_MC shape: {Mic_dD_MC.shape}")
print(f"FF_d13C_MC_EDGAR shape: {FF_d13C_MC_EDGAR.shape}")
print(f"CH4 shape: {CH4.shape}")
print(f"dD_AnnAvg_matrix shape: {dD_AnnAvg_matrix.shape}")
print(f"d13C_glob_iterations shape: {d13C_glob_iterations.shape}")

# Define years - use d13C years (1999-2022, 24 years) but take last 23 to match compiled data
# The compiled data has 23 rows from SumSource (CH4[4:28] gives 24, so 24-1=23)
year = FF_d13C_data.iloc[28:28+23,0].values  # 23 years from 1999 to 2021
# Define result matrices for 3-source model
d13C_Source_compiled = pd.DataFrame()
d13C_Source_RedCl_compiled = pd.DataFrame()
d13C_Source_IncOH_compiled = pd.DataFrame()
dD_Source_compiled = pd.DataFrame()
dD_Source_RedCl_compiled = pd.DataFrame()
dD_Source_IncOH_compiled = pd.DataFrame()

# 3-source Results: Shape (years, 1000)
BB_compiled = np.zeros((len(SumSource), 1000))  # Biomass Burning
FF_compiled = np.zeros((len(SumSource), 1000))  # Fossil Fuel
Mic_compiled = np.zeros((len(SumSource), 1000))  # Microbial

# Sensitivity tests for 3-source model
BB_RedCl_compiled = np.zeros((len(SumSource), 1000))
FF_RedCl_compiled = np.zeros((len(SumSource), 1000))
Mic_RedCl_compiled = np.zeros((len(SumSource), 1000))

BB_IncOH_compiled = np.zeros((len(SumSource), 1000))
FF_IncOH_compiled = np.zeros((len(SumSource), 1000))
Mic_IncOH_compiled = np.zeros((len(SumSource), 1000))

BB_BBdrop_compiled = np.zeros((len(SumSource), 1000))
FF_BBdrop_compiled = np.zeros((len(SumSource), 1000))
Mic_BBdrop_compiled = np.zeros((len(SumSource), 1000))


# ============================================================================
# Begin 3-Source, Dual-Isotope Monte Carlo Mass Balance Analysis
# ============================================================================
# Mathematical Framework:
# For each year t and MC iteration k, solve the system:
#
#   x1 + x2 + x3 = Total_Source
#   x1·δ13C_BB + x2·δ13C_FF + x3·δ13C_Mic = Total_Source·δ13C_source
#   x1·δD_BB + x2·δD_FF + x3·δD_Mic = Total_Source·δD_source
#
# where x1, x2, x3 are emissions from Biomass Burning, Fossil Fuel, Microbial
# ============================================================================

print("Starting 3-source Monte Carlo mass balance analysis...")
print(f"Data shapes: SumSource={len(SumSource)}, FF d13C MC={FF_d13C_MC_EDGAR.shape}, " +
      f"Mic d13C MC={Mic_d13C_MC.shape}, Mic dD MC={Mic_dD_MC.shape}")

MB_DEBUG = "--debug" in sys.argv
first_nonphysical = None  # (k, year, x, A, B)

for k in range(0, 1000):
    if (k + 1) % 100 == 0:
        print(f"  Processing iteration {k + 1}/1000...")
    
    # Sample atmospheric source signatures from MC iterations
    # These are calculated from the box model for each MC iteration
    d13C_atm_MC = (
        d13C_glob_iterations[:, k]
        if k < d13C_glob_iterations.shape[1]
        else d13C_glob_iterations[:, -1]
    )  # Shape: (len(CH4),)
    dD_atm_MC = (
        dD_AnnAvg_matrix[:, k]
        if k < dD_AnnAvg_matrix.shape[1]
        else dD_AnnAvg_matrix[:, -1]
    )  # Shape: (len(CH4),)
    
    # ----------- d13C SOURCE SIGNATURE CALCULATIONS ---------
    # Calculate atmospheric CO composition and derive source signatures
    f13_atm = delta_to_fraction_d13C(d13C_atm_MC)
    n13C = f13_atm * CH4 * PT
    
    # Calculate fractionation alphas for different sensitivity scenarios
    alpha_13C_base = 1 / Sink_13C  # Scalar
    alpha_13C_RedCl = 1 / Sink_13C_RedCl  # Array
    alpha_13C_IncOH = 1 / Sink_13C_OH_inc  # Array
    
    # Initialize arrays to store total-source heavy-isotope fractions (linear quantity)
    d13C_source_base = np.zeros(len(CH4) - 1)
    d13C_source_RedCl = np.zeros(len(CH4) - 1)
    d13C_source_IncOH = np.zeros(len(CH4) - 1)
    
    # Calculate source signatures using box model inversion
    for j in range(len(CH4) - 1):
        # Base case: standard Cl sink (alpha_13C_base is scalar)
        d13C_source_base[j] = (n13C[j + 1] - n13C[j] + n13C[j] * alpha_13C_base / Lifetime) / SumSource[j]
        
        # Reduced Cl scenario (alpha_13C_RedCl is array)
        d13C_source_RedCl[j] = (n13C[j + 1] - n13C[j] + n13C[j] * alpha_13C_RedCl[j] / Lifetime) / SumSource[j]
        
        # Increased OH scenario (alpha_13C_IncOH is array)
        d13C_source_IncOH[j] = (n13C[j + 1] - n13C[j] + n13C[j] * alpha_13C_IncOH[j] / Lifetime_OH_inc[j]) / SumSource_OH_inc[j]
    
    # ----------- dD SOURCE SIGNATURE CALCULATIONS ---------
    # Calculate atmospheric dD composition and derive source signatures
    fD_atm = delta_to_fraction_dD(dD_atm_MC)
    ndD = fD_atm * CH4 * PT
    
    # Calculate fractionation alphas for dD (simplified: use same sink fractionation)
    alpha_dD_base = 1 / Sink_D  # Scalar
    alpha_dD_RedCl = 1 / Sink_D_RedCl  # Array
    alpha_dD_IncOH = 1 / Sink_D_OH_inc  # Array
    
    # Initialize arrays for total-source heavy-isotope fractions
    dD_source_base = np.zeros(len(CH4) - 1)
    dD_source_RedCl = np.zeros(len(CH4) - 1)
    dD_source_IncOH = np.zeros(len(CH4) - 1)
    
    # Calculate source signatures using box model inversion
    for j in range(len(CH4) - 1):
        # Base case: standard Cl sink (alpha_dD_base is scalar)
        dD_source_base[j] = (ndD[j + 1] - ndD[j] + ndD[j] * alpha_dD_base / Lifetime) / SumSource[j]
        
        # Reduced Cl scenario (alpha_dD_RedCl is array)
        dD_source_RedCl[j] = (ndD[j + 1] - ndD[j] + ndD[j] * alpha_dD_RedCl[j] / Lifetime) / SumSource[j]
        
        # Increased OH scenario (alpha_dD_IncOH is array)
        dD_source_IncOH[j] = (ndD[j + 1] - ndD[j] + ndD[j] * alpha_dD_IncOH[j] / Lifetime_OH_inc[j]) / SumSource_OH_inc[j]
    
    # ---------- SOURCE SIGNATURE MONTE CARLO SAMPLING ----------
    # Sample source signatures with including isotopic measurement uncertainty
    # Draw new random samples for source uncertainties
    RandomGauss_FF_d13C = np.random.normal(0, 1)
    RandomGauss_BB_d13C = np.random.normal(0, 1)
    RandomGauss_FF_dD = np.random.normal(0, 1)
    RandomGauss_BB_dD = np.random.normal(0, 1)
    
    # Apply MC uncertainty to FF and BB source signatures
    # FF: use MC iteration k if available, otherwise use mean
    ff_d13C_MC_iter = np.array(FF_d13C_MC_EDGAR.iloc[:, k]) if k < FF_d13C_MC_EDGAR.shape[1] else ff_d13C + RandomGauss_FF_d13C * ff_d13C_U
    ff_dD_MC_iter = np.array(FF_dD_MC_EDGAR.iloc[:, k]) if k < FF_dD_MC_EDGAR.shape[1] else ff_dd + RandomGauss_FF_dD * ff_dd_U
    
    # BB: constant per year with MC sampling of uncertainty
    bb_d13C_MC_iter = bb_d13C + RandomGauss_BB_d13C * bb_d13C_U
    bb_dD_MC_iter = bb_dd + RandomGauss_BB_dD * bb_dd_U
    
    # Mic: use MC column k
    mic_d13C_MC_iter = np.array(Mic_d13C_MC.iloc[:, k])
    mic_dD_MC_iter = np.array(Mic_dD_MC.iloc[:, k])

    # Convert endmember deltas (‰) to heavy-isotope fractions for linear mass balance
    f13_bb = delta_to_fraction_d13C(bb_d13C_MC_iter)
    f13_ff = delta_to_fraction_d13C(ff_d13C_MC_iter)
    f13_mic = delta_to_fraction_d13C(mic_d13C_MC_iter)

    fD_bb = delta_to_fraction_dD(bb_dD_MC_iter)
    fD_ff = delta_to_fraction_dD(ff_dD_MC_iter)
    fD_mic = delta_to_fraction_dD(mic_dD_MC_iter)
    
    # =====================================================================
    # SCENARIO 1: BASE CASE (Standard Cl Sink)
    # =====================================================================
    for j in range(len(SumSource)):
        try:
            # Build 3×3 coefficient matrix A:
            # Row 0: Mass balance coefficients (all ones)
            # Row 1: δ13C coefficients (isotopic signatures from each source)
            # Row 2: δD coefficients (isotopic signatures from each source)
            A = np.array([
                [1.0, 1.0, 1.0],                                  # Mass: BB + FF + Mic = Total
                [bb_d13C_MC_iter[j], ff_d13C_MC_iter[j], mic_d13C_MC_iter[j]],  # δ13C balance
                [bb_dD_MC_iter[j], ff_dD_MC_iter[j], mic_dD_MC_iter[j]]         # δD balance
            ])
            
            # Unify isotope definition: use heavy-isotope fractions (linear for mass balance)
            A[1, :] = [f13_bb[j], f13_ff[j], f13_mic[j]]
            A[2, :] = [fD_bb[j], fD_ff[j], fD_mic[j]]

            # Build emissions vector B:
            # Element 0: Total source strength (Tg yr⁻¹)
            # Element 1: Total isotopic mass for d13C (Tg·‰ yr⁻¹)
            # Element 2: Total isotopic mass for dD (Tg·‰ yr⁻¹)
            A[1, :] = [f13_bb[j], f13_ff[j], f13_mic[j]]
            A[2, :] = [fD_bb[j], fD_ff[j], fD_mic[j]]

            B = np.array([
                SumSource[j],
                SumSource[j] * d13C_source_base[j],
                SumSource[j] * dD_source_base[j]
            ])
            
            # Solve Ax = B for x = [BB, FF, Mic] emissions
            x = np.linalg.solve(A, B)
            BB_compiled[j, k] = x[0]
            FF_compiled[j, k] = x[1]
            Mic_compiled[j, k] = x[2]

            if MB_DEBUG and first_nonphysical is None:
                if np.any(~np.isfinite(x)) or np.any(x < 0):
                    year_dbg = int(1999 + j)
                    first_nonphysical = (k, year_dbg, x.copy(), A.copy(), B.copy())
                    break
            
        except (np.linalg.LinAlgError, ValueError):
            # Matrix is singular or solution invalid -> store NaN
            BB_compiled[j, k] = np.nan
            FF_compiled[j, k] = np.nan
            Mic_compiled[j, k] = np.nan

    if MB_DEBUG and first_nonphysical is not None:
        break
    
    # =====================================================================
    # SCENARIO 2: REDUCED Cl SINK
    # =====================================================================
    for j in range(len(SumSource)):
        try:
            A = np.array([
                [1.0, 1.0, 1.0],
                [bb_d13C_MC_iter[j], ff_d13C_MC_iter[j], mic_d13C_MC_iter[j]],
                [bb_dD_MC_iter[j], ff_dD_MC_iter[j], mic_dD_MC_iter[j]]
            ])
            
            A[1, :] = [f13_bb[j], f13_ff[j], f13_mic[j]]
            A[2, :] = [fD_bb[j], fD_ff[j], fD_mic[j]]

            B = np.array([
                SumSource[j],
                SumSource[j] * d13C_source_RedCl[j],
                SumSource[j] * dD_source_RedCl[j]
            ])
            
            x = np.linalg.solve(A, B)
            BB_RedCl_compiled[j, k] = x[0]
            FF_RedCl_compiled[j, k] = x[1]
            Mic_RedCl_compiled[j, k] = x[2]
            
        except (np.linalg.LinAlgError, ValueError):
            BB_RedCl_compiled[j, k] = np.nan
            FF_RedCl_compiled[j, k] = np.nan
            Mic_RedCl_compiled[j, k] = np.nan
    
    # =====================================================================
    # SCENARIO 3: INCREASED OH SINK
    # =====================================================================
    for j in range(len(SumSource_OH_inc)):
        try:
            A = np.array([
                [1.0, 1.0, 1.0],
                [bb_d13C_MC_iter[j], ff_d13C_MC_iter[j], mic_d13C_MC_iter[j]],
                [bb_dD_MC_iter[j], ff_dD_MC_iter[j], mic_dD_MC_iter[j]]
            ])
            
            A[1, :] = [f13_bb[j], f13_ff[j], f13_mic[j]]
            A[2, :] = [fD_bb[j], fD_ff[j], fD_mic[j]]

            B = np.array([
                SumSource_OH_inc[j],
                SumSource_OH_inc[j] * d13C_source_IncOH[j],
                SumSource_OH_inc[j] * dD_source_IncOH[j]
            ])
            
            x = np.linalg.solve(A, B)
            BB_IncOH_compiled[j, k] = x[0]
            FF_IncOH_compiled[j, k] = x[1]
            Mic_IncOH_compiled[j, k] = x[2]
            
        except (np.linalg.LinAlgError, ValueError):
            BB_IncOH_compiled[j, k] = np.nan
            FF_IncOH_compiled[j, k] = np.nan
            Mic_IncOH_compiled[j, k] = np.nan
    
    # =====================================================================
    # SCENARIO 4: REDUCED BIOMASS BURNING
    # =====================================================================
    for j in range(len(SumSource)):
        try:
            # Same sink as base case but with different forcing/context
            A = np.array([
                [1.0, 1.0, 1.0],
                [bb_d13C_MC_iter[j], ff_d13C_MC_iter[j], mic_d13C_MC_iter[j]],
                [bb_dD_MC_iter[j], ff_dD_MC_iter[j], mic_dD_MC_iter[j]]
            ])
            
            B = np.array([
                SumSource[j],
                SumSource[j] * d13C_source_base[j],
                SumSource[j] * dD_source_base[j]
            ])
            
            x = np.linalg.solve(A, B)
            BB_BBdrop_compiled[j, k] = x[0]
            FF_BBdrop_compiled[j, k] = x[1]
            Mic_BBdrop_compiled[j, k] = x[2]
            
        except (np.linalg.LinAlgError, ValueError):
            BB_BBdrop_compiled[j, k] = np.nan
            FF_BBdrop_compiled[j, k] = np.nan
            Mic_BBdrop_compiled[j, k] = np.nan

print("Monte Carlo analysis complete!")

if MB_DEBUG and first_nonphysical is not None:
    k_dbg, year_dbg, x_dbg, A_dbg, B_dbg = first_nonphysical
    print("\n[MB_DEBUG] First non-physical solution detected")
    print(f"[MB_DEBUG] iteration k={k_dbg}, year={year_dbg}")
    print(f"[MB_DEBUG] solution x=[BB, FF, Mic]={x_dbg}")
    print(f"[MB_DEBUG] A=\n{A_dbg}")
    print(f"[MB_DEBUG] B={B_dbg}")
    f13_src = float(B_dbg[1] / B_dbg[0])
    fD_src = float(B_dbg[2] / B_dbg[0])
    f13_min, f13_max = float(np.min(A_dbg[1, :])), float(np.max(A_dbg[1, :]))
    fD_min, fD_max = float(np.min(A_dbg[2, :])), float(np.max(A_dbg[2, :]))
    print(f"[MB_DEBUG] implied f13_source={f13_src:.10f} (endmembers range {f13_min:.10f}..{f13_max:.10f})")
    print(f"[MB_DEBUG] implied fD_source={fD_src:.10f} (endmembers range {fD_min:.10f}..{fD_max:.10f})")
    sys.exit(0)

    
#%% Post-processing: Compute statistics and smooth results

# Calculate statistics for 3-source model (mean and std across 1000 iterations)
print("Computing statistics from Monte Carlo results...")

# Base case (standard sink)
BB_mean = np.mean(BB_compiled, axis=1)
BB_std = np.std(BB_compiled, axis=1)
FF_mean = np.mean(FF_compiled, axis=1)
FF_std = np.std(FF_compiled, axis=1)
Mic_mean = np.mean(Mic_compiled, axis=1)
Mic_std = np.std(Mic_compiled, axis=1)

# Red Cl scenario
BB_RedCl_mean = np.mean(BB_RedCl_compiled, axis=1)
BB_RedCl_std = np.std(BB_RedCl_compiled, axis=1)
FF_RedCl_mean = np.mean(FF_RedCl_compiled, axis=1)
FF_RedCl_std = np.std(FF_RedCl_compiled, axis=1)
Mic_RedCl_mean = np.mean(Mic_RedCl_compiled, axis=1)
Mic_RedCl_std = np.std(Mic_RedCl_compiled, axis=1)

# IncOH scenario
BB_IncOH_mean = np.mean(BB_IncOH_compiled, axis=1)
BB_IncOH_std = np.std(BB_IncOH_compiled, axis=1)
FF_IncOH_mean = np.mean(FF_IncOH_compiled, axis=1)
FF_IncOH_std = np.std(FF_IncOH_compiled, axis=1)
Mic_IncOH_mean = np.mean(Mic_IncOH_compiled, axis=1)
Mic_IncOH_std = np.std(Mic_IncOH_compiled, axis=1)

# BBdrop scenario
BB_BBdrop_mean = np.mean(BB_BBdrop_compiled, axis=1)
BB_BBdrop_std = np.std(BB_BBdrop_compiled, axis=1)
FF_BBdrop_mean = np.mean(FF_BBdrop_compiled, axis=1)
FF_BBdrop_std = np.std(FF_BBdrop_compiled, axis=1)
Mic_BBdrop_mean = np.mean(Mic_BBdrop_compiled, axis=1)
Mic_BBdrop_std = np.std(Mic_BBdrop_compiled, axis=1)

# Optional: Apply 5-year smoothing to results
def smooth_5year_1D(data):
    """Apply 5-year smoothing to 1D array"""
    smoothed = np.zeros(len(data))
    for i in range(len(data)):
        if i < 2:
            # Use fewer points at beginning
            smoothed[i] = np.nanmean(data[:i+3])
        elif i > len(data) - 3:
            # Use fewer points at end
            smoothed[i] = np.nanmean(data[i-2:])
        else:
            # Full 5-year window
            smoothed[i] = np.nanmean(data[i-2:i+3])
    return smoothed

# Smooth the mean values (optional, apply if desired)
# BB_mean_smooth = smooth_5year_1D(BB_mean)
# FF_mean_smooth = smooth_5year_1D(FF_mean)
# Mic_mean_smooth = smooth_5year_1D(Mic_mean)

print("Statistics computed successfully!")

#%% Save results to output files

print("Saving MC results to output files...")

# Get year vector for output
years_output = np.arange(1999, 1999 + len(SumSource))

# Function to save results
def save_source_results(filename, BB_array, FF_array, Mic_array, years_vec, scenario_name="Base"):
    """Save 3-source results to CSV file (all iterations)"""
    df_list = []
    for i, year in enumerate(years_vec):
        df_list.append({
            'Year': year,
            'Source': 'BB_mean', 'Value': BB_array[i],
        })
    # This is a simplified version - you can expand for full MC iterations if needed

def save_summary_statistics(filename_prefix, years_vec, scenario_name="Base"):
    """Save summary statistics for all scenarios"""
    results = pd.DataFrame({
        'Year': years_vec,
        'BB_mean': BB_mean, 'BB_std': BB_std,
        'FF_mean': FF_mean, 'FF_std': FF_std,
        'Mic_mean': Mic_mean, 'Mic_std': Mic_std,
    })
    results.to_csv(OUT_DIR / f'{filename_prefix}_base_{scenario_name}.csv', index=False)
    
    results_redcl = pd.DataFrame({
        'Year': years_vec,
        'BB_mean': BB_RedCl_mean, 'BB_std': BB_RedCl_std,
        'FF_mean': FF_RedCl_mean, 'FF_std': FF_RedCl_std,
        'Mic_mean': Mic_RedCl_mean, 'Mic_std': Mic_RedCl_std,
    })
    results_redcl.to_csv(OUT_DIR / f'{filename_prefix}_RedCl_{scenario_name}.csv', index=False)
    
    results_incoh = pd.DataFrame({
        'Year': years_vec,
        'BB_mean': BB_IncOH_mean, 'BB_std': BB_IncOH_std,
        'FF_mean': FF_IncOH_mean, 'FF_std': FF_IncOH_std,
        'Mic_mean': Mic_IncOH_mean, 'Mic_std': Mic_IncOH_std,
    })
    results_incoh.to_csv(OUT_DIR / f'{filename_prefix}_IncOH_{scenario_name}.csv', index=False)
    
    results_bbdrop = pd.DataFrame({
        'Year': years_vec,
        'BB_mean': BB_BBdrop_mean, 'BB_std': BB_BBdrop_std,
        'FF_mean': FF_BBdrop_mean, 'FF_std': FF_BBdrop_std,
        'Mic_mean': Mic_BBdrop_mean, 'Mic_std': Mic_BBdrop_std,
    })
    results_bbdrop.to_csv(OUT_DIR / f'{filename_prefix}_BBdrop_{scenario_name}.csv', index=False)

# Save summary statistics
save_summary_statistics('ThreeSource_MassBalance_MC', years_output, 'DualIsotope')

# Save detailed MC results (all 1000 iterations) - optional
# Convert to DataFrames for easy CSV export
BB_compiled_df = pd.DataFrame(BB_compiled, index=years_output, columns=[f'Iter_{i}' for i in range(1000)])
FF_compiled_df = pd.DataFrame(FF_compiled, index=years_output, columns=[f'Iter_{i}' for i in range(1000)])
Mic_compiled_df = pd.DataFrame(Mic_compiled, index=years_output, columns=[f'Iter_{i}' for i in range(1000)])

BB_compiled_df.to_csv(OUT_DIR / 'BB_3source_MC_alliterations.csv')
FF_compiled_df.to_csv(OUT_DIR / 'FF_3source_MC_alliterations.csv')
Mic_compiled_df.to_csv(OUT_DIR / 'Mic_3source_MC_alliterations.csv')

print("Results saved to Output/ directory")
print(f"  - ThreeSource_MassBalance_MC_base_DualIsotope.csv")
print(f"  - ThreeSource_MassBalance_MC_RedCl_DualIsotope.csv")
print(f"  - ThreeSource_MassBalance_MC_IncOH_DualIsotope.csv")
print(f"  - ThreeSource_MassBalance_MC_BBdrop_DualIsotope.csv")
print(f"  - BB_3source_MC_alliterations.csv")
print(f"  - FF_3source_MC_alliterations.csv")
print(f"  - Mic_3source_MC_alliterations.csv")
#%% Visualization: Create spaghetti plots for 3-source model

print("Creating visualization plots...")

# Create figure with subplots for all three sources
plt.rc('font', size=12) 
fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=300, sharex=True)

# Row 1: Biomass Burning
axes[0, 0].plot(years_output, BB_compiled, linewidth=0.5, alpha=0.3, color='red')
axes[0, 0].plot(years_output, BB_mean, 'r-', linewidth=2.5, label='Mean')
axes[0, 0].fill_between(years_output, BB_mean - BB_std, BB_mean + BB_std, alpha=0.2, color='red')
axes[0, 0].set_ylabel('BB Emissions (Tg yr$^{-1}$)')
axes[0, 0].set_title('Biomass Burning - All MC Iterations')
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].legend()

axes[0, 1].plot(years_output, BB_mean, 'r-', linewidth=2.5, label='Base', marker='o')
axes[0, 1].plot(years_output, BB_RedCl_mean, 'b--', linewidth=2, label='RedCl', marker='s')
axes[0, 1].plot(years_output, BB_IncOH_mean, 'g--', linewidth=2, label='IncOH', marker='^')
axes[0, 1].set_ylabel('BB Emissions (Tg yr$^{-1}$)')
axes[0, 1].set_title('Biomass Burning - Scenario Comparison')
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].legend()

# Row 2: Fossil Fuel
axes[1, 0].plot(years_output, FF_compiled, linewidth=0.5, alpha=0.3, color='blue')
axes[1, 0].plot(years_output, FF_mean, 'b-', linewidth=2.5, label='Mean')
axes[1, 0].fill_between(years_output, FF_mean - FF_std, FF_mean + FF_std, alpha=0.2, color='blue')
axes[1, 0].set_ylabel('FF Emissions (Tg yr$^{-1}$)')
axes[1, 0].set_title('Fossil Fuel - All MC Iterations')
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].legend()

axes[1, 1].plot(years_output, FF_mean, 'b-', linewidth=2.5, label='Base', marker='o')
axes[1, 1].plot(years_output, FF_RedCl_mean, 'r--', linewidth=2, label='RedCl', marker='s')
axes[1, 1].plot(years_output, FF_IncOH_mean, 'g--', linewidth=2, label='IncOH', marker='^')
axes[1, 1].set_ylabel('FF Emissions (Tg yr$^{-1}$)')
axes[1, 1].set_title('Fossil Fuel - Scenario Comparison')
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].legend()

# Row 3: Microbial
axes[2, 0].plot(years_output, Mic_compiled, linewidth=0.5, alpha=0.3, color='green')
axes[2, 0].plot(years_output, Mic_mean, 'g-', linewidth=2.5, label='Mean')
axes[2, 0].fill_between(years_output, Mic_mean - Mic_std, Mic_mean + Mic_std, alpha=0.2, color='green')
axes[2, 0].set_ylabel('Mic Emissions (Tg yr$^{-1}$)')
axes[2, 0].set_xlabel('Year')
axes[2, 0].set_title('Microbial - All MC Iterations')
axes[2, 0].grid(True, alpha=0.3)
axes[2, 0].legend()

axes[2, 1].plot(years_output, Mic_mean, 'g-', linewidth=2.5, label='Base', marker='o')
axes[2, 1].plot(years_output, Mic_RedCl_mean, 'r--', linewidth=2, label='RedCl', marker='s')
axes[2, 1].plot(years_output, Mic_IncOH_mean, 'b--', linewidth=2, label='IncOH', marker='^')
axes[2, 1].set_ylabel('Mic Emissions (Tg yr$^{-1}$)')
axes[2, 1].set_xlabel('Year')
axes[2, 1].set_title('Microbial - Scenario Comparison')
axes[2, 1].grid(True, alpha=0.3)
axes[2, 1].legend()

plt.tight_layout()
plt.savefig(OUT_DIR / 'ThreeSource_MassBalance_MC_allScenarios.png', dpi=300, bbox_inches='tight')
#plt.show() here we do not show the plot to avoid blocking execution in some environments

print(f"Plot saved: {OUT_DIR / 'ThreeSource_MassBalance_MC_allScenarios.png'}")
print("\n" + "="*70)
print("3-SOURCE MONTE CARLO MASS BALANCE MODEL ANALYSIS COMPLETE")
print("="*70)
print(f"Processed {len(years_output)} years with 1000 Monte Carlo iterations")
print(f"Sources: Biomass Burning (BB), Fossil Fuel (FF), Microbial (Mic)")
print(f"Isotopes: δ¹³C and δD")
print(f"Scenarios: Base, RedCl (Reduced Cl Sink), IncOH (Increased OH), BBdrop (Reduced BB)")
print("="*70)
# # 3 year smooth: 
# # Smooth FF_compiled before calculating stats (3-year smoothing)
# moving_average = np.zeros((FF_compiled.shape[0] - 2, FF_compiled.shape[1]))
# for col in range(FF_compiled.shape[1]):
#     for k in range(FF_compiled.shape[0] - 2):
#         moving_average[k, col] = (
#             FF_compiled.iloc[k, col] +
#             FF_compiled.iloc[k+1, col] +
#             FF_compiled.iloc[k+2, col]
#         ) / 3
# # first and last points use shorter averages
# mov_start = ((FF_compiled.iloc[0, :] + FF_compiled.iloc[1, :]) / 2).to_numpy().reshape(1, -1)
# mov_end   = ((FF_compiled.iloc[-2, :] + FF_compiled.iloc[-1, :]) / 2).to_numpy().reshape(1, -1)

# FF_compiled_smoothed  = np.vstack([mov_start, moving_average, mov_end])
# FF_compiled_smoothedRtest = FF_compiled_smoothed - np.mean(FF_compiled_smoothed[6:9, :], axis=0)

# # Smooth FF_compiled_R before calculating stats (3-year smoothing)
# moving_averageR = np.zeros((FF_compiledR.shape[0] - 2, FF_compiledR.shape[1]))
# for col in range(FF_compiledR.shape[1]):
#     for k in range(FF_compiledR.shape[0] - 2):
#         moving_averageR[k, col] = (
#             FF_compiledR.iloc[k, col] +
#             FF_compiledR.iloc[k+1, col] +
#             FF_compiledR.iloc[k+2, col]
#         ) / 3
# # first and last points use shorter averages
# mov_startR = ((FF_compiledR.iloc[0, :] + FF_compiledR.iloc[1, :]) / 2).to_numpy().reshape(1, -1)
# mov_endR   = ((FF_compiledR.iloc[-2, :] + FF_compiledR.iloc[-1, :]) / 2).to_numpy().reshape(1, -1)

# FF_compiledR_smoothed  = np.vstack([mov_startR, moving_averageR, mov_endR])

# # Smooth Mic_compiled (3-year smoothing)
# moving_average = np.zeros((Mic_compiled.shape[0] - 2, Mic_compiled.shape[1]))
# for col in range(Mic_compiled.shape[1]):
#     for k in range(Mic_compiled.shape[0] - 2):
#         moving_average[k, col] = (
#             Mic_compiled.iloc[k, col] +
#             Mic_compiled.iloc[k+1, col] +
#             Mic_compiled.iloc[k+2, col]
#         ) / 3
# mov_start = ((Mic_compiled.iloc[0, :] + Mic_compiled.iloc[1, :]) / 2).to_numpy().reshape(1, -1)
# mov_end   = ((Mic_compiled.iloc[-2, :] + Mic_compiled.iloc[-1, :]) / 2).to_numpy().reshape(1, -1)

# Mic_compiled_smoothed  = np.vstack([mov_start, moving_average, mov_end])
# Mic_compiled_smoothedRtest = Mic_compiled_smoothed - np.mean(Mic_compiled_smoothed[6:9, :], axis=0)

# # Smooth Mic_compiled (3-year smoothing)
# moving_averageR = np.zeros((Mic_compiledR.shape[0] - 2, Mic_compiledR.shape[1]))
# for col in range(Mic_compiledR.shape[1]):
#     for k in range(Mic_compiledR.shape[0] - 2):
#         moving_averageR[k, col] = (
#             Mic_compiledR.iloc[k, col] +
#             Mic_compiledR.iloc[k+1, col] +
#             Mic_compiledR.iloc[k+2, col]
#         ) / 3
# mov_startR = ((Mic_compiledR.iloc[0, :] + Mic_compiledR.iloc[1, :]) / 2).to_numpy().reshape(1, -1)
# mov_endR   = ((Mic_compiledR.iloc[-2, :] + Mic_compiledR.iloc[-1, :]) / 2).to_numpy().reshape(1, -1)

# Mic_compiledR_smoothed  = np.vstack([mov_startR, moving_averageR, mov_endR])


# 5 year smooth
# Smooth FF_compiled before calculating stats
moving_average = np.zeros((FF_compiled.shape[0] - 4, FF_compiled.shape[1]))
for col in range(FF_compiled.shape[1]):
    for k in range(FF_compiled.shape[0] - 4):
        moving_average[k, col] = (
            FF_compiled[k, col] +
            FF_compiled[k+1, col] +
            FF_compiled[k+2, col] +
            FF_compiled[k+3, col] +
            FF_compiled[k+4, col]
        ) / 5
# first two and last two points use shorter averages
mov_start = np.vstack([
    (FF_compiled[0, :] + FF_compiled[1, :] + FF_compiled[2, :]) / 3,
    (FF_compiled[0, :] + FF_compiled[1, :] + FF_compiled[2, :] + FF_compiled[3, :]) / 4
])
mov_end = np.vstack([
    (FF_compiled[-4, :] + FF_compiled[-3, :] + FF_compiled[-2, :] + FF_compiled[-1, :]) / 4,
    (FF_compiled[-3, :] + FF_compiled[-2, :] + FF_compiled[-1, :]) / 3
])
FF_compiled_smoothed = np.vstack([mov_start, moving_average, mov_end])
FF_compiledR_smoothed = FF_compiled_smoothed - np.mean(FF_compiled_smoothed[6:9,:],axis=0) # Relative to 2005 to 2007 average
FF_compiledR_smoothed_2013 = FF_compiled_smoothed - FF_compiled_smoothed[14,:] # relative to 2013

# Smooth Mic_compiled
moving_average = np.zeros((Mic_compiled.shape[0] - 4, Mic_compiled.shape[1]))
for col in range(Mic_compiled.shape[1]):
    for k in range(Mic_compiled.shape[0] - 4):
        moving_average[k, col] = (
            Mic_compiled[k, col] +
            Mic_compiled[k+1, col] +
            Mic_compiled[k+2, col] +
            Mic_compiled[k+3, col] +
            Mic_compiled[k+4, col]
        ) / 5
mov_start = np.vstack([
    (Mic_compiled[0, :] + Mic_compiled[1, :] + Mic_compiled[2, :]) / 3,
    (Mic_compiled[0, :] + Mic_compiled[1, :] + Mic_compiled[2, :] + Mic_compiled[3, :]) / 4
])
mov_end = np.vstack([
    (Mic_compiled[-4, :] + Mic_compiled[-3, :] + Mic_compiled[-2, :] + Mic_compiled[-1, :]) / 4,
    (Mic_compiled[-3, :] + Mic_compiled[-2, :] + Mic_compiled[-1, :]) / 3
])
Mic_compiled_smoothed = np.vstack([mov_start, moving_average, mov_end])
Mic_compiledR_smoothed = Mic_compiled_smoothed - np.mean(Mic_compiled_smoothed[6:9,:],axis=0) # Relative to 2005 to 2007 average
Mic_compiledR_smoothed_2013 = Mic_compiled_smoothed - Mic_compiled_smoothed[14,:] # relative to 2013


# Calculate statistics
FF_mean = FF_compiled_smoothed.mean(axis=1)
FF_std = FF_compiled_smoothed.std(axis=1)
Mic_mean = Mic_compiled_smoothed.mean(axis=1)
Mic_std = Mic_compiled_smoothed.std(axis=1)
d13C_source_meanMC = d13C_Source_compiled.mean(axis=1).to_numpy()
Mic_RedCl_mean = Mic_RedCl_compiled.mean(axis=1)
FF_RedCl_mean = FF_RedCl_compiled.mean(axis=1)
Mic_IncOH_mean = Mic_IncOH_compiled.mean(axis=1)
FF_IncOH_mean = FF_IncOH_compiled.mean(axis=1)
Mic_BBdrop_mean = Mic_BBdrop_compiled.mean(axis=1)
FF_BBdrop_mean = FF_BBdrop_compiled.mean(axis=1)
# Calculate relative 
FF_meanR = FF_compiledR_smoothed.mean(axis=1)
FF_stdR = FF_compiledR_smoothed.std(axis=1)
Mic_meanR = Mic_compiledR_smoothed.mean(axis=1)
Mic_stdR = Mic_compiledR_smoothed.std(axis=1)
# Calculate relative to 2013 
FF_meanR_2013 = FF_compiledR_smoothed_2013.mean(axis=1)
FF_stdR_2013 = FF_compiledR_smoothed_2013.std(axis=1)
Mic_meanR_2013 = Mic_compiledR_smoothed_2013.mean(axis=1)
Mic_stdR_2013 = Mic_compiledR_smoothed_2013.std(axis=1)
# Combine into results matrix
results = np.stack((year, FF_mean, Mic_mean), axis=-1)

# plot spaghetti diagram for d13C
plt.rc('font', size=14) 
fig, axes = plt.subplots(2, 2, figsize=(10, 10), dpi=300, sharex=True, gridspec_kw={'hspace': 0.35, 'wspace': 0.3})

# top row title
fig.text(0.5, 0.90, r'${\delta}^{13}C$-derived Fossil Fuel Emissions', ha='center', fontsize=16)
# top row
axes[0, 0].plot(year, FF_compiled_smoothed - FF_compiled_smoothed[0, :])
axes[0, 0].set_ylabel('Relative FF Emissions Tg yr$^{-1}$')
axes[0, 0].set_xlim(1998, 2022)
axes[0, 0].set_ylim(-120, 120)
axes[0, 0].grid(False)
axes[0, 1].plot(year, FF_compiled_smoothed)
axes[0, 1].set_ylabel('FF Emissions Tg yr$^{-1}$')
axes[0, 1].set_xlim(1998, 2022)
axes[0, 1].set_ylim(10, 250)
axes[0, 1].grid(False)

# bottom row title
fig.text(0.5, 0.48, r'${\delta}^{13}C$-derived Microbial Emissions', ha='center', fontsize=16)
# bottom row
axes[1, 0].plot(year, Mic_compiled_smoothed - Mic_compiled_smoothed[0, :])
axes[1, 0].set_ylabel('Relative Mic Emissions Tg yr$^{-1}$')
axes[1, 0].set_xlim(1998, 2022)
axes[1, 0].set_ylim(-90, 200)
axes[1, 0].set_xlabel('Year')
axes[1, 0].grid(False)
axes[1, 1].plot(year, Mic_compiled_smoothed)
axes[1, 1].set_ylabel('Mic Emissions Tg yr$^{-1}$')
axes[1, 1].set_xlim(1998, 2022)
axes[1, 1].set_ylim(290, 580)
axes[1, 1].set_xlabel('Year')
axes[1, 1].grid(False)
plt.tight_layout(rect=[0, 0.05, 1, 0.90])
if sys.stdin is not None and sys.stdin.isatty():
    plt.show()
else:
    plt.close(fig)


# # Save spaghetti emissions 
# df_spaghetti = pd.DataFrame(FF_compiled_smoothed)
# # save to Excel
# df_spaghetti.to_excel('Output/FF_compiled_smoothed_d13C.xlsx', index=False)


#%% Plot histogram of difference (2020 to 2022 - 2005 to 2007)

# Calculate Deltas (2020 to 2022 - 2005 to 2007)
Mic_Delta_MC = Mic_compiled[-3:, :].mean(axis=0) - Mic_compiled[6:8, :].mean(axis=0)
FF_Delta_MC = FF_compiled[-3:, :].mean(axis=0) - FF_compiled[6:8, :].mean(axis=0)
Mic_Delta_MC_RedCl = Mic_RedCl_compiled[-3:, :].mean(axis=0) - Mic_RedCl_compiled[6:8, :].mean(axis=0)
FF_Delta_MC_RedCl = FF_RedCl_compiled[-3:, :].mean(axis=0) - FF_RedCl_compiled[6:8, :].mean(axis=0)
Mic_Delta_MC_IncOH = Mic_IncOH_compiled[-3:, :].mean(axis=0) - Mic_IncOH_compiled[6:8, :].mean(axis=0)
FF_Delta_MC_IncOH = FF_IncOH_compiled[-3:, :].mean(axis=0) - FF_IncOH_compiled[6:8, :].mean(axis=0)
Mic_Delta_MC_BBdrop = Mic_BBdrop_compiled[-3:, :].mean(axis=0) - Mic_BBdrop_compiled[6:8, :].mean(axis=0)
FF_Delta_MC_BBdrop = FF_BBdrop_compiled[-3:, :].mean(axis=0) - FF_BBdrop_compiled[6:8, :].mean(axis=0)
# Calculate Deltas  (2020 - 2000 to 2009)
Mic_Delta_MC_2020 = Mic_compiled[-3, :] - Mic_compiled[1:11, :].mean(axis=0)
FF_Delta_MC_2020 = FF_compiled[-3, :] - FF_compiled[1:11, :].mean(axis=0)


# Save the data
# Prepare the data for the first set of plots (absolute emissions)
df_absolute = pd.DataFrame({'Year': year,'FF_mean': FF_mean,'FF_std': FF_std,'Mic_mean': Mic_mean,'Mic_std': Mic_std,})
# Save the first DataFrame to a CSV file
df_absolute.to_csv(OUT_DIR / 'Results_d13C-MassBalance_Cantrell.csv', index=False)
# Prepare the data for the second set of plots (relative emissions)
df_relative = pd.DataFrame({'Year': year,'FF_meanR': FF_meanR,'FF_stdR': FF_stdR,'Mic_meanR': Mic_meanR,'Mic_stdR': Mic_stdR,})
# Save the second DataFrame to a CSV file
df_relative.to_csv(OUT_DIR / 'Results_Rd13C-MassBalance_Cantrell.csv', index=False)

# Save histogram data
matrix = np.column_stack((Mic_Delta_MC, FF_Delta_MC, Mic_Delta_MC_2020, FF_Delta_MC_2020))
np.savetxt(OUT_DIR / "d13C_histogram_Cantrell.csv", matrix, delimiter=",")


# # Prepare sensitivity test results and export 
# df_sensitivity = pd.DataFrame({'Year': year,'FF Reduced BB': FF_BBdrop_mean_mov,'Mic Reduced BB': Mic_BBdrop_mean_mov,'FF Reduced Cl': FF_RedCl_mean_mov,'Mic Reduced Cl': Mic_RedCl_mean_mov,})
# df_sensitivity.to_csv('Output/Results_d13C-BBClSensitivity.csv', index=False)  
# # Prepare FF source sig test results and export 
# df_FFsrcsig = pd.DataFrame({'Year': year,'FF, Stable Mic': FFS_StableMic_mean_mov,'FF, Half Mic': FFS_Half_mean_mov,'d13C FF, Stable Mic': FF_SrcSig_mean_mov,'d13C FF, Half Mic': FF_SrcSig_Half_mean_mov,})
# df_FFsrcsig.to_csv('Output/Results_d13C-FF_Sensitivity.csv', index=False) 
if sys.stdin is not None and sys.stdin.isatty():
    input("Press Enter to continue...")
