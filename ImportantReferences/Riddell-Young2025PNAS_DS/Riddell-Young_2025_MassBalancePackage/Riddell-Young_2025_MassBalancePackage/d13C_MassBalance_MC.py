#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 15:07:05 2024

@author: ryoung
"""

# Code for Monte Carlo analysis of onebox delta13C-CH4 mass balance
locals().clear()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

# Load INSTAAR DEI d13C-CH4
C13data = pd.read_excel('data/ch4c13_nh_sh_mean.xlsx').to_numpy()
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
CH4data = pd.read_excel('data/GML_CH4_AnnualMean.xlsx').to_numpy()
CH4data = CH4data[11:,:]
CH4 = CH4data[4:,1] # Starting in 1999
CH4year = CH4data[4:,0]

# Load global annual mean DEI iterations
d13C_glob_iterations_data = np.loadtxt("data/d13C_dei_compiled.txt")
d13C_glob_iterations = d13C_glob_iterations_data[:,1:]

# Load Carbon tracker methane
data2 = pd.read_excel('data/CarbonTracker_CH4.xlsx')
yearss = data2.iloc[:,0].values
micCT = data2.iloc[:,7].values
ffCT = data2.iloc[:,3].values
bbCT = data2.iloc[:,9].values #Prior emissions (GFED4)

# Load annual mean Mic and BB source signatures
BB_d13C_data = pd.read_csv('Output/BB_d13C_annual.csv', delimiter=',', header = None) # Calculated from Luo C3 C4 map and CTCH4 BB emissions
Mic_d13C_data = pd.read_csv('Output/Mic_d13C_annual.csv', delimiter=',', header = None)
Mic_d13C_MC_trends = pd.read_csv('Output/Mic_d13C_MC.csv', delimiter=',', header = None) 
Mic_d13C_MC = Mic_d13C_MC_trends.iloc[:,1:]

# Load FF source sigantures
FF_d13C_data = pd.read_csv('Output/FF_d13C_GlobUnc.csv', delimiter=',')
FF_d13C_MC_CTCH4_data = pd.read_csv('Output/FF_d13C_GlobMC_CTCH4.csv', delimiter=',', header = None)
FF_d13C_MC_EDGAR_data = pd.read_csv('Output/FF_d13C_GlobMC_EDGAR.csv', delimiter=',')
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

# Calculate sink KIE
OH_KIE = 1.0054    # Saueressig et al,. 2001 is 1.0039. Best estimate in Ben Li's opinion (Cantrell is 1.0054)
Cl_KIE = 1.066     # Also Saueressig
Strat_KIE = 1.003  # Also Saueressig, Lassey et al., 2007
Soil_KIE = 1.0201  # Average of snover and quay, tyler, and reeburgh
# Sink strengths: carbontracker documentation (Thanwerdas)
OH_Sink = .835 #.859
OH_Sink_Than = .899 #.899 .859
Cl_Sink = .035 #.008
Cl_Sink_Than = .006 #.006 #.008
Strat_Sink = .07 #.072
Strat_Sink_Than = .03 #.03 #.072
Soil_Sink = .06 #.062
Soil_Sink_Than = .065 #.065 #.062
# Net KIE
Sink = OH_KIE*OH_Sink + Cl_KIE*Cl_Sink + Strat_KIE*Strat_Sink + Soil_KIE*Soil_Sink  #Thanwerdas 2022 uses this number 

# Calculate change in KIE if Cl sink decreases
Cl_Sink_red = np.concatenate((np.full(6, Cl_Sink), np.linspace(Cl_Sink, .011, 18)))
Sink_Rem = Cl_Sink - Cl_Sink_red
OH_Sink_RedCl = OH_Sink + Sink_Rem * OH_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Strat_Sink_RedCl = Strat_Sink + Sink_Rem * Strat_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Soil_Sink_RedCl = Soil_Sink + Sink_Rem * Soil_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Check = OH_Sink_RedCl + Strat_Sink_RedCl + Soil_Sink_RedCl + Cl_Sink_red
Sink_RedCl = OH_KIE*OH_Sink_RedCl + Cl_KIE*Cl_Sink_red + Strat_KIE*Strat_Sink_RedCl + Soil_KIE*Soil_Sink_RedCl

# Calculate change in KIE if Cl sink increases (Van Herpen et al., 2023)
Cl_Sink_inc_frac = np.array([1, 1, 1, 1, 1, 1, 1, 0.97, 0.94, 0.58, 0.70, 1.13, 0.79, 1.09, 0.80, 1.48, 1.86, 2.08, 1.62, 1.85, 1.85, 1.85, 1.85, 1.85])
Cl_Sink_inc = Cl_Sink_inc_frac*Cl_Sink
Sink_Cl_inc = (OH_KIE*OH_Sink + Cl_KIE*Cl_Sink_inc + Strat_KIE*Strat_Sink + Soil_KIE*Soil_Sink) / (OH_Sink + Cl_Sink_inc + Soil_Sink + Strat_Sink)

# Create increasing OH following Olaf et al., trends
OH_ann_change = OH_Sink*.003
OH_Sink_inc = np.concatenate((np.full(6, OH_Sink), np.linspace(OH_Sink, OH_Sink + OH_ann_change * (18 - 1), 18)))
# Calculate net KIE 
Sink_OH_inc = (OH_KIE*OH_Sink_inc + Cl_KIE*Cl_Sink + Strat_KIE*Strat_Sink + Soil_KIE*Soil_Sink) / (OH_Sink_inc + Cl_Sink + Strat_Sink + Soil_Sink)


# Create sink fractionation scenario based on varying OH (no chemistry)
# Constants for model
C13Std = 0.011113                   # standard for 13C-12C isotopic ratio (IUPAC 2024). Old: .0112020968 
Watm = 28.96           # Molecular weight of atmosphere (g/mole)
Matm = 5.15 * 10**21   # Mass of atmosphere (g)
Lifetime = 9
Lifetime_OH_inc = Lifetime / (OH_Sink_inc + Cl_Sink + Strat_Sink + Soil_Sink)
Lifetime_Cl_inc = Lifetime / (OH_Sink + Cl_Sink_inc + Strat_Sink + Soil_Sink)
PT = 2.815 # Conversion factor for ppb to Tg using the molar mass of the atmosphere. 
# Conversion from delta notation to 13C/(13C + 12C)
def Ratio2(y):
    return ((y / 1000 + 1) * C13Std) / ((y / 1000 + 1) * C13Std + 1)


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
Cl_Tg_double = np.linspace(19, 38, 24) 
Strat_Tg = 39
Soil_Tg = 33
OH_Tg = SumSink - Cl_Tg - Strat_Tg - Soil_Tg
OH_Tg_doubleCl = SumSink - Cl_Tg_double - Strat_Tg - Soil_Tg
OH_Tg_StableSource = SumSink_StableSource - Cl_Tg - Strat_Tg - Soil_Tg


# Calculate KIEs
# C13 first
C13KIE_OnlyOH = (OH_KIE*OH_Tg + Cl_KIE*Cl_Tg + Strat_KIE*Strat_Tg + Soil_KIE*Soil_Tg) / SumSink
C13KIE_DoubleCl = (OH_KIE*OH_Tg_doubleCl + Cl_KIE*Cl_Tg_double + Strat_KIE*Strat_Tg + Soil_KIE*Soil_Tg) / SumSink
C13KIE_StableSource = (OH_KIE*OH_Tg_StableSource + Cl_KIE*Cl_Tg + Strat_KIE*Strat_Tg + Soil_KIE*Soil_Tg) / SumSink_StableSource


#%% Begin Monte Carlo Box model analysis

# Define constants and uncertainties
# Microbial
mic_d13C =  Mic_d13C_data.iloc[:, 1].mean() # Currently, no trend in microbial d13C-CH4. Could add this to test sensitivity
mic_d13C_U =  Mic_d13C_data.iloc[:, 2].mean() 
# Fossil
ff_d13C = FF_d13C_data.iloc[28:,1] #Starting in 1999
ff_d13C_U = FF_d13C_data.iloc[28:,2] 
# Biomass burning
bb_d13C = BB_d13C_data.iloc[:,1]
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

# Define years
year = FF_d13C_data.iloc[28:,0]
# Define result matrices
d13C_Source_compiled = pd.DataFrame()
d13C_Source_RedCl_compiled = pd.DataFrame()
d13C_Source_IncOH_compiled = pd.DataFrame()
FF_compiled = pd.DataFrame()
Mic_compiled = pd.DataFrame()
FF_compiledR = pd.DataFrame()
Mic_compiledR = pd.DataFrame()
FF_RedCl_compiled = pd.DataFrame()
Mic_RedCl_compiled = pd.DataFrame()
FF_IncOH_compiled = pd.DataFrame()
Mic_IncOH_compiled = pd.DataFrame()
FF_BBdrop_compiled = pd.DataFrame()
Mic_BBdrop_compiled = pd.DataFrame()


# Begin MC analysis
for k in range(0, 1000):
    # First, define random gaussian numbers for each type of uncertainty
    RandomGauss1 = np.random.normal(0, 1, len(d13C_glob_iterations)) #Uncertainty for atmopsheriric value
    RandomGauss2 = np.random.normal(0, 1)
    RandomGauss3 = np.random.normal(0, 1)
    RandomGauss4 = np.random.normal(0, 1)
    # Run box model
    d13C_atm_MC = d13C_glob_iterations[:,k] # Add uncertainty to atmospheric d13C
    d13C_RD = Ratio2(d13C_atm_MC)
    n13C = d13C_RD*CH4*PT
    n13Cc = d13C_RD*CH4*PT
    # Calculate alphas for different Cl Scenarios
    alpha_13C = 1 / Sink
    alpha_13C_OH_inc = 1 / Sink_OH_inc 
    alpha_13C_Cl_inc = 1 / Sink_Cl_inc 
    alpha_13C_RedCl = 1 / Sink_RedCl
    alpha_13C_OnlyOH = 1 / C13KIE_OnlyOH 
    alpha_13C_DoubleCl = 1 / C13KIE_DoubleCl
    # More constants
    n13Cb = d13C_RD[0]*CH4[0]*PT
    R13Cb = Ratio2(-54.37) # The average of the 1999 and 2000 d13C_Source from the inverison
    sumR13C = np.zeros(len(CH4) - 1)
    sumR13C_RedCl = np.zeros(len(CH4) - 1)
    sumR13C_IncOH = np.zeros(len(CH4) - 1)
    sum13Cn = np.zeros(len(CH4) - 1)
    for j in range(len(CH4) - 1):
        # Calculation of source 13C:12C ratio for umoothed data
        R13C = (n13C[j + 1] - n13C[j] + n13C[j] * alpha_13C / (Lifetime)) / SumSource[j]
        sumR13C[j] = R13C
        # Calculation of source 13C:12C ratio for scenario where Cl sink proportion decreases
        R13C = (n13Cc[j + 1] - n13Cc[j] + n13Cc[j] * alpha_13C_RedCl[j] / (Lifetime)) / SumSource[j]
        sumR13C_RedCl[j] = R13C
        # Calculation of source 13C:12C ratio for scenario where OH sink increases
        R13C = (n13Cc[j + 1] - n13Cc[j] + n13Cc[j] * alpha_13C_OH_inc[j] / (Lifetime_OH_inc[j])) / SumSource_OH_inc[j]
        sumR13C_IncOH[j] = R13C
        # Calculate atmospheric d13C assuming stable source signature
        n13Cb = n13Cb + R13Cb*SumSource[j] - n13Cb*(alpha_13C)/Lifetime;
        sum13Cn[j] = n13Cb/PT;
    # calculate delta13C of source
    d13C_source = ((sumR13C - C13Std + sumR13C * C13Std) / ((C13Std - sumR13C * C13Std) / 1000))
    d13C_Source_compiled[f'Iteration_{k}'] = d13C_source
    # calculate delta13C of source for reduced cl scenario
    d13C_source_RedCl = ((sumR13C_RedCl - C13Std + sumR13C_RedCl * C13Std) / ((C13Std - sumR13C_RedCl * C13Std) / 1000))
    d13C_Source_RedCl_compiled[f'Iteration_{k}'] = d13C_source_RedCl
    # calculate delta13C of source for increased cl scenario
    d13C_source_IncOH = ((sumR13C_IncOH - C13Std + sumR13C_IncOH * C13Std) / ((C13Std - sumR13C_IncOH * C13Std) / 1000))
    d13C_Source_IncOH_compiled[f'Iteration_{k}'] = d13C_source_IncOH
    # Calculate atmospheric delta13C
    delta13C_CH4 = (((sum13Cn / (CH4[1:] - sum13Cn)) / C13Std) - 1) * 1000
    
    # Calculate mass balance with temporally varying FF signature
    # calulcate new source sigs for MC analysis
    ff_d13C_MC = ff_d13C + RandomGauss2 * ff_d13C_U
    bb_d13C_MC = mean_bb + RandomGauss4 * np.mean(bb_d13C_U) 
    FFS_ffvary = (SumSource * d13C_source - Mic_d13C_MC.iloc[:,k].values * (SumSource - BB) - bb_d13C_MC * BB) / (FF_d13C_MC_EDGAR.iloc[:,k].values - Mic_d13C_MC.iloc[:,k].values)
    MicS_ffvary = SumSource  - BB - FFS_ffvary
    # Calculate emissions relative to 1999
    FFS_ffvaryR = FFS_ffvary - FFS_ffvary[6:9].mean() #Relative to 2005 - 2007 average
    MicS_ffvaryR = MicS_ffvary - MicS_ffvary[6:9].mean() #Relative to 2005 - 2007 average
    # Calculate emissions with drop in Cl sink proportion
    FFS_RedCl = (SumSource * d13C_source_RedCl - Mic_d13C_MC.iloc[:,k].values * (SumSource - BB) - bb_d13C_MC * BB) / (FF_d13C_MC_EDGAR.iloc[:,k].values - Mic_d13C_MC.iloc[:,k].values)
    MicS_RedCl = SumSource  - BB - FFS_RedCl
    # Calculate emissions with increase in OH
    FFS_IncOH = (SumSource_OH_inc * d13C_source_IncOH - Mic_d13C_MC.iloc[:,k].values * (SumSource_OH_inc - BB) - bb_d13C_MC * BB) / (FF_d13C_MC_EDGAR.iloc[:,k].values - Mic_d13C_MC.iloc[:,k].values)
    MicS_IncOH = SumSource_OH_inc  - BB - FFS_IncOH
    # Calculate emissions assuming drop in biomass burning emissions
    FFS_BBdrop = (SumSource * d13C_source - Mic_d13C_MC.iloc[:,k].values * (SumSource - BBneg) - bb_d13C_MC * BBneg) / (FF_d13C_MC_EDGAR.iloc[:,k].values - Mic_d13C_MC.iloc[:,k].values)
    MicS_BBdrop = SumSource  - BBneg - FFS_BBdrop
  
    # Compile in results matrix
    FF_compiled[f'Iteration_{k}'] = FFS_ffvary
    Mic_compiled[f'Iteration_{k}'] = MicS_ffvary
    FF_compiledR[f'Iteration_{k}'] = FFS_ffvaryR
    Mic_compiledR[f'Iteration_{k}'] = MicS_ffvaryR
    FF_RedCl_compiled[f'Iteration_{k}'] = FFS_RedCl
    Mic_RedCl_compiled[f'Iteration_{k}'] = MicS_RedCl
    FF_IncOH_compiled[f'Iteration_{k}'] = FFS_IncOH
    Mic_IncOH_compiled[f'Iteration_{k}'] = MicS_IncOH
    FF_BBdrop_compiled[f'Iteration_{k}'] = FFS_BBdrop
    Mic_BBdrop_compiled[f'Iteration_{k}'] = MicS_BBdrop

    
#%% Smooth and calculate means/stdevs

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
            FF_compiled.iloc[k, col] +
            FF_compiled.iloc[k+1, col] +
            FF_compiled.iloc[k+2, col] +
            FF_compiled.iloc[k+3, col] +
            FF_compiled.iloc[k+4, col]
        ) / 5
# first two and last two points use shorter averages
mov_start = np.vstack([
    ((FF_compiled.iloc[0, :] + FF_compiled.iloc[1, :] + FF_compiled.iloc[2, :]) / 3).to_numpy(),
    ((FF_compiled.iloc[0, :] + FF_compiled.iloc[1, :] + FF_compiled.iloc[2, :] + FF_compiled.iloc[3, :]) / 4).to_numpy()
])
mov_end = np.vstack([
    ((FF_compiled.iloc[-4, :] + FF_compiled.iloc[-3, :] + FF_compiled.iloc[-2, :] + FF_compiled.iloc[-1, :]) / 4).to_numpy(),
    ((FF_compiled.iloc[-3, :] + FF_compiled.iloc[-2, :] + FF_compiled.iloc[-1, :]) / 3).to_numpy()
])
FF_compiled_smoothed = np.vstack([mov_start, moving_average, mov_end])
FF_compiledR_smoothed = FF_compiled_smoothed - np.mean(FF_compiled_smoothed[6:9,:],axis=0) # Relative to 2005 to 2007 average
FF_compiledR_smoothed_2013 = FF_compiled_smoothed - FF_compiled_smoothed[14,:] # relative to 2013

# Smooth FF_compiled before calculating stats
moving_averageR = np.zeros((FF_compiledR.shape[0] - 4, FF_compiledR.shape[1]))
for col in range(FF_compiledR.shape[1]):
    for k in range(FF_compiledR.shape[0] - 4):
        moving_averageR[k, col] = (
            FF_compiledR.iloc[k, col] +
            FF_compiledR.iloc[k+1, col] +
            FF_compiledR.iloc[k+2, col] +
            FF_compiledR.iloc[k+3, col] +
            FF_compiledR.iloc[k+4, col]
        ) / 5
# first two and last two points use shorter averages
mov_startR = np.vstack([
    ((FF_compiledR.iloc[0, :] + FF_compiledR.iloc[1, :] + FF_compiledR.iloc[2, :]) / 3).to_numpy(),
    ((FF_compiledR.iloc[0, :] + FF_compiledR.iloc[1, :] + FF_compiledR.iloc[2, :] + FF_compiledR.iloc[3, :]) / 4).to_numpy()
])
mov_endR = np.vstack([
    ((FF_compiledR.iloc[-4, :] + FF_compiledR.iloc[-3, :] + FF_compiledR.iloc[-2, :] + FF_compiledR.iloc[-1, :]) / 4).to_numpy(),
    ((FF_compiledR.iloc[-3, :] + FF_compiledR.iloc[-2, :] + FF_compiledR.iloc[-1, :]) / 3).to_numpy()
])
FF_compiled_smoothedRtest = np.vstack([mov_startR, moving_averageR, mov_endR])

# Smooth Mic_compiled
moving_average = np.zeros((Mic_compiled.shape[0] - 4, Mic_compiled.shape[1]))
for col in range(Mic_compiled.shape[1]):
    for k in range(Mic_compiled.shape[0] - 4):
        moving_average[k, col] = (
            Mic_compiled.iloc[k, col] +
            Mic_compiled.iloc[k+1, col] +
            Mic_compiled.iloc[k+2, col] +
            Mic_compiled.iloc[k+3, col] +
            Mic_compiled.iloc[k+4, col]
        ) / 5
mov_start = np.vstack([
    ((Mic_compiled.iloc[0, :] + Mic_compiled.iloc[1, :] + Mic_compiled.iloc[2, :]) / 3).to_numpy(),
    ((Mic_compiled.iloc[0, :] + Mic_compiled.iloc[1, :] + Mic_compiled.iloc[2, :] + Mic_compiled.iloc[3, :]) / 4).to_numpy()
])
mov_end = np.vstack([
    ((Mic_compiled.iloc[-4, :] + Mic_compiled.iloc[-3, :] + Mic_compiled.iloc[-2, :] + Mic_compiled.iloc[-1, :]) / 4).to_numpy(),
    ((Mic_compiled.iloc[-3, :] + Mic_compiled.iloc[-2, :] + Mic_compiled.iloc[-1, :]) / 3).to_numpy()
])
Mic_compiled_smoothed = np.vstack([mov_start, moving_average, mov_end])
Mic_compiledR_smoothed = Mic_compiled_smoothed - np.mean(Mic_compiled_smoothed[6:9,:],axis=0) # Relative to 2005 to 2007 average
Mic_compiledR_smoothed_2013 = Mic_compiled_smoothed - Mic_compiled_smoothed[14,:] # relative to 2013

# Smooth MicR_compiled
moving_averageR = np.zeros((Mic_compiledR.shape[0] - 4, Mic_compiledR.shape[1]))
for col in range(Mic_compiledR.shape[1]):
    for k in range(Mic_compiledR.shape[0] - 4):
        moving_averageR[k, col] = (
            Mic_compiledR.iloc[k, col] +
            Mic_compiledR.iloc[k+1, col] +
            Mic_compiledR.iloc[k+2, col] +
            Mic_compiledR.iloc[k+3, col] +
            Mic_compiledR.iloc[k+4, col]
        ) / 5
mov_startR = np.vstack([
    ((Mic_compiledR.iloc[0, :] + Mic_compiledR.iloc[1, :] + Mic_compiledR.iloc[2, :]) / 3).to_numpy(),
    ((Mic_compiledR.iloc[0, :] + Mic_compiledR.iloc[1, :] + Mic_compiledR.iloc[2, :] + Mic_compiledR.iloc[3, :]) / 4).to_numpy()
])
mov_endR = np.vstack([
    ((Mic_compiledR.iloc[-4, :] + Mic_compiledR.iloc[-3, :] + Mic_compiledR.iloc[-2, :] + Mic_compiledR.iloc[-1, :]) / 4).to_numpy(),
    ((Mic_compiledR.iloc[-3, :] + Mic_compiledR.iloc[-2, :] + Mic_compiledR.iloc[-1, :]) / 3).to_numpy()
])
Mic_compiledR_smoothedRtest = np.vstack([mov_startR, moving_averageR, mov_endR])


# Calculate statistics
FF_mean = FF_compiled_smoothed.mean(axis=1)
FF_std = FF_compiled_smoothed.std(axis=1)
Mic_mean = Mic_compiled_smoothed.mean(axis=1)
Mic_std = Mic_compiled_smoothed.std(axis=1)
d13C_source_meanMC = d13C_Source_compiled.mean(axis=1).to_numpy()
Mic_RedCl_mean = Mic_RedCl_compiled.mean(axis=1).to_numpy()
FF_RedCl_mean = FF_RedCl_compiled.mean(axis=1).to_numpy()
Mic_IncOH_mean = Mic_IncOH_compiled.mean(axis=1).to_numpy()
FF_IncOH_mean = FF_IncOH_compiled.mean(axis=1).to_numpy()
Mic_BBdrop_mean = Mic_BBdrop_compiled.mean(axis=1).to_numpy()
FF_BBdrop_mean = FF_BBdrop_compiled.mean(axis=1).to_numpy()
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
axes[0, 0].plot(years[:-1], FF_compiled_smoothed - FF_compiled_smoothed[0, :])
axes[0, 0].set_ylabel('Relative FF Emissions Tg yr$^{-1}$')
axes[0, 0].set_xlim(1998, 2022)
axes[0, 0].set_ylim(-120, 120)
axes[0, 0].grid(False)
axes[0, 1].plot(years[:-1], FF_compiled_smoothed)
axes[0, 1].set_ylabel('FF Emissions Tg yr$^{-1}$')
axes[0, 1].set_xlim(1998, 2022)
axes[0, 1].set_ylim(10, 250)
axes[0, 1].grid(False)

# bottom row title
fig.text(0.5, 0.48, r'${\delta}^{13}C$-derived Microbial Emissions', ha='center', fontsize=16)
# bottom row
axes[1, 0].plot(years[:-1], Mic_compiled_smoothed - Mic_compiled_smoothed[0, :])
axes[1, 0].set_ylabel('Relative Mic Emissions Tg yr$^{-1}$')
axes[1, 0].set_xlim(1998, 2022)
axes[1, 0].set_ylim(-90, 200)
axes[1, 0].set_xlabel('Year')
axes[1, 0].grid(False)
axes[1, 1].plot(years[:-1], Mic_compiled_smoothed)
axes[1, 1].set_ylabel('Mic Emissions Tg yr$^{-1}$')
axes[1, 1].set_xlim(1998, 2022)
axes[1, 1].set_ylim(290, 580)
axes[1, 1].set_xlabel('Year')
axes[1, 1].grid(False)
plt.tight_layout(rect=[0, 0.05, 1, 0.90])
plt.show()


# # Save spaghetti emissions 
# df_spaghetti = pd.DataFrame(FF_compiled_smoothed)
# # save to Excel
# df_spaghetti.to_excel('Output/FF_compiled_smoothed_d13C.xlsx', index=False)


#%% Plot histogram of difference (2020 to 2022 - 2005 to 2007)

# Calculate Deltas (2020 to 2022 - 2005 to 2007)
Mic_Delta_MC = Mic_compiled.iloc[-3:, :].mean(axis=0) - Mic_compiled.iloc[6:8, :].mean(axis=0)
FF_Delta_MC = FF_compiled.iloc[-3:, :].mean(axis=0) - FF_compiled.iloc[6:8, :].mean(axis=0)
Mic_Delta_MC_RedCl = Mic_RedCl_compiled.iloc[-3:, :].mean(axis=0) - Mic_RedCl_compiled.iloc[6:8, :].mean(axis=0)
FF_Delta_MC_RedCl = FF_RedCl_compiled.iloc[-3:, :].mean(axis=0) - FF_RedCl_compiled.iloc[6:8, :].mean(axis=0)
Mic_Delta_MC_IncOH = Mic_IncOH_compiled.iloc[-3:, :].mean(axis=0) - Mic_IncOH_compiled.iloc[6:8, :].mean(axis=0)
FF_Delta_MC_IncOH = FF_IncOH_compiled.iloc[-3:, :].mean(axis=0) - FF_IncOH_compiled.iloc[6:8, :].mean(axis=0)
Mic_Delta_MC_BBdrop = Mic_BBdrop_compiled.iloc[-3:, :].mean(axis=0) - Mic_BBdrop_compiled.iloc[6:8, :].mean(axis=0)
FF_Delta_MC_BBdrop = FF_BBdrop_compiled.iloc[-3:, :].mean(axis=0) - FF_BBdrop_compiled.iloc[6:8, :].mean(axis=0)
# Calculate Deltas  (2020 - 2000 to 2009)
Mic_Delta_MC_2020 = Mic_compiled.iloc[-3, :] - Mic_compiled.iloc[1:11, :].mean(axis=0)
FF_Delta_MC_2020 = FF_compiled.iloc[-3, :] - FF_compiled.iloc[1:11, :].mean(axis=0)


# Save the data
# Prepare the data for the first set of plots (absolute emissions)
df_absolute = pd.DataFrame({'Year': year,'FF_mean': FF_mean,'FF_std': FF_std,'Mic_mean': Mic_mean,'Mic_std': Mic_std,})
# Save the first DataFrame to a CSV file
df_absolute.to_csv('Output/Results_d13C-MassBalance_Cantrell.csv', index=False)
# Prepare the data for the second set of plots (relative emissions)
df_relative = pd.DataFrame({'Year': year,'FF_meanR': FF_meanR,'FF_stdR': FF_stdR,'Mic_meanR': Mic_meanR,'Mic_stdR': Mic_stdR,})
# Save the second DataFrame to a CSV file
df_relative.to_csv('Output/Results_Rd13C-MassBalance_Cantrell.csv', index=False)

# Save histogram data
matrix = np.column_stack((Mic_Delta_MC, FF_Delta_MC, Mic_Delta_MC_2020, FF_Delta_MC_2020))
np.savetxt("Output/d13C_histogram_Cantrell.csv", matrix, delimiter=",")


# # Prepare sensitivity test results and export 
# df_sensitivity = pd.DataFrame({'Year': year,'FF Reduced BB': FF_BBdrop_mean_mov,'Mic Reduced BB': Mic_BBdrop_mean_mov,'FF Reduced Cl': FF_RedCl_mean_mov,'Mic Reduced Cl': Mic_RedCl_mean_mov,})
# df_sensitivity.to_csv('Output/Results_d13C-BBClSensitivity.csv', index=False)  
# # Prepare FF source sig test results and export 
# df_FFsrcsig = pd.DataFrame({'Year': year,'FF, Stable Mic': FFS_StableMic_mean_mov,'FF, Half Mic': FFS_Half_mean_mov,'d13C FF, Stable Mic': FF_SrcSig_mean_mov,'d13C FF, Half Mic': FF_SrcSig_Half_mean_mov,})
# df_FFsrcsig.to_csv('Output/Results_d13C-FF_Sensitivity.csv', index=False) 
    