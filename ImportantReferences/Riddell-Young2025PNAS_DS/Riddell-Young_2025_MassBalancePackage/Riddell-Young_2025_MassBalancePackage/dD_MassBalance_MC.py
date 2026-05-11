#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 15:03:37 2024

@author: ryoung
"""

# Code for Monte Carlo analysis of onebox deltaD-CH4 mass balance
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
        # The same is true for tests using the global mean FF source signature using CT-CH4 posterior flux-weighting.


#%% First, load data

# Load updated annual mean DEI
glob_ann_dD_data = np.loadtxt("../Riddell-Young_2025_dD_GlobMean/output/GlobMean_dD_dei_UmezawaCal_noBUDS.csv", delimiter=",", dtype=np.float64, skiprows=1)
glob_ann_dD = glob_ann_dD_data[:,1] - 0.5
glob_ann_dD_unc = glob_ann_dD_data[:,2]
glob_ann_dD_uncR = glob_ann_dD_data[:,4]
glob_ann_dD_years = glob_ann_dD_data[:,0]

# # Load all dD DEI iterations
dD_AnnAvg_matrix = pd.read_excel('../Riddell-Young_2025_dD_GlobMean/output/GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx', header=None).to_numpy()

# Load GML global annual means for CH4
CH4data = pd.read_excel('data/GML_CH4_AnnualMean.xlsx').to_numpy()
CH4data = CH4data[21:,:]
CH4 = CH4data[:,1] # Starting in 1995

# Load Carbon tracker methane
data2 = pd.read_excel('data/CarbonTracker_CH4.xlsx')
years = data2.iloc[:,0].values
bbCT = data2.iloc[:,9].values #Prior emissions (GFED4)

# Load annual mean source signatures
BB_dD_data = pd.read_csv('Output/BB_dD_annual.csv', delimiter=',', header = None)
Mic_dD_data = pd.read_csv('Output/Mic_dD_AnnGlob.csv', delimiter=',', header = None)
Mic_dD_MC_trends = pd.read_csv('Output/Mic_dD_MC.csv', delimiter=',', header = None) 
Mic_dD_MC = Mic_dD_MC_trends.iloc[6:,1:]

# Load FF source sigantures
FF_dD_data = pd.read_csv('Output/FF_dD_GlobUnc.csv', delimiter=',')
FF_dD_MC_CTCH4_data = pd.read_csv('Output/FF_dD_GlobMC_CTCH4.csv', delimiter=',', header = None)
FF_dD_MC_EDGAR_data = pd.read_csv('Output/FF_dD_GlobMC_EDGAR.csv', delimiter=',')
FF_dD_MC_EDGAR = FF_dD_MC_EDGAR_data.iloc[34:,1:]
FF_dD_MC_CTCH4 = FF_dD_MC_CTCH4_data.iloc[7:,1:]
FF_dD_MC_CTCH4 = pd.concat([FF_dD_MC_CTCH4, FF_dD_MC_CTCH4.iloc[[-1]]], ignore_index=True)
# Calculate CTCH4 statistics
FF_dD_mean_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
FF_dD_mean_EDGAR = FF_dD_MC_EDGAR.iloc[:, 1:].mean(axis=1).to_numpy()
FF_dD_std_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()
R_FF_dD_mean_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
R_FF_dD_std_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy() 


#%% Load model components

# Calculate sink KIE
OH_KIE = 1.294     #Saueressig et al,. 2001. Best estimate in Ben Li's opinion (Whitehill-Joelson Avg = 1.327)
Cl_KIE = 1.52     #Also Saueressig
Strat_KIE = 1.179 # Dyonisius et al., 2020; Beck et al., 2018 
Soil_KIE = 1.083  #Snover and Quay. Only number available.
# Sink strengths (carbontracker documentation)
OH_Sink = .835
Cl_Sink = .035
Strat_Sink = .07
Soil_Sink = .06
# Net KIE
Sink = OH_KIE*OH_Sink + Cl_KIE*Cl_Sink + Strat_KIE*Strat_Sink + Soil_KIE*Soil_Sink  #Thanwerdas 2022 uses this number

# Calculate change in KIE if Cl sink decreases
Cl_Sink_red = np.linspace(Cl_Sink, .011, 18)
Sink_Rem = Cl_Sink - Cl_Sink_red
OH_Sink_RedCl = OH_Sink + Sink_Rem * OH_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Strat_Sink_RedCl = Strat_Sink + Sink_Rem * Strat_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Soil_Sink_RedCl = Soil_Sink + Sink_Rem * Soil_Sink / (OH_Sink + Strat_Sink + Soil_Sink)
Check = OH_Sink_RedCl + Strat_Sink_RedCl + Soil_Sink_RedCl + Cl_Sink_red
Sink_RedCl = OH_KIE*OH_Sink_RedCl + Cl_KIE*Cl_Sink_red + Strat_KIE*Strat_Sink_RedCl + Soil_KIE*Soil_Sink_RedCl

# Calculate change in KIE if Cl sink increases (Van Herpen et al., 2023)
Cl_Sink_inc_frac = np.array([1, 0.97, 0.94, 0.58, 0.70, 1.13, 0.79, 1.09, 0.80, 1.48, 1.86, 2.08, 1.62, 1.85, 1.85, 1.85, 1.85, 1.85])
Cl_Sink_inc = Cl_Sink_inc_frac*Cl_Sink
Sink_IncCl = (OH_KIE*OH_Sink + Cl_KIE*Cl_Sink_inc + Strat_KIE*Strat_Sink + Soil_KIE*Soil_Sink) / (OH_Sink + Cl_Sink_inc + Soil_Sink + Strat_Sink)

# Create increasing OH following Olaf et al., trends
OH_ann_change = OH_Sink*.003
OH_Sink_inc = np.linspace(OH_Sink, OH_Sink + OH_ann_change * (18 - 1), 18)
# Calculate net KIE 
Sink_OH_inc = (OH_KIE*OH_Sink_inc + Cl_KIE*Cl_Sink + Strat_KIE*Strat_Sink + Soil_KIE*Soil_Sink) / (OH_Sink_inc + Cl_Sink + Strat_Sink + Soil_Sink)

# Constants for model
DStd = 0.00015576      # standard for D-H isotopic ratio
Watm = 28.96           # Molecular weight of atmosphere (g/mole)
Matm = 5.15 * 10**21   # Mass of atmosphere (g)
Lifetime = 9
Lifetime_OH_inc = Lifetime / (OH_Sink_inc + Cl_Sink + Strat_Sink + Soil_Sink)
Lifetime_Cl_inc = Lifetime / (OH_Sink + Cl_Sink_inc + Strat_Sink + Soil_Sink)
Dstart = -94.94
PT = 2.815 # Conversion factor for ppb to Tg using the molar mass of the atmosphere. 

# Conversion from delta notation to 13C/(13C + 12C)
def Ratio2(y):
    return ((y / 1000 + 1) * DStd) / ((y / 1000 + 1) * DStd + 1)


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
    
    
#%% Begin Monte Carlo Box model analysis

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
# Define years
year = FF_dD_data.iloc[34:,0]
# Define result matrices
dD_Source_compiled = pd.DataFrame()
dD_Source_RedCl_compiled = pd.DataFrame()
dD_Source_IncOH_compiled = pd.DataFrame()
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
    RandomGauss2 = np.random.normal(0, 1)
    RandomGauss3 = np.random.normal(0, 1)
    RandomGauss4 = np.random.normal(0, 1)
    # Run box model
    dD_atm_MC = dD_AnnAvg_matrix[:,(k+1)] # Add uncertainty to atmospheric dD
    dD_RD = Ratio2(dD_atm_MC)
    nD = dD_RD*CH4*PT
    nDc = dD_RD*CH4*PT
    nDd = dD_RD*CH4*PT
    alpha_D = 1 / Sink
    alpha_D_RedCl = 1 / Sink_RedCl
    alpha_D_IncCl = 1 / Sink_IncCl
    alpha_D_OH_inc = 1 / Sink_OH_inc
    RDb = Ratio2(-267.5) # The first value of dD_Source from the inverison
    sumRD = np.zeros(len(CH4) - 1)
    sumRD_RedCl = np.zeros(len(CH4) - 1)
    sumRD_IncOH = np.zeros(len(CH4) - 1)
    for j in range(len(CH4) - 1):
        # Calculation of source D:H ratio for unsmoothed data
        RD = (nD[j + 1] - nD[j] + nD[j] * alpha_D / (Lifetime)) / SumSource[j]
        sumRD[j] = RD
        # Calculation of source D:H ratio for scenario where Cl sink proportion decreases
        RD = (nDc[j + 1] - nDc[j] + nDc[j] * alpha_D_RedCl[j] / (Lifetime)) / SumSource[j]
        sumRD_RedCl[j] = RD
        # Calculation of source D:H ratio for scenario where OH sink increases
        RD = (nDd[j + 1] - nDd[j] + nDd[j] * alpha_D_OH_inc[j] / (Lifetime_OH_inc[j])) / SumSource_OH_inc[j]
        sumRD_IncOH[j] = RD
    # calculate deltaD of source
    dD_source = ((sumRD - DStd + sumRD * DStd) / ((DStd - sumRD * DStd) / 1000))
    dD_Source_compiled[f'Iteration_{k}'] = dD_source
    # calculate deltaD of source for reduced cl scenario
    dD_source_RedCl = ((sumRD_RedCl - DStd + sumRD_RedCl * DStd) / ((DStd - sumRD_RedCl * DStd) / 1000))
    dD_Source_RedCl_compiled[f'Iteration_{k}'] = dD_source_RedCl
    # calculate deltaD of source for increased cl scenario
    dD_source_IncOH = ((sumRD_IncOH - DStd + sumRD_IncOH * DStd) / ((DStd - sumRD_IncOH * DStd) / 1000))
    dD_Source_IncOH_compiled[f'Iteration_{k}'] = dD_source_IncOH
    
    # Calculate mass balance with temporally varying FF signature
    bb_dd_MC = mean_bb + RandomGauss4 * mean_bb_U
    # Calculate Mic and FF emissions
    FFS_ffvary = (SumSource * dD_source - Mic_dD_MC.iloc[:,k].values * (SumSource - BB) - bb_dd_MC * BB) / (FF_dD_MC_EDGAR.iloc[:,k].values - Mic_dD_MC.iloc[:,k].values)
    MicS_ffvary = SumSource  - BB - FFS_ffvary
    # Calculate emissions relative to 1998
    FFS_ffvaryR = FFS_ffvary - FFS_ffvary[0:3].mean() #Relative to 2005 - 2007 average
    MicS_ffvaryR = MicS_ffvary - MicS_ffvary[0:3].mean() #Relative to 2005 - 2007 average
    # Calculate emissions with drop in Cl sink proportion
    FFS_RedCl = (SumSource * dD_source_RedCl - Mic_dD_MC.iloc[:,k].values * (SumSource - BB) - bb_dd_MC * BB) / (FF_dD_MC_EDGAR.iloc[:,k].values - Mic_dD_MC.iloc[:,k].values)
    MicS_RedCl = SumSource  - BB - FFS_RedCl
    # Calculate emissions with increase in OH sink
    FFS_IncOH = (SumSource_OH_inc * dD_source_IncOH - Mic_dD_MC.iloc[:,k].values * (SumSource_OH_inc - BB) - bb_dd_MC * BB) / (FF_dD_MC_EDGAR.iloc[:,k].values - Mic_dD_MC.iloc[:,k].values)
    MicS_IncOH = SumSource_OH_inc  - BB - FFS_IncOH
    # Calculate emissions assuming drop in biomass burning emissions
    FFS_BBdrop = (SumSource * dD_source - Mic_dD_MC.iloc[:,k].values * (SumSource - BBneg) - bb_dd_MC * BBneg) / (FF_dD_MC_EDGAR.iloc[:,k].values - Mic_dD_MC.iloc[:,k].values)
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

# Convert to numpy    
Mic_compiled = Mic_compiled.to_numpy()
Mic_compiledR = Mic_compiledR.to_numpy()
FF_compiled = FF_compiled.to_numpy()
FF_compiledR = FF_compiledR.to_numpy()
Mic_RedCl_compiled = Mic_RedCl_compiled.to_numpy()
FF_RedCl_compiled = FF_RedCl_compiled.to_numpy()
Mic_IncOH_compiled = Mic_IncOH_compiled.to_numpy()
FF_IncOH_compiled = FF_IncOH_compiled.to_numpy()
Mic_BBdrop_compiled = Mic_BBdrop_compiled.to_numpy()
FF_BBdrop_compiled = FF_BBdrop_compiled.to_numpy()


#%% Plot histogram of difference (2020 to 2022 - 2005 to 2007)

# Calculate Deltas
Mic_Delta_MC = Mic_compiled[-3:, :].mean(axis=0) - Mic_compiled[:3, :].mean(axis=0)
FF_Delta_MC = FF_compiled[-3:, :].mean(axis=0) - FF_compiled[:3, :].mean(axis=0)
Mic_Delta_MC_RedCl = Mic_RedCl_compiled[-3:, :].mean(axis=0) - Mic_RedCl_compiled[:3, :].mean(axis=0)
FF_Delta_MC_RedCl = FF_RedCl_compiled[-3:, :].mean(axis=0) - FF_RedCl_compiled[:3, :].mean(axis=0)
Mic_Delta_MC_IncOH = Mic_IncOH_compiled[-3:, :].mean(axis=0) - Mic_IncOH_compiled[:3, :].mean(axis=0)
FF_Delta_MC_Incv = FF_IncOH_compiled[-3:, :].mean(axis=0) - FF_IncOH_compiled[:3, :].mean(axis=0)
Mic_Delta_MC_BBdrop = Mic_BBdrop_compiled[-3:, :].mean(axis=0) - Mic_BBdrop_compiled[:3, :].mean(axis=0)
FF_Delta_MC_BBdrop = FF_BBdrop_compiled[-3:, :].mean(axis=0) - FF_BBdrop_compiled[:3, :].mean(axis=0)
# Calculate delta from 2020
Mic_Delta_MC_2020 = Mic_compiled[-3, :] - Mic_compiled[:3, :].mean(axis=0)
FF_Delta_MC_2020 = FF_compiled[-3, :] - FF_compiled[:3, :].mean(axis=0)


#%% Smooth results

# Create 5 year smoothing function
def smooth5(data):
    smoothed = np.zeros_like(data)
    for col in range(data.shape[1]):
        # first two points
        smoothed[0, col] = np.mean(data[0:3, col])
        smoothed[1, col] = np.mean(data[0:4, col])
        # centered 5-point moving average
        for k in range(2, data.shape[0] - 2):
            smoothed[k, col] = np.mean(data[k-2:k+3, col])
        # last two points
        smoothed[-2, col] = np.mean(data[-4:, col])
        smoothed[-1, col] = np.mean(data[-3:, col])
    return smoothed

# # create double 3-year smoothing function
# def smooth3_double(data):
#     def smooth3_once(arr):
#         smoothed = np.zeros_like(arr)
#         for col in range(arr.shape[1]):
#             # first point (average of first 2)
#             smoothed[0, col] = np.mean(arr[0:2, col])
#             # middle points (3-year moving average)
#             for k in range(1, arr.shape[0] - 1):
#                 smoothed[k, col] = np.mean(arr[k-1:k+2, col])
#             # last point (average of last 2)
#             smoothed[-1, col] = np.mean(arr[-2:, col])
#         return smoothed
    
#     # apply twice
#     return smooth3_once(smooth3_once(data))

# Smooth the FF_compiled record
FF_compiled_smoothed = smooth5(FF_compiled)
FF_compiledR_smoothed = FF_compiled_smoothed - np.mean(FF_compiled_smoothed[0:3,:],axis=0)
FF_compiledR_smoothed_2013 = FF_compiled_smoothed - FF_compiled_smoothed[8,:] # Relative to 2013
# Smooth the FF_compiledR record
FF_compiled_smoothedRtest = smooth5(FF_compiledR)
# Smooth the Mic_compiled record
Mic_compiled_smoothed = smooth5(Mic_compiled)
Mic_compiledR_smoothed = Mic_compiled_smoothed - np.mean(Mic_compiled_smoothed[0:3,:],axis=0)
Mic_compiledR_smoothed_2013 = Mic_compiled_smoothed - Mic_compiled_smoothed[8,:] # Relative to 2013
# Smooth the Mic_compiledR record
Mic_compiled_smoothedRtest = smooth5(Mic_compiledR)

#Smooth sensitivity tests
# Reduced Chlorine sink Sensitivity tests
FF_RedCl_compiled_smoothed = smooth5(FF_RedCl_compiled)
Mic_RedCl_compiled_smoothed = smooth5(Mic_RedCl_compiled)

# Reduced Biomass Burning Sensitivity tests
FF_BBdrop_compiled_smoothed = smooth5(FF_BBdrop_compiled)
Mic_BBdrop_compiled_smoothed = smooth5(Mic_BBdrop_compiled)

# Increased OH sink Sensitivity tests
FF_IncOH_compiled_smoothed = smooth5(FF_IncOH_compiled)
Mic_IncOH_compiled_smoothed = smooth5(Mic_IncOH_compiled)

# Calculate statistics
FF_mean = FF_compiled_smoothed.mean(axis=1) 
FF_std = FF_compiled_smoothed.std(axis=1) 
Mic_mean = Mic_compiled_smoothed.mean(axis=1) 
Mic_std = Mic_compiled_smoothed.std(axis=1) 
dD_source_meanMC = dD_Source_compiled.mean(axis=1) 
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
# Calculate sensitivity statistics
FF_RedCl_mean = FF_RedCl_compiled_smoothed.mean(axis=1)  
Mic_RedCl_mean = Mic_RedCl_compiled_smoothed.mean(axis=1) 
FF_BBdrop_mean = FF_BBdrop_compiled_smoothed.mean(axis=1)  
Mic_BBdrop_mean = Mic_BBdrop_compiled_smoothed.mean(axis=1) 
FF_IncOH_mean = FF_IncOH_compiled_smoothed.mean(axis=1)  
Mic_IncOH_mean = Mic_IncOH_compiled_smoothed.mean(axis=1) 


# plot spaghetti diagram for dD
plt.rc('font', size=14) 
fig, axes = plt.subplots(2, 2, figsize=(10, 10), dpi=300, sharex=True, gridspec_kw={'hspace': 0.35, 'wspace': 0.3})

years_dD = years[7:] + 1

# top row title
fig.text(0.5, 0.90, '${\delta}D$-derived Fossil Fuel Emissions', ha='center', fontsize=16)
# top row
axes[0, 0].plot(years_dD, FF_compiledR_smoothed - FF_compiledR_smoothed[0, :])
axes[0, 0].set_ylabel('Relative FF Emissions Tg yr$^{-1}$')
axes[0, 0].set_xlim(1998, 2022)
axes[0, 0].set_ylim(-130, 100)
axes[0, 0].grid(False)
axes[0, 1].plot(years_dD, FF_compiled_smoothed)
axes[0, 1].set_ylabel('FF Emissions Tg yr$^{-1}$')
axes[0, 1].set_xlim(1998, 2022)
axes[0, 1].set_ylim(0, 230)
axes[0, 1].grid(False)

# bottom row title
fig.text(0.5, 0.48, '${\delta}D$-derived Microbial Emissions', ha='center', fontsize=16)
# bottom row
axes[1, 0].plot(years_dD, Mic_compiledR_smoothed - Mic_compiledR_smoothed[0, :])
axes[1, 0].set_ylabel('Relative Mic Emissions Tg yr$^{-1}$')
axes[1, 0].set_xlim(1998, 2022)
axes[1, 0].set_ylim(-80, 220)
axes[1, 0].set_xlabel('Year')
axes[1, 0].grid(False)
axes[1, 1].plot(years_dD, Mic_compiled_smoothed)
axes[1, 1].set_ylabel('Mic Emissions Tg yr$^{-1}$')
axes[1, 1].set_xlim(1998, 2022)
axes[1, 1].set_ylim(310, 610)
axes[1, 1].set_xlabel('Year')
axes[1, 1].grid(False)
plt.tight_layout(rect=[0, 0.05, 1, 0.90])
plt.show()

# Save spaghetti emissions 
df_spaghetti = pd.DataFrame(FF_compiled_smoothed)
# save to Excel
df_spaghetti.to_excel('Output/FF_compiled_smoothed_dD.xlsx', index=False)

    
#%% Save the data
# Prepare the data for the first set of plots (absolute emissions)
df_absolute = pd.DataFrame({'Year': year,'FF_mean': FF_mean,'FF_std': FF_std,'Mic_mean': Mic_mean,'Mic_std': Mic_std,})
# Save the first DataFrame to a CSV file
df_absolute.to_csv('Output/Results_dD-MassBalance_UmezawaCal_noBUDS.csv', index=False)
# Prepare the data for the second set of plots (relative emissions)
df_relative = pd.DataFrame({'Year': year,'FF_meanR': FF_meanR,'FF_stdR': FF_stdR,'Mic_meanR': Mic_meanR,'Mic_stdR': Mic_stdR,})
# Save the second DataFrame to a CSV file
df_relative.to_csv('Output/Results_RdD-MassBalance_UmezawaCal_noBUDS.csv', index=False)

# Save histogram data
matrix = np.column_stack((Mic_Delta_MC, FF_Delta_MC, Mic_Delta_MC_2020, FF_Delta_MC_2020))
np.savetxt("Output/dD_histogram_UmezawaCal_noBUDS.csv", matrix, delimiter=",")


# # Prepare sensitivity test results and export 
# df_sensitivity = pd.DataFrame({'Year': year,'FF Reduced BB': FF_BBdrop_mean_mov,'Mic Reduced BB': Mic_BBdrop_mean_mov,'FF Reduced Cl': FF_RedCl_mean_mov,'Mic Reduced Cl': Mic_RedCl_mean_mov,})
# df_sensitivity.to_csv('Output/Results_dD-BBClSensitivity.csv', index=False)  


    
    