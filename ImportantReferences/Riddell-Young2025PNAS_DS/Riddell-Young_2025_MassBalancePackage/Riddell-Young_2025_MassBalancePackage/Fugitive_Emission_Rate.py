#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 10:58:52 2024

@author: ryoung
"""

# This code calculates the fugutive emission rate of fossil emissions based on the mass balance data and bottom up estimates, 
# extending MC error propagation from main mass balance codes.

locals().clear()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#%% Load posterior data

# Load EDGAR dD mass balance
# Load dD absolute emissions
df_absolute = pd.read_csv('Output/Results_dD-MassBalance_UmezawaCal_noBUDS.csv')
# Separate out each variable from the absolute emissions dataset
year_dD = df_absolute['Year']
FF_mean_dD = df_absolute['FF_mean']
FF_std_dD = df_absolute['FF_std']
Mic_mean_dD = df_absolute['Mic_mean']
Mic_std_dD = df_absolute['Mic_std']
# Load dD relative emissions
df_relative_dD = pd.read_csv('Output/Results_RdD-MassBalance_UmezawaCal_noBUDS.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD = df_relative_dD['FF_meanR']
FF_stdR_dD = df_relative_dD['FF_stdR']
#FF_stdR_dD = FF_stdR_dD - FF_stdR_dD[5] # Make relative to 2010
Mic_meanR_dD = df_relative_dD['Mic_meanR']
Mic_stdR_dD = df_relative_dD['Mic_stdR']

# Load the spaghetti output 
FF_dD_MC = pd.read_excel('Output/FF_compiled_smoothed_dD.xlsx').to_numpy()
FF_d13C_MC = pd.read_excel('Output/FF_compiled_smoothed_d13C.xlsx').to_numpy()

# Load EDGAR d13C mass balance
# Load d13C absolute emissions
df_absolute = pd.read_csv('Output/Results_d13C-MassBalance_Cantrell.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C = df_absolute['Year']
FF_mean_d13C = df_absolute['FF_mean']
FF_std_d13C = df_absolute['FF_std']
Mic_mean_d13C = df_absolute['Mic_mean']
Mic_std_d13C = df_absolute['Mic_std']
# Load d13C relative emissions
df_relative = pd.read_csv('Output/Results_Rd13C-MassBalance_Cantrell.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C = df_relative['FF_meanR']
FF_stdR_d13C = df_relative['FF_stdR']
#FF_stdR_d13C = FF_stdR_d13C - FF_stdR_d13C[11] # Make relative to 2010
Mic_meanR_d13C = df_relative['Mic_meanR']
Mic_stdR_d13C = df_relative['Mic_stdR']

# Load Carbon tracker methane
data2 = pd.read_excel('data/CarbonTracker_CH4.xlsx')
totalCT = data2.iloc[:,15].values
yearss = data2.iloc[:,0].values
micCT = data2.iloc[:,7].values
ffCT = data2.iloc[:,3].values
bbCT = data2.iloc[:,11].values

# Load Schwietzke et al., 2016 box model results
SchwietzkeData = pd.read_excel('data/Schwietzke_2016.xlsx').to_numpy()
SchwietzkeYear = SchwietzkeData[:,0]
SchwietzkeFF = SchwietzkeData[:,1]
SchwietzkeFFStd = SchwietzkeData[:,2]
SchwietzkeMic = SchwietzkeData[:,3]
SchwietzkeMicStd = SchwietzkeData[:,4]
# Load Schwietzke et al., 2016 FER
Schwietzke_FER_Data = pd.read_excel('data/Schwietzke_2016_FER.xlsx').to_numpy()
SchwietzkeYear_FER = Schwietzke_FER_Data[:,0]
SchwietzkeFER = Schwietzke_FER_Data[:,1]

# Load methane intensity estimates
CH4_intensity = pd.read_excel('data/CH4_intensity.xlsx')
CH4int_year = CH4_intensity.iloc[:, 0].to_numpy()
CH4int_Lu = CH4_intensity.iloc[:, 1].to_numpy()
CH4int_Mcdonald = CH4_intensity.iloc[:, 2].to_numpy()


#%% Load prior/bottom-up data

# Load EDGAR emissions
Country_ONG_emis = pd.read_csv('data/EDGAR8_ONG.csv') #Ignores shale gas contributions country-wide for simplicity
Country_Coal_emis = pd.read_csv('data/EDGAR8_Coal.csv') #Ignores shale gas contributions country-wide for simplicity
glob_ONG_emis = Country_ONG_emis.iloc[:, 1:].sum(axis=0).values.reshape(-1, 1)
glob_Coal_emis = Country_Coal_emis.iloc[:, 1:].sum(axis=0).values.reshape(-1, 1)
glob_Coal_emis = glob_Coal_emis.ravel()/1000 # Convert from gigagrams to teragrams
glob_ONG_emis = glob_ONG_emis.ravel()/1000 # Convert from gigagrams to teragrams
year_EDGAR = np.arange(1970, 2023) 

# Load EIA production data
EIA_data = pd.read_excel('data/EIA_ONG_prod.xlsx')
EIA_year = EIA_data.iloc[:,0].values
oil_prod = EIA_data.iloc[:,1].values #million barrels per day
NG_prod = EIA_data.iloc[:,2].values #billion cubic feet


#%% Calculations

# Calculate emissions 
oil_ch4_conv = 5000
oil_emis = oil_prod/oil_ch4_conv

# Define constants
# Define geologcial emissions
geo_emis = 5 # Hmiel et al., 2020
# Define NG CH4 content
NG_content = 0.86 # Schwietzke et al., 2016
# Tg Ch4 per billion cubic foot
Tg_per_bcf = 1000000000*1.1981*16.04/1000000000000
# Calculate total NG CH4 production 
NG_CH4_prod_bcf = NG_prod * NG_content
# Convert from billion cubic feet to teragrams
NG_CH4_prod = NG_CH4_prod_bcf * Tg_per_bcf

# Combine data into dataframe
# Create dataframes for each pairing
df_CT = pd.DataFrame({'Year': yearss, 'FF_mean_dD': ffCT})
df_Schwietzke = pd.DataFrame({'Year': SchwietzkeYear, 'FF_mean_dD': SchwietzkeFF})
df_EDGAR = pd.DataFrame({'Year': year_EDGAR, 'glob_Coal_emis': glob_Coal_emis})
df_EIA = pd.DataFrame({'Year': EIA_year, 'oil_emis': oil_emis, 'NG Prod': NG_CH4_prod})

FER_dD_list = []
FER_d13C_list = []
FER_dD_list_StableCoal = []
FER_d13C_list_StableCoal = []

for i in range(1000):
    # Combine data into dataframe
    # Create dataframes for each pairing
    df_dD = pd.DataFrame({'Year': year_dD, 'FF_mean_dD': FF_dD_MC[:,i]})
    df_d13C = pd.DataFrame({'Year': year_d13C, 'FF_mean_dD': FF_d13C_MC[:,i]})
    
    # for dD, merge all three dataframes on 'Year', preserving all years using an outer join
    merged_dD = pd.merge(df_dD, df_EDGAR, on='Year', how='outer')
    merged_dD = pd.merge(merged_dD, df_EIA, on='Year', how='outer')
    # Sort by Year to ensure the data is in chronological order
    merged_dD = merged_dD.sort_values('Year')
    merged_dD = merged_dD.dropna()

    # Repeat for d13C
    merged_d13C = pd.merge(df_d13C, df_EDGAR, on='Year', how='outer')
    merged_d13C = pd.merge(merged_d13C, df_EIA, on='Year', how='outer')
    # Sort by Year to ensure the data is in chronological order
    merged_d13C = merged_d13C.sort_values('Year')
    merged_d13C = merged_d13C.dropna()
    
    # Calculate FER using dD and d13C fossil emissions
    # Calculate NG emissions
    merged_dD['NG Emissions'] = merged_dD.iloc[:, 1] - merged_dD.iloc[:, 2] - merged_dD.iloc[:, 3] - geo_emis
    merged_dD['FER'] = merged_dD['NG Emissions'] / merged_dD['NG Prod'] * 100
    merged_d13C['NG Emissions'] = merged_d13C.iloc[:, 1] - merged_d13C.iloc[:, 2] - merged_d13C.iloc[:, 3] - geo_emis
    merged_d13C['FER'] = merged_d13C['NG Emissions'] / merged_d13C['NG Prod'] * 100
    
    # Calculate assuming stable coal
    merged_dD['NG Emissions: Stable Coal'] = merged_dD.iloc[:, 1] - merged_dD.iloc[5, 2] - merged_dD.iloc[:, 3] - geo_emis # Stable at 2010
    merged_dD['FER: Stable Coal'] = merged_dD['NG Emissions: Stable Coal'] / merged_dD['NG Prod'] * 100
    merged_d13C['NG Emissions: Stable Coal'] = merged_d13C.iloc[:, 1] - merged_d13C.iloc[11, 2] - merged_d13C.iloc[:, 3] - geo_emis # Stable at 2010
    merged_d13C['FER: Stable Coal'] = merged_d13C['NG Emissions: Stable Coal'] / merged_d13C['NG Prod'] * 100
    
    # Compile
    FER_dD_list.append(merged_dD['FER'].values)
    FER_d13C_list.append(merged_d13C['FER'].values)
    FER_dD_list_StableCoal.append(merged_dD['FER: Stable Coal'].values)
    FER_d13C_list_StableCoal.append(merged_d13C['FER: Stable Coal'].values)
    
# Calculate statistics and means
# Convert lists to 2D arrays: rows = years, columns = simulations
FER_dD_array = np.column_stack(FER_dD_list)
FER_d13C_array = np.column_stack(FER_d13C_list)
FER_dD_array_StableCoal = np.column_stack(FER_dD_list_StableCoal)
FER_d13C_array_StableCoal = np.column_stack(FER_d13C_list_StableCoal)

# Calculate ecah relative to 2010
FER_dD_arrayR = FER_dD_array - FER_dD_array[5,:]
FER_d13C_arrayR = FER_d13C_array - FER_d13C_array[11,:]
FER_dD_arrayR_StableCoal = FER_dD_array_StableCoal - FER_dD_array_StableCoal[5,:]
FER_d13C_arrayR_StableCoal = FER_d13C_array_StableCoal - FER_d13C_array_StableCoal[11,:]

# Calculate averages and standard deviations
# First total values
FER_dD_mean = np.mean(FER_dD_array, axis=1)
FER_dD_stdev =  np.std(FER_dD_array, axis=1)
FER_d13C_mean = np.mean(FER_d13C_array, axis=1)
FER_d13C_stdev =  np.std(FER_d13C_array, axis=1)
# Stable coal means
FER_dD_mean_StableCoal = np.mean(FER_dD_array_StableCoal, axis=1)
FER_d13C_mean_StableCoal = np.mean(FER_d13C_array_StableCoal, axis=1)

# Next relative
FER_dD_meanR = np.mean(FER_dD_arrayR, axis=1)
FER_dD_stdevR =  np.std(FER_dD_arrayR, axis=1)
FER_d13C_meanR = np.mean(FER_d13C_arrayR, axis=1)
FER_d13C_stdevR =  np.std(FER_d13C_arrayR, axis=1)
# Stable coal means
FER_dD_meanR_StableCoal =  np.mean(FER_dD_arrayR_StableCoal, axis=1)
FER_d13C_meanR_StableCoal =  np.mean(FER_d13C_arrayR_StableCoal, axis=1)

# Repeat for CarbonTracker Methane
merged_CT = pd.merge(df_CT, df_EDGAR, on='Year', how='outer')
merged_CT = pd.merge(merged_CT, df_EIA, on='Year', how='outer')
# Sort by Year to ensure the data is in chronological order
merged_CT = merged_CT.sort_values('Year')
merged_CT = merged_CT.dropna()

# Repeat for Schwietzke et al., 2016
merged_Schwietzke = pd.merge(df_Schwietzke, df_EDGAR, on='Year', how='outer')
merged_Schwietzke = pd.merge(merged_Schwietzke, df_EIA, on='Year', how='outer')
# Sort by Year to ensure the data is in chronological order
merged_Schwietzke = merged_Schwietzke.sort_values('Year')
merged_Schwietzke = merged_Schwietzke.dropna()

# Calculate FER using dD and d13C fossil emissions
merged_CT['NG Emissions'] = merged_CT.iloc[:, 1] - merged_CT.iloc[:, 2] - merged_CT.iloc[:, 3] - geo_emis
merged_CT['FER'] = merged_CT['NG Emissions'] / merged_CT['NG Prod'] * 100
merged_Schwietzke['NG Emissions'] = merged_Schwietzke.iloc[:, 1] - merged_Schwietzke.iloc[:, 2] - merged_Schwietzke.iloc[:, 3] - geo_emis
merged_Schwietzke['FER'] = merged_Schwietzke['NG Emissions'] / merged_Schwietzke['NG Prod'] * 100

#%% Plot it

# Create a subplot with two rows and one column (2x1 layout)
fig, ax = plt.subplots(2, 1, figsize=(7, 12), dpi=300)  
#fig.subplots_adjust(hspace=0.4)  # Adjust space between subplots
fig.set_facecolor('w')

# First subplot: Absolute FER values
# plot existing studies
ax[0].plot(CH4int_year, CH4int_Mcdonald, color='black', marker='*', linestyle='-', label='US (McDonald, 2023)')
ax[0].plot(CH4int_year, CH4int_Lu, color='cyan', marker='*', linestyle='-', label='US (Lu, 2023)')
# plot d13C
ax[0].plot(merged_d13C['Year'], FER_d13C_mean, color='firebrick', label=r"${\delta}^{13}\mathrm{C}$ derived", linewidth = 3)
ax[0].plot(merged_d13C['Year'], FER_d13C_mean_StableCoal, color='firebrick', label=r"${\delta}^{13}\mathrm{C}$-derived, stable coal", linestyle='--', linewidth = 2)
ax[0].fill_between(merged_d13C['Year'], (FER_d13C_mean + 2*FER_d13C_stdev), (FER_d13C_mean - 2*FER_d13C_stdev), color='firebrick', alpha=0.3)
# plot dD
ax[0].plot(merged_dD['Year'], FER_dD_mean, color='blue', label=r"${\delta}D$-derived", linewidth = 3)
ax[0].plot(merged_dD['Year'], FER_dD_mean_StableCoal, color='blue', label=r"${\delta}D$-derived, stable coal", linestyle='--', linewidth = 2)
ax[0].fill_between(merged_dD['Year'], (FER_dD_mean + 2*FER_dD_stdev), (FER_dD_mean - 2*FER_dD_stdev), color='blue', alpha=0.3)
ax[0].set_ylabel(r'CH$_4$ Intensity (%)', fontsize=15)
ax[0].grid(False)

# Second subplot: FER Difference from Baseline
ax[1].plot(CH4int_year, CH4int_Mcdonald - CH4int_Mcdonald[0], color='black', marker='*', linestyle='-', label='US (McDonald, 2023)')
ax[1].plot(CH4int_year, CH4int_Lu - CH4int_Lu[0], color='cyan', marker='*', linestyle='-', label='US (Lu, 2023)')
# plot d13C
ax[1].plot(merged_d13C['Year'], FER_d13C_meanR, color='firebrick', label=r"${\delta}^{13}\mathrm{C}$-Derived", linewidth = 3)
ax[1].plot(merged_d13C['Year'], FER_d13C_meanR_StableCoal, color='firebrick', label=r"${\delta}^{13}\mathrm{C}$-Derived: Stable Coal", linestyle='--', linewidth = 2)
ax[1].fill_between(merged_d13C['Year'],  (FER_d13C_meanR + 2*FER_d13C_stdevR), (FER_d13C_meanR - 2*FER_d13C_stdevR), color='firebrick', alpha=0.3)
# plot dD
ax[1].plot(merged_dD['Year'], FER_dD_meanR, color='blue', label=r"${\delta}D$-Derived", linewidth = 3)
ax[1].plot(merged_dD['Year'], FER_dD_meanR_StableCoal, color='blue', label=r"${\delta}D$-Derived: Stable Coal", linestyle='--', linewidth = 2)
ax[1].fill_between(merged_dD['Year'], (FER_dD_meanR + 2*FER_dD_stdevR), (FER_dD_meanR - 2*FER_dD_stdevR), color='blue', alpha=0.3)
ax[1].axhline(0, color='black', linestyle=':', linewidth=1)
ax[1].set_ylabel(r'CH$_4$ Intensity (%) Relative to 2010', fontsize=15)
ax[1].legend(fontsize=11, loc='upper right')
ax[1].grid(False)


#%% PLot figure 4
import matplotlib.pyplot as plt

# set font size ~9 pt for PNAS small figure
plt.rcParams.update({'font.size': 10})

def cm2inch(x): return x / 2.54

# small format: 9 × 6 cm
fig, axs = plt.subplots(2, 1, figsize=(cm2inch(18), cm2inch(22)), dpi=200, sharex=False)
fig.set_facecolor('w')

# unpack axes
axA, axB = axs

# --- subplot A: absolute FER values ---
axA.plot(CH4int_year, CH4int_Mcdonald, color='black', marker='*', linestyle='-', label='US (McDonald, 2023)')
axA.plot(CH4int_year, CH4int_Lu, color='cyan', marker='*', linestyle='-', label='US (Lu, 2023)')

# δ13C
axA.plot(merged_d13C['Year'], FER_d13C_mean, color='firebrick', label=r"${\delta}^{13}\mathrm{C}$-derived", linewidth=2)
axA.plot(merged_d13C['Year'], FER_d13C_mean_StableCoal, color='firebrick',
         linestyle='--', linewidth=2, label=r"${\delta}^{13}\mathrm{C}$-derived, stable coal")
axA.fill_between(merged_d13C['Year'],
                 FER_d13C_mean + 2*FER_d13C_stdev,
                 FER_d13C_mean - 2*FER_d13C_stdev,
                 color='firebrick', alpha=0.3)

# δD
axA.plot(merged_dD['Year'], FER_dD_mean, color='blue', label=r"${\delta}D$-derived", linewidth=2)
axA.plot(merged_dD['Year'], FER_dD_mean_StableCoal, color='blue',
         linestyle='--', linewidth=2, label=r"${\delta}D$-derived, stable coal")
axA.fill_between(merged_dD['Year'],
                 FER_dD_mean + 2*FER_dD_stdev,
                 FER_dD_mean - 2*FER_dD_stdev,
                 color='blue', alpha=0.3)

axA.set_ylabel(r'CH$_4$ Intensity (%)')
axA.grid(False)

# --- subplot B: FER difference from baseline ---
axB.plot(CH4int_year, CH4int_Mcdonald - CH4int_Mcdonald[0],
         color='black', marker='*', linestyle='-', label='US (McDonald, 2023)')
axB.plot(CH4int_year, CH4int_Lu - CH4int_Lu[0],
         color='cyan', marker='*', linestyle='-', label='US (Lu, 2023)')

# δ13C
axB.plot(merged_d13C['Year'], FER_d13C_meanR, color='firebrick',
         label=r"${\delta}^{13}\mathrm{C}$-derived", linewidth=2)
axB.plot(merged_d13C['Year'], FER_d13C_meanR_StableCoal, color='firebrick',
         linestyle='--', linewidth=2, label=r"${\delta}^{13}\mathrm{C}$-derived, stable coal")
axB.fill_between(merged_d13C['Year'],
                 FER_d13C_meanR + 2*FER_d13C_stdevR,
                 FER_d13C_meanR - 2*FER_d13C_stdevR,
                 color='firebrick', alpha=0.3)

# δD
axB.plot(merged_dD['Year'], FER_dD_meanR, color='blue',
         label=r"${\delta}D$-derived", linewidth=2)
axB.plot(merged_dD['Year'], FER_dD_meanR_StableCoal, color='blue',
         linestyle='--', linewidth=2, label=r"${\delta}D$-derived, stable coal")
axB.fill_between(merged_dD['Year'],
                 FER_dD_meanR + 2*FER_dD_stdevR,
                 FER_dD_meanR - 2*FER_dD_stdevR,
                 color='blue', alpha=0.3)

axB.axhline(0, color='black', linestyle=':', linewidth=1)
axB.set_ylabel(r'CH$_4$ Intensity (%) Rel. to 2010')
axB.legend(fontsize=10, loc='upper right')
axB.grid(False)

# panel labels
axA.text(0.02, 0.95, 'A', transform=axA.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='left')
axB.text(0.02, 0.95, 'B', transform=axB.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='left')

plt.tight_layout()
plt.savefig("Output/Fig4.pdf", dpi=200, bbox_inches="tight")
plt.show()
plt.close()

