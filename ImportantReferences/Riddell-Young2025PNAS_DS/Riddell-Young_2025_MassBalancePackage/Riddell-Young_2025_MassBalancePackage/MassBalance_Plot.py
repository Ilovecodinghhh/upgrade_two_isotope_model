#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:26:56 2024

@author: ryoung
"""

# This code plots the results of the dD and d13C mass balance together
# Code for Monte Carlo analysis of onebox deltaD-CH4 mass balance
locals().clear()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%% Load all ensemenle model outputs


# Load dD first

# Load EDGAR dD mass balance with Umezawa calibration and no buds (BASELINE SCENARIO)
# Load dD absolute emissions
df_absolute = pd.read_csv('Output/Results_dD-MassBalance_UmezawaCal_noBUDS.csv')
# Separate out each variable from the absolute emissions dataset
year_dD = df_absolute['Year']
FF_mean_dD = df_absolute['FF_mean']
FF_std_dD = df_absolute['FF_std']
Mic_mean_dD = df_absolute['Mic_mean']
Mic_std_dD = df_absolute['Mic_std']
# Load dD relative emissions
df_relative = pd.read_csv('Output/Results_RdD-MassBalance_UmezawaCal_noBUDS.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD = df_relative['FF_meanR']
FF_stdR_dD = df_relative['FF_stdR']
Mic_meanR_dD = df_relative['Mic_meanR']
Mic_stdR_dD = df_relative['Mic_stdR']

# Load EDGAR dD mass balance with Dasgupta calibration and no buds
# Load dD absolute emissions
df_absolute_Dasgupta = pd.read_csv('Output/Results_dD-MassBalance_DasguptaCal_noBUDS.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_Dasgupta = df_absolute_Dasgupta['Year']
FF_mean_dD_Dasgupta = df_absolute_Dasgupta['FF_mean']
FF_std_dD_Dasgupta = df_absolute_Dasgupta['FF_std']
Mic_mean_dD_Dasgupta = df_absolute_Dasgupta['Mic_mean']
Mic_std_dD_Dasgupta = df_absolute_Dasgupta['Mic_std']
# Load dD relative emissions
df_relative_Dasgupta = pd.read_csv('Output/Results_RdD-MassBalance_DasguptaCal_noBUDS.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_Dasgupta = df_relative_Dasgupta['FF_meanR']
FF_stdR_dD_Dasgupta = df_relative_Dasgupta['FF_stdR']
Mic_meanR_dD_Dasgupta = df_relative_Dasgupta['Mic_meanR']
Mic_stdR_dD_Dasgupta = df_relative_Dasgupta['Mic_stdR']

# Load EDGAR dD mass balance with Umezawa  calibration and stable mic and ff src sigs
# Load dD absolute emissions
df_absolute_StableMicFF = pd.read_csv('Output/Results_dD-MassBalance_StableMicFF.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_StableMicFF = df_absolute_StableMicFF['Year']
FF_mean_dD_StableMicFF = df_absolute_StableMicFF['FF_mean']
FF_std_dD_StableMicFF = df_absolute_StableMicFF['FF_std']
Mic_mean_dD_StableMicFF = df_absolute_StableMicFF['Mic_mean']
Mic_std_dD_StableMicFF = df_absolute_StableMicFF['Mic_std']
# Load dD relative emissions
df_relative_StableMicFF = pd.read_csv('Output/Results_RdD-MassBalance_StableMicFF.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_StableMicFF = df_relative_StableMicFF['FF_meanR']
FF_stdR_dD_StableMicFF = df_relative_StableMicFF['FF_stdR']
Mic_meanR_dD_StableMicFF = df_relative_StableMicFF['Mic_meanR']
Mic_stdR_dD_StableMicFF = df_relative_StableMicFF['Mic_stdR']

# Load EDGAR dD mass balance with Umezawa calibration and stable mic src sigs
# Load dD absolute emissions
df_absolute_StableMic = pd.read_csv('Output/Results_dD-MassBalance_StableMic.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_StableMic = df_absolute_StableMic['Year']
FF_mean_dD_StableMic = df_absolute_StableMic['FF_mean']
FF_std_dD_StableMic = df_absolute_StableMic['FF_std']
Mic_mean_dD_StableMic = df_absolute_StableMic['Mic_mean']
Mic_std_dD_StableMic = df_absolute_StableMic['Mic_std']
# Load dD relative emissions
df_relative_StableMic = pd.read_csv('Output/Results_RdD-MassBalance_StableMic.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_StableMic = df_relative_StableMic['FF_meanR']
FF_stdR_dD_StableMic = df_relative_StableMic['FF_stdR']
Mic_meanR_dD_StableMic = df_relative_StableMic['Mic_meanR']
Mic_stdR_dD_StableMic = df_relative_StableMic['Mic_stdR']
 
# Load EDGAR dD mass balance with Umezawa calibration
# Load dD absolute emissions
df_absolute_StableFF = pd.read_csv('Output/Results_dD-MassBalance_StableFF.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_StableFF = df_absolute_StableFF['Year']
FF_mean_dD_StableFF = df_absolute_StableFF['FF_mean']
FF_std_dD_StableFF = df_absolute_StableFF['FF_std']
Mic_mean_dD_StableFF = df_absolute_StableFF['Mic_mean']
Mic_std_dD_StableFF = df_absolute_StableFF['Mic_std']
# Load dD relative emissions
df_relative_StableFF = pd.read_csv('Output/Results_RdD-MassBalance_StableFF.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_StableFF = df_relative_StableFF['FF_meanR']
FF_stdR_dD_StableFF = df_relative_StableFF['FF_stdR']
Mic_meanR_dD_StableFF = df_relative_StableFF['Mic_meanR']
Mic_stdR_dD_StableFF = df_relative_StableFF['Mic_stdR']

# Load CTCH4 dD mass balance
# Load dD absolute emissions
df_absolute_CTCH4 = pd.read_csv('Output/Results_dD-MassBalance_FF-CTCH4.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_CTCH4 = df_absolute_CTCH4['Year']
FF_mean_dD_CTCH4 = df_absolute_CTCH4['FF_mean']
FF_std_dD_CTCH4 = df_absolute_CTCH4['FF_std']
Mic_mean_dD_CTCH4 = df_absolute_CTCH4['Mic_mean']
Mic_std_dD_CTCH4 = df_absolute_CTCH4['Mic_std']
# Load dD relative emissions
df_relative_CTCH4 = pd.read_csv('Output/Results_RdD-MassBalance_FF-CTCH4.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_CTCH4 = df_relative_CTCH4['FF_meanR']
FF_stdR_dD_CTCH4 = df_relative_CTCH4['FF_stdR']
Mic_meanR_dD_CTCH4 = df_relative_CTCH4['Mic_meanR']
Mic_stdR_dD_CTCH4 = df_relative_CTCH4['Mic_stdR']


# Now load d13C results

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
Mic_meanR_d13C = df_relative['Mic_meanR']
Mic_stdR_d13C = df_relative['Mic_stdR']

# Load EDGAR d13C mass balance with stable mic and ff src sigs
# Load d13C absolute emissions
df_absolute_StableMicFF = pd.read_csv('Output/Results_d13C-MassBalance_StableMicFF.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_StableMicFF = df_absolute_StableMicFF['Year']
FF_mean_d13C_StableMicFF = df_absolute_StableMicFF['FF_mean']
FF_std_d13C_StableMicFF = df_absolute_StableMicFF['FF_std']
Mic_mean_d13C_StableMicFF = df_absolute_StableMicFF['Mic_mean']
Mic_std_d13C_StableMicFF = df_absolute_StableMicFF['Mic_std']
# Load d13C relative emissions
df_relative_StableMicFF = pd.read_csv('Output/Results_Rd13C-MassBalance_StableMicFF.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_StableMicFF = df_relative_StableMicFF['FF_meanR']
FF_stdR_d13C_StableMicFF = df_relative_StableMicFF['FF_stdR']
Mic_meanR_d13C_StableMicFF = df_relative_StableMicFF['Mic_meanR']
Mic_stdR_d13C_StableMicFF = df_relative_StableMicFF['Mic_stdR']

# Load EDGAR d13C mass balance with stable mic src sigs
# Load d13C absolute emissions
df_absolute_StableMic = pd.read_csv('Output/Results_d13C-MassBalance_StableMic.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_StableMic = df_absolute_StableMic['Year']
FF_mean_d13C_StableMic = df_absolute_StableMic['FF_mean']
FF_std_d13C_StableMic = df_absolute_StableMic['FF_std']
Mic_mean_d13C_StableMic = df_absolute_StableMic['Mic_mean']
Mic_std_d13C_StableMic = df_absolute_StableMic['Mic_std']
# Load d13C relative emissions
df_relative_StableMic = pd.read_csv('Output/Results_Rd13C-MassBalance_StableMic.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_StableMic = df_relative_StableMic['FF_meanR']
FF_stdR_d13C_StableMic = df_relative_StableMic['FF_stdR']
Mic_meanR_d13C_StableMic = df_relative_StableMic['Mic_meanR']
Mic_stdR_d13C_StableMic = df_relative_StableMic['Mic_stdR']

# Load EDGAR d13C mass balance with stable ff src sigs
# Load d13C absolute emissions
df_absolute_StableFF = pd.read_csv('Output/Results_d13C-MassBalance_StableFF.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_StableFF = df_absolute_StableFF['Year']
FF_mean_d13C_StableFF = df_absolute_StableFF['FF_mean']
FF_std_d13C_StableFF = df_absolute_StableFF['FF_std']
Mic_mean_d13C_StableFF = df_absolute_StableFF['Mic_mean']
Mic_std_d13C_StableFF = df_absolute_StableFF['Mic_std']
# Load d13C relative emissions
df_relative_StableFF = pd.read_csv('Output/Results_Rd13C-MassBalance_StableFF.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_StableFF = df_relative_StableFF['FF_meanR']
FF_stdR_d13C_StableFF = df_relative_StableFF['FF_stdR']
Mic_meanR_d13C_StableFF = df_relative_StableFF['Mic_meanR']
Mic_stdR_d13C_StableFF = df_relative_StableFF['Mic_stdR']

# Load FF CTCH4 mass balance
# Load d13C absolute emissions
df_absolute_CTCH4 = pd.read_csv('Output/Results_d13C-MassBalance_FF-CTCH4.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_CTCH4 = df_absolute_CTCH4['Year']
FF_mean_d13C_CTCH4 = df_absolute_CTCH4['FF_mean']
FF_std_d13C_CTCH4 = df_absolute_CTCH4['FF_std']
Mic_mean_d13C_CTCH4 = df_absolute_CTCH4['Mic_mean']
Mic_std_d13C_CTCH4 = df_absolute_CTCH4['Mic_std']
# Load d13C relative emissions
df_relative_CTCH4 = pd.read_csv('Output/Results_Rd13C-MassBalance_FF-CTCH4.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_CTCH4 = df_relative_CTCH4['FF_meanR']
FF_stdR_d13C_CTCH4 = df_relative_CTCH4['FF_stdR']
Mic_meanR_d13C_CTCH4 = df_relative_CTCH4['Mic_meanR']
Mic_stdR_d13C_CTCH4 = df_relative_CTCH4['Mic_stdR']

# Load EDGAR Saueressig Mass balance
# Load d13C absolute emissions
df_absolute_Saueressig = pd.read_csv('Output/Results_d13C-MassBalance_Saueressig.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_Saueressig = df_absolute_Saueressig['Year']
FF_mean_d13C_Saueressig = df_absolute_Saueressig['FF_mean']
FF_std_d13C_Saueressig = df_absolute_Saueressig['FF_std']
Mic_mean_d13C_Saueressig = df_absolute_Saueressig['Mic_mean']
Mic_std_d13C_Saueressig = df_absolute_Saueressig['Mic_std']
# Load d13C relative emissions
df_relative_Saueressig = pd.read_csv('Output/Results_Rd13C-MassBalance_Saueressig.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_Saueressig = df_relative_Saueressig['FF_meanR']
FF_stdR_d13C_Saueressig = df_relative_Saueressig['FF_stdR']
Mic_meanR_d13C_Saueressig = df_relative_Saueressig['Mic_meanR']
Mic_stdR_d13C_Saueressig = df_relative_Saueressig['Mic_stdR']

# Load Carbon tracker methane
data2 = pd.read_excel('data/CarbonTracker_CH4.xlsx')
totalCT = data2.iloc[:,15].values
yearss = data2.iloc[:,0].values
micCT = data2.iloc[:,7].values
ffCT = data2.iloc[:,3].values
bbCT = data2.iloc[:,11].values

# Define Geologic emission strength
geo = 5 


#%% Plot figure 2

import matplotlib.pyplot as plt
import numpy as np

# set font size ~10 pt (PNAS guideline 6–12 pt)
plt.rcParams.update({'font.size': 10})

# helper to convert cm → inch
def cm2inch(x): return x / 2.54

# create a 2x2 subplot grid, 11x11 cm
fig, axs = plt.subplots(2, 2, figsize=(cm2inch(18), cm2inch(22)), dpi=200)  
axs = axs.flatten()

# panel A: fossil emissions
axs[0].fill_between(year_dD, FF_mean_dD + 2*FF_std_dD - geo,
                    FF_mean_dD - 2*FF_std_dD - geo, color='blue', alpha=0.3)
axs[0].plot(year_dD, FF_mean_dD - geo, color='blue', linewidth=1.5)
axs[0].plot(year_d13C, FF_mean_d13C - geo, color='firebrick', linewidth=1.5)
axs[0].fill_between(year_d13C, FF_mean_d13C + 2*FF_std_d13C - geo,
                    FF_mean_d13C - 2*FF_std_d13C - geo, color='firebrick', alpha=0.3)
axs[0].set_ylabel('FF Emissions (Tg yr$^{-1}$)')
axs[0].set_title('FF Emissions',fontsize=12)
axs[0].set_ylim(40, 300)

# panel B: microbial emissions
axs[1].fill_between(year_dD, Mic_mean_dD + 2*Mic_std_dD,
                    Mic_mean_dD - 2*Mic_std_dD, color='blue', alpha=0.3)
axs[1].plot(year_dD, Mic_mean_dD, color='blue', linewidth=1.5, label='${\delta}D$ Mass Balance')
axs[1].plot(year_d13C, Mic_mean_d13C, color='firebrick', linewidth=1.5, label=r"${\delta}^{13}C$ Mass Balance")
axs[1].fill_between(year_d13C, Mic_mean_d13C + 2*Mic_std_d13C,
                    Mic_mean_d13C - 2*Mic_std_d13C, color='firebrick', alpha=0.3)
axs[1].legend(loc='upper center', fontsize=10)
axs[1].set_ylabel('Mic Emissions (Tg yr$^{-1}$)')
axs[1].set_title('Mic Emissions',fontsize=12)
axs[1].set_ylim(300, 560)

# panel C: relative fossil emissions
axs[2].fill_between(year_dD, FF_meanR_dD + 2*FF_stdR_dD,
                    FF_meanR_dD - 2*FF_stdR_dD, color='blue', alpha=0.3)
axs[2].plot(year_dD, FF_meanR_dD, color='blue', linewidth=1.5)
axs[2].plot(year_d13C, FF_meanR_d13C, color='firebrick', linewidth=1.5)
axs[2].fill_between(year_d13C, FF_meanR_d13C + 2*FF_stdR_d13C,
                    FF_meanR_d13C - 2*FF_stdR_d13C, color='firebrick', alpha=0.3)
axs[2].axhline(0, color='black', linestyle='--', linewidth=0.8)
axs[2].yaxis.set_major_locator(plt.MultipleLocator(20))
axs[2].grid(axis='y', linestyle='--', alpha=0.6)
axs[2].set_ylabel('Δ FF Emissions (Tg yr$^{-1}$)')
axs[2].set_title('Relative FF Emissions',fontsize=12)
axs[2].set_ylim(-40, 120)

# panel D: relative microbial emissions
axs[3].fill_between(year_dD, Mic_meanR_dD + 2*Mic_stdR_dD,
                    Mic_meanR_dD - 2*Mic_stdR_dD, color='blue', alpha=0.3)
axs[3].plot(year_dD, Mic_meanR_dD, color='blue', linewidth=1.5, label='${\delta}D$ Mass Balance')
axs[3].plot(year_d13C, Mic_meanR_d13C, color='firebrick', linewidth=1.5, label=r"${\delta}^{13}C$ Mass Balance")
axs[3].fill_between(year_d13C, Mic_meanR_d13C + 2*Mic_stdR_d13C,
                    Mic_meanR_d13C - 2*Mic_stdR_d13C, color='firebrick', alpha=0.3)
axs[3].axhline(0, color='black', linestyle='--', linewidth=0.8)
axs[3].yaxis.set_major_locator(plt.MultipleLocator(20))
axs[3].grid(axis='y', linestyle='--', alpha=0.6)
axs[3].set_ylabel('Δ Mic Emissions (Tg yr$^{-1}$)')
axs[3].set_title('Relative Mic Emissions',fontsize=12)
axs[3].set_ylim(-40, 120)

# panel labels
axs[0].text(0.02, 0.98, 'A', transform=axs[0].transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')
axs[1].text(0.02, 0.98, 'B', transform=axs[1].transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')
axs[2].text(0.02, 0.98, 'C', transform=axs[2].transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')
axs[3].text(0.02, 0.98, 'D', transform=axs[3].transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')

plt.tight_layout()

# adjust subplot area inside fixed canvas
plt.subplots_adjust(top=0.94, bottom=0.08)  
# top < 1.0 leaves space at the top
# bottom > 0.0 leaves space at the bottom

# save WITHOUT bbox_inches="tight" so canvas stays 18 × 22 cm
plt.savefig("Output/fig2.pdf", dpi=200)

plt.show()

# # tidy layout
# plt.tight_layout()

# # save as pdf (PNAS-compliant)
# plt.savefig("Output/fig2.pdf", dpi=200, bbox_inches="none")

# # show in spyder plots pane
# plt.show()


#%% Plot the data and sensitivity tests

# Plot combined net and relative results for Fossil and Microbial emissions
fig, axs = plt.subplots(2, 2, figsize=(12, 9), dpi=500, sharex=True)
axs = axs.flatten()
plt.rcParams.update({'font.size': 14})

# Fossil - Net
axs[0].fill_between(year_dD, FF_mean_dD + 2 * FF_std_dD - geo, FF_mean_dD - 2 * FF_std_dD - geo, color='blue', alpha=0.3)
axs[0].plot(year_dD, FF_mean_dD_StableMic - geo, color='indigo', linewidth=2, label='Stable Mic')
axs[0].plot(year_dD, FF_mean_dD_StableFF - geo, color='deeppink', linewidth=2, label='Stable FF')
axs[0].plot(year_dD_CTCH4, FF_mean_dD_CTCH4 - geo, color='blue', linewidth=2, label='dD-CTCH4', linestyle='dashed')
axs[0].plot(year_dD_Dasgupta, FF_mean_dD_Dasgupta - geo, color='blue', linewidth=2, linestyle='dotted')
axs[0].plot(year_dD, FF_mean_dD - geo, color='blue', linewidth=3)
axs[0].plot(year_d13C, FF_mean_d13C - geo, color='firebrick', linewidth=3)
axs[0].fill_between(year_d13C, FF_mean_d13C + 2 * FF_std_d13C - geo, FF_mean_d13C - 2 * FF_std_d13C - geo, color='firebrick', alpha=0.3)
axs[0].plot(year_d13C, FF_mean_d13C_StableMic - geo, color='grey', linewidth=2, label='Stable Mic')
axs[0].plot(year_d13C, FF_mean_d13C_StableFF - geo, color='navy', linewidth=2, label='Stable FF')
axs[0].plot(year_d13C_CTCH4[1:], FF_mean_d13C_CTCH4[1:] - geo, color='firebrick', linewidth=2, linestyle='dashed')
axs[0].set_ylabel('FF Emissions (Tg yr$^{-1}$)')
axs[0].set_title('Fossil Emissions')
axs[0].set_ylim(40, 290)

# Fossil - Relative
axs[2].fill_between(year_dD, FF_meanR_dD + 2 * FF_stdR_dD, FF_meanR_dD - 2 * FF_stdR_dD, color='blue', alpha=0.3)
axs[2].plot(year_dD, FF_meanR_dD_StableMic, color='indigo', linewidth=2, label='Stable Mic')
axs[2].plot(year_dD, FF_meanR_dD_StableFF, color='deeppink', linewidth=2, label='Stable FF')
axs[2].plot(year_dD_CTCH4, FF_meanR_dD_CTCH4, color='blue', linewidth=2, linestyle='dashed')
axs[2].plot(year_dD_Dasgupta, FF_meanR_dD_Dasgupta, color='blue', linewidth=2, linestyle='dotted')
axs[2].plot(year_dD, FF_meanR_dD, color='blue', linewidth=3)
axs[2].plot(year_d13C, FF_meanR_d13C, color='firebrick', linewidth=3)
axs[2].fill_between(year_d13C, FF_meanR_d13C + 2 * FF_stdR_d13C, FF_meanR_d13C - 2 * FF_stdR_d13C, color='firebrick', alpha=0.3)
axs[2].plot(year_d13C, FF_meanR_d13C_StableMic, color='grey', linewidth=2, label='Stable Mic')
axs[2].plot(year_d13C, FF_meanR_d13C_StableFF, color='navy', linewidth=2, label='Stable FF')
axs[2].plot(year_d13C_CTCH4[1:], FF_meanR_d13C_CTCH4[1:], color='firebrick', linewidth=2, linestyle='dashed')
axs[2].axhline(0, color='black', linestyle='--')
axs[2].yaxis.set_major_locator(plt.MultipleLocator(10))
axs[2].grid(axis='y', linestyle='--', alpha=0.6)
axs[2].set_ylabel('Δ FF Emissions (Tg yr$^{-1}$)')
axs[2].set_title('Relative Fossil Emissions')
axs[2].set_ylim(-40, 120)

# Microbial - Net
axs[1].fill_between(year_dD, Mic_mean_dD + 2 * Mic_std_dD, Mic_mean_dD - 2 * Mic_std_dD, color='blue', alpha=0.3)
axs[1].plot(year_dD, Mic_mean_dD_StableMic, color='indigo', linewidth=2, label='${\delta}D$ Stable Mic')
axs[1].plot(year_dD, Mic_mean_dD_StableFF, color='deeppink', linewidth=2, label='${\delta}D$ Stable FF')
axs[1].plot(year_dD_CTCH4, Mic_mean_dD_CTCH4, color='blue', linewidth=2, label='dD-CTCH4', linestyle='dashed')
axs[1].plot(year_dD_Dasgupta, Mic_mean_dD_Dasgupta, color='blue', linewidth=2, linestyle='dotted', label='${\delta}D$ Dasgupta Mass Balance')
axs[1].plot(year_dD, Mic_mean_dD, color='blue', linewidth=3, label='${\delta}D$ Mass Balance')
axs[1].plot(year_d13C, Mic_mean_d13C, color='firebrick', linewidth=3, label=r"${\delta}^{13}C$ Mass Balance")
axs[1].fill_between(year_d13C, Mic_mean_d13C + 2 * Mic_std_d13C, Mic_mean_d13C - 2 * Mic_std_d13C, color='firebrick', alpha=0.3)
axs[1].plot(year_d13C, Mic_mean_d13C_StableMic, color='grey', linewidth=2, label=r"${\delta}^{13}C$ Stable Mic")
axs[1].plot(year_d13C, Mic_mean_d13C_StableFF, color='navy', linewidth=2, label=r"${\delta}^{13}C$ Stable FF")
axs[1].plot(year_d13C_CTCH4[1:], Mic_mean_d13C_CTCH4[1:], color='firebrick', linewidth=2, linestyle='dashed')
axs[1].set_ylabel('Mic Emissions (Tg yr$^{-1}$)')
axs[1].set_title('Microbial Emissions')
axs[1].set_ylim(310, 560)

# Microbial - Relative
axs[3].fill_between(year_dD, Mic_meanR_dD + 2 * Mic_stdR_dD, Mic_meanR_dD - 2 * Mic_stdR_dD, color='blue', alpha=0.3)
axs[3].plot(year_dD, Mic_meanR_dD_StableMic, color='indigo', linewidth=2, label='${\delta}D$ Stable Mic')
axs[3].plot(year_dD, Mic_meanR_dD_StableFF, color='deeppink', linewidth=2, label='${\delta}D$ Stable FF')
axs[3].plot(year_dD_CTCH4, Mic_meanR_dD_CTCH4, color='blue', linewidth=2, linestyle='dashed', label='${\delta}D$ w/FF CTCH4')
axs[3].plot(year_dD_Dasgupta, Mic_meanR_dD_Dasgupta, color='blue', linewidth=2, linestyle='dotted', label='${\delta}D$ w/Dasgupta Offsets')
axs[3].plot(year_dD, Mic_meanR_dD, color='blue', linewidth=3, label='${\delta}D$ Mass Balance')
axs[3].plot(year_d13C, Mic_meanR_d13C, color='firebrick', linewidth=3, label=r"${\delta}^{13}C$ Mass Balance")
axs[3].fill_between(year_d13C, Mic_meanR_d13C + 2 * Mic_stdR_d13C, Mic_meanR_d13C - 2 * Mic_stdR_d13C, color='firebrick', alpha=0.3)
axs[3].plot(year_d13C, Mic_meanR_d13C_StableMic, color='grey', linewidth=2, label=r"${\delta}^{13}C$ Stable Mic")
axs[3].plot(year_d13C, Mic_meanR_d13C_StableFF, color='navy', linewidth=2, label=r"${\delta}^{13}C$ Stable FF")
axs[3].plot(year_d13C_CTCH4[1:], Mic_meanR_d13C_CTCH4[1:], color='firebrick', linewidth=2, linestyle='dashed', label=r"${\delta}^{13}C$ w/FF CT-CH4")
axs[3].axhline(0, color='black', linestyle='--')
axs[3].set_ylabel('Δ Mic Emissions (Tg yr$^{-1}$)')
axs[3].legend(loc='upper left', ncol=2, fontsize=10)
axs[3].yaxis.set_major_locator(plt.MultipleLocator(10))
axs[3].grid(axis='y', linestyle='--', alpha=0.6)
axs[3].set_title('Relative Microbial Emissions')
axs[3].set_ylim(-40, 120)

plt.tight_layout()
plt.show()


#%% Histogram plot comparison to previous studies

# Load data
dD_histogram = pd.read_csv('Output/dD_histogram_UmezawaCal_noBUDS.csv')
Mic_Delta_MC_2020_dD = dD_histogram.iloc[:,2]
FF_Delta_MC_2020_dD = dD_histogram.iloc[:,3]
Mic_Delta_MC_total_dD = dD_histogram.iloc[:,0]
FF_Delta_MC_total_dD = dD_histogram.iloc[:,1]
d13C_histogram = pd.read_csv('Output/d13C_histogram_Cantrell.csv')
Mic_Delta_MC_2020_d13C = d13C_histogram.iloc[:,2]
FF_Delta_MC_2020_d13C = d13C_histogram.iloc[:,3]
Mic_Delta_MC_total_d13C = d13C_histogram.iloc[:,0]
FF_Delta_MC_total_d13C = d13C_histogram.iloc[:,1]
GCP_ensembles = pd.read_csv('data/GCP_budget_deltas.csv') # From GCP, Saunois 2020
GCP_FF = GCP_ensembles.iloc[:,1]
GCP_Mic = GCP_ensembles.iloc[:,2]

# Calculate carbontracker difference
# First, 2020 minus 2000 to 2009 average
Mic_Delta_CT_2020 = micCT[-2] - np.mean(micCT[3:13])
FF_Delta_CT_2020 = ffCT[-2] - np.mean(ffCT[3:13])
# Now, 2020 to 2021 average minus 2005 to 2007 average 
Mic_Delta_CT_total = np.mean(micCT[-2:]) - np.mean(micCT[8:11])
FF_Delta_CT_total = np.mean(ffCT[-2:]) - np.mean(ffCT[8:11])

bin_width = 3
mic_bins = np.arange(min(GCP_Mic), max(Mic_Delta_MC_2020_dD) + bin_width, bin_width)
ff_bins = np.arange(min(FF_Delta_MC_2020_dD), max(GCP_FF) + bin_width, bin_width)

fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(7, 10), dpi=300)

# microbial deltas
closest_bin_mic = min(mic_bins, key=lambda x: abs(x - Mic_Delta_CT_2020))
ax1.axvspan(closest_bin_mic, closest_bin_mic + bin_width, color='black', alpha=0.5, label="CT-CH4")
ax1.hist(Mic_Delta_MC_2020_dD, bins=mic_bins, color='maroon', edgecolor='black', alpha=0.7, label="${\delta}D$")
ax1.hist(Mic_Delta_MC_2020_d13C, bins=mic_bins, color='blue', edgecolor='black', alpha=0.7, label=r"$\delta^{13}C$")
ax1.axvline(0, color='black', linestyle='--', linewidth=1.5)
ax1.set_xlabel("2020, 2000-2009 Microbial Emission Change (Tg/yr)")
ax1.set_ylabel("Mass Balance MC Frequency", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(axis='y', linestyle='--', alpha=0.6)
ax1_twin = ax1.twinx()
ax1_twin.hist(GCP_Mic, bins=mic_bins, color='orange', edgecolor='black', alpha=0.7, label="GCP")
ax1_twin.set_ylabel("GCP Ensemble Frequency", color='orange')
ax1_twin.tick_params(axis='y', labelcolor='orange')
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax1_twin.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc='upper right')

# fossil fuel deltas
closest_bin_ff = min(ff_bins, key=lambda x: abs(x - FF_Delta_CT_2020))
ax3.axvspan(closest_bin_ff, closest_bin_ff + bin_width, color='black', alpha=0.5, label="CT-CH4")
ax3.hist(FF_Delta_MC_2020_d13C, bins=ff_bins, color='blue', edgecolor='black', alpha=0.7, label=r"$\delta^{13}C$")
ax3.hist(FF_Delta_MC_2020_dD, bins=ff_bins, color='maroon', edgecolor='black', alpha=0.7, label="${\delta}D$")
ax3.axvline(0, color='black', linestyle='--', linewidth=1.5)
ax3.set_xlabel("2020, 2000-2009 FF Emission Change (Tg/yr)")
ax3.set_ylabel("Mass Balance MC Frequency", color='black')
ax3.tick_params(axis='y', labelcolor='black')
ax3.grid(axis='y', linestyle='--', alpha=0.6)
ax3_twin = ax3.twinx()
ax3_twin.hist(GCP_FF, bins=ff_bins, color='orange', edgecolor='black', alpha=0.7, label="GCP")
ax3_twin.set_ylabel("GCP Ensemble Frequency", color='orange')
ax3_twin.tick_params(axis='y', labelcolor='orange')

plt.tight_layout()
plt.show()


#%% Plot sensitivity test results

# Load results first
# First dD
# Load OH increase sensitivity test
# Load dD absolute emissions
df_absolute_OH_inc = pd.read_csv('Output/Results_dD-MassBalance_OH_inc.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_OH_inc = df_absolute_OH_inc['Year']
FF_mean_dD_OH_inc = df_absolute_OH_inc['FF_mean']
FF_std_dD_OH_inc = df_absolute_OH_inc['FF_std']
Mic_mean_dD_OH_inc = df_absolute_OH_inc['Mic_mean']
Mic_std_dD_OH_inc = df_absolute_OH_inc['Mic_std']
# Load dD relative emissions
df_relative_OH_inc = pd.read_csv('Output/Results_RdD-MassBalance_OH_inc.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_OH_inc = df_relative_OH_inc['FF_meanR']
FF_stdR_dD_OH_inc = df_relative_OH_inc['FF_stdR']
Mic_meanR_dD_OH_inc = df_relative_OH_inc['Mic_meanR']
Mic_stdR_dD_OH_inc = df_relative_OH_inc['Mic_stdR']

# Load BB decrease sensitivity test
# Load dD absolute emissions
df_absolute_redBB = pd.read_csv('Output/Results_dD-MassBalance_redBB.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_redBB = df_absolute_redBB['Year']
FF_mean_dD_redBB = df_absolute_redBB['FF_mean']
FF_std_dD_redBB = df_absolute_redBB['FF_std']
Mic_mean_dD_redBB = df_absolute_redBB['Mic_mean']
Mic_std_dD_redBB = df_absolute_redBB['Mic_std']
# Load dD relative emissions
df_relative_redBB = pd.read_csv('Output/Results_RdD-MassBalance_redBB.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_redBB = df_relative_redBB['FF_meanR']
FF_stdR_dD_redBB = df_relative_redBB['FF_stdR']
Mic_meanR_dD_redBB = df_relative_redBB['Mic_meanR']
Mic_stdR_dD_redBB = df_relative_redBB['Mic_stdR']

# Load Cl decrease sensitivity test
# Load dD absolute emissions
df_absolute_redCl = pd.read_csv('Output/Results_dD-MassBalance_redCl.csv')
# Separate out each variable from the absolute emissions dataset
year_dD_redCl = df_absolute_redCl['Year']
FF_mean_dD_redCl = df_absolute_redCl['FF_mean']
FF_std_dD_redCl = df_absolute_redCl['FF_std']
Mic_mean_dD_redCl = df_absolute_redCl['Mic_mean']
Mic_std_dD_redCl = df_absolute_redCl['Mic_std']
# Load dD relative emissions
df_relative_redCl = pd.read_csv('Output/Results_RdD-MassBalance_redCl.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_dD_redCl = df_relative_redCl['FF_meanR']
FF_stdR_dD_redCl = df_relative_redCl['FF_stdR']
Mic_meanR_dD_redCl = df_relative_redCl['Mic_meanR']
Mic_stdR_dD_redCl = df_relative_redCl['Mic_stdR']


# # Next d13C
# Load OH increase sensitivity test
# Load d13C absolute emissions
df_absolute_OH_inc = pd.read_csv('Output/Results_d13C-MassBalance_OH_inc.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_OH_inc = df_absolute_OH_inc['Year']
FF_mean_d13C_OH_inc = df_absolute_OH_inc['FF_mean']
FF_std_d13C_OH_inc = df_absolute_OH_inc['FF_std']
Mic_mean_d13C_OH_inc = df_absolute_OH_inc['Mic_mean']
Mic_std_d13C_OH_inc = df_absolute_OH_inc['Mic_std']
# Load d13C relative emissions
df_relative_OH_inc = pd.read_csv('Output/Results_Rd13C-MassBalance_OH_inc.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_OH_inc = df_relative_OH_inc['FF_meanR']
FF_stdR_d13C_OH_inc = df_relative_OH_inc['FF_stdR']
Mic_meanR_d13C_OH_inc = df_relative_OH_inc['Mic_meanR']
Mic_stdR_d13C_OH_inc = df_relative_OH_inc['Mic_stdR']

# Load BB decrease sensitivity test
# Load d13C absolute emissions
df_absolute_redBB = pd.read_csv('Output/Results_d13C-MassBalance_redBB.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_redBB = df_absolute_redBB['Year']
FF_mean_d13C_redBB = df_absolute_redBB['FF_mean']
FF_std_d13C_redBB = df_absolute_redBB['FF_std']
Mic_mean_d13C_redBB = df_absolute_redBB['Mic_mean']
Mic_std_d13C_redBB = df_absolute_redBB['Mic_std']
# Load d13C relative emissions
df_relative_redBB = pd.read_csv('Output/Results_Rd13C-MassBalance_redBB.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_redBB = df_relative_redBB['FF_meanR']
FF_stdR_d13C_redBB = df_relative_redBB['FF_stdR']
Mic_meanR_d13C_redBB = df_relative_redBB['Mic_meanR']
Mic_stdR_d13C_redBB = df_relative_redBB['Mic_stdR']

# Load Cl decrease sensitivity test
# Load d13C absolute emissions
df_absolute_redCl = pd.read_csv('Output/Results_d13C-MassBalance_redCl.csv')
# Separate out each variable from the absolute emissions dataset
year_d13C_redCl = df_absolute_redCl['Year']
FF_mean_d13C_redCl = df_absolute_redCl['FF_mean']
FF_std_d13C_redCl = df_absolute_redCl['FF_std']
Mic_mean_d13C_redCl = df_absolute_redCl['Mic_mean']
Mic_std_d13C_redCl = df_absolute_redCl['Mic_std']
# Load d13C relative emissions
df_relative_redCl = pd.read_csv('Output/Results_Rd13C-MassBalance_redCl.csv')
# Separate out each variable from the relative emissions dataset
FF_meanR_d13C_redCl = df_relative_redCl['FF_meanR']
FF_stdR_d13C_redCl = df_relative_redCl['FF_stdR']
Mic_meanR_d13C_redCl = df_relative_redCl['Mic_meanR']
Mic_stdR_d13C_redCl = df_relative_redCl['Mic_stdR']


# Load FF dD source sigantures
FF_dD_data = pd.read_csv('Output/FF_dD_GlobUnc.csv', delimiter=',')
FF_dD_MC_CTCH4 = pd.read_csv('Output/FF_dD_GlobMC_CTCH4.csv', delimiter=',', header = None)
FF_dD_MC_EDGAR_data = pd.read_csv('Output/FF_dD_GlobMC_EDGAR.csv', delimiter=',')
FF_dD_MC_EDGAR = FF_dD_MC_EDGAR_data.iloc[34:,1:]
# Calculate EDGAR statistics
rFF_dD_MC_EDGAR = FF_dD_MC_EDGAR-FF_dD_MC_EDGAR.iloc[0,:]
rFF_dD_mean_EDGAR = rFF_dD_MC_EDGAR.iloc[:, 1:].mean(axis=1).to_numpy()
rFF_dD_std_EDGAR = rFF_dD_MC_EDGAR.iloc[:, 1:].std(axis=1).to_numpy()
# Calculate CTCH4 statistics
FF_dD_mean_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
FF_dD_mean_EDGAR = FF_dD_MC_EDGAR.iloc[:, 1:].mean(axis=1).to_numpy()
FF_dD_std_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()
R_FF_dD_mean_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
R_FF_dD_std_CTCH4 = FF_dD_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()
# Load FF d13C source sigantures
FF_d13C_data = pd.read_csv('Output/FF_d13C_GlobUnc.csv', delimiter=',')
FF_d13C_MC_CTCH4 = pd.read_csv('Output/FF_d13C_GlobMC_CTCH4.csv', delimiter=',', header = None)
FF_d13C_MC_EDGAR_data = pd.read_csv('Output/FF_d13C_GlobMC_EDGAR.csv', delimiter=',')
FF_d13C_MC_EDGAR = FF_d13C_MC_EDGAR_data.iloc[28:,1:] # starting in 1999
# Calculate EDGAR statistics
rFF_d13C_MC_EDGAR = FF_d13C_MC_EDGAR-FF_d13C_MC_EDGAR.iloc[0,:]
rFF_d13C_mean_EDGAR = rFF_d13C_MC_EDGAR.iloc[:, 1:].mean(axis=1).to_numpy()
rFF_d13C_std_EDGAR = rFF_d13C_MC_EDGAR.iloc[:, 1:].std(axis=1).to_numpy()
# Calculate CTCH4 statistics
FF_d13C_mean_EDGAR = FF_d13C_MC_EDGAR.iloc[:, 1:].mean(axis=1).to_numpy()
FF_d13C_std_EDGAR = FF_d13C_MC_EDGAR.iloc[:, 1:].std(axis=1).to_numpy()
FF_d13C_mean_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
FF_d13C_std_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()
R_FF_d13C_mean_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].mean(axis=1).to_numpy()
R_FF_d13C_std_CTCH4 = FF_d13C_MC_CTCH4.iloc[:, 1:].std(axis=1).to_numpy()


# First plot Cl and Biomass Burning sensitivity tests
# Plot relative results for Fossil
plt.rcParams.update({'font.size': 16})  # Adjust this value as needed
# Create a figure with 2 subplots
fig, axs = plt.subplots(2, 1, figsize=(8, 12), dpi=300)
axs[0].plot(year_dD, FF_meanR_dD, color='black', linewidth=3, label="${\delta}D$ Mass balance")
axs[0].plot(year_d13C, FF_meanR_d13C, color='black', linewidth=3, label=r"$\delta^{13}C$ Mass Balance", linestyle='dashed')
axs[0].fill_between(year_dD, FF_meanR_dD + 2 * FF_stdR_dD, FF_meanR_dD - 2 * FF_stdR_dD, color='black', alpha=0.3)
axs[0].plot(year_dD_redBB, FF_meanR_dD_redBB, color='chocolate', linewidth=3, label='BB decrease to 0')
axs[0].plot(year_dD_redCl, FF_meanR_dD_redCl, color='cyan', linewidth=3, label='reduced Cl Sink')
axs[0].plot(year_dD_OH_inc, FF_meanR_dD_OH_inc, color='magenta', linewidth=3, label='0.3%/yr OH increase')
axs[0].fill_between(year_d13C, FF_meanR_d13C + 2 * FF_stdR_d13C, FF_meanR_d13C - 2 * FF_stdR_d13C, color='black', alpha=0.3)
axs[0].plot(year_d13C_redBB, FF_meanR_d13C_redBB, color='chocolate', linewidth=3, linestyle='dashed')
axs[0].plot(year_d13C_redCl, FF_meanR_d13C_redCl, color='cyan', linewidth=3, linestyle='dashed')
axs[0].plot(year_d13C, FF_meanR_d13C_OH_inc, color='magenta', linewidth=3, linestyle='dashed')
axs[0].yaxis.set_major_locator(plt.MultipleLocator(10))  # Set grid interval to 10
axs[0].grid(axis='y', linestyle='--', alpha=0.6)
axs[0].set_ylabel('Δ FF Emissions (Tg yr$^{-1}$)')
axs[0].set_title('Relative FF Emissions')
axs[0].legend(loc='lower left', fontsize=10)

# Plot relative results for microbial
axs[1].plot(year_dD, Mic_meanR_dD, color='black', linewidth=3)
axs[1].fill_between(year_dD, Mic_meanR_dD + 2 * Mic_stdR_dD, Mic_meanR_dD - 2 * FF_stdR_dD, color='black', alpha=0.3)
axs[1].plot(year_d13C, Mic_meanR_d13C, color='black', linewidth=3, linestyle='dashed')
axs[1].plot(year_dD_redBB, Mic_meanR_dD_redBB, color='chocolate', linewidth=3)
axs[1].plot(year_dD_redBB, Mic_meanR_dD_redCl, color='cyan', linewidth=3)
axs[1].plot(year_dD_OH_inc, Mic_meanR_dD_OH_inc, color='magenta', linewidth=3)
axs[1].fill_between(year_d13C, Mic_meanR_d13C + 2 * Mic_stdR_d13C, Mic_meanR_d13C - 2 * Mic_stdR_d13C, color='black', alpha=0.3)
axs[1].plot(year_d13C_redBB, Mic_meanR_d13C_redBB, color='chocolate', linewidth=3, linestyle='dashed')
axs[1].plot(year_d13C_redCl, Mic_meanR_d13C_redCl, color='cyan', linewidth=3, linestyle='dashed')
axs[1].plot(year_d13C, Mic_meanR_d13C_OH_inc, color='magenta', linewidth=3, linestyle='dashed')
axs[1].yaxis.set_major_locator(plt.MultipleLocator(10))  # Set grid interval to 10
axs[1].grid(axis='y', linestyle='--', alpha=0.6)
axs[1].set_ylabel('Δ Mic Emissions (Tg yr$^{-1}$)')
axs[1].set_title('Relative Mic Emissions')

axs[0].text(0.01, 0.97, 'A', transform=axs[0].transAxes,
             fontsize=25, fontweight='bold', va='top', ha='left')

axs[1].text(0.01, 0.97, 'B', transform=axs[1].transAxes,
             fontsize=25, fontweight='bold', va='top', ha='left')



#%% Now plot fossil source signature results

# Load Thanwerdas fossil source signature
import xarray as xr

# Load the NetCDF file
ds = xr.open_dataset("data/Thanwerdas_2024_posterior.nc")

# Extract variables while preserving their dimensions
time = ds['time'].values
scenario = ds['scenario'].values
sign = ds['sign'].values
cat = ds['cat'].values
region = ds['region'].values

# Extract global mean FF signature starting in 2009
FF_d13C = sign[0, 2, -1, 1:]
AGW_d13C = sign[0, 0, -1, 1:]
WET_d13C = np.mean(sign, axis=0)[0, -1, 1:]
timeB = time[1:]
# Calculate relative to 2005
FF_d13C_rel = FF_d13C - FF_d13C[0]
AGW_d13C_rel = AGW_d13C - AGW_d13C[0]

# Load calculated AGW source signature using Chang et al.,2019 and 
AGW_Calc = pd.read_excel('data/AGW_d13C_Calc.xlsx', header=None)
AGW_Calc_years = AGW_Calc.iloc[:,0].values
AGW_Calc_d13C = AGW_Calc.iloc[:,1].values
AGW_Calc_d13C_rel = AGW_Calc_d13C - AGW_Calc_d13C[0]

# Plot comparison between Thanwerdas and our source signature trends

# Plot relative results for Fossil
plt.rcParams.update({'font.size': 16})  # Adjust this value as needed
# Create a figure with 2 subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
# Plot change in AGW d13C
axs[1, 0].plot(AGW_Calc_years, AGW_Calc_d13C_rel, color='blue', linewidth=3, label="This study")
axs[1, 0].plot(timeB, AGW_d13C_rel, color='firebrick', linewidth=3, linestyle = 'solid', label="Thanwerdas et al., 2024")
axs[1, 0].fill_between(AGW_Calc_years, AGW_Calc_d13C_rel + 0.5 * AGW_Calc_d13C_rel, AGW_Calc_d13C_rel - 0.5 * AGW_Calc_d13C_rel, color='blue', alpha=0.3)
axs[1, 0].grid(axis='y', linestyle='--', alpha=0.6)
axs[1, 0].legend(loc='upper right', fontsize=14)
axs[1, 0].set_ylabel(r"Δ $\delta^{13}C$  AGW")
# Plot change in fossil d13C
axs[1, 1].plot(year_d13C, (rFF_d13C_mean_EDGAR), color='blue', linewidth=3, linestyle = 'solid')
axs[1, 1].plot(timeB, FF_d13C_rel, color='firebrick', linewidth=3, linestyle = 'solid')
axs[1, 1].fill_between(year_d13C, rFF_d13C_mean_EDGAR + 2 * rFF_d13C_std_EDGAR, rFF_d13C_mean_EDGAR - 2 * rFF_d13C_std_EDGAR, color='blue', alpha=0.3)
axs[1, 1].yaxis.set_major_locator(plt.MultipleLocator(.5))  # Set grid interval to 10
axs[1, 1].grid(axis='y', linestyle='--', alpha=0.6)
axs[1, 1].set_ylabel(r"Δ $\delta^{13}C$  FF")
# Plot total AGW d13C
axs[0, 0].plot(AGW_Calc_years, AGW_Calc_d13C, color='blue', linewidth=3, label="This study")
axs[0, 0].plot(timeB, AGW_d13C, color='firebrick', linewidth=3, linestyle = 'solid', label="Thanwerdas et al., 2024")
axs[0, 0].fill_between(AGW_Calc_years, AGW_Calc_d13C + 2.6, AGW_Calc_d13C - 2.6, color='blue', alpha=0.3)
axs[0, 0].grid(axis='y', linestyle='--', alpha=0.6)
axs[0, 0].set_ylabel(r"$\delta^{13}C$  AGW")
# Plot total fossil d13C
axs[0, 1].plot(year_d13C, (FF_d13C_mean_EDGAR), color='blue', linewidth=3, linestyle = 'solid')
axs[0, 1].plot(timeB, FF_d13C, color='firebrick', linewidth=3, linestyle = 'solid')
axs[0, 1].fill_between(year_d13C, FF_d13C_mean_EDGAR + 2 * FF_d13C_std_EDGAR, FF_d13C_mean_EDGAR - 2 * FF_d13C_std_EDGAR, color='blue', alpha=0.3)
axs[0, 1].grid(axis='y', linestyle='--', alpha=0.6)
axs[0, 1].set_ylabel(r"$\delta^{13}C$  FF")

axs[0, 0].text(0.01, 0.95, 'A', transform=axs[0, 0].transAxes,
             fontsize=20, fontweight='bold', va='top', ha='left')

axs[0, 1].text(0.01, 0.95, 'B', transform=axs[0,1].transAxes,
             fontsize=20, fontweight='bold', va='top', ha='left')

axs[1, 0].text(0.01, 0.95, 'C', transform=axs[1,0].transAxes,
             fontsize=20, fontweight='bold', va='top', ha='left')

axs[1, 1].text(0.01, 0.95, 'D', transform=axs[1,1].transAxes,
             fontsize=20, fontweight='bold', va='top', ha='left')

plt.tight_layout()  # automatic but you can add a pad: pad=1.5, e.g.
plt.subplots_adjust(wspace=0.3)  # increase spacing between left/right plots
plt.show()

