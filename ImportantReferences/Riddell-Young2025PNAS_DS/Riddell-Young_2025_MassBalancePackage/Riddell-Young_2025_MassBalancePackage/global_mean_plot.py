#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:07:29 2024

@author: ryoung
"""

# This code makes nice plots of all CH4, C13 and dD data in one stacked figure starting in 1998

locals().clear()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
d13C_glob = Glob_annual_avg.iloc[:, 1].values
d13C_years = Glob_annual_avg.iloc[:, 0].values

# Load d13C-CH4 uncertainties
C13data_ann = pd.read_excel('data/ch4c13_glob_mean_unc.xlsx').to_numpy()
C13_glob_unc = C13data_ann[:,2]

# Load GML global annual means for CH4
CH4data = pd.read_excel('data/GML_CH4_AnnualMean.xlsx').to_numpy()
CH4data = CH4data[11:,:]
CH4 = CH4data[3:,1] # Starting in 1998
CH4year = CH4data[3:,0]

# Load monthly CH4 data 
CH4monthlydata = pd.read_excel('data/ch4_mm_gl.xlsx').to_numpy()
CH4monthly = CH4monthlydata[:,1] # Starting in 1998
CH4monthlyyear = CH4monthlydata[:,0]

# dD data
dDdata = pd.read_csv('../Riddell-Young_2025_dD_GlobMean/output/GlobMean_dD_DEI_UmezawaCal_noBUDS.csv')
dDann = dDdata.iloc[:,1]
dDyear = dDdata.iloc[:,0]
dDannR = dDdata.iloc[:,3]
dDann_unc = dDdata.iloc[:,2]
dDannR_unc = dDdata.iloc[:,4]

# Weekly dD data
dDdata_weekly = pd.read_csv('../Riddell-Young_2025_dD_GlobMean/output/HemMean_dD_dei_UmezawaCal_noBUDS.csv')
dDdata_weekly = dDdata_weekly[~np.isnan(dDdata_weekly).any(axis=1)]
dDann_weekly = dDdata_weekly.iloc[:-1,1].to_numpy()
dDNH_weekly = dDdata_weekly.iloc[:-1,2].to_numpy()
dDSH_weekly = dDdata_weekly.iloc[:-1,3].to_numpy()
dDyear_weekly = dDdata_weekly.iloc[:-1,0].to_numpy()
dDyear_weekly = dDyear_weekly[~np.isnan(dDann_weekly)]
dDann_weekly = dDann_weekly[~np.isnan(dDann_weekly)]
dDNH_weekly = dDNH_weekly[~np.isnan(dDNH_weekly)]
dDSH_weekly = dDSH_weekly[~np.isnan(dDSH_weekly)]

# Make dD data monthly
dD_monthly = dDann_weekly[:-2].reshape(-1, 4).mean(axis=1)
dDyear_monthly = dDyear_weekly[:-2].reshape(-1, 4).mean(axis=1)


#%% PNAS-style Figure 1

import matplotlib.pyplot as plt

# set global font size ~10 pt
plt.rc('font', size=10)

# helper to convert cm to inch
def cm2inch(x): return x/2.54

# make figure: 3 stacked subplots, large format
fig, axs = plt.subplots(
    3, 1, sharex=True,
    figsize=(cm2inch(18), cm2inch(22)),
    dpi=200
)

# unpack axes for clarity
ax1, ax2, ax3 = axs

# --- ch4 data ---
ax1.plot((CH4year+0.5), CH4, label="CH$_4$ Annual Mean", color='black')
ax1.plot(CH4monthlyyear, CH4monthly, label="CH$_4$ Monthly Mean", color='black', alpha=0.5)
ax1.set_ylabel("CH$_4$ (ppb)")
ax1.legend(loc="center left", fontsize=10)
ax1.set_xlim(1998, 2024)
ax1.set_ylim(1750, 1930)
ax1.tick_params(labelbottom=False)

# --- d13c data ---
ax2.plot((d13C_years[1:]+0.5), d13C_glob[1:], label=r"${\delta}^{13}\mathrm{C}$ Annual Mean", color='red')
ax2.plot((d13C_years[1:]+0.5), (d13C_glob[1:] + 2*C13_glob_unc),
         color='red', linestyle='dashed', linewidth=0.5)
ax2.plot((d13C_years[1:]+0.5), (d13C_glob[1:] - 2*C13_glob_unc),
         color='red', linestyle='dashed', linewidth=0.5)
ax2.plot(glob_dates[24:], glob_mean[24:], label=r"${\delta}^{13}\mathrm{C}$ Monthly Mean",
         color='red', alpha=0.5)
ax2.set_ylabel(r"${\delta}^{13}\mathrm{C-CH_4}$ (‰)")
ax2.legend(loc="lower left", fontsize=10)
ax2.set_xlim(1998, 2024)
ax2.tick_params(labelbottom=False)

# --- dd data ---
ax3.plot((dDyear+0.5), dDann, label='${\delta}D$ Annual Mean', color='blue')
ax3.plot((dDyear+0.5), (dDann + 2*dDann_unc),
         color='blue', linestyle='dashed', linewidth=0.5)
ax3.plot((dDyear+0.5), (dDann - 2*dDann_unc),
         color='blue', linestyle='dashed', linewidth=0.5)
ax3.plot(dDyear_monthly, dD_monthly, label='${\delta}D$ Monthly Mean', color='blue', alpha=0.5)
ax3.set_ylabel('${\delta}D$-CH$_4$ (‰)')
ax3.legend(loc="lower left", fontsize=10)
ax3.set_xlim(1998, 2024)

# panel labels
ax1.text(0.01, 0.97, 'A', transform=ax1.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='left')
ax2.text(0.01, 0.97, 'B', transform=ax2.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='left')
ax3.text(0.01, 0.97, 'C', transform=ax3.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='left')

# tidy up
for ax in axs:
    ax.grid(False)

plt.tight_layout(rect=[0, 0, 1, 0.96])

# save and show
plt.savefig("Output/fig1.pdf", dpi=200, bbox_inches="tight")
plt.show()
plt.close()



#%% Plot up interhemispheric difference 

plt.gcf().set_facecolor('w')
plt.figure(dpi=500)
plt.plot(dDyear_weekly, (dDNH_weekly - dDSH_weekly),color='blue')
plt.ylabel('NH - SH ${\delta}D$-CH$_4$ (‰)')
plt.xlim(2005,2024)

