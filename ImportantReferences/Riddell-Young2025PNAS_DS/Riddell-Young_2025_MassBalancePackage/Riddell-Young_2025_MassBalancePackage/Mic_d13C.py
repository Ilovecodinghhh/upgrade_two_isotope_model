#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:36:38 2024

@author: ryoung
"""

# This code calculates the global mean and uncertainty d13C of microbial emissions

locals().clear()
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd
import xarray as xr


#%% Load data 

# Load C3/C4 map, 2001 to 2019 (Luo et al., 2023)
C3C4_data = nc.Dataset('data/C4_distribution_NUS_v2.2.nc', 'r')
C4_maps = C3C4_data.variables['C4_area'][:]
C4_maps = np.nan_to_num(C4_maps)
# Compress to a 1x1 grid
C4_maps_compressed = C4_maps.reshape(19, 360, 2, 180, 2)
C3C4 = C4_maps_compressed.mean(axis=(2, 4))
C3C4 = np.transpose(C3C4, (0, 2, 1))
# Expand array to include 1998, 1999, 2000, 2020 and 2021
C4exp = np.concatenate((np.repeat(C3C4[0:1, :180, :380], 3, axis=0), C3C4), axis=0)
C4exp = np.concatenate((C4exp, np.repeat(C4exp[-1:], 2, axis=0)), axis=0)
# Flip C4 map
C4exp = np.flip(C4exp, axis=1)
# Create C3 map
C3exp = 100 - C4exp

# Load old still et al., 2003 map
C4_still = pd.read_excel('data/Still2003_C4.xlsx', header=None)
C4_still = np.flipud(C4_still)
C4_still = np.where(C4_still < 0, 0, C4_still)
# Deine C3 matrix
C3_still = 100 - C4_still

# Load posterior microbial fluxes
flux_data = nc.Dataset('data/CTCH4_2023_flux3x2.nc', 'r')
time_start = flux_data.variables['time_start'][:]
time_mid = flux_data.variables['time_mid'][:]
time_end = flux_data.variables['time_end'][:]
# Convert to decimal date
time_end_dcml = []
time_start_dcml = []
time_mid_dcml = []
for i in range(len(time_start)):
    time_end_dcml.append(time_end[i, 0] + (time_end[i, 1] - 1) / 12 + (time_end[i, 2] - 1) / 365.25)
    time_start_dcml.append(time_start[i, 0] + (time_start[i, 1] - 1) / 12 + (time_start[i, 2] - 1) / 365.25)
    time_mid_dcml.append(time_mid[i, 0] + (time_mid[i, 1] - 1) / 12 + (time_mid[i, 2] - 1) / 365.25)
# Load emissions data
MicA = flux_data.variables['microbial_flux'][:]
# expand to 1x1 grid
Mic_exp = np.repeat(MicA, repeats=2, axis=1)
Mic_exp = np.repeat(Mic_exp, repeats=3, axis=2)/6 
# Convert from kg/s to Tg/month
# Constants
kpt = 1000000000  # Kg per Tg
spm = 2.628e+6  # Sec per Month
Mic = Mic_exp / kpt * spm
# Flip on latitudinal axis
Mic = np.flip(Mic, axis=1)
# Calculate the annual sum and weight each year by its respective total
Mic_ann = Mic_exp.reshape(24, 12, 180, 360).sum(axis=1)
Mic_ann_sum = Mic_ann.sum(axis=(1, 2))
Mic_ann_sum = Mic_ann_sum[:, np.newaxis, np.newaxis] # Prep ann sum for the multiplication

# Load CTCH4 Prior emissions
Priors = xr.open_dataset('data/prior_monthly_emission_kg_lei.nc')
# Rice 
flux_Rice = Priors['flux_rice']
Rice_ann = flux_Rice.sum(dim='month')
Rice_ann = Rice_ann.values.astype(np.float64)
Rice_ann = np.transpose(Rice_ann[:, :, -24:], (2, 1, 0))
# Ruminant
flux_Ruminant = Priors['flux_ruminant']
Ruminant_ann = flux_Ruminant.sum(dim='month')
Ruminant_ann = Ruminant_ann.values.astype(np.float64)
Ruminant_ann = np.transpose(Ruminant_ann[:, :, -24:], (2, 1, 0))
# Wetland
flux_Wetland = Priors['flux_wetland']
Wetland_ann = flux_Wetland.sum(dim='month')
Wetland_ann = Wetland_ann.values.astype(np.float64)
Wetland_ann = np.transpose(Wetland_ann[:, :, -24:], (2, 1, 0))
# Termite
flux_Termite = Priors['flux_termite']
Termite_ann = flux_Termite.sum(dim='month')
Termite_ann = Termite_ann.values.astype(np.float64)
Termite_ann = np.transpose(Termite_ann[:, :, -24:], (2, 1, 0))
# Landfill
flux_Landfill = Priors['flux_waste_landfill']
Landfill_ann = flux_Landfill.sum(dim='month')
Landfill_ann = Landfill_ann.values.astype(np.float64)
Landfill_ann = np.transpose(Landfill_ann[:, :, -24:], (2, 1, 0))
# Wild Animal
flux_WildAnimal = Priors['flux_wild_animals']
WildAnimal_ann = flux_WildAnimal.sum(dim='month')
WildAnimal_ann = WildAnimal_ann.values.astype(np.float64)
WildAnimal_ann = np.transpose(WildAnimal_ann[:, :, -24:], (2, 1, 0))

# Load temporal wetland map
wetland_d13C_nc = nc.Dataset('data/isotem_wetland_d13C-CH4.nc4', 'r') # Oh et al., 2022
wetland_d13C = wetland_d13C_nc.variables['wetland_d13C-CH4'][:]
wetland_d13C_ann = np.mean(wetland_d13C, axis=0)  # Resulting in a shape of (33, 720, 360)
wetland_d13C_ann = np.swapaxes(wetland_d13C_ann, 1, 2)
wetland_d13C_ann = wetland_d13C_ann[:, ::-1, :]
# Replace NaN with annual averge
wetland_d13C_ann = np.ma.filled(wetland_d13C_ann, fill_value=np.nan)
wetland_d13C_annavg = np.nanmean(wetland_d13C_ann, axis=(1, 2)).reshape(33, 1)
for i in range(wetland_d13C_ann.shape[0]):
    wetland_d13C_ann[i, np.isnan(wetland_d13C_ann[i])] = wetland_d13C_annavg[i]
wetland_d13C_ann[1, np.isnan(wetland_d13C_ann[1])] = wetland_d13C_annavg[1]
# Change years 
# Assuming wetland_d13C_final has dimensions 33x180x360
years = np.arange(1984, 2017)  # Array of years from 1984 to 2016
start_year = 1998
end_year = 2016
# Indices for years from 1998 to 2016
indices = (years >= start_year) & (years <= end_year)
trimmed_array = wetland_d13C_ann[indices]
# Repeat the year 2016 five more times
wetland_d13C_trimmed = np.concatenate([trimmed_array, np.repeat(trimmed_array[-1:, :, :], 5, axis=0)], axis=0)
# convert to 24x180x360
wetland_d13C_final = wetland_d13C_trimmed.reshape(24, 360 // 2, 2, 720 // 2, 2).mean(axis=(2, 4))

# save 1x1 degree maps
years = np.arange(1998, 2022)  # 1998–2021

for i, year in enumerate(years):
    arr2d = np.flipud(wetland_d13C_final[i, :, :])  # flip so row 0 = 90N
    filename = f"Output/Wetland_d13C_1x1_{year}.txt"
    np.savetxt(filename, arr2d, delimiter=',', fmt='%+.3f')


#%% Weighting Calculations

# Sum all prior emissions
prior_mic_total = np.sum(Rice_ann + Ruminant_ann + Wetland_ann + Termite_ann + Landfill_ann + WildAnimal_ann)
# Calculate total prior mic
prior_mic = (Rice_ann + Ruminant_ann + Wetland_ann + Termite_ann + Landfill_ann + WildAnimal_ann) / spm
# Create weighting scheme 
Posterior_to_prior = Mic_ann/prior_mic
Posterior_to_prior = np.where(Posterior_to_prior == '--', 1, Posterior_to_prior).astype(float)
Posterior_to_prior = np.where(Posterior_to_prior == 0, 1, Posterior_to_prior)

# Weight priors based on posterior ratios
Rice_ann_posterior = Rice_ann*Posterior_to_prior
Ruminant_ann_posterior = Ruminant_ann*Posterior_to_prior
Wetland_ann_posteriorB = Wetland_ann*Posterior_to_prior
Landfill_ann_posterior = Landfill_ann*Posterior_to_prior
Termite_ann_posterior = Termite_ann*Posterior_to_prior
WildAnimal_ann_posterior = WildAnimal_ann*Posterior_to_prior
prior_mic_posterior = prior_mic*Posterior_to_prior
# Add unassigned values to wetlands
Wetland_ann_posterior = Wetland_ann_posteriorB + Mic_ann*spm - Rice_ann_posterior - Ruminant_ann_posterior - Wetland_ann_posteriorB - Termite_ann_posterior - Landfill_ann_posterior - WildAnimal_ann_posterior

# Weight subcategories based on annual sums using posterior weighting
Rice_ann_weighted = Rice_ann_posterior / np.sum(Rice_ann_posterior, axis=(1, 2), keepdims=True) 
Ruminant_ann_weighted = Ruminant_ann_posterior / np.sum(Ruminant_ann_posterior, axis=(1, 2), keepdims=True) 
Wetland_ann_weighted = Wetland_ann_posterior / np.sum(Wetland_ann_posterior, axis=(1, 2), keepdims=True) 
Landfill_ann_weighted = Landfill_ann_posterior / np.sum(Landfill_ann_posterior, axis=(1, 2), keepdims=True)
Termite_ann_weighted = Termite_ann_posterior / np.sum(Termite_ann_posterior, axis=(1, 2), keepdims=True)
WildAnimal_ann_weighted = WildAnimal_ann_posterior / np.sum(WildAnimal_ann_posterior, axis=(1, 2), keepdims=True) 

# Calculate posterior mic annual totals.
Rice_Post_AnnTotal = np.sum(Rice_ann_posterior, axis=(1, 2))
Ruminant_Post_AnnTotal = np.sum(Ruminant_ann_posterior, axis=(1, 2))
Wetland_Post_AnnTotal = np.sum(Wetland_ann_posterior, axis=(1, 2))
Landfill_Post_AnnTotal = np.sum(Landfill_ann_posterior, axis=(1, 2))
Termite_Post_AnnTotal = np.sum(Termite_ann_posterior, axis=(1, 2))
WildAnimal_Post_AnnTotal = np.sum(WildAnimal_ann_posterior, axis=(1, 2))
Mic_Posterior_AnnTotal = np.sum(Mic_ann, axis=(1, 2))

# Calulate posterior fractions
Rice_Post_Frac = (Rice_Post_AnnTotal / spm) / Mic_Posterior_AnnTotal
Ruminant_Post_Frac = (Ruminant_Post_AnnTotal / spm) / Mic_Posterior_AnnTotal
Wetland_Post_Frac = (Wetland_Post_AnnTotal / spm) / Mic_Posterior_AnnTotal
Landfill_Post_Frac = (Landfill_Post_AnnTotal / spm) / Mic_Posterior_AnnTotal
Termite_Post_Frac = (Termite_Post_AnnTotal / spm) / Mic_Posterior_AnnTotal
WildAnimal_Post_Frac = (WildAnimal_Post_AnnTotal / spm) / Mic_Posterior_AnnTotal


#%% Now do runimants with MC analysis for error propagation

# Calculate grid level emissions, posterior
C3emissions_posterior = C3exp * Ruminant_ann_posterior
C4emissions_posterior = C4exp * Ruminant_ann_posterior
# Repeat for prior
C3emissions = C3_still * Ruminant_ann
C4emissions = C4_still * Ruminant_ann

# Sum to get emissions by year
C3total = C3emissions.sum(axis=(1, 2))
C4total = C4emissions.sum(axis=(1, 2))
# Repeat for posterior
C3total_posterior = C3emissions_posterior.sum(axis=(1, 2))
C4total_posterior = C4emissions_posterior.sum(axis=(1, 2))

# Get fraction for each year
C3f = C3total/(C3total+C4total)
C4f = C4total/(C3total+C4total)
# Repeat for _posterior
C3f_posterior = C3total_posterior/(C3total_posterior+C4total_posterior)
C4f_posterior = C4total_posterior/(C3total_posterior+C4total_posterior)

# Put into array with years
# Define years
years = np.linspace(1998, 2021, 24)
C3C4summary = np.stack((years, C3total, C4total, C3f, C4f), axis=1)
C3C4summary_posterior = np.stack((years, C3total_posterior, C4total_posterior, C3f_posterior, C4f_posterior), axis=1)

C3_Rum_d13C = -66.63846154	
C3_Rum_d13C_std = 3.385886999
C4_Rum_d13C = -54.95714286
C4_Rum_d13C_std = 3.430188111

# Initialize an array to store the results
result_array_Rum_posterior = np.zeros((24, 1000))
result_array_Rum = np.zeros((24, 1000))

# Perform the Monte Carlo analysis
for i in range(1000):
    # Generate random multipliers with specified standard deviations
    C3_Rum_d13C_multiplier = C3_Rum_d13C + np.random.normal(1, C3_Rum_d13C_std, size=24)
    C4_Rum_d13C_multiplier = C4_Rum_d13C + np.random.normal(1, C4_Rum_d13C_std, size=24)
    # Multiply the columns by these random numbers
    multiplied_c3_Rum = C3C4summary[:, 3] * C3_Rum_d13C_multiplier
    multiplied_c4_Rum = C3C4summary[:, 4] * C4_Rum_d13C_multiplier
    # Repeat for posterior
    multiplied_c3_Rum_posterior = C3C4summary_posterior[:, 3] * C3_Rum_d13C_multiplier
    multiplied_c4_Rum_posterior = C3C4summary_posterior[:, 4] * C4_Rum_d13C_multiplier
    # Calculate the sum of these multiplied rows and store in the result array
    result_array_Rum[:,i] = (multiplied_c3_Rum + multiplied_c4_Rum)
    result_array_Rum_posterior[:,i] = (multiplied_c3_Rum_posterior + multiplied_c4_Rum_posterior)

# Calculate statistics
Rum_means = result_array_Rum.mean(axis=1)
stdev_mean = np.std(Rum_means)
Rum_std_devs = result_array_Rum.std(axis=1)
# Repeat for posterior
Rum_means_posterior = result_array_Rum_posterior.mean(axis=1)
Rum_std_devs_posterior = result_array_Rum_posterior.std(axis=1)


#%% Repeat for wild animals

# Calculate grid level emissions, posterior
C3emissions_posterior = C3exp * WildAnimal_ann_posterior
C4emissions_posterior = C4exp * WildAnimal_ann_posterior
# Repeat for prior
C3emissions = C3_still * WildAnimal_ann
C4emissions = C4_still * WildAnimal_ann

# Sum to get emissions by year
C3total = C3emissions.sum(axis=(1, 2))
C4total = C4emissions.sum(axis=(1, 2))
# Repeat for posterior
C3total_posterior = C3emissions_posterior.sum(axis=(1, 2))
C4total_posterior = C4emissions_posterior.sum(axis=(1, 2))

# Get fraction for each year
C3f = C3total/(C3total+C4total)
C4f = C4total/(C3total+C4total)
# Repeat for _posterior
C3f_posterior = C3total_posterior/(C3total_posterior+C4total_posterior)
C4f_posterior = C4total_posterior/(C3total_posterior+C4total_posterior)

# Put into array with years
# Define years
years = np.linspace(1998, 2021, 24)
C3C4summary = np.stack((years, C3total, C4total, C3f, C4f), axis=1)
C3C4summary_posterior = np.stack((years, C3total_posterior, C4total_posterior, C3f_posterior, C4f_posterior), axis=1)

# Used same number for ruminants
C3_WildAn_d13C = -66.63846154	
C3_WildAn_d13C_std = 3.385886999
C4_WildAn_d13C = -54.95714286
C4_WildAn_d13C_std = 3.430188111

# Initialize an array to store the results
result_array_WildAn_posterior = np.zeros((24, 1000))
result_array_WildAn = np.zeros((24, 1000))

# Perform the Monte Carlo analysis
for i in range(1000):
    # Generate random multipliers with specified standard deviations
    C3_WildAn_d13C_multiplier = C3_WildAn_d13C + np.random.normal(1, C3_WildAn_d13C_std, size=24)
    C4_WildAn_d13C_multiplier = C4_WildAn_d13C + np.random.normal(1, C4_WildAn_d13C_std, size=24)
    # Multiply the columns by these random numbers
    multiplied_c3_WildAn = C3C4summary[:, 3] * C3_WildAn_d13C_multiplier
    multiplied_c4_WildAn = C3C4summary[:, 4] * C4_WildAn_d13C_multiplier
    # Repeat for posterior
    multiplied_c3_WildAn_posterior = C3C4summary_posterior[:, 3] * C3_WildAn_d13C_multiplier
    multiplied_c4_WildAn_posterior = C3C4summary_posterior[:, 4] * C4_WildAn_d13C_multiplier
    # Calculate the sum of these multiplied rows and store in the result array
    result_array_WildAn[:,i] = (multiplied_c3_WildAn + multiplied_c4_WildAn)
    result_array_WildAn_posterior[:,i] = (multiplied_c3_WildAn_posterior + multiplied_c4_WildAn_posterior)

# Calculate statistics
WildAn_means = result_array_WildAn.mean(axis=1)
stdev_mean = np.std(WildAn_means)
WildAn_std_devs = result_array_WildAn.std(axis=1)
# Repeat for still
WildAn_means_posterior = result_array_WildAn_posterior.mean(axis=1)
WildAn_std_devs_posterior = result_array_WildAn_posterior.std(axis=1)


#%% Check matrix orientation

from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize

# Plot BB annual emissions gridded
lower_left_lat = -90
lower_left_lon = -180
upper_right_lat = 90
upper_right_lon = 180
vmin = -75  # Minimum value
vmax = -30  # Maximum value
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(wetland_d13C_final[-1,:,:], cmap= 'gist_earth', interpolation='nearest', norm=Normalize(vmin=vmin, vmax=vmax),origin="upper")
plt.set_cmap('gist_earth')
plt.set_cmap('gist_earth')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('Wetland ${\delta}13C$-CH$_4$ (‰)', fontsize=8)
m.drawcoastlines(linewidth=.4,color='k')


#%% Claculate global mean microbial in a mass balance framework

# Load ruminant and wetland data
Ruminants = pd.read_excel('data/Chang_2019_ruminants.xlsx').values
year = Ruminants[37:-3,0]
Rum_d13C = Ruminants[37:-3,1]
Rum_d13C_unc = 1.45

Wetlands = pd.read_excel('data/Oh_2022_Wetlands.xlsx').values
yearW = Wetlands[:,0]
Wetland_d13C = Wetlands[:,1]
Wetland_d13C_unc = 0.7

# Define rice, termite, and waste/landfill source sigs and uncertainties from Riddell-Young 2025 source signature inventory
Waste_d13C = -54.8 #-54.6
Waste_d13C_unc = 4.4
Rice_d13C = -59.9 #-59.4
Rice_d13C_unc = 4.5
Termite_d13C = -65.2
Termite_d13C_unc = 7.6
WildAn_d13C = np.mean(WildAn_means_posterior)
WildAn_d13C_unc = np.mean(WildAn_std_devs_posterior)

# Define Suess Effect trend
Trend = -.024 # Annual trend in source signature due to suess effect. 
Trend_Unc = 0.005 #If we want to add uncertainty here
# Apply Suess Effect trend to remaining sources
Waste_d13C_Suess = Waste_d13C + (np.arange(24) - 12) * Trend
Rice_d13C_Suess = Rice_d13C + (np.arange(24) - 12) * Trend
Termite_d13C_Suess = Termite_d13C + (np.arange(24) - 12) * Trend
WildAn_d13C_Suess = WildAn_d13C + (np.arange(24) - 12) * Trend

Mic_d13C_MC = []
for i in range(1000):
    # Add random uncertainty to source signatures
    Waste_d13C_MC = Waste_d13C_Suess + Waste_d13C_unc*np.random.normal(0,1)
    Rice_d13C_MC = Rice_d13C_Suess + Rice_d13C_unc*np.random.normal(0,1)
    Termite_d13C_MC = Termite_d13C_Suess + Termite_d13C_unc*np.random.normal(0,1)
    WildAn_d13C_MC = WildAn_d13C_Suess + WildAn_d13C_unc*np.random.normal(0,1)
    Wetland_d13C_MC = Wetland_d13C + Wetland_d13C_unc*np.random.normal(0,1)
    Rum_d13C_MC = Rum_d13C + Rum_d13C_unc*np.random.normal(0,1)
    
    # Add random uncertainty to subcategory proportions 
    perturbations = np.random.normal(loc=1, scale=0.1, size=6)
    # Combine all subcategories into matrix
    SubCatProp = np.array([Rice_Post_Frac,Ruminant_Post_Frac,Wetland_Post_Frac,Landfill_Post_Frac,Termite_Post_Frac,WildAnimal_Post_Frac]).T
    SubCatProp_MC = SubCatProp*perturbations
    SubCatProp_MC /= SubCatProp_MC.sum(axis=1, keepdims=True) # Rescale each row to 1

    #Calculate weighted mean
    Mic_d13C_MCrun = SubCatProp_MC[:,0]*Rice_d13C_MC + SubCatProp_MC[:,1]*Rum_d13C_MC + SubCatProp_MC[:,2]*Wetland_d13C_MC + SubCatProp_MC[:,3]*Waste_d13C_MC + SubCatProp_MC[:,4]*Termite_d13C_MC + SubCatProp_MC[:,5]*WildAn_d13C_MC
    Mic_d13C_MC.append(Mic_d13C_MCrun)

Mic_d13C_MC = np.array(Mic_d13C_MC).T

# Calculate statistics
mean_d13C_mic  = np.mean(Mic_d13C_MC,axis=1)
meanmean = np.mean(mean_d13C_mic)
stdev_d13C_mic = np.std(Mic_d13C_MC,axis=1)
meanstd = np.mean(stdev_d13C_mic)


# Plot results
plt.rcParams.update({'font.size': 14})  # Adjust this value as needed
plt.figure(figsize=(6, 10), dpi=1000)
plt.plot(year, Rum_d13C, color='blue', label='Ruminants', linewidth=1.5)
plt.plot(year, Wetland_d13C, color='green', label='Wetlands', linewidth=1.5)
plt.plot(year, Waste_d13C_Suess, color='red', label='Waste', linewidth=1.5)
plt.plot(year, Rice_d13C_Suess, color='cyan', label='Rice', linewidth=1.5)
plt.plot(year, Termite_d13C_Suess , color='orange', label='Termites', linewidth=1.5)
plt.plot(year, WildAn_d13C_Suess , color='gray', label='Wild Animals', linewidth=1.5)
plt.plot(year, mean_d13C_mic, color='black', label='All Mic', linewidth=4)
plt.xlabel('Year')
#plt.title('Microbial ${\delta}13C$-CH$_4$ (‰)')
plt.ylabel('Global Annual Mean ${\delta}13C$-CH$_4$ (‰)')
plt.legend(loc='upper left', ncol=2)
plt.ylim(-66, -50)  
plt.show()

# Compile and save as .csv
Mic_globmean = np.stack((year, mean_d13C_mic, stdev_d13C_mic), axis=1)
file_path = 'Output/Mic_d13C_annual.csv'
np.savetxt(file_path, Mic_globmean, delimiter=',')


#%% Plot summary figure

import numpy as np
import matplotlib.pyplot as plt

# Sample data
categories = ['Wetlands', 'Rice', 'Ruminants', 'Wild Animals', 'Termites', 'Landfills', 'Total']
SrcSigs = np.array([Wetland_d13C[-1], Rice_d13C_Suess[-1], Rum_d13C[-1], Termite_d13C_Suess[-1], WildAn_d13C_Suess[-1], Waste_d13C_Suess[-1], mean_d13C_mic[-1]])  # Below the line
Proportions = np.array([np.mean(Wetland_Post_Frac), np.mean(Rice_Post_Frac), np.mean(Ruminant_Post_Frac), 
                        np.mean(WildAnimal_Post_Frac), np.mean(Termite_Post_Frac), np.mean(Landfill_Post_Frac), 1])  # Above the line
uncertainty_SrcSigs = np.array([0.7, Rice_d13C_unc, 1.45, WildAn_d13C_unc, Termite_d13C_unc, Waste_d13C_unc,stdev_d13C_mic[-1]])  # Error bars for below bars
uncertainty_Prop = np.array([np.mean(Wetland_Post_Frac)*0.2, np.mean(Rice_Post_Frac)*0.2, np.mean(Ruminant_Post_Frac)*0.2, 
                        np.mean(WildAnimal_Post_Frac)*0.2, np.mean(Termite_Post_Frac)*0.2, np.mean(Landfill_Post_Frac)*0.2, 0]) # Error bars for above bars

x = np.arange(len(categories))  # X positions for bars
bar_width = 0.4  # Width of the bars
colors = ['royalblue', 'green', 'purple', 'orange', 'red', 'cyan', 'gray']  # Unique colors per category

plt.rcParams.update({'font.size': 12}) 
fig, ax1 = plt.subplots(figsize=(6, 6), dpi=1000)  # Increase figure height
ax2 = ax1.twinx()  # Create a twin y-axis for separate limits

# Plot bars above the line (solid) on ax1
for i, color in enumerate(colors):
    ax1.bar(x[i], Proportions[i], yerr=uncertainty_Prop[i], capsize=5, 
            width=bar_width, color=color, label=f'Category {categories[i]}' if i == 0 else "")

# Plot bars below the line (dashed edge) on ax2
for i, color in enumerate(colors):
    ax2.bar(x[i], SrcSigs[i] - (-47), yerr=uncertainty_SrcSigs[i], capsize=5, 
            width=bar_width, color=color, edgecolor='black', linestyle='dashed', hatch='//', bottom=-47)

# Set separate y-axis limits
ax1.set_ylim(-1, 1)  
ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_ylim(-75, -19)
ax2.set_yticks([-50, -60, -70])

# Adjust positioning to make 0 to 1 take up half of the figure
ax1.set_position([0.15, 0.55, 0.7, 0.4])  # Move ax1 higher and shrink height
ax2.set_position([0.15, 0.05, 0.7, 0.4])  # Move ax2 lower and shrink height

# Formatting
ax1.axhline(0, color='black', linewidth=1)  # Reference line at y=0
ax1.set_xticks(x)
ax1.set_xticklabels(categories, rotation=20)
ax1.set_ylabel("Proportions", color='black')
ax2.set_ylabel(r"${\delta}^{13}\mathrm{C-CH_4}$ (‰)", color='black')
ax1.set_xlabel("Microbial Emission Category")
ax1.set_title("Microbial ${\delta}^{13}\mathrm{C-CH_4}$  for 2022")
ax1.yaxis.set_label_coords(-0.14, 0.75)  # Moves ax1 label higher
ax2.yaxis.set_label_coords(1.14, 0.25)  # Moves ax2 label lower

# Custom legend
handles = [plt.Rectangle((0,0),1,1, color=color) for color in colors]
# ax1.legend(handles, categories, title="Categories", loc='upper right')

plt.show()


#%% Create MC array of mic scenarios with varying trend and varying mean value

# First calculate relative trend in mean
Mic_trend = mean_d13C_mic - np.mean(mean_d13C_mic)
# Calculate mean mean and mean standard deviation
Mic_stdev_mean = np.mean(stdev_d13C_mic)
Mic_mean_mean = np.mean(mean_d13C_mic)


# Initialize an empty array to store results (24 rows, 1000 cols)
Mic_d13C_MC_trend = np.zeros((24, 1000))
rMic_d13C_MC_trend = np.zeros((24, 1000))
for i in range(1000):
    # Generate random numbers
    rn = np.random.normal(loc=0, scale=0.25)  # Single random number
    rn2 = np.random.normal(loc=0, scale=Mic_stdev_mean)  # 24-element vector
    
    # Create MC mic trends
    Mic_mean_MC = Mic_trend * rn + Mic_trend + Mic_mean_mean + rn2  # Shape (24,)
    rMic_mean_MC = Mic_trend * rn + Mic_trend # Shape (24,)
    
    # Store result in the i-th column
    Mic_d13C_MC_trend[:, i] = Mic_mean_MC
    rMic_d13C_MC_trend[:, i] = rMic_mean_MC
    
    
#%% Plot and save data

plt.rc('font', size=14) 
fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=1000, sharex=True)  # 2 rows, 1 column
# First subplot: Relative Mic δD-CH4
plt.suptitle('Microbial ${\delta}13C$ Simulations')
axes[0].plot(years, rMic_d13C_MC_trend)
axes[0].set_ylabel('Relative Mic ${\delta}13C$-CH$_4$  (‰)')
axes[0].grid(False)
# Second subplot: Mic δD-CH4
axes[1].plot(years, Mic_d13C_MC_trend)
axes[1].set_ylabel('Mic ${\delta}13C$-CH$_4$  (‰)')
axes[1].grid(False)
axes[1].set_xlabel('Year')  # Shared x-axis label
plt.tight_layout()  # Adjust layout for better spacing
plt.show()


# Compile and save as .csv
Mic_MC = np.column_stack((year, Mic_d13C_MC_trend))
file_path = 'Output/Mic_d13C_MC.csv'
np.savetxt(file_path, Mic_MC, delimiter=',')
















