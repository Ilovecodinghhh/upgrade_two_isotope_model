#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:10:55 2024

@author: ryoung
"""

# This code calculates the global mean d13C of fire emissions

locals().clear()
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize
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

# Load CarbonTracker Emissions, 1998 to 2021
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
pyroA = flux_data.variables['pyrogenic_flux'][:] # Kg/s
# Expand to 1x1 grid
BB_exp = np.repeat(pyroA, repeats=2, axis=1) # Kg/s
BB_exp = np.repeat(BB_exp, repeats=3, axis=2)/6  # Kg/s with smaller grid cells
# Sum to get annual emissions
BB_ann = BB_exp.reshape(24, 12, 180, 360).mean(axis=1) # Kg/s, average per month

# Load GFED Emissions
data_list = []
# Loop through the years from 2002 to 2020
for year in range(2002, 2021):
    # Construct the file name for the current year
    filename = f'data/GFED5_Beta/GFED5_Beta_monthly_{year}.nc'
    # Open the NetCDF file and append the CH4 data array to the list
    with xr.open_dataset(filename) as ds:
        data_list.append(ds['CH4'])
compiled_GFED5 = xr.concat(data_list, dim='time')
compiled_GFED5 = compiled_GFED5.values
# Compress to 1x1 degree map
GFED5_compressed = compiled_GFED5.reshape(228, 180, 4, 360, 4).sum(axis=(2, 4))
# Sum to get annual emissions
GFED5_ann = GFED5_compressed.reshape(19, 12, 180, 360).sum(axis=1)
GFED5_ann = np.flip(GFED5_ann, axis=1)
# Expand to encompass 1998 to 2021
GFED5_ann = np.concatenate((np.repeat(GFED5_ann[0:1, :180, :380], 4, axis=0), GFED5_ann), axis=0)
GFED5_ann = np.concatenate((GFED5_ann, np.repeat(GFED5_ann[-1:], 1, axis=0)), axis=0)

# Load carbontracker prior emissions
Priors = xr.open_dataset('data/prior_monthly_emission_kg_lei.nc')
# Rice 
flux_BB = Priors['flux_biomass_burning']
BB_prior = flux_BB.sum(dim='month')
BB_prior = BB_prior.values.astype(np.float64)
BB_prior_ann = np.transpose(BB_prior[:, :, -24:], (2, 1, 0))
# Sum data annally and weight
BB_prior_ann_sum = BB_prior_ann.sum(axis=(1, 2))
BB_prior_ann_sum = BB_prior_ann_sum[:, np.newaxis, np.newaxis] # Prep ann sum for the multiplication
BB_prior_ann_weight = BB_prior_ann  / BB_prior_ann_sum

# Sum annual and compare
years = np.arange(1998, 2022)


#%% Create BB d13C map 

# Excluding crop residue, weighted mean to heavily weight Cerling et al., 1997
C4_d13C = -12.7
C4_d13C_std = 4.6

C3_d13C = -26.8
C3_d13C_std = 2.9


# Create map for each year and save as txt files
for i, year in enumerate(years):
    # calculation for this year
    BB_d13C_map = C4exp[i, :, :] * C4_d13C / 100 + (100 - C4exp[i, :, :]) * C3_d13C / 100
    
    # flip so row 0 = 90N
    BB_d13C_map_flipped = np.flipud(BB_d13C_map)

    # filename
    filename = f"Output/BB_d13C_1x1_{year}.txt"

    # save as comma-separated text file
    np.savetxt(filename, BB_d13C_map_flipped, delimiter=',', fmt='%+.3f')

# Create source signature map for the 15th year of Luo as a midpoint. 
BB_d13C_map = C4exp[11,:,:]*C4_d13C/100 + (100-C4exp[11,:,:])*C3_d13C/100

# plot C4 proportions
lower_left_lat = -90
lower_left_lon = -180
upper_right_lat = 90
upper_right_lon = 180
vmin = -35  # Minimum value
vmax = -12  # Maximum value
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(BB_d13C_map, cmap= 'inferno', interpolation='nearest', norm=Normalize(vmin=vmin, vmax=vmax))
plt.set_cmap('inferno')
plt.set_cmap('inferno')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('BB ${\delta}13C$-CH$_4$ (‰)', fontsize=8)
m.drawcoastlines(linewidth=.4,color='k')


#%% Calculate Total C3 and C4 BB fluxes for each year

# This multiplies the sum of CTCH4 fire emissions for each grid cell for each year by the corresponding years C3 C4 fraction

# Calculate grid level emissions
C3emissions = C3exp * BB_ann
C4emissions = C4exp * BB_ann
# Repeat for Still
C3emissions_still = C3_still * BB_ann
C4emissions_still = C4_still * BB_ann
# Repeat for GFED
C3emissions_GFED5 = C3exp * GFED5_ann
C4emissions_GFED5 = C4exp * GFED5_ann
# Repeat for prior
C3emissions_prior = C3exp * BB_prior_ann
C4emissions_prior = C4exp * BB_prior_ann

# Sum to get emissions by year
C3total = C3emissions.sum(axis=(1, 2))
C4total = C4emissions.sum(axis=(1, 2))
# Repeat for still
C3total_still = C3emissions_still.sum(axis=(1, 2))
C4total_still = C4emissions_still.sum(axis=(1, 2))
# Repeat for GFED
C3total_GFED5 = C3emissions_GFED5.sum(axis=(1, 2))
C4total_GFED5 = C4emissions_GFED5.sum(axis=(1, 2))
# Repeat for prior
C3total_prior = C3emissions_prior.sum(axis=(1, 2))
C4total_prior = C4emissions_prior.sum(axis=(1, 2))

# Get fraction for each year
C3f = C3total/(C3total+C4total)
C4f = C4total/(C3total+C4total)
# Repeat for still
C3f_still = C3total_still/(C3total_still+C4total_still)
C4f_still = C4total_still/(C3total_still+C4total_still)
# Repeat for GFED
C3f_GFED5 = C3total_GFED5/(C3total_GFED5+C4total_GFED5)
C4f_GFED5 = C4total_GFED5/(C3total_GFED5+C4total_GFED5)
# Repeat for prior
C3f_prior = C3total_prior/(C3total_prior+C4total_prior)
C4f_prior = C4total_prior/(C3total_prior+C4total_prior)

# Put into array with years
# Define years
years = np.linspace(1998, 2021, 24)
C3C4summary = np.stack((years, C3total, C4total, C3f, C4f), axis=1)
C3C4summary_still = np.stack((years, C3total_still, C4total_still, C3f_still, C4f_still), axis=1)
C3C4summary_GFED5 = np.stack((years, C3total_GFED5, C4total_GFED5, C3f_GFED5, C4f_GFED5), axis=1)
C3C4summary_prior = np.stack((years, C3total_prior, C4total_prior, C3f_prior, C4f_prior), axis=1)

#%% Load pyrogenic d13C-CH4 data and perform MC analysis to calculate global mean for each year

# Initialize an array to store the results
result_array = np.zeros((24, 1000))
result_array_still = np.zeros((24, 1000))
result_array_GFED5 = np.zeros((24, 1000))
result_array_prior = np.zeros((24, 1000))

# Perform the Monte Carlo analysis
for i in range(1000):
    # Generate random multipliers with specified standard deviations
    C3_d13C_multiplier = C3_d13C + np.random.normal(1, C3_d13C_std, size=24)
    C4_d13C_multiplier = C4_d13C + np.random.normal(1, C4_d13C_std, size=24)
    
    # Multiply the columns by these random numbers
    multiplied_c3 = C3C4summary[:, 3] * C3_d13C_multiplier
    multiplied_c4 = C3C4summary[:, 4] * C4_d13C_multiplier
    # Repeat for still
    multiplied_c3_still = C3C4summary_still[:, 3] * C3_d13C_multiplier
    multiplied_c4_still = C3C4summary_still[:, 4] * C4_d13C_multiplier
    # Repeat for GFED
    multiplied_c3_GFED5 = C3C4summary_GFED5[:, 3] * C3_d13C_multiplier
    multiplied_c4_GFED5 = C3C4summary_GFED5[:, 4] * C4_d13C_multiplier
    # Repeat for prior
    multiplied_c3_prior = C3C4summary_prior[:, 3] * C3_d13C_multiplier
    multiplied_c4_prior = C3C4summary_prior[:, 4] * C4_d13C_multiplier
    
    # Calculate the sum of these multiplied rows and store in the result array
    result_array[:,i] = (multiplied_c3 + multiplied_c4)
    result_array_still[:,i] = (multiplied_c3_still + multiplied_c4_still)
    result_array_GFED5[:,i] = (multiplied_c3_GFED5 + multiplied_c4_GFED5)
    result_array_prior[:,i] = (multiplied_c3_prior + multiplied_c4_prior)

# Calculate statistics: Luo and CTCH4 emissions
means = result_array.mean(axis=1)
stdev_mean = np.std(means)
std_devs = result_array.std(axis=1)
mean_all = result_array.mean(axis=(0, 1))
std_all = result_array.std(axis=(0, 1))
# Repeat for still
means_still = result_array_still.mean(axis=1)
std_devs_still = result_array_still.std(axis=1)
# Repeat for GFED5
means_GFED5 = result_array_GFED5.mean(axis=1)
std_devs_GFED5 = result_array_GFED5.std(axis=1)
# Repeat for prior
means_prior = result_array_prior.mean(axis=1)
std_devs_prior = result_array_prior.std(axis=1)

# Compile and save as .csv
LUO_CTCH4 = np.stack((years, means, std_devs), axis=1)
LUO_CTCH4 = np.vstack((np.column_stack((years, means, std_devs)), [2030, mean_all, std_all]))
file_path = 'Output/BB_d13C_annual.csv'
np.savetxt(file_path, LUO_CTCH4, delimiter=',')

#%% Plot results

plt.rcParams.update({'font.size': 14})  # Adjust this value as needed
plt.figure(figsize=(10, 6), dpi=1000)

# Plotting the data
plt.plot(years, means, color='blue', label='Luo & CTCH4', linewidth=3)
# plt.fill_between(years, means - std_devs, means + std_devs, color='blue', alpha=0.2)
plt.plot(years, means_still, color='red', label='Still & CTCH4')
#plt.plot(years, means_GFED5, color='orange', label='Luo & GFED5')
plt.plot(years, means_prior, color='orange', label='Luo & GFED4')
plt.errorbar((years[-1]+2), np.mean(mean_all), yerr=np.mean(std_all), fmt='s', color='blue', markersize=8, capsize=5, capthick=1.5, elinewidth=1.5,label='Mean ± 1$\sigma$')
plt.legend()
plt.xlabel('Year')
plt.ylabel('Pyrogenic ${\delta}13C$-CH$_4$ (‰)')
plt.show()


#%% Calculate trend statistics

import numpy as np, statsmodels.api as sm 

print(sm.OLS(means, sm.add_constant(np.arange(len(means)))).fit().summary())

# Final Conclusion
# There is no statistically significant linear trend in your time series.
# The slope is negative but not significant (p = 0.170), meaning the trend could be due to random noise.
# If you still suspect a trend, you might try non-parametric trend tests (e.g., Mann-Kendall test) or smoothing techniques to reveal any hidden patterns.


