#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:59:49 2024

@author: ryoung
"""

# Code for developing global average microbial dD-CH4
locals().clear()
import tifffile
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import xarray as xr
from matplotlib.colors import Normalize

plt.clf()
plt.close('all')

#%% load data

# Load prior netcdf file
Priors = xr.open_dataset('data/prior_monthly_emission_kg_lei.nc')
flux_Rice = Priors['flux_rice']
Rice_ann = flux_Rice.sum(dim='month')
Rice_ann = Rice_ann.values.astype(np.float64)
Rice_ann = np.transpose(Rice_ann[:, :, -24:], (2, 1, 0))

flux_Ruminant = Priors['flux_ruminant']
Ruminant_ann = flux_Ruminant.sum(dim='month')
Ruminant_ann = Ruminant_ann.values.astype(np.float64)
Ruminant_ann = np.transpose(Ruminant_ann[:, :, -24:], (2, 1, 0))

flux_Wetland = Priors['flux_wetland']
Wetland_ann = flux_Wetland.sum(dim='month')
Wetland_ann = Wetland_ann.values.astype(np.float64)
Wetland_ann = np.transpose(Wetland_ann[:, :, -24:], (2, 1, 0))

flux_Termite = Priors['flux_termite']
Termite_ann = flux_Termite.sum(dim='month')
Termite_ann = Termite_ann.values.astype(np.float64)
Termite_ann = np.transpose(Termite_ann[:, :, -24:], (2, 1, 0))

flux_Landfill = Priors['flux_waste_landfill']
Landfill_ann = flux_Landfill.sum(dim='month')
Landfill_ann = Landfill_ann.values.astype(np.float64)
Landfill_ann = np.transpose(Landfill_ann[:, :, -24:], (2, 1, 0))

flux_WildAnimal = Priors['flux_wild_animals']
WildAnimal_ann = flux_WildAnimal.sum(dim='month')
WildAnimal_ann = WildAnimal_ann.values.astype(np.float64)
WildAnimal_ann = np.transpose(WildAnimal_ann[:, :, -24:], (2, 1, 0))

# Load posterior microbial fluxes
flux_data = nc.Dataset('data/CTCH4_2023_flux3x2.nc', 'r')
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


#%% Visualize the latitudinal contribution of different subcategories

# Define latitude
latitude = np.linspace(-90, 90, 180).reshape(180, 1)

# Normalize by latitude
# Rice
Rice_lat = Rice_ann.sum(axis=(0, 2)).reshape(180, 1)
Rice_total = Rice_lat.sum()
Rice_lat_rel = Rice_lat / Rice_total
# Ruminant
Ruminant_lat = Ruminant_ann.sum(axis=(0, 2)).reshape(180, 1)
Ruminant_total = Ruminant_lat.sum()
Ruminant_lat_rel = Ruminant_lat / Ruminant_total
# Wetland
Wetland_lat = Wetland_ann.sum(axis=(0, 2)).reshape(180, 1)
Wetland_total = Wetland_lat.sum()
Wetland_lat_rel = Wetland_lat / Wetland_total
# WildAnimal
WildAnimal_lat = WildAnimal_ann.sum(axis=(0, 2)).reshape(180, 1)
WildAnimal_total = WildAnimal_lat.sum()
WildAnimal_lat_rel = WildAnimal_lat / WildAnimal_total
# Landfill
Landfill_lat = Landfill_ann.sum(axis=(0, 2)).reshape(180, 1)
Landfill_total = Landfill_lat.sum()
Landfill_lat_rel = Landfill_lat / Landfill_total
# Termite
Termite_lat = Termite_ann.sum(axis=(0, 2)).reshape(180, 1)
Termite_total = Termite_lat.sum()
Termite_lat_rel = Termite_lat / Termite_total

# Plot
plt.figure(figsize=(6, 10), dpi=1000)
plt.plot(Rice_lat_rel, latitude, label='Rice')
plt.plot(Ruminant_lat_rel, latitude, label='Ruminant')
plt.plot(Wetland_lat_rel, latitude, label='Wetland')
plt.plot(Landfill_lat_rel, latitude, label='Landfill')
plt.plot(Termite_lat_rel, latitude, label='Termite')
plt.plot(WildAnimal_lat_rel, latitude, label='Wild Animal')

# Labels and title
plt.xlabel('Density')
plt.ylabel('Latitude')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.grid(False)

# Show the plot
plt.show()


#%% Load the dD-H2O map
with tifffile.TiffFile('data/d2h_MA.tif') as tif:
    # Access metadata
    metadata = tif.pages[0].tags

    # Access the image data
    image = tif.pages[0].asarray()
    ocean_value = image[78,1]
    new_rows = np.full((77, 4320), -3.4e38) # 77 is the number of 5' rows in the arctic that are missing
    image = np.vstack((new_rows, image)) #Now we have a full 180degree by 360 degree matrix at 1 degree resolution
# Mask ocean values
masked_data = np.isclose(image, ocean_value)
masked_dataflip = np.flipud(masked_data)

# Define the special value
Special = 0#-3.4e10
# Calculate the size of each block for averaging
block_size = 12
# Calculate the dimensions of the new matrix
new_rows = image.shape[0] // block_size
new_cols = image.shape[1] // block_size
# Reshape the data matrix into a 4D array with dimensions (new_rows, block_size, new_cols, block_size)
reshaped_image = image.reshape(new_rows, block_size, new_cols, block_size)
# Initialize the averaged matrix with zeros
averaged_image = np.zeros((new_rows, new_cols))
# Iterate over each block
for i in range(new_rows):
    for j in range(new_cols):
        block = reshaped_image[i, :, j, :]
        # Flatten the block to 1D array for easier processing
        block_flat = block.flatten()
        # Filter out special values
        valid_values = block_flat[block_flat >= -1000]
        if len(valid_values) > 0:
            # If there are valid values in the block, calculate their average
            averaged_image[i, j] = np.mean(valid_values)
        else:
            # If all values in the block are special values, set the average to the special value
            averaged_image[i, j] = Special
# Create a mask for special values
image_shrink = np.isclose(averaged_image, Special)
# Flip image 
averaged_image = np.flipud(averaged_image)

# Create latitudinal weighting factor to calculate hemispheric and global means
latitudes = np.linspace(-89.5, 89.5, 180)
weighting_factor = np.cos(np.radians(latitudes))
weighting_factor /= np.sum(weighting_factor)


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


#%% Determine global mean mic dD and uncertainty

# Create empty matrices to append results onto
Rice_ann_dD_MC = []
Wetland_ann_dD_MC = []
Ruminant_ann_dD_MC = []
Landfill_ann_dD_MC = []
Termite_ann_dD_MC = []
WildAnimal_ann_dD_MC = []

# Run MC analysis
Mic_dD_MC = []
for k in range(0,1000):
    #Add random error
    random_slope_wetland = np.random.normal(0, 1)
    random_intercept_wetland = np.random.normal(0, 1)
    random_slope_rice = np.random.normal(0, 1)
    random_intercept_rice = np.random.normal(0, 1)
    random_slope_landfill = np.random.normal(0, 1)
    random_intercept_landfill = np.random.normal(0, 1)
    #Calculate map for this iteration
    wetland_dD_map = (0.6088+random_slope_wetland*.072)*averaged_image-(285.7+random_intercept_wetland*6.9)
    rice_dD_map = (.5520+random_slope_rice*.191)*averaged_image-(309.7+random_intercept_rice*5.6)
    landfill_dD_map = (0.757+random_slope_landfill*.436)*averaged_image-(245.8+random_intercept_landfill*23.7)
    if np.abs((0.757 + random_slope_landfill * 0.336)) > 5:
        print(f"HULK SLOPE: {random_slope_landfill:.2f}")
    if np.abs((245.8 + random_intercept_landfill * 23.7)) > 500:
        print(f"HULK INTERCEPT: {random_intercept_landfill:.2f}")
    
    # This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
   
    # Calculate latitudinal average dD for wetland
    mask_wetland = wetland_dD_map < -1000
    # Apply the mask to replace values below -1000 with NaN (Not a Number)
    masked_data_wetland = np.where(mask_wetland, np.nan, wetland_dD_map)
    # Calculate the average of each row, excluding NaN values
    arr_wetland = np.nanmean(masked_data_wetland, axis=1)
    # Function to replace NaN values with the average of nearest non-NaN values
    def replace_nan_with_avg(arr_wetland):
        for i in range(len(arr_wetland)):
            if np.isnan(arr_wetland[i]):
                # Find nearest valid value above (backward)
                idx_above = i - 1
                while idx_above >= 0 and np.isnan(arr_wetland[idx_above]):
                    idx_above -= 1
                # Find nearest valid value below (forward)
                idx_below = i + 1
                while idx_below < len(arr_wetland) and np.isnan(arr_wetland[idx_below]):
                    idx_below += 1
                # Calculate average of nearest non-NaN values
                if idx_above >= 0 and idx_below < len(arr_wetland):
                    arr_wetland[i] = np.mean([arr_wetland[idx_above], arr_wetland[idx_below]])
                elif idx_above >= 0:
                    arr_wetland[i] = arr_wetland[idx_above]
                elif idx_below < len(arr_wetland):
                    arr_wetland[i] = arr_wetland[idx_below]
        return arr_wetland
    # Replace NaN values in the array
    arr_replaced_wetland = replace_nan_with_avg(arr_wetland)
    # Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
    wetland_dD_map[mask_wetland] = np.tile(arr_replaced_wetland, (360, 1)).T[mask_wetland]
    
    # Calculate latitudinal average dD for rice
    mask_rice = rice_dD_map < -1000
    # Apply the mask to replace values below -1000 with NaN (Not a Number)
    masked_data_rice = np.where(mask_rice, np.nan, rice_dD_map)
    # Calculate the average of each row, excluding NaN values
    arr_rice = np.nanmean(masked_data_rice, axis=1)
    # Replace NaN values in the array
    arr_replaced_rice = replace_nan_with_avg(arr_rice)
    # Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
    rice_dD_map[mask_rice] = np.tile(arr_replaced_rice, (360, 1)).T[mask_rice]
    
    # Calculate latitudinal average dD for landfill
    mask_landfill = landfill_dD_map < -1000
    # Apply the mask to replace values below -1000 with NaN (Not a Number)
    masked_data_landfill = np.where(mask_landfill, np.nan, landfill_dD_map)
    # Calculate the average of each row, excluding NaN values
    arr_landfill = np.nanmean(masked_data_landfill, axis=1)
    # Replace NaN values in the array
    arr_replaced_landfill = replace_nan_with_avg(arr_landfill)
    # Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
    landfill_dD_map[mask_landfill] = np.tile(arr_replaced_landfill, (360, 1)).T[mask_landfill]

    # Multiply by subcategory map and then sum each year and append onto matrix
    Rice_ann_dD_it = np.sum((Rice_ann_weighted*wetland_dD_map), axis=(1, 2))
    Ruminant_ann_dD_it = np.sum((Ruminant_ann_weighted*wetland_dD_map), axis=(1, 2))
    Wetland_ann_dD_it = np.sum((Wetland_ann_weighted*wetland_dD_map), axis=(1, 2))
    Landfill_ann_dD_it = np.sum((Landfill_ann_weighted*landfill_dD_map), axis=(1, 2))
    Termite_ann_dD_it = np.sum((Termite_ann_weighted*wetland_dD_map), axis=(1, 2))
    WildAnimal_ann_dD_it = np.sum((WildAnimal_ann_weighted*wetland_dD_map), axis=(1, 2))
    
    # Compile
    Rice_ann_dD_MC.append(np.sum((Rice_ann_weighted*wetland_dD_map), axis=(1, 2)))
    Ruminant_ann_dD_MC.append(np.sum((Ruminant_ann_weighted*wetland_dD_map), axis=(1, 2))) 
    Wetland_ann_dD_MC.append(np.sum((Wetland_ann_weighted*wetland_dD_map), axis=(1, 2)))
    Landfill_ann_dD_MC.append(np.sum((Landfill_ann_weighted*landfill_dD_map), axis=(1, 2)))
    Termite_ann_dD_MC.append(np.sum((Termite_ann_weighted*wetland_dD_map), axis=(1, 2)))
    WildAnimal_ann_dD_MC.append(np.sum((WildAnimal_ann_weighted*wetland_dD_map), axis=(1, 2)))
    
    # Add random uncertainty to subcategory proportions 
    perturbations = np.random.normal(loc=1, scale=0.1, size=6)
    # Combine all subcategories into matrix
    SubCatProp = np.array([Rice_Post_Frac,Ruminant_Post_Frac,Wetland_Post_Frac,Landfill_Post_Frac,Termite_Post_Frac,WildAnimal_Post_Frac]).T
    SubCatProp_MC = SubCatProp*perturbations
    SubCatProp_MC /= SubCatProp_MC.sum(axis=1, keepdims=True) # Rescale each row to 1

    #Calculate weighted mean
    Mic_dD_MCrun = SubCatProp_MC[:,0]*Rice_ann_dD_it + SubCatProp_MC[:,1]*Ruminant_ann_dD_it + SubCatProp_MC[:,2]*Wetland_ann_dD_it + SubCatProp_MC[:,3]*Landfill_ann_dD_it + SubCatProp_MC[:,4]*Termite_ann_dD_it + SubCatProp_MC[:,5]*WildAnimal_ann_dD_it
    Mic_dD_MC.append(Mic_dD_MCrun)

# Calculate total mean and uncertainty
#Mic_dD_MC = np.vstack(Mic_dD_MC)
Mic_dD_MC = np.array(Mic_dD_MC).T
# Calculate statistics
mean_dD_mic  = np.mean(Mic_dD_MC,axis=1)
meanmean = np.mean(mean_dD_mic)
stdev_dD_mic = np.std(Mic_dD_MC,axis=1)
meanstd = np.mean(stdev_dD_mic)

# Convert the matrix to an array
Rice_ann_dD_MC = np.array(Rice_ann_dD_MC)
Ruminant_ann_dD_MC = np.array(Ruminant_ann_dD_MC)
Wetland_ann_dD_MC = np.array(Wetland_ann_dD_MC)
Landfill_ann_dD_MC = np.array(Landfill_ann_dD_MC)
Termite_ann_dD_MC = np.array(Termite_ann_dD_MC)
WildAnimal_ann_dD_MC = np.array(WildAnimal_ann_dD_MC)

# Calculate mean and standard deviation for each subcategory
Rice_ann_dD_std = np.std(Rice_ann_dD_MC, axis=0)
Rice_ann_dD_mean = np.mean(Rice_ann_dD_MC, axis=0)
Ruminant_ann_dD_std = np.std(Ruminant_ann_dD_MC, axis=0)
Ruminant_ann_dD_mean = np.mean(Ruminant_ann_dD_MC, axis=0)
Wetland_ann_dD_std = np.std(Wetland_ann_dD_MC, axis=0)
Wetland_ann_dD_mean = np.mean(Wetland_ann_dD_MC, axis=0)
Landfill_ann_dD_std = np.std(Landfill_ann_dD_MC, axis=0)
Landfill_ann_dD_mean = np.mean(Landfill_ann_dD_MC, axis=0)
Termite_ann_dD_std = np.std(Termite_ann_dD_MC, axis=0)
Termite_ann_dD_mean = np.mean(Termite_ann_dD_MC, axis=0)
WildAnimal_ann_dD_std = np.std(WildAnimal_ann_dD_MC, axis=0)
WildAnimal_ann_dD_mean = np.mean(WildAnimal_ann_dD_MC, axis=0)

# Calculate global mean microbial over time
dD_Mic_total = (Rice_Post_Frac * Rice_ann_dD_mean + Ruminant_Post_Frac * Ruminant_ann_dD_mean +
                Landfill_Post_Frac * Landfill_ann_dD_mean + Wetland_Post_Frac * Wetland_ann_dD_mean +
                Termite_Post_Frac * Termite_ann_dD_mean + WildAnimal_Post_Frac * WildAnimal_ann_dD_mean)
# Calculate weighted mean uncertainty
# Stack the arrays along a new axis to create 24x6 matrices
stds = np.column_stack([Rice_ann_dD_std, Ruminant_ann_dD_std, Wetland_ann_dD_std, Landfill_ann_dD_std, Termite_ann_dD_std, WildAnimal_ann_dD_std])
means = np.column_stack([Rice_ann_dD_mean, Ruminant_ann_dD_mean, Wetland_ann_dD_mean, Landfill_ann_dD_mean, Termite_ann_dD_mean, WildAnimal_ann_dD_mean])
weights = np.column_stack([Rice_Post_Frac, Ruminant_Post_Frac, Wetland_Post_Frac, Landfill_Post_Frac, Termite_Post_Frac, WildAnimal_Post_Frac])
# Calculate the weighted mean uncertainty
dD_Mic_total_unc = np.sqrt(np.sum(weights * stds**2, axis=1) / np.sum(weights, axis=1))


#%% plot the results

# Create an array of years from 1998 to 2021

plt.rc('font', size=14) 
years = np.arange(1998, 2022).reshape(24, 1)

# Plot
plt.figure(figsize=(6, 10), dpi=1000)
plt.plot(years, Ruminant_ann_dD_mean, color='blue', label='Ruminants', linewidth=1.5)
plt.plot(years, Wetland_ann_dD_mean, color='green', label='Wetlands', linewidth=1.5)
plt.plot(years, Landfill_ann_dD_mean, color='red', label='Waste', linewidth=1.5)
plt.plot(years, Rice_ann_dD_mean, color='cyan', label='Rice', linewidth=1.5)
plt.plot(years, Rice_ann_dD_mean, color='cyan', label='Rice', linewidth=1.5)
plt.plot(years, Termite_ann_dD_mean , color='orange', label='Termites', linewidth=1.5)
plt.plot(years, WildAnimal_ann_dD_mean , color='gray', label='Wild Animals', linewidth=1.5)
plt.plot(years, mean_dD_mic, color='black', label='All Mic', linewidth=4)

# Labels and title
plt.xlabel('Year')
plt.ylabel('Global Annual Mean ${\delta}D$-CH$_4$ (‰)')
plt.title('Microbial ${\delta}D$-CH$_4$ (‰)')
plt.grid(False)
plt.ylim(-340, -265)

# Show the plot
plt.show()

#%% Save global annual

#Create result matrix
dD_Mic_results = np.column_stack([years, mean_dD_mic,dD_Mic_total_unc])
# Specify the file path
file_path = 'Output/Mic_dD_AnnGlob.csv'
# Save the array to a CSV file
np.savetxt(file_path, dD_Mic_results , delimiter=',')

#%% Main Figure plot

import numpy as np
import matplotlib.pyplot as plt

# Sample data
categories = ['Wetlands', 'Rice', 'Ruminants', 'Wild Animals', 'Termites', 'Landfills', 'Total']
SrcSigs = np.array([Wetland_ann_dD_mean[-1], Rice_ann_dD_mean[-1], Ruminant_ann_dD_mean[-1], Termite_ann_dD_mean[-1], WildAnimal_ann_dD_mean[-1], Landfill_ann_dD_mean[-1], mean_dD_mic[-1]])  # Below the line
Proportions = np.array([np.mean(Wetland_Post_Frac), np.mean(Rice_Post_Frac), np.mean(Ruminant_Post_Frac), 
                        np.mean(WildAnimal_Post_Frac), np.mean(Termite_Post_Frac), np.mean(Landfill_Post_Frac), 1])  # Above the line
uncertainty_SrcSigs = np.array([Wetland_ann_dD_std[-1], Rice_ann_dD_std[-1], Ruminant_ann_dD_std[-1], Termite_ann_dD_std[-1], WildAnimal_ann_dD_std[-1], Landfill_ann_dD_std[-1], stdev_dD_mic[-1]])  # Below the line
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
    ax2.bar(x[i], SrcSigs[i] - (-240), yerr=uncertainty_SrcSigs[i], capsize=5, 
            width=bar_width, color=color, edgecolor='black', linestyle='dashed', hatch='//', bottom=-240)


# Set separate y-axis limits
ax1.set_ylim(-1, 1)  
ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_ylim(-340, -140)
ax2.set_yticks([-250, -270, -290, -310, -330])

# Adjust positioning to make 0 to 1 take up half of the figure
ax1.set_position([0.15, 0.55, 0.7, 0.4])  # Move ax1 higher and shrink height
ax2.set_position([0.15, 0.05, 0.7, 0.4])  # Move ax2 lower and shrink height

# Formatting
ax1.axhline(0, color='black', linewidth=1)  # Reference line at y=0
ax1.set_xticks(x)
ax1.set_xticklabels(categories, rotation=20)
ax1.set_ylabel("Proportions", color='black')
ax2.set_ylabel("${\delta}D$-CH$_4$ (‰)", color='black')
ax1.set_xlabel("Microbial Emission Category")
ax1.set_title("Microbial ${\delta}D$-CH$_4$ for 2022")
ax1.yaxis.set_label_coords(-0.14, 0.75)  # Moves ax1 label higher
ax2.yaxis.set_label_coords(1.14, 0.25)  # Moves ax2 label lower

# Custom legend
handles = [plt.Rectangle((0,0),1,1, color=color) for color in colors]
# ax1.legend(handles, categories, title="Categories", loc='upper right')

plt.show()


#%% Calculate trend statistics

import numpy as np, statsmodels.api as sm 
print(sm.OLS(Wetland_ann_dD_mean, sm.add_constant(np.arange(len(Wetland_ann_dD_mean)))).fit().summary())
print(sm.OLS(mean_dD_mic, sm.add_constant(np.arange(len(dD_Mic_total)))).fit().summary())


#%% Create MC array of mic scenarios with varying trend and varying mean value

# First calculate relative trend in mean
Mic_trend = mean_dD_mic - np.mean(mean_dD_mic)
# Calculate mean mean and mean standard deviation
Mic_stdev_mean = np.mean(stdev_dD_mic)
Mic_mean_mean = np.mean(mean_dD_mic)


# Initialize an empty array to store results (24 rows, 1000 cols)
Mic_dD_MC_trend = np.zeros((24, 1000))
rMic_dD_MC_trend = np.zeros((24, 1000))
for i in range(1000):
    # Generate random numbers
    rn = np.random.normal(loc=0, scale=0.25)  # Single random number
    rn2 = np.random.normal(loc=0, scale=Mic_stdev_mean)  # 24-element vector
    
    # Create MC mic trends
    Mic_mean_MC = Mic_trend * rn + Mic_trend + Mic_mean_mean + rn2  # Shape (24,)
    rMic_mean_MC = Mic_trend * rn + Mic_trend # Shape (24,)
    
    # Store result in the i-th column
    Mic_dD_MC_trend[:, i] = Mic_mean_MC
    rMic_dD_MC_trend[:, i] = rMic_mean_MC
    
    
#%% Plot and save data

plt.rc('font', size=14) 
fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=1000, sharex=True)  # 2 rows, 1 column
# First subplot: Relative Mic δD-CH4
plt.suptitle('Microbial ${\delta}D$ Simulations')
axes[0].plot(years, rMic_dD_MC_trend)
axes[0].set_ylabel('Relative Mic ${\delta}D$-CH$_4$ (‰)')
axes[0].grid(False)
# Second subplot: Mic δD-CH4
axes[1].plot(years, Mic_dD_MC_trend)
axes[1].set_ylabel('Mic ${\delta}D$-CH$_4$ (‰)')
axes[1].grid(False)
axes[1].set_xlabel('Year')  # Shared x-axis label
plt.tight_layout()  # Adjust layout for better spacing
plt.show()


# Compile and save as .csv
Mic_MC = np.column_stack((years, Mic_dD_MC_trend))
file_path = 'Output/Mic_dD_MC.csv'
np.savetxt(file_path, Mic_MC, delimiter=',')



