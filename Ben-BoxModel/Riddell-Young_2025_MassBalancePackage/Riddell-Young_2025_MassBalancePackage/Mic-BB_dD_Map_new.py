#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:40:00 2024

@author: ryoung
"""
# This code creates 1x1 degree maps of BB and microbial dD-CH4 source signatures. 
# It also calculates global mean and uncertainty for BB dD.
# The code can be modified to generate microbial subcategory maps with new regressions if needed.

locals().clear()
import tifffile
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import netCDF4 as nc
import xarray as xr
plt.clf()
plt.close('all')

# Open the TIFF file
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
Special = -3.4e38
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


# Create latitudinal weighting factor to calculate hemispheric and global means
latitudes = np.linspace(-89.5, 89.5, 180)
weighting_factor = np.cos(np.radians(latitudes))
weighting_factor /= np.sum(weighting_factor)

#%% Calculate the expected dD-CH4 of microbial emissions based on Douglas et al., 2021

# Calculate map
dD = .6088*averaged_image-285.7
dDflip = np.flipud(dD)


# This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
# Calculate latitudinal average dD 
mask = dDflip < -1000
# Apply the mask to replace values below -1000 with NaN (Not a Number)
masked_data = np.where(mask, np.nan, dDflip)
# Calculate the average of each row, excluding NaN values
arr = np.nanmean(masked_data, axis=1)
# Function to replace NaN values with the average of nearest non-NaN values
def replace_nan_with_avg(arr):
    for i in range(len(arr)):
        if np.isnan(arr[i]):
            # Find nearest valid value above (backward)
            idx_above = i - 1
            while idx_above >= 0 and np.isnan(arr[idx_above]):
                idx_above -= 1
            # Find nearest valid value below (forward)
            idx_below = i + 1
            while idx_below < len(arr) and np.isnan(arr[idx_below]):
                idx_below += 1
            # Calculate average of nearest non-NaN values
            if idx_above >= 0 and idx_below < len(arr):
                arr[i] = np.mean([arr[idx_above], arr[idx_below]])
            elif idx_above >= 0:
                arr[i] = arr[idx_above]
            elif idx_below < len(arr):
                arr[i] = arr[idx_below]
    return arr
# Replace NaN values in the array
arr_replaced = replace_nan_with_avg(arr)
# Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
dDflip[mask] = np.tile(arr_replaced, (360, 1)).T[mask]

# Plot the grid for microbial
# Set the minimum and maximum values for the color bar scale
vmin = -500  # Minimum value
vmax = -250  # Maximum value
lower_left_lat = -90
lower_left_lon = -180
upper_right_lat = 90
upper_right_lon = 180
# Set the minimum and maximum values for the color bar scale that is the ocean
vminB = -100  # Minimum value
vmaxB = 0  # Maximum value
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(dDflip, cmap= 'ocean', interpolation='nearest', norm=Normalize(vmin=vmin, vmax=vmax))
plt.set_cmap('ocean')
plt.set_cmap('ocean')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('Mic ${\delta}D$-CH$_4$ (‰)', fontsize=8)
# Uncomment to cover oceans
#m.imshow(np.where(masked_dataflip, 0, np.nan), cmap='gray', norm=Normalize(vmin=vminB, vmax=vmaxB))
m.drawcoastlines(linewidth=.4,color='k')


#%% MC error propagation to get grid level uncertainty for microbial emissions

# TO CHANGE WHICH SLOPE YOU CALCULATE FOR, CHANGE LINES 165, 221 and 222 
# Define Wetland regression characteristics (Douglas et al., 2021)
Slope_wetland = .6088
Slope_Unc_wetland = .072
Intercept_wetland = 285.7
Incercept_Unc_wetland = 6.9
# Define Rice regression characteristics
Slope_rice = .552
Slope_Unc_rice  = .191
Intercept_rice  = 309.7
Incercept_Unc_rice  = 5.6
# Define Lake regression characteristics
Slope_lake = .4215
Slope_Unc_lake = .069
Intercept_lake = 294
Incercept_Unc_lake = 7
# Define Landfill regression characteristics
Slope_landfill = .757
Slope_Unc_landfill = .436
Intercept_landfill = 245.8
Incercept_Unc_landfill = 23.7
# Define All Mic regression characteristics
Slope_AllMic = .519
Slope_Unc_AllMic = .047
Intercept_AllMic = 278.3
Incercept_Unc_AllMic = 4.1

# Calculate map
dD = Slope_landfill*averaged_image-Intercept_landfill
dDflip = np.flipud(dD)


# This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
# Calculate latitudinal average dD 
mask = dDflip < -1000
# Apply the mask to replace values below -1000 with NaN (Not a Number)
masked_data = np.where(mask, np.nan, dDflip)
# Calculate the average of each row, excluding NaN values
arr = np.nanmean(masked_data, axis=1)
# Function to replace NaN values with the average of nearest non-NaN values
def replace_nan_with_avg(arr):
    for i in range(len(arr)):
        if np.isnan(arr[i]):
            # Find nearest valid value above (backward)
            idx_above = i - 1
            while idx_above >= 0 and np.isnan(arr[idx_above]):
                idx_above -= 1
            # Find nearest valid value below (forward)
            idx_below = i + 1
            while idx_below < len(arr) and np.isnan(arr[idx_below]):
                idx_below += 1
            # Calculate average of nearest non-NaN values
            if idx_above >= 0 and idx_below < len(arr):
                arr[i] = np.mean([arr[idx_above], arr[idx_below]])
            elif idx_above >= 0:
                arr[i] = arr[idx_above]
            elif idx_below < len(arr):
                arr[i] = arr[idx_below]
    return arr
# Replace NaN values in the array
arr_replaced = replace_nan_with_avg(arr)
# Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
dDflip[mask] = np.tile(arr_replaced, (360, 1)).T[mask]

# # This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
# # Calculate latitudinal average dD 
# mask = dDflip < -1000
# # Apply the mask to replace values below -1000 with NaN (Not a Number)
# masked_data = np.where(mask, np.nan, dDflip)
# # Calculate the average of each row, excluding NaN values
# arr = np.nanmean(masked_data, axis=1)
# # Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
# dDflip[mask] = np.tile(arr_replaced, (360, 1)).T[mask]

# Initialize an empty list to store matrices
matrices = []
# Run MC analysis
for k in range(0,1000):
    #Add random error
    random_gaussian1 = np.random.normal(0, 1)
    random_gaussian2 = np.random.normal(0, 1)
    dDmic_Unc = (Slope_landfill+random_gaussian1*Slope_Unc_landfill)*averaged_image-(Intercept_landfill+random_gaussian2*Incercept_Unc_landfill)
    dDmicflip_Unc = np.flipud(dDmic_Unc)
    # Append to the list
    matrices.append(dDmicflip_Unc)
    
# Convert the list of matrices to a 3D NumPy array
matrices_array = np.array(matrices)
# Calculate standard deviation for each cell
std_dev_matrix_mic = np.std(matrices_array, axis=0)

# This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
# Calculate latitudinal average dD 
mask = std_dev_matrix_mic > 1000
# Apply the mask to replace values below -1000 with NaN (Not a Number)
masked_data = np.where(mask, np.nan, std_dev_matrix_mic)
# Calculate the average of each row, excluding NaN values
arr = np.nanmean(masked_data, axis=1)
# Function to replace NaN values with the average of nearest non-NaN values
def replace_nan_with_avg(arr):
    for i in range(len(arr)):
        if np.isnan(arr[i]):
            # Find nearest valid value above (backward)
            idx_above = i - 1
            while idx_above >= 0 and np.isnan(arr[idx_above]):
                idx_above -= 1
            # Find nearest valid value below (forward)
            idx_below = i + 1
            while idx_below < len(arr) and np.isnan(arr[idx_below]):
                idx_below += 1
            # Calculate average of nearest non-NaN values
            if idx_above >= 0 and idx_below < len(arr):
                arr[i] = np.mean([arr[idx_above], arr[idx_below]])
            elif idx_above >= 0:
                arr[i] = arr[idx_above]
            elif idx_below < len(arr):
                arr[i] = arr[idx_below]
    return arr
# Replace NaN values in the array
arr_replaced = replace_nan_with_avg(arr)
# Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
std_dev_matrix_mic[mask] = np.tile(arr_replaced, (360, 1)).T[mask]


# # This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
# # mask values > 1000
# mask_Unc = std_dev_matrix_mic > 1000
# # replace with nan
# masked_data = np.where(mask_Unc, np.nan, std_dev_matrix_mic)
# # row means (ignoring nans), shape (180,)
# row_means = np.nanmean(masked_data, axis=1)
# # broadcast row means to full 180x360 grid
# row_means_full = np.tile(row_means[:, np.newaxis], (1, std_dev_matrix_mic.shape[1]))
# # fill masked values with row mean
# std_dev_matrix_mic[mask_Unc] = row_means_full[mask_Unc]
# # replace rows that are still all nan with average of nearest rows above/below
# for i in range(std_dev_matrix_mic.shape[0]):
#     if np.all(np.isnan(std_dev_matrix_mic[i])):
#         if 0 < i < std_dev_matrix_mic.shape[0] - 1:
#             std_dev_matrix_mic[i] = np.nanmean([std_dev_matrix_mic[i-1], std_dev_matrix_mic[i+1]], axis=0)
# Flip it
std_dev_matrix_mic_90N = np.flipud(std_dev_matrix_mic)
dD_90N = np.flipud(dDflip)


# Plot test
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(dD_90N, cmap= 'ocean', interpolation='nearest', norm=Normalize(vmin=-500, vmax=-250), origin="upper")
plt.set_cmap('ocean')
plt.set_cmap('ocean')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('Mic ${\delta}D$-CH$_4$ (‰)', fontsize=8)
# Uncomment to cover oceans
# m.imshow(np.where(masked_dataflip, 0, np.nan), cmap='gray', norm=Normalize(vmin=vminB, vmax=vmaxB))
m.drawcoastlines(linewidth=.4,color='k')

# Export!
# Specify the file path
file_path = 'Output/Landfill_dD_1x1.txt'
file_path2 = 'Output/Landfill_dD_1x1_Unc.txt'
# Save the array to a CSV file
np.savetxt(file_path, dD_90N, delimiter=',')
np.savetxt(file_path2, std_dev_matrix_mic_90N, delimiter=',')


#%% Calculate the dD-CH4 of biomass burning emissions based on Umezawa et al., 2011 regression

# Define regression (Umezawa et al., 2011)
dDBB = 1.16*averaged_image-177
dDBBflip = np.flipud(dDBB)

# This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
# Calculate latitudinal average dD 
mask = dDBBflip < -1000
# Apply the mask to replace values below -1000 with NaN (Not a Number)
masked_data = np.where(mask, np.nan, dDBBflip)
# Calculate the average of each row, excluding NaN values
arr = np.nanmean(masked_data, axis=1)
# Function to replace NaN values with the average of nearest non-NaN values
def replace_nan_with_avg(arr):
    for i in range(len(arr)):
        if np.isnan(arr[i]):
            # Find nearest valid value above (backward)
            idx_above = i - 1
            while idx_above >= 0 and np.isnan(arr[idx_above]):
                idx_above -= 1
            # Find nearest valid value below (forward)
            idx_below = i + 1
            while idx_below < len(arr) and np.isnan(arr[idx_below]):
                idx_below += 1
            # Calculate average of nearest non-NaN values
            if idx_above >= 0 and idx_below < len(arr):
                arr[i] = np.mean([arr[idx_above], arr[idx_below]])
            elif idx_above >= 0:
                arr[i] = arr[idx_above]
            elif idx_below < len(arr):
                arr[i] = arr[idx_below]
    return arr
# Replace NaN values in the array
arr_replaced = replace_nan_with_avg(arr)
# Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
dDBBflip[mask] = np.tile(arr_replaced, (360, 1)).T[mask]

# Plot the grid for BB
# Set the minimum and maximum values for the color bar scale
vmin = -500  # Minimum value
vmax = -150  # Maximum value
lower_left_lat = -90
lower_left_lon = -180
upper_right_lat = 90
upper_right_lon = 180
# Set the minimum and maximum values for the color bar scale that is the ocean
vminB = -100  # Minimum value
vmaxB = 0  # Maximum value
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(dDBBflip, cmap= 'inferno', interpolation='nearest', norm=Normalize(vmin=vmin, vmax=vmax))
plt.set_cmap('inferno')
plt.set_cmap('inferno')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('BB ${\delta}D$-CH$_4$ (‰)', fontsize=8)
# Uncomment to cover oceans
#m.imshow(np.where(masked_dataflip, 0, np.nan), cmap='gray', norm=Normalize(vmin=vminB, vmax=vmaxB))
m.drawcoastlines(linewidth=.4,color='k')

#%% MC error propagation to get grid level uncertainty

# Define regression uncertainties (Umezawa et al., 2011)
Slope_Unc = 0.09
Incercept_Unc = 6.5

# Initialize an empty list to store matrices
matrices = []
for k in range(0,1000):
    #Add random error
    random_gaussian1 = np.random.normal(0, 1)
    random_gaussian2 = np.random.normal(0, 1)
    dDBB_Unc = (1.16+random_gaussian1*Slope_Unc)*averaged_image-(177+random_gaussian2*Incercept_Unc)
    dDBBflip_Unc = np.flipud(dDBB_Unc)
    # Append to the list
    matrices.append(dDBBflip_Unc)
    
# Convert the list of matrices to a 3D NumPy array
matrices_array = np.array(matrices)
# Calculate standard deviation for each cell
std_dev_matrix = np.std(matrices_array, axis=0)

# This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
# Calculate latitudinal average dD 
mask = std_dev_matrix > 1000
# Apply the mask to replace values below -1000 with NaN (Not a Number)
masked_data = np.where(mask, np.nan, std_dev_matrix)
# Calculate the average of each row, excluding NaN values
arr = np.nanmean(masked_data, axis=1)
# Function to replace NaN values with the average of nearest non-NaN values
def replace_nan_with_avg(arr):
    for i in range(len(arr)):
        if np.isnan(arr[i]):
            # Find nearest valid value above (backward)
            idx_above = i - 1
            while idx_above >= 0 and np.isnan(arr[idx_above]):
                idx_above -= 1
            # Find nearest valid value below (forward)
            idx_below = i + 1
            while idx_below < len(arr) and np.isnan(arr[idx_below]):
                idx_below += 1
            # Calculate average of nearest non-NaN values
            if idx_above >= 0 and idx_below < len(arr):
                arr[i] = np.mean([arr[idx_above], arr[idx_below]])
            elif idx_above >= 0:
                arr[i] = arr[idx_above]
            elif idx_below < len(arr):
                arr[i] = arr[idx_below]
    return arr
# Replace NaN values in the array
arr_replaced = replace_nan_with_avg(arr)
# Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
std_dev_matrix[mask] = np.tile(arr_replaced, (360, 1)).T[mask]

# Flip so 90N is row 1
std_dev_matrix_90N = np.flipud(std_dev_matrix)
dDBB_90N = np.flipud(dDBBflip)

# Plot standard deviation
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(std_dev_matrix, cmap= 'inferno', interpolation='nearest', norm=Normalize(vmin=0, vmax=30))
plt.set_cmap('inferno')
plt.set_cmap('inferno')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('BB ${\delta}D$-CH$_4$ (‰)', fontsize=8)
#Ucomment to cover oceans
#m.imshow(np.where(masked_dataflip, 0, np.nan), cmap='gray', norm=Normalize(vmin=vminB, vmax=vmaxB))
m.drawcoastlines(linewidth=.4,color='k')

#%% Do a MC average to estiamte global mean and uncertainty
    
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
# Sum data annally and weight
GFED5_ann_sum = GFED5_ann.sum(axis=(1, 2))
GFED5_ann_sum = GFED5_ann_sum[:, np.newaxis, np.newaxis] # Prep ann sum for the multiplication
GFED5_ann_weight = GFED5_ann  / GFED5_ann_sum


# Load CTCH4 Prior emissions
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


# Load emissions data
# First load CT-CH4 data to get weighted matrix of mean annual Mic emissions
# Load Monthly emissions (Carbon Tracker Methane)
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
pyroA = flux_data.variables['pyrogenic_flux'][:]
# expand to 1x1 grid
BB_exp = np.repeat(pyroA, repeats=2, axis=1)
BB_exp = np.repeat(BB_exp, repeats=3, axis=2)/6 
# Convert from kg/s to Tg/month
# Constants
kpt = 1000000000  # Kg per Tg
spm = 2.628e+6  # Sec per Month
BB = BB_exp / kpt * spm
# Flip on latitudinal axis
BB = np.flip(BB, axis=1)
# Calculate the annual sum and weight each year by its respective total
BB_ann = BB_exp.reshape(24, 12, 180, 360).sum(axis=1)
BB_ann_sum = BB_ann.sum(axis=(1, 2))
BB_ann_sum = BB_ann_sum[:, np.newaxis, np.newaxis] # Prep ann sum for the multiplication
BB_ann_weight = BB_ann  / BB_ann_sum
BB_ann_weight_flipped = BB_ann_weight[:, ::-1, :] # Flip it
# Calculate average annual emissions
BB_sum = (BB[:, :, :].sum(axis=(0)))/24
BB_total = BB_sum[:, :].sum(axis=(0,1))
# Create weighted matrix
BB_weight = BB_sum / BB_total 
BB_weight = np.flipud(BB_weight)

# Run MC analysis
Results = []
Results_ann = []
Results_ann_GFED5 = []
Results_ann_prior = []
for k in range(0,1000):
    #Add random error
    random_gaussian1 = np.random.normal(0, 1)
    random_gaussian2 = np.random.normal(0, 1)
    dDBB_Map_MC = (1.16+random_gaussian1*Slope_Unc)*averaged_image-(177+random_gaussian2*Incercept_Unc) # Regression from umezawa et al., 2011
    dDBB_Map_MC_flip = np.flipud(dDBB_Map_MC)
    # This next portion of code provides latitudinal averages to ocean cells. This ensures that any cell with out an estiamte of dDmic that is multiplied by a flux is multiplied by a reasonable value
    # Calculate latitudinal average dD 
    mask_Unc = dDBB_Map_MC_flip < -1000
    # Apply the mask to replace values below -1000 with NaN (Not a Number)
    masked_data = np.where(mask_Unc, np.nan, dDBB_Map_MC_flip)
    # Calculate the average of each row, excluding NaN values
    arr = np.nanmean(masked_data, axis=1)
    # Use broadcasting to replace values in arr_180x360 with values from arr_180_rows
    dDBB_Map_MC_flip[mask_Unc] = np.tile(arr_replaced, (360, 1)).T[mask_Unc]
    dDBB_Map_MC_weight = dDBB_Map_MC_flip*BB_weight
    dDBB_globMC = dDBB_Map_MC_weight [:, :].sum(axis=(0,1))
    # Append to the list
    Results.append(dDBB_globMC)
    # Do the same but calculate annual averages
    dDBB_Map_MC_flip[mask_Unc] = np.tile(arr_replaced, (360, 1)).T[mask_Unc]
    dDBB_Map_MC_flip_ann = dDBB_Map_MC_flip[np.newaxis, :, :] # Prep dDBB_Map_MC_flipt for the multiplication
    dDBB_Map_MC_weight_ann = dDBB_Map_MC_flip_ann*BB_ann_weight
    dDBB_globMC_ann = dDBB_Map_MC_weight_ann[:, :].sum(axis=(1,2))
    # Append to the list
    Results_ann.append(dDBB_globMC_ann)
    # Repeat for GFED emissions
    dDGFED5_Map_MC_weight_ann = dDBB_Map_MC_flip_ann*GFED5_ann_weight
    dDGFED5_globMC_ann = dDGFED5_Map_MC_weight_ann[:, :].sum(axis=(1,2))
    Results_ann_GFED5.append(dDGFED5_globMC_ann)
    # Repeat for BB CTCH4 prior (GFED4)
    dDprior_Map_MC_weight_ann = dDBB_Map_MC_flip_ann*BB_prior_ann_weight
    dDprior_globMC_ann = dDprior_Map_MC_weight_ann[:, :].sum(axis=(1,2))
    Results_ann_prior.append(dDprior_globMC_ann)
    

Results_ann = np.stack(Results_ann, axis=1)
Results_ann_GFED5 = np.stack(Results_ann_GFED5, axis=1)
Results_ann_prior = np.stack(Results_ann_prior, axis=1)
# Calculate averages 
Glob_BB_stdev = np.std(Results)
Glob_BB_mean = np.mean(Results)
# Calculate annual averages 
ann_means = Results_ann.mean(axis=1)
ann_stds = Results_ann.std(axis=1)
ann_means_GFED5 = Results_ann_GFED5.mean(axis=1)
ann_stds_GFED5 = Results_ann_GFED5.std(axis=1)
ann_means_prior = Results_ann_prior.mean(axis=1)
ann_stds_prior = Results_ann_prior.std(axis=1)
mean_all = Results_ann.mean(axis=(0, 1))
std_all = Results_ann.std(axis=(0, 1))

#%% Export annual averages and plot

years = np.linspace(1998, 2021, 24)
BB_dD = np.vstack((np.column_stack((years, ann_means, ann_stds)), [2030, mean_all, std_all]))
file_path = 'Output/BB_dD_annual.csv'
np.savetxt(file_path, BB_dD, delimiter=',')

# Plotting the data
plt.rcParams.update({'font.size': 14})  # Adjust this value as needed
plt.figure(figsize=(10, 6), dpi=1000)
plt.plot(years, ann_means, color='blue',linewidth=3,label='CTCH4')
#plt.plot(years[4:-1], ann_means_GFED5[4:-1], color='orange',label='GFED5')
plt.plot(years, ann_means_prior, color='orange',label='GFED4/prior')
plt.errorbar((years[-1]+2), np.mean(mean_all), yerr=np.mean(std_all), fmt='s', color='blue', markersize=8, capsize=5, capthick=1.5, elinewidth=1.5,label='Mean ± 1$\sigma$')
plt.xlabel('Year')
plt.legend()
plt.ylabel('Pyrogenic ${\delta}D$-CH$_4$ (‰)')
plt.show()


#%% Export!

# Specify the file path
file_path = 'Output/BB_dD_1x1.txt'
file_path_Unc = 'Output/BB_dD_1x1_Unc.txt'
# Save the array to a CSV file
np.savetxt(file_path, dDBB_90N, delimiter=',')
np.savetxt(file_path_Unc, std_dev_matrix_90N, delimiter=',')


#%% Perform statistical analysis 

import numpy as np, statsmodels.api as sm 

# Exclude 2021 due to anamolous fires
print(sm.OLS(ann_means[:-1], sm.add_constant(np.arange(len(ann_means[:-1])))).fit().summary())

