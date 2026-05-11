#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:15:23 2024

@author: ryoung
"""

# Code to convert all files to netcdf files

locals().clear()
import numpy as np
from mpl_toolkits.basemap import Basemap
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from netCDF4 import Dataset
plt.clf()
plt.close('all')


#%% Load data

#Fossil
# First ONG dD
ONG_dD_stderr = pd.read_csv('Output/ONGStdErr_dD_1x1.csv', header = None);
ONG_dD_stderr = np.flipud(ONG_dD_stderr.values)
ONG_dD_stdev = pd.read_csv('Output/ONGStdDev_dD_1x1.csv', header = None);
ONG_dD_stdev = np.flipud(ONG_dD_stdev.values)
ONG_dD_N = pd.read_csv('Output/ONGN_dD_1x1.csv', header = None);
ONG_dD_N = np.flipud(ONG_dD_N.values)
# THen Coal dD
Coal_dD_stderr = pd.read_csv('Output/CoalStdErr_dD_1x1.csv', header = None);
Coal_dD_stderr = np.flipud(Coal_dD_stderr.values)
Coal_dD_stdev = pd.read_csv('Output/CoalStdDev_dD_1x1.csv', header = None);
Coal_dD_stdev = np.flipud(Coal_dD_stdev.values)
Coal_dD_N = pd.read_csv('Output/CoalN_dD_1x1.csv', header = None);
Coal_dD_N = np.flipud(Coal_dD_N.values)
Coal_dD = pd.read_csv('Output/Coal_dD_1x1.csv', header = None);
Coal_dD = np.flipud(Coal_dD.values)
# Then ONG d13C
ONG_d13C_stderr = pd.read_csv('Output/ONGStdErr_d13C_1x1.csv', header = None);
ONG_d13C_stderr = np.flipud(ONG_d13C_stderr.values)
ONG_d13C_stdev = pd.read_csv('Output/ONGStdDev_d13C_1x1.csv', header = None);
ONG_d13C_stdev = np.flipud(ONG_d13C_stdev.values)
ONG_d13C_N = pd.read_csv('Output/ONGN_d13C_1x1.csv', header = None);
ONG_d13C_N = np.flipud(ONG_d13C_N.values)
# Then Coal d13C
Coal_d13C_stderr = pd.read_csv('Output/CoalStdErr_d13C_1x1.csv', header = None);
Coal_d13C_stderr = np.flipud(Coal_d13C_stderr.values)
Coal_d13C_stdev = pd.read_csv('Output/CoalStdDev_d13C_1x1.csv', header = None);
Coal_d13C_stdev = np.flipud(Coal_d13C_stdev.values)
Coal_d13C_N = pd.read_csv('Output/CoalN_d13C_1x1.csv', header = None);
Coal_d13C_N = np.flipud(Coal_d13C_N.values)
Coal_d13C = pd.read_csv('Output/Coal_d13C_1x1.csv', header = None);
Coal_d13C = np.flipud(Coal_d13C.values)

# Biomass Burning
BB_dD = pd.read_csv('Output/BB_dD_1x1.csv', header = None);
BB_dD  = BB_dD.values
BB_dD_Unc = pd.read_csv('Output/BB_dD_1x1_Unc.csv', header = None);
BB_dD_Unc  = BB_dD_Unc.values

# Microbial
# First wetland
Wetland_dD = pd.read_csv('Output/Wetland_dD_1x1.csv', header = None);
Wetland_dD  = Wetland_dD.values
Wetland_dD_Unc = pd.read_csv('Output/Wetland_dD_1x1_Unc.csv', header = None);
Wetland_dD_Unc  = Wetland_dD_Unc.values
# Then Lake
Lake_dD = pd.read_csv('Output/Lake_dD_1x1.csv', header = None);
Lake_dD  = Lake_dD.values
Lake_dD_Unc = pd.read_csv('Output/Lake_dD_1x1_Unc.csv', header = None);
Lake_dD_Unc  = Lake_dD_Unc.values
# Then Rice
Rice_dD = pd.read_csv('Output/Rice_dD_1x1.csv', header = None);
Rice_dD  = Rice_dD.values
Rice_dD_Unc = pd.read_csv('Output/Rice_dD_1x1_Unc.csv', header = None);
Rice_dD_Unc  = Rice_dD_Unc.values
# Then Landfill
Landfill_dD = pd.read_csv('Output/Landfill_dD_1x1.csv', header = None);
Landfill_dD  = Landfill_dD.values
Landfill_dD_Unc = pd.read_csv('Output/Landfill_dD_1x1_Unc.csv', header = None);
Landfill_dD_Unc  = Landfill_dD_Unc.values
# Then All Microbial
AllMic_dD = pd.read_csv('Output/AllMic_dD_1x1.csv', header = None);
AllMic_dD  = AllMic_dD.values
AllMic_dD_Unc = pd.read_csv('Output/AllMic_dD_1x1_Unc.csv', header = None);
AllMic_dD_Unc  = AllMic_dD_Unc.values


#%% Plot to test dataset

from mpl_toolkits.basemap import Basemap
# Set the minimum and maximum values for the color bar scale
lower_left_lat = -90
lower_left_lon = -180
upper_right_lat = 90
upper_right_lon = 180
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
# plot ONG dD stderr 
plt.figure(dpi=1000)
m.imshow(Coal_d13C, cmap= 'YlGnBu', interpolation='nearest')
plt.set_cmap('YlGnBu')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
plt.title('ONG ${\delta}D$-CH$_4$ (‰) Uncertainty', fontsize=8)
# Uncomment to cover oceans
#m.imshow(np.where(masked_dataflip, 0, np.nan), cmap='gray', norm=Normalize(vmin=vminB, vmax=vmaxB))
m.drawcoastlines(linewidth=.4,color='k')
m.drawcountries(linewidth=.4,color='k')


#%% Now create large netcdf with temporal ONG maps for dD

# Define latitude and longitude ranges
latitudes = np.linspace(-89.5, 89.5, 180)
longitudes = np.linspace(-179.5, 179.5, 360)
years = np.arange(1970, 2022)  # Years from 1970 to 2021

# Load each CSV file as a NumPy array and flip the matrix vertically
def load_and_process_data(years, prefix='Output/ONG_dD_1x1_', suffix='.csv'):
    data_array = []
    for year in years:
        filename = f'{prefix}{year}{suffix}'
        df = pd.read_csv(filename, header=None)
        array_data = df.to_numpy()
        flipped_array = np.flipud(array_data)
        data_array.append(flipped_array)
    return np.array(data_array)

# Create a NetCDF file with a time dimension
def create_netcdf_with_time_dimension(filename, data_array, description, units):
    with Dataset(filename, 'w', format='NETCDF4') as ds:
        # Create dimensions
        ds.createDimension('lat', len(latitudes))
        ds.createDimension('lon', len(longitudes))
        ds.createDimension('time', len(years))

        # Create latitude, longitude, and time variables
        latitudes_var = ds.createVariable('lat', 'f4', ('lat',))
        longitudes_var = ds.createVariable('lon', 'f4', ('lon',))
        time_var = ds.createVariable('time', 'i4', ('time',))

        # Assign attributes to lat, lon, and time variables
        latitudes_var.units = 'degrees north'
        latitudes_var.long_name = 'Latitude'
        longitudes_var.units = 'degrees east'
        longitudes_var.long_name = 'Longitude'
        time_var.units = 'years'
        time_var.long_name = 'Time'

        # Write data to lat, lon, and time variables
        latitudes_var[:] = latitudes
        longitudes_var[:] = longitudes
        time_var[:] = years

        # Create the main data variable
        data_var = ds.createVariable('ONG_dD', 'f4', ('time', 'lat', 'lon'))
        data_var.units = units
        data_var.long_name = description

        # Write the combined data array to the NetCDF variable
        data_var[:, :, :] = data_array

# Example usage
data_array_dD = load_and_process_data(years)
output_file = "Output/netcdf/ONG_dD_1970-2021.nc"
create_netcdf_with_time_dimension(output_file, data_array_dD, "The country-specific dD-CH4 of ONG CH4 emissions over time", "permil VSMOW")


#%% Now create large netcdf with temporal ONG maps for d13C

# Define latitude and longitude ranges
latitudes = np.linspace(-89.5, 89.5, 180)
longitudes = np.linspace(-179.5, 179.5, 360)
years = np.arange(1970, 2022)  # Years from 1970 to 2021

# Load each CSV file as a NumPy array and flip the matrix vertically
def load_and_process_data(years, prefix='Output/ONG_d13C_1x1_', suffix='.csv'):
    data_array = []
    for year in years:
        filename = f'{prefix}{year}{suffix}'
        df = pd.read_csv(filename, header=None)
        array_data = df.to_numpy()
        flipped_array = np.flipud(array_data)
        data_array.append(flipped_array)
    return np.array(data_array)

# Create a NetCDF file with a time dimension
def create_netcdf_with_time_dimension(filename, data_array, description, units):
    with Dataset(filename, 'w', format='NETCDF4') as ds:
        # Create dimensions
        ds.createDimension('lat', len(latitudes))
        ds.createDimension('lon', len(longitudes))
        ds.createDimension('time', len(years))

        # Create latitude, longitude, and time variables
        latitudes_var = ds.createVariable('lat', 'f4', ('lat',))
        longitudes_var = ds.createVariable('lon', 'f4', ('lon',))
        time_var = ds.createVariable('time', 'i4', ('time',))

        # Assign attributes to lat, lon, and time variables
        latitudes_var.units = 'degrees north'
        latitudes_var.long_name = 'Latitude'
        longitudes_var.units = 'degrees east'
        longitudes_var.long_name = 'Longitude'
        time_var.units = 'years'
        time_var.long_name = 'Time'

        # Write data to lat, lon, and time variables
        latitudes_var[:] = latitudes
        longitudes_var[:] = longitudes
        time_var[:] = years

        # Create the main data variable
        data_var = ds.createVariable('ONG_d13C', 'f4', ('time', 'lat', 'lon'))
        data_var.units = units
        data_var.long_name = description

        # Write the combined data array to the NetCDF variable
        data_var[:, :, :] = data_array

# Example usage
data_array_d13C = load_and_process_data(years)
output_file = "Output/netcdf/ONG_d13C_1970-2021.nc"
create_netcdf_with_time_dimension(output_file, data_array_d13C, "The country-specific d13C-CH4 of ONG CH4 emissions over time", "permil VPDB")


#%% Create netcdf for dD Coal and ONG

# Define latitude and longitude ranges
latitudes = np.linspace(-89.5, 89.5, 180)
longitudes = np.linspace(-179.5, 179.5, 360)

# Function to create a NetCDF file with multiple variables
def create_netcdf_with_multiple_vars(filename, data_dict, descriptions, units):
    with Dataset(filename, 'w', format='NETCDF4') as ds:
        # Create dimensions
        ds.createDimension('lat', len(latitudes))
        ds.createDimension('lon', len(longitudes))

        # Create latitude and longitude variables
        latitudes_var = ds.createVariable('lat', 'f4', ('lat',))
        longitudes_var = ds.createVariable('lon', 'f4', ('lon',))

        # Assign attributes to lat and lon variables
        latitudes_var.units = 'degrees north'
        latitudes_var.long_name = 'Latitude'
        longitudes_var.units = 'degrees east'
        longitudes_var.long_name = 'Longitude'

        # Write data to lat and lon variables
        latitudes_var[:] = latitudes
        longitudes_var[:] = longitudes

        # Create and write each data variable
        for variable_name, data in data_dict.items():
            var = ds.createVariable(variable_name, 'f4', ('lat', 'lon'))
            var.units = units[variable_name]
            var.long_name = descriptions[variable_name]
            var[:, :] = data


# Coal dD first
data_dict_dD_Coal = {
    "Coal_dD-CH4": Coal_dD,
    "stdev_Coal_dD-CH4": Coal_dD_stdev,
    "N_Coal_dD-CH4": Coal_dD_N,
    "stderr_Coal_dD-CH4": ONG_dD_stderr} 
descriptions_dD_Coal = {
    "Coal_dD-CH4": "The country-specific dD-CH4 of coal CH4 emissions",
    "stdev_Coal_dD-CH4": "The country-specific standard deviation of the dD-CH4 of Coal CH4 emissions",
    "N_Coal_dD-CH4": "The number of measurements of coal emission dD-CH4 for each country",
    "stderr_Coal_dD-CH4": "The country-specific standard error of the dD-CH4 of Coal CH4 emissions"} 
units_dD_Coal = {
    "Coal_dD-CH4": "permil VSMOW",
    "stdev_Coal_dD-CH4": "permil VSMOW",
    "N_Coal_dD-CH4": "N",
    "stderr_Coal_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables
output_file_dD_Coal = "Output/netcdf/Coal_dD_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_Coal, data_dict_dD_Coal, descriptions_dD_Coal, units_dD_Coal)

# ONG dD
data_dict_dD_ONG = {
    
    "stdev_ONG_dD-CH4": ONG_dD_stdev,
    "N_ONG_dD-CH4": ONG_dD_N,
    "stderr_ONG_dD-CH4": ONG_dD_stderr} 
descriptions_dD_ONG = {
    "stdev_ONG_dD-CH4": "The country-specific standard deviation of the dD-CH4 of ONG CH4 emissions",
    "N_ONG_dD-CH4": "The number of measurements of ONG emission dD-CH4 for each country",
    "stderr_ONG_dD-CH4": "The country-specific standard error of the dD-CH4 of ONG CH4 emissions"} 
units_dD_ONG = {
    "stdev_ONG_dD-CH4": "permil VSMOW",
    "N_ONG_dD-CH4": "N",
    "stderr_ONG_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables
output_file_dD_ONG = "Output/netcdf/ONG_dD_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_ONG, data_dict_dD_ONG, descriptions_dD_ONG, units_dD_ONG)

# Coal d13C
data_dict_d13C_Coal = {
    "Coal_d13C-CH4": Coal_d13C,
    "stdev_Coal_d13C-CH4": Coal_d13C_stdev,
    "N_Coal_d13C-CH4": Coal_d13C_N,
    "stderr_Coal_d13C-CH4": ONG_d13C_stderr}
descriptions_d13C_Coal = {
    "Coal_d13C-CH4": "The country-specific d13C-CH4 of coal CH4 emissions",
    "stdev_Coal_d13C-CH4": "The country-specific standard deviation of the d13C-CH4 of Coal CH4 emissions",
    "N_Coal_d13C-CH4": "The number of measurements of coal emission d13C-CH4 for each country",
    "stderr_Coal_d13C-CH4": "The country-specific standard error of the d13C-CH4 of Coal CH4 emissions"}
units_d13C_Coal = {
    "Coal_d13C-CH4": "permil VSMOW",
    "stdev_Coal_d13C-CH4": "permil VPDB",
    "N_Coal_d13C-CH4": "N",
    "stderr_Coal_d13C-CH4": "permil VPDB"}
# Create NetCDF file with multiple variables for Coal d13C
output_file_d13C_Coal = "Output/netcdf/Coal_d13C_1x1.nc"
create_netcdf_with_multiple_vars(output_file_d13C_Coal, data_dict_d13C_Coal, descriptions_d13C_Coal, units_d13C_Coal)

# ONG d13C
data_dict_d13C_ONG = {
    "stdev_ONG_d13C-CH4": ONG_d13C_stdev,
    "N_ONG_d13C-CH4": ONG_d13C_N,
    "stderr_ONG_d13C-CH4": ONG_d13C_stderr}
descriptions_d13C_ONG = {
    "stdev_ONG_d13C-CH4": "The country-specific standard deviation of the d13C-CH4 of ONG CH4 emissions",
    "N_ONG_d13C-CH4": "The number of measurements of ONG emission d13C-CH4 for each country",
    "stderr_ONG_d13C-CH4": "The country-specific standard error of the d13C-CH4 of ONG CH4 emissions"}
units_d13C_ONG = {
    "stdev_ONG_d13C-CH4": "permil VPDB",
    "N_ONG_d13C-CH4": "N",
    "stderr_ONG_d13C-CH4": "permil VPDB"}
# Create NetCDF file with multiple variables for ONG d13C
output_file_d13C_ONG = "Output/netcdf/ONG_d13C_1x1.nc"
create_netcdf_with_multiple_vars(output_file_d13C_ONG, data_dict_d13C_ONG, descriptions_d13C_ONG, units_d13C_ONG)

# BB dD
data_dict_dD_BB = {
    "BB_dD-CH4": BB_dD,
    "stdev_BB_dD-CH4": BB_dD_Unc}
descriptions_dD_BB = {
    "BB_dD-CH4": "The dD-CH4 of BB CH4 emissions",
    "stdev_BB_dD-CH4": "The standard deviation of the dD-CH4 of BB CH4 emissions"}
units_dD_BB = {
    "BB_dD-CH4": "permil VSMOW",
    "stdev_BB_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables for ONG d13C
output_file_dD_BB = "Output/netcdf/dD_BB_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_BB, data_dict_dD_BB, descriptions_dD_BB, units_dD_BB)

# Wetland dD
data_dict_dD_Wetland = {
    "Wetland_dD-CH4": Wetland_dD,
    "stdev_Wetland_dD-CH4": Wetland_dD_Unc}
descriptions_dD_Wetland = {
    "Wetland_dD-CH4": "The dD-CH4 of Wetland CH4 emissions",
    "stdev_Wetland_dD-CH4": "The standard deviation of the dD-CH4 of Wetland CH4 emissions"}
units_dD_Wetland = {
    "Wetland_dD-CH4": "permil VSMOW",
    "stdev_Wetland_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables for Wetland dD
output_file_dD_Wetland = "Output/netcdf/dD_Wetland_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_Wetland, data_dict_dD_Wetland, descriptions_dD_Wetland, units_dD_Wetland)

# Lake dD
data_dict_dD_Lake = {
    "Lake_dD-CH4": Lake_dD,
    "stdev_Lake_dD-CH4": Lake_dD_Unc}
descriptions_dD_Lake = {
    "Lake_dD-CH4": "The dD-CH4 of Lake CH4 emissions",
    "stdev_Lake_dD-CH4": "The standard deviation of the dD-CH4 of Lake CH4 emissions"}
units_dD_Lake = {
    "Lake_dD-CH4": "permil VSMOW",
    "stdev_Lake_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables for Lake dD
output_file_dD_Lake = "Output/netcdf/dD_Lake_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_Lake, data_dict_dD_Lake, descriptions_dD_Lake, units_dD_Lake)

# Landfill dD
data_dict_dD_Landfill = {
    "Landfill_dD-CH4": Landfill_dD,
    "stdev_Landfill_dD-CH4": Landfill_dD_Unc}
descriptions_dD_Landfill = {
    "Landfill_dD-CH4": "The dD-CH4 of Landfill CH4 emissions",
    "stdev_Landfill_dD-CH4": "The standard deviation of the dD-CH4 of Landfill CH4 emissions"}
units_dD_Landfill = {
    "Landfill_dD-CH4": "permil VSMOW",
    "stdev_Landfill_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables for Landfill dD
output_file_dD_Landfill = "Output/netcdf/dD_Landfill_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_Landfill, data_dict_dD_Landfill, descriptions_dD_Landfill, units_dD_Landfill)

# Rice dD
data_dict_dD_Rice = {
    "Rice_dD-CH4": Rice_dD,
    "stdev_Rice_dD-CH4": Rice_dD_Unc}
descriptions_dD_Rice = {
    "Rice_dD-CH4": "The dD-CH4 of Rice CH4 emissions",
    "stdev_Rice_dD-CH4": "The standard deviation of the dD-CH4 of Rice CH4 emissions"}
units_dD_Rice = {
    "Rice_dD-CH4": "permil VSMOW",
    "stdev_Rice_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables for Rice dD
output_file_dD_Rice = "Output/netcdf/dD_Rice_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_Rice, data_dict_dD_Rice, descriptions_dD_Rice, units_dD_Rice)

# AllMic dD
data_dict_dD_AllMic = {
    "AllMic_dD-CH4": AllMic_dD,
    "stdev_AllMic_dD-CH4": AllMic_dD_Unc}
descriptions_dD_AllMic = {
    "AllMic_dD-CH4": "The dD-CH4 of AllMic CH4 emissions",
    "stdev_AllMic_dD-CH4": "The standard deviation of the dD-CH4 of AllMic CH4 emissions"}
units_dD_AllMic = {
    "AllMic_dD-CH4": "permil VSMOW",
    "stdev_AllMic_dD-CH4": "permil VSMOW"}
# Create NetCDF file with multiple variables for AllMic dD
output_file_dD_AllMic = "Output/netcdf/dD_AllMic_1x1.nc"
create_netcdf_with_multiple_vars(output_file_dD_AllMic, data_dict_dD_AllMic, descriptions_dD_AllMic, units_dD_AllMic)








