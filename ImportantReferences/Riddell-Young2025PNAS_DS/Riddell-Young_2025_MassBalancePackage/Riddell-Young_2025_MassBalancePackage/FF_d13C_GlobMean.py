#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 10:44:08 2025

@author: ryoung
"""

# Code for developing global FF emission maps
locals().clear()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import netCDF4 as nc
import os
from mpl_toolkits.basemap import Basemap
import xarray as xr

plt.clf()
plt.close('all')

#%% First load data

# load new d13C data (already summarized)
coal_d13C = pd.read_csv('data/coal_d13C.csv')
ONG_d13C  = pd.read_csv('data/ONG_d13C.csv')  # does not include shale gas

# load EDGAR emissions
Country_ONG_emis  = pd.read_csv('data/EDGAR8_ONG.csv')  # ignores shale gas contributions country-wide
Country_Coal_emis = pd.read_csv('data/EDGAR8_Coal.csv') # ignores shale gas contributions country-wide

# load temporal trends in ONG d13C for china, us and canada
US_ONG_d13C_data   = pd.read_csv('data/US_ONG_trends.csv')
US_ONG_d13C        = US_ONG_d13C_data.iloc[:53, 2]
CanadaChina_ONG_d13C = pd.read_csv('data/China_Canada_ONG_trends.csv')
China_ONG_d13C     = CanadaChina_ONG_d13C.iloc[:, 1]
Canada_ONG_d13C    = CanadaChina_ONG_d13C.iloc[:, 3]

# calculate historical average emissions
ONG_average = Country_ONG_emis.iloc[:, 1:].mean(axis=1)
Country_ONG_emis_mean  = pd.DataFrame({'COUNTRY': Country_ONG_emis.iloc[:, 0], 'ONG_Avg': ONG_average})
Coal_average = Country_Coal_emis.iloc[:, 1:].mean(axis=1)
Country_Coal_emis_mean = pd.DataFrame({'COUNTRY': Country_Coal_emis.iloc[:, 0], 'Coal_Avg': Coal_average})

# combine country name with mean ch4 emissions for each country
country_CH4_emis = pd.merge(Country_ONG_emis_mean, Country_Coal_emis_mean, on='COUNTRY', how='outer')

# calculate global total ch4 emissions per type
glob_coal_emis = Country_Coal_emis.iloc[:, 53].sum()
glob_ONG_emis  = Country_ONG_emis.iloc[:, 53].sum()

#%% use new preprocessed isotope data

# determine global means
d13C_ONG_glob   = ONG_d13C["mean"].mean()
d13C_coal_glob  = coal_d13C["mean"].mean()

# normalize names
ONG_d13C['COUNTRY']  = ONG_d13C['COUNTRY'].str.lower().replace({"russia": "russian federation"})
coal_d13C['COUNTRY'] = coal_d13C['COUNTRY'].str.lower().replace({"russia": "russian federation"})
country_CH4_emis['COUNTRY'] = country_CH4_emis['COUNTRY'].str.lower()

# rename columns to match old variable names
ONG_d13C = ONG_d13C.rename(columns={
    'mean': 'd13C_CH4_MEAN',
    'std': 'Weighted_Std_ONG',
    'n': 'd13C_CH4_N'
})
coal_d13C = coal_d13C.rename(columns={
    'mean': 'd13C_CH4_MEAN',
    'std': 'Weighted_Std_coal',
    'n': 'd13C_CH4_N'
})

# build countrymean equivalents (include stds too)
ONG_d13C_countrymean  = ONG_d13C[['COUNTRY', 'd13C_CH4_MEAN', 'Weighted_Std_ONG', 'd13C_CH4_N']]
coal_d13C_countrymean = coal_d13C[['COUNTRY', 'd13C_CH4_MEAN', 'Weighted_Std_coal', 'd13C_CH4_N']]

# merge data into one dataframe
merged_df = pd.merge(country_CH4_emis, ONG_d13C_countrymean, on='COUNTRY', how='left')
merged_df = pd.merge(
    merged_df,
    coal_d13C_countrymean,
    on='COUNTRY', how='left',
    suffixes=('_x', '_y')
)

# replace missing emissions with 0
merged_df['ONG_Avg'].fillna(0, inplace=True)
merged_df['Coal_Avg'].fillna(0, inplace=True)

# calculate global mean FF for countries with no data
d13C_FF_glob = (
    d13C_ONG_glob * merged_df['ONG_Avg'].sum() +
    d13C_coal_glob * merged_df['Coal_Avg'].sum()
) / (merged_df['ONG_Avg'].sum() + merged_df['Coal_Avg'].sum())

# replace mean ong values for china, us, canada with 2022 trend data
merged_df.loc[merged_df['COUNTRY'] == 'china', 'd13C_CH4_MEAN_x'] = China_ONG_d13C.iat[-1]
merged_df.loc[merged_df['COUNTRY'] == 'united states', 'd13C_CH4_MEAN_x'] = US_ONG_d13C.iat[-1]
merged_df.loc[merged_df['COUNTRY'] == 'canada', 'd13C_CH4_MEAN_x'] = Canada_ONG_d13C.iat[-1]

# calculate country level fossil mean (weighted)
merged_df['FF_d13C_mean'] = (
    merged_df['ONG_Avg'] * merged_df['d13C_CH4_MEAN_x'] +
    merged_df['Coal_Avg'] * merged_df['d13C_CH4_MEAN_y']
) / (merged_df['ONG_Avg'] + merged_df['Coal_Avg'])

# calculate % of emissions with no isotope data
ONG_nodata = merged_df['ONG_Avg'][merged_df['d13C_CH4_MEAN_x'].isna()].sum()
Coal_nodata = merged_df['Coal_Avg'][merged_df['d13C_CH4_MEAN_y'].isna()].sum()
ONG_nodata_percent = ONG_nodata / merged_df['ONG_Avg'].sum()
Coal_nodata_percent = Coal_nodata / merged_df['Coal_Avg'].sum()

# replace nans with global averages
merged_df['d13C_CH4_MEAN_x'].fillna(d13C_ONG_glob, inplace=True)
merged_df['d13C_CH4_MEAN_y'].fillna(d13C_coal_glob, inplace=True)

#%% Calculate country specific standard deviation and standard error (using preprocessed data)

# global stdevs from preprocessed data
ONG_stdev_glob  = ONG_d13C["Weighted_Std_ONG"].mean()
coal_stdev_glob = coal_d13C["Weighted_Std_coal"].mean()

# carry over merged_df and ensure stds are present
merged_df_std = merged_df.copy()

# fill missing stdevs with global averages
merged_df_std['Weighted_Std_ONG'].fillna(ONG_stdev_glob, inplace=True)
merged_df_std['Weighted_Std_coal'].fillna(coal_stdev_glob, inplace=True)

# calculate country level fossil stdev using isotope mass balance
merged_df_std['FF_d13C_stdev'] = (
    merged_df_std['ONG_Avg'] * merged_df_std['Weighted_Std_ONG'] +
    merged_df_std['Coal_Avg'] * merged_df_std['Weighted_Std_coal']
) / (merged_df_std['ONG_Avg'] + merged_df_std['Coal_Avg'])

FF_StdDev_Mean = merged_df_std['FF_d13C_stdev'].mean(skipna=True)

# calculate standard errors
merged_df_std['ONG_StdErr']  = np.where(
    merged_df_std['d13C_CH4_N_x'] > 0,
    merged_df_std['Weighted_Std_ONG'] / np.sqrt(merged_df_std['d13C_CH4_N_x']),
    np.nan
)
merged_df_std['Coal_StdErr'] = np.where(
    merged_df_std['d13C_CH4_N_y'] > 0,
    merged_df_std['Weighted_Std_coal'] / np.sqrt(merged_df_std['d13C_CH4_N_y']),
    np.nan
)
merged_df_std['FF_StdErr'] = np.where(
    (merged_df_std['d13C_CH4_N_x'].fillna(0) + merged_df_std['d13C_CH4_N_y'].fillna(0)) > 0,
    merged_df_std['FF_d13C_stdev'] / np.sqrt(
        merged_df_std['d13C_CH4_N_x'].fillna(0) + merged_df_std['d13C_CH4_N_y'].fillna(0)
    ),
    np.nan
)

# compute global averages from real data only
ONG_StdErr_Mean  = merged_df_std['ONG_StdErr'].mean(skipna=True)
Coal_StdErr_Mean = merged_df_std['Coal_StdErr'].mean(skipna=True)
FF_StdErr_Mean   = merged_df_std['FF_StdErr'].mean(skipna=True)

# replace NaNs or infs with global means
merged_df_std['ONG_StdErr'].replace([np.inf, -np.inf], np.nan, inplace=True)
merged_df_std['Coal_StdErr'].replace([np.inf, -np.inf], np.nan, inplace=True)
merged_df_std['FF_StdErr'].replace([np.inf, -np.inf], np.nan, inplace=True)

merged_df_std['ONG_StdErr'].fillna(ONG_StdErr_Mean, inplace=True)
merged_df_std['Coal_StdErr'].fillna(Coal_StdErr_Mean, inplace=True)
merged_df_std['FF_StdErr'].fillna(FF_StdErr_Mean, inplace=True)


#%% Calculate global mean ONG and Coal d13C for each year using EDGAR emissions

# First sum coal and fossil emissions for each year and then weight each country
Sum_Coal = Country_Coal_emis.iloc[:, 1:].sum()
Sum_ONG = Country_ONG_emis.iloc[:, 1:].sum()
# Then weight each country's d13C based on emissions
Country_Coal_weight = Country_Coal_emis.iloc[:,1:]/Sum_Coal
Country_Coal_emis['COUNTRY'] = Country_Coal_emis['COUNTRY'].str.lower() #Convert to lower case
Country_Coal_weight = pd.concat([Country_Coal_emis['COUNTRY'], Country_Coal_weight], axis=1)
Country_ONG_weight = Country_ONG_emis.iloc[:,1:]/Sum_ONG
Country_ONG_emis['COUNTRY'] = Country_ONG_emis['COUNTRY'].str.lower() #Convert to lower case
Country_ONG_weight = pd.concat([Country_ONG_emis['COUNTRY'], Country_ONG_weight], axis=1)

# Now, add in temporal component for US, China, and Canada
# Extract country level mean for ONG and coal
coal_d13C_All = merged_df[['COUNTRY', 'd13C_CH4_MEAN_y']]
ONG_d13C_All = merged_df[['COUNTRY', 'd13C_CH4_MEAN_x']]
# Expand dataframe
repeated_values = np.tile(ONG_d13C_All['d13C_CH4_MEAN_x'].values.reshape(-1, 1), (1, 53))
repeated_df = pd.DataFrame(repeated_values, columns=[f'{1969+i}' for i in range(1, 54)])
ONG_d13C_All_new = pd.concat([ONG_d13C_All[['COUNTRY']], repeated_df], axis=1)

# Now add temporal components to 3 countries
for i, row in ONG_d13C_All_new.iterrows():
    if row['COUNTRY'] == 'china':
        ONG_d13C_All_new.iloc[i, 1:] = China_ONG_d13C.values
    elif row['COUNTRY'] == 'united states':
        ONG_d13C_All_new.iloc[i, 1:] = US_ONG_d13C.values
    elif row['COUNTRY'] == 'canada':
        ONG_d13C_All_new.iloc[i, 1:] = Canada_ONG_d13C.values
    else:
        ONG_d13C_All_new.iloc[i, 1:] = ONG_d13C_All.loc[i, 'd13C_CH4_MEAN_x']

# Merge with country mean ONG and Coal d13C
Country_coal_weight_merge = pd.merge(coal_d13C_All, Country_Coal_weight, on='COUNTRY', how='left')
Country_ONG_weight_merge = pd.merge(ONG_d13C_All_new, Country_ONG_weight, on='COUNTRY', how='left')
Country_coal_weight_merge.fillna(0, inplace=True)
Country_ONG_weight_merge.fillna(0, inplace=True)

# Calculate weighted average for each year for each country
Country_coal_weight_merge.iloc[:, 2:] = Country_coal_weight_merge.iloc[:, 2:].mul(Country_coal_weight_merge['d13C_CH4_MEAN_y'], axis=0)
coal_d13C_AnnAvg = Country_coal_weight_merge.iloc[:, 2:].sum()
for i in range(1, 54):
    Country_ONG_weight_merge.iloc[:, i + 53] = Country_ONG_weight_merge.iloc[:, i] * Country_ONG_weight_merge.iloc[:, i + 53]
ONG_d13C_AnnAvg = Country_ONG_weight_merge.iloc[:, 54:].sum()


#%% Calculate global mean Fossil for each year

# First, calculate total emissions by country each year
Total_emissions = pd.merge(Country_ONG_emis, Country_Coal_emis, on='COUNTRY', how='left')
Total_emissions = pd.merge(coal_d13C_All, Total_emissions, on='COUNTRY', how='left')
Total_emissions = pd.merge(ONG_d13C_All_new, Total_emissions, on='COUNTRY', how='left')
Total_emissions['d13C_CH4_MEAN_y'].fillna(d13C_coal_glob, inplace=True)
Total_emissions.fillna(0, inplace=True)
Emissions_df = pd.DataFrame()
Country_FF_d13C = pd.DataFrame()
for i in range(2, 55):
    column_name1 = f'{1968+i}_x'
    column_name2 = f'{1968+i}_y'  # Column55 corresponds to Column2, Column56 to Column3, etc.
    column_name3 = f'{1968+i}'
    Emissions_df[f'{1968+i}'] = Total_emissions[column_name1] + Total_emissions[column_name2]
    Country_FF_d13C[f'{1968+i}'] = (Total_emissions[column_name1]*Total_emissions[column_name3]
                                    + Total_emissions[column_name2]*Total_emissions['d13C_CH4_MEAN_y']) / (Total_emissions[column_name1] + Total_emissions[column_name2])
    
# Add the 'Country' column to the result DataFrame
Emissions_df.insert(0, 'COUNTRY', Total_emissions['COUNTRY'])
Country_FF_d13C.insert(0, 'COUNTRY', Total_emissions['COUNTRY'])

# Now calculate global average trend over time for d13C
# First sum coal and fossil emissions for each year and then weight each country
Sum_FF = Emissions_df.iloc[:, 1:].sum()
Country_FF_weight = Emissions_df.iloc[:, 1:] / Sum_FF
Emissions_df['COUNTRY'] = Emissions_df['COUNTRY'].str.lower()  # Convert to lower case
Country_FF_weight = pd.concat([Emissions_df['COUNTRY'], Country_FF_weight], axis=1)
# Calculate weighted average for each year for each country
d13C_FF_AnnAvgA = Country_FF_weight.copy()  # Make a copy of df1 to store the result
d13C_FF_AnnAvgA.iloc[:, 1:] = Country_FF_weight.iloc[:, 1:].values * Country_FF_d13C.iloc[:, 1:].values
# Sum each column (excluding the first column)
sums = d13C_FF_AnnAvgA.iloc[:, 1:].sum()

# Create a new DataFrame with the results
d13C_FF_AnnAvg = pd.DataFrame({
    'Column': sums.index,
    'Sum': sums.values
}).reset_index(drop=True)
# For any time with NaN for a country, insert global average
for i in range(0, 53):
    Country_FF_d13C[f'{1970+i}'].fillna(d13C_FF_AnnAvg.iloc[i, 1], inplace=True)


#%% Monte Carlo analysis incoroporating uncertainty in source signatures

# First coal
coal_d13C_AnnAvg_MC = pd.DataFrame()
R_coal_d13C_AnnAvg_MC = pd.DataFrame()
for i in range(0, 1000):
    # First, create and reset dataframe of coal with production weights, mean coal sig and mean coal stdev
    Country_coal_std_merge = pd.merge(coal_d13C_All, merged_df_std[['COUNTRY', 'Coal_StdErr']], on='COUNTRY')
    Country_coal_std_merge_weight = pd.merge(Country_coal_std_merge, Country_Coal_weight, on='COUNTRY', how='left')
    Country_coal_std_merge_weight.fillna(0, inplace=True)
    Country_coal_std_merge_weightB = Country_coal_std_merge_weight
    # Now perform calculations to get global mean coal of each iteration
    random_gaussians = np.random.normal(0, 1, len(Country_coal_std_merge_weight))
    Country_coal_std_merge_weight['MC'] = Country_coal_std_merge_weight['d13C_CH4_MEAN_y'] + (random_gaussians * Country_coal_std_merge_weight['Coal_StdErr'])
    Country_coal_std_merge_weightB.iloc[:, 3:-1] = Country_coal_std_merge_weight.iloc[:, 3:-1].mul(Country_coal_std_merge_weight['MC'], axis=0)
    coal_d13C_AnnAvg_MC[f'{i}'] = Country_coal_std_merge_weightB.iloc[:, 3:-1].sum()
    # Calculate relative
    R_coal_d13C_AnnAvg_MC[f'{i}'] = Country_coal_std_merge_weightB.iloc[:, 3:-1].sum() - Country_coal_std_merge_weightB.iloc[:, 31].sum() #Relative to 1998

# Calculate statistics
coal_d13C_mean = coal_d13C_AnnAvg_MC.mean(axis=1).to_numpy()
coal_d13C_std = coal_d13C_AnnAvg_MC.std(axis=1).to_numpy()
R_coal_d13C_mean = R_coal_d13C_AnnAvg_MC.mean(axis=1).to_numpy()
R_coal_d13C_std = R_coal_d13C_AnnAvg_MC.std(axis=1).to_numpy()


#%% Next ONG
# A potentially important point here is that there is no temporal variability in the d13C stdev of temporally varying countries 

ONG_d13C_AnnAvg_MC = pd.DataFrame()
R_ONG_d13C_AnnAvg_MC = pd.DataFrame()
for i in range(0, 1000):
    # First, create and reset dataframe of coal with production weights, mean coal sig and mean coal stdev
    Country_ONG_std_merge = pd.merge(merged_df_std[['COUNTRY', 'ONG_StdErr']], ONG_d13C_All_new, on='COUNTRY')
    Country_ONG_std_merge_weight = pd.merge(Country_ONG_std_merge, Country_ONG_weight, on='COUNTRY', how='left')
    Country_ONG_std_merge_weight.fillna(0, inplace=True)
    Country_ONG_std_merge_weightB = Country_ONG_std_merge_weight
    # Now perform calculations to get global mean coal of each iteration
    random_gaussians = np.random.normal(0, 1, len(Country_ONG_std_merge_weight))
    for j in range(1, 54):
        Country_ONG_std_merge_weight[f'{j} MC'] = Country_ONG_std_merge_weight.iloc[:, j + 1] + (random_gaussians * Country_ONG_std_merge_weight['ONG_StdErr'])
        Country_ONG_std_merge_weightB.iloc[:, j + 54] = Country_ONG_std_merge_weight.iloc[:, j + 54] * Country_ONG_std_merge_weight.iloc[:, j + 107]
    ONG_d13C_AnnAvg_MC[f'{i}'] = Country_ONG_std_merge_weightB.iloc[:, 55:108].sum()
    # Calculate relative
    R_ONG_d13C_AnnAvg_MC[f'{i}'] = Country_ONG_std_merge_weightB.iloc[:, 55:108].sum() - Country_ONG_std_merge_weightB.iloc[:, 83].sum() #Relative to 1998

# Calculate statistics
ONG_d13C_mean = ONG_d13C_AnnAvg_MC.mean(axis=1).to_numpy()
ONG_d13C_std = ONG_d13C_AnnAvg_MC.std(axis=1).to_numpy()
R_ONG_d13C_mean = R_ONG_d13C_AnnAvg_MC.mean(axis=1).to_numpy()
R_ONG_d13C_std = R_ONG_d13C_AnnAvg_MC.std(axis=1).to_numpy()


#%% Now total fossil

FF_d13C_AnnAvg_MC_EDGAR = pd.DataFrame()
R_FF_d13C_AnnAvg_MC_EDGAR = pd.DataFrame()
for k in range(0, 1000):
    random_gaussian_ONG = np.random.normal(0, 1, len(Country_ONG_std_merge_weight))
    random_gaussian_Coal = np.random.normal(0, 1, len(Country_coal_std_merge_weight))
    Emissions_df_std = pd.DataFrame()
    Country_FF_d13C_std = pd.DataFrame()
    for i in range(2, 55):
        column_name1 = f'{1968+i}_x'
        column_name2 = f'{1968+i}_y'  # Column55 corresponds to Column2, Column56 to Column3, etc.
        column_name3 = f'{1968+i}'
        Emissions_df_std[f'{1968+i}'] = Total_emissions[column_name1] + Total_emissions[column_name2]
        Country_FF_d13C_std[f'{1968+i}'] = (Total_emissions[column_name1] * (Total_emissions[column_name3] + (random_gaussian_ONG * Country_ONG_std_merge_weight['ONG_StdErr']))
                                      + Total_emissions[column_name2] * (Total_emissions['d13C_CH4_MEAN_y'] + (random_gaussian_Coal * Country_coal_std_merge_weight['Coal_StdErr']))) / (Total_emissions[column_name1] + Total_emissions[column_name2])
    # Add the 'Country' column to the result DataFrame
    Emissions_df_std.insert(0, 'COUNTRY', Total_emissions['COUNTRY'])
    Country_FF_d13C_std.insert(0, 'COUNTRY', Total_emissions['COUNTRY'])
    # Now calculate global average trend over time for d13C
    # First sum coal and fossil emissions for each year and then weight each country
    Sum_FF = Emissions_df_std.iloc[:, 1:].sum()
    Country_FF_weight = Emissions_df_std.iloc[:,1:]/Sum_FF
    Emissions_df_std['COUNTRY'] = Emissions_df_std['COUNTRY'].str.lower() #Convert to lower case
    Country_FF_weight = pd.concat([Emissions_df_std['COUNTRY'],Country_FF_weight],axis=1) 
    # Calculate weighted average for each year for each country
    d13C_FF_AnnAvgA_std = Country_FF_weight.copy()  # Make a copy of df1 to store the result
    d13C_FF_AnnAvgA_std.iloc[:, 1:] = Country_FF_weight.iloc[:, 1:].values * Country_FF_d13C_std.iloc[:, 1:].values
    # Sum each column (excluding the first column)
    sums_std = d13C_FF_AnnAvgA_std.iloc[:, 1:].sum()
    # Create a new DataFrame with the results
    FF_d13C_AnnAvg_MC_EDGAR[f'{k}'] = pd.DataFrame({'Sum': sums_std.values}).reset_index(drop=True)
    # Calculate relative
    R_FF_d13C_AnnAvg_MC_EDGAR[f'{k}'] = pd.DataFrame({'Sum': sums_std.values}).reset_index(drop=True) - sums_std.iloc[28]
    
# Calculate statistics
FF_d13C_mean = FF_d13C_AnnAvg_MC_EDGAR.mean(axis=1).to_numpy()
FF_d13C_std = FF_d13C_AnnAvg_MC_EDGAR.std(axis=1).to_numpy()
R_FF_d13C_mean = R_FF_d13C_AnnAvg_MC_EDGAR.mean(axis=1).to_numpy()
R_FF_d13C_std = R_FF_d13C_AnnAvg_MC_EDGAR.std(axis=1).to_numpy()


#%% plot it

years = np.arange(1970, 2023)

# Simplified figure plot 
plt.rcParams.update({'font.size': 14})  # Adjust this value as needed
# Create a figure with 3 subplots
plt.figure(figsize=(10, 6), dpi=1000)
# Plot DataFrame 3 in the third subplot
plt.plot(years, FF_d13C_mean, color='green')
plt.plot(years, FF_d13C_mean + 2 * FF_d13C_std, color='green', linestyle=':', linewidth=3)
plt.plot(years, FF_d13C_mean - 2 * FF_d13C_std, color='green', linestyle=':', linewidth=3, label=r'2$\sigma$ total uncertainty')
plt.plot(years, R_FF_d13C_mean + FF_d13C_mean[28], label='Global Mean Fossil', color='black', linewidth=4)
plt.plot(years, R_FF_d13C_mean + 2 * R_FF_d13C_std+ FF_d13C_mean[28], color='black', linestyle=':', linewidth=3)
plt.plot(years, R_FF_d13C_mean - 2 * R_FF_d13C_std+ FF_d13C_mean[28], color='black', linestyle=':', linewidth=3, label=r'2$\sigma$ relative uncertainty')
plt.ylabel(r"${\delta}^{13}\mathrm{C-CH_4}$ (‰)")
plt.legend(loc='lower right')
plt.title(r"Fossil Fuel ${\delta}^{13}\mathrm{C-CH_4}$ (‰)")
plt.xlim(1998, 2022)
plt.ylim(-45,-43)

# Save data
ONG_d13C = np.column_stack((years, ONG_d13C_mean, ONG_d13C_std))
np.savetxt('Output/ONG_d13C_GlobUnc.csv', ONG_d13C, delimiter = ',', fmt = '%+.3f')
Coal_d13C = np.column_stack((years, coal_d13C_mean, coal_d13C_std))
np.savetxt('Output/Coal_d13C_GlobUnc.csv', Coal_d13C, delimiter = ',', fmt = '%+.3f')
FF_d13C = np.column_stack((years, FF_d13C_mean, FF_d13C_std))
np.savetxt('Output/FF_d13C_GlobUnc.csv', FF_d13C, delimiter = ',', fmt = '%+.3f')
R_FF_d13C = np.column_stack((years, R_FF_d13C_mean, R_FF_d13C_std))
np.savetxt('Output/R_FF_d13C_GlobUnc.csv', R_FF_d13C, delimiter = ',', fmt = '%+.3f')
# Save MC analysis for total fossil
FF_d13C_MC_EDGAR = np.column_stack((years, FF_d13C_AnnAvg_MC_EDGAR))
np.savetxt('Output/FF_d13C_GlobMC_EDGAR.csv', FF_d13C_MC_EDGAR, delimiter = ',', fmt = '%+.3f')


#%% Calculate global average with FF map and CT-CH4 emissions

# Load CT-CH4 fluxes
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
FossilA = flux_data.variables['fossil_flux'][:]
# Expand to 1x1 grid
FF_exp = np.repeat(FossilA, repeats=2, axis=1)
FF_exp = np.repeat(FF_exp, repeats=3, axis=2) / 6 
# Convert from kg/s to Tg/month
kpt = 1000000000  # Kg per Tg
spm = 2.628e+6  # Sec per Month
FF = FF_exp / kpt * spm
# Sum annual emissions and create weighted array
FF_ann = FF.reshape(24, 12, 180, 360).sum(axis=1)


# Load prior netcdf file
Priors = xr.open_dataset('data/prior_monthly_emission_kg_lei.nc')
flux_Coal = Priors['flux_coal']
Coal_ann = flux_Coal.sum(dim='month')
Coal_ann = Coal_ann.values.astype(np.float64)
Coal_ann = np.transpose(Coal_ann[:, :, -24:], (2, 1, 0))

flux_ONG = Priors['flux_oil_gas']
ONG_ann = flux_ONG.sum(dim='month')
ONG_ann = ONG_ann.values.astype(np.float64)
ONG_ann = np.transpose(ONG_ann[:, :, -24:], (2, 1, 0))

flux_Geo = Priors['flux_geologic_seep']
Geo_ann = flux_Geo.sum(dim='month')
Geo_ann = Geo_ann.values.astype(np.float64)
Geo_ann = np.transpose(Geo_ann[:, :, -24:], (2, 1, 0))

flux_Industry = Priors['flux_other_industry']
Industry_ann = flux_Industry.sum(dim='month')
Industry_ann = Industry_ann.values.astype(np.float64)
Industry_ann = np.transpose(Industry_ann[:, :, -24:], (2, 1, 0))

# Calculate posterior emissions minus prior non fossil fuel thermogenic emissions
# Sum all prior emissions
prior_fossil_total = np.sum(Coal_ann + ONG_ann + Geo_ann + Industry_ann)
# Calculate total prior mic
spm = 2.628e+6  # Sec per Month
prior_FF = (Coal_ann + ONG_ann + Geo_ann + Industry_ann) / spm
# Create weighting scheme 
Posterior_to_prior = FF_ann/prior_FF
Posterior_to_prior = np.where(Posterior_to_prior == '--', 1, Posterior_to_prior).astype(float)
Posterior_to_prior = np.where(Posterior_to_prior == 0, 1, Posterior_to_prior)

# Weight priors based on posterior ratios
Coal_ann_posterior = Coal_ann*Posterior_to_prior
ONG_ann_posteriorB = ONG_ann*Posterior_to_prior
Geo_ann_posterior = Geo_ann*Posterior_to_prior
Industry_ann_posterior = Industry_ann*Posterior_to_prior
prior_FF_posterior = prior_FF*Posterior_to_prior
# Add unassigned values to ONG
ONG_ann_posterior = ONG_ann_posteriorB + FF_ann*spm - Coal_ann_posterior - Geo_ann_posterior - ONG_ann_posteriorB - Industry_ann_posterior

# Calculate posterior FF minus adjusted Geo and other industry
FF_CoalONG_posterior = FF_ann*spm - Geo_ann_posterior - Industry_ann_posterior
FF_ann_weight = FF_CoalONG_posterior / FF_CoalONG_posterior.sum(axis=(1, 2), keepdims=True)  # Create annual weighted array

# plot test map
lower_left_lat = -90
lower_left_lon = -180
upper_right_lat = 90
upper_right_lon = 180

FFtest = FF_ann.mean(axis=(0))*spm

m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
            urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(FFtest, cmap='viridis', interpolation='nearest')
cbar = plt.colorbar(shrink=.5)
cbar.ax.tick_params(labelsize=8)
cbar.set_label('FF ${\delta}^{13}C$-CH$_4$ (‰)', fontsize=8)
m.drawcoastlines(linewidth=.4, color='k')
m.drawcountries(linewidth=.4, color='k')

# Load FF_d13C_Uncertainty
start_year = 1998
FF_d13C_Unc = pd.read_csv('Output/FFStdErr_d13C_1x1.csv', header=None)
FF_d13C_Unc = FF_d13C_Unc.values
unique_values = np.unique(FF_d13C_Unc)
# Load 1x1 FF_d13C maps
directory = 'Output/'
compiled_FF_d13C = np.zeros((24, 180, 360))
for i in range(24):
    filename = f'FF_d13C_1x1_{start_year + i}.csv'
    file_path = os.path.join(directory, filename)
    FF_d13C = np.loadtxt(file_path, delimiter=',')
    compiled_FF_d13C[i, :, :] = FF_d13C
compiled_FF_d13C_flip = np.flip(compiled_FF_d13C, axis=1)


#%% Now run MC error analysis using 1x1 d13C map and uncertainties 

# Define number of runs
numMC = 1000
# Initialize an array to store the results
FF_d13C_AnnAvg_MC = np.zeros((24, numMC))
rFF_d13C_AnnAvg_MC = np.zeros((24, numMC))
# Monte Carlo simulation
for i in range(numMC):
    # Generate random Gaussian numbers for each unique value
    gaussian_map = {value: np.random.normal(0, 1) for value in unique_values}
    # Create the new array with Gaussian numbers
    gaussian_array = np.vectorize(lambda x: gaussian_map[x])(FF_d13C_Unc)
    MC_Unc = FF_d13C_Unc * gaussian_array
    MC_Unc = np.tile(MC_Unc, (24, 1, 1))
    MC_Unc_flip = np.flip(MC_Unc, axis=1)
    FF_d13C_rand = compiled_FF_d13C_flip + MC_Unc_flip
    # Multiply source sig map by emission weighting
    FF_d13C_weight = FF_d13C_rand * FF_ann_weight
    FF_d13C_glob = FF_d13C_weight.sum(axis=(1, 2))
    # Store the result in the array
    FF_d13C_AnnAvg_MC[:, i] = FF_d13C_glob
    rFF_d13C_AnnAvg_MC[:, i] = FF_d13C_glob - FF_d13C_glob[0]
    
# Calculate statistics
FF_d13C_AnnAvg = FF_d13C_AnnAvg_MC.mean(axis=1)
FF_d13C_AnnAvg_std = FF_d13C_AnnAvg_MC.std(axis=1)
rFF_d13C_AnnAvg = rFF_d13C_AnnAvg_MC.mean(axis=1)
rFF_d13C_AnnAvg_std = rFF_d13C_AnnAvg_MC.std(axis=1)

# Save MC analysis for edgar
FF_d13C_MC_CTCH4 = np.column_stack((years[28:-1], FF_d13C_AnnAvg_MC))
np.savetxt('Output/FF_d13C_GlobMC_CTCH4.csv', FF_d13C_MC_CTCH4, delimiter = ',', fmt = '%+.3f')


#%% Plot carbontracker vs EDGAR

years = np.arange(1970, 2023)

# Plot the absolute results
plt.rcParams.update({'font.size': 20})  # Adjust this value as needed
# Create a figure with 3 subplots
fig, axs = plt.subplots(3, 1, figsize=(10, 15), dpi=500)
# Plot DataFrame 1 in the first subplot
axs[0].plot(years, coal_d13C_mean, label='Coal avg', color='blue')
axs[0].plot(years, coal_d13C_mean[28] + R_coal_d13C_mean + 2 * R_coal_d13C_std, color='blue', linestyle='--')
axs[0].plot(years, coal_d13C_mean[28] + R_coal_d13C_mean - 2 * R_coal_d13C_std, color='blue', linestyle='--')
axs[0].plot(years, coal_d13C_mean + 2 * coal_d13C_std, color='blue', linestyle=':')
axs[0].plot(years, coal_d13C_mean - 2 * coal_d13C_std, color='blue', linestyle=':')
axs[0].set_ylabel(r"${\delta}^{13}\mathrm{C-CH_4}$ (‰)")
axs[0].legend(loc='lower right')
axs[0].set_title('Coal Data')
# Plot DataFrame 2 in the second subplot
axs[1].plot(years, ONG_d13C_mean, label='ONG avg', color='orange')
axs[1].plot(years, ONG_d13C_mean[28] + R_ONG_d13C_mean + 2*R_ONG_d13C_std, color='orange', linestyle='--')
axs[1].plot(years, ONG_d13C_mean[28] + R_ONG_d13C_mean - 2*R_ONG_d13C_std, color='orange', linestyle='--')
axs[1].plot(years, ONG_d13C_mean + 2 * ONG_d13C_std, color='orange', linestyle=':')
axs[1].plot(years, ONG_d13C_mean - 2 * ONG_d13C_std, color='orange', linestyle=':')
axs[1].set_ylabel(r"${\delta}^{13}\mathrm{C-CH_4}$ (‰)")
axs[1].legend(loc='lower right')
axs[1].set_title('ONG Data')
# Plot DataFrame 3 in the third subplot
axs[2].plot(years, FF_d13C_mean, label='EDGAR', color='green')
axs[2].plot(years, FF_d13C_mean + 2 * FF_d13C_std, color='green', linestyle=':')
axs[2].plot(years, FF_d13C_mean - 2 * FF_d13C_std, color='green', linestyle=':')
axs[2].plot(years, FF_d13C_mean[28] + R_FF_d13C_mean + 2 * R_FF_d13C_std, color='green', linestyle='--')
axs[2].plot(years, FF_d13C_mean[28] + R_FF_d13C_mean - 2 * R_FF_d13C_std, color='green', linestyle='--')
axs[2].plot(years[28:-1], FF_d13C_AnnAvg, label='CTCH4', color='black', linewidth=5)
axs[2].plot(years[28:-1], FF_d13C_AnnAvg[0] + rFF_d13C_AnnAvg + 2 * rFF_d13C_AnnAvg_std, color='black', linestyle='--')
axs[2].plot(years[28:-1], FF_d13C_AnnAvg[0] + rFF_d13C_AnnAvg - 2 * rFF_d13C_AnnAvg_std, color='black', linestyle='--')
axs[2].plot(years[28:-1], FF_d13C_AnnAvg + 2 * FF_d13C_AnnAvg_std, color='black', linestyle=':')
axs[2].plot(years[28:-1], FF_d13C_AnnAvg - 2 * FF_d13C_AnnAvg_std, color='black', linestyle=':')
axs[2].set_ylim(-46.5,-42.7)
axs[2].set_xlabel('Year')
axs[2].set_ylabel('${\delta}^{13}C$-CH$_4$ (‰)')
axs[2].legend(loc='lower left')
axs[2].set_title('Total Fossil Data')
# Set common x-axis ticks and labels
for ax in axs:
    ax.set_xticks([1970, 1980, 1990, 2000, 2010, 2020])
plt.tight_layout()  # Adjust subplot parameters to give some padding
plt.show()



