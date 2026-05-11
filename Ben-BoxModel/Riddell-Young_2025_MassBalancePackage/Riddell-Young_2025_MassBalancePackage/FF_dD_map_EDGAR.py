#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 12:52:49 2025

@author: ryoung
"""

# code for developing global ff emission maps
locals().clear()
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
import pandas as pd
import geopandas as gpd
import cartopy.io.shapereader as shpreader
import shapely.geometry
from shapely.geometry import Point

plt.clf()
plt.close('all')


#%% First load data

# load new dD data (already summarized)
coal_dD = pd.read_csv('data/coal_dD.csv')
ONG_dD  = pd.read_csv('data/ONG_dD.csv')  # does not include shale gas

# load EDGAR emissions
Country_ONG_emis  = pd.read_csv('data/EDGAR8_ONG.csv')  # ignores shale gas contributions country-wide
Country_Coal_emis = pd.read_csv('data/EDGAR8_Coal.csv') # ignores shale gas contributions country-wide

# load temporal trends in ONG dD for china, us and canada
US_ONG_dD_data   = pd.read_csv('data/US_ONG_trends.csv')
US_ONG_dD        = US_ONG_dD_data.iloc[:53, 1]
CanadaChina_ONG_dD = pd.read_csv('data/China_Canada_ONG_trends.csv')
China_ONG_dD     = CanadaChina_ONG_dD.iloc[:, 2]

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
dD_ONG_glob   = ONG_dD["mean"].mean()
dD_coal_glob  = coal_dD["mean"].mean()

# normalize names
ONG_dD['COUNTRY']  = ONG_dD['COUNTRY'].str.lower().replace({"russia": "russian federation"})
coal_dD['COUNTRY'] = coal_dD['COUNTRY'].str.lower().replace({"russia": "russian federation"})
country_CH4_emis['COUNTRY'] = country_CH4_emis['COUNTRY'].str.lower()

# rename columns to match old variable names
ONG_dD = ONG_dD.rename(columns={
    'mean': 'dD_CH4_MEAN',
    'std': 'Weighted_Std_ONG',
    'n': 'dD_CH4_N'
})
coal_dD = coal_dD.rename(columns={
    'mean': 'dD_CH4_MEAN',
    'std': 'Weighted_Std_coal',
    'n': 'dD_CH4_N'
})

# build countrymean equivalents (include stds too)
ONG_dD_countrymean  = ONG_dD[['COUNTRY', 'dD_CH4_MEAN', 'Weighted_Std_ONG', 'dD_CH4_N']]
coal_dD_countrymean = coal_dD[['COUNTRY', 'dD_CH4_MEAN', 'Weighted_Std_coal', 'dD_CH4_N']]

# merge data into one dataframe
merged_df = pd.merge(country_CH4_emis, ONG_dD_countrymean, on='COUNTRY', how='left')
merged_df = pd.merge(
    merged_df,
    coal_dD_countrymean,
    on='COUNTRY', how='left',
    suffixes=('_x', '_y')
)

# replace nans with global averages
merged_df['dD_CH4_MEAN_x'].fillna(dD_ONG_glob, inplace=True)
merged_df['dD_CH4_MEAN_y'].fillna(dD_coal_glob, inplace=True)

# replace missing emissions with 0
merged_df['ONG_Avg'].fillna(0, inplace=True)
merged_df['Coal_Avg'].fillna(0, inplace=True)

# calculate global mean FF for countries with no data
dD_FF_glob = (
    dD_ONG_glob * merged_df['ONG_Avg'].sum() +
    dD_coal_glob * merged_df['Coal_Avg'].sum()
) / (merged_df['ONG_Avg'].sum() + merged_df['Coal_Avg'].sum())

# replace mean ong values for china, us, canada with 2022 trend data
merged_df.loc[merged_df['COUNTRY'] == 'china', 'dD_CH4_MEAN_x'] = China_ONG_dD.iat[-1]
merged_df.loc[merged_df['COUNTRY'] == 'united states', 'dD_CH4_MEAN_x'] = US_ONG_dD.iat[-1]

# calculate country level fossil mean (weighted)
merged_df['FF_dD_mean'] = (
    merged_df['ONG_Avg'] * merged_df['dD_CH4_MEAN_x'] +
    merged_df['Coal_Avg'] * merged_df['dD_CH4_MEAN_y']
) / (merged_df['ONG_Avg'] + merged_df['Coal_Avg'])

# calculate % of emissions with no isotope data
ONG_nodata = merged_df['ONG_Avg'][merged_df['Weighted_Std_ONG'].isna()].sum()
Coal_nodata = merged_df['Coal_Avg'][merged_df['Weighted_Std_coal'].isna()].sum()
ONG_nodata_percent = ONG_nodata / merged_df['ONG_Avg'].sum()
Coal_nodata_percent = Coal_nodata / merged_df['Coal_Avg'].sum()


#%% Calculate country specific standard deviation and standard error (using preprocessed data)

# global stdevs from preprocessed data
ONG_stdev_glob  = ONG_dD["Weighted_Std_ONG"].mean()
coal_stdev_glob = coal_dD["Weighted_Std_coal"].mean()

# carry over merged_df and ensure stds are present
merged_df_std = merged_df.copy()

# fill missing stdevs with global averages
merged_df_std['Weighted_Std_ONG'].fillna(ONG_stdev_glob, inplace=True)
merged_df_std['Weighted_Std_coal'].fillna(coal_stdev_glob, inplace=True)

# calculate country level fossil stdev using isotope mass balance
merged_df_std['FF_dD_stdev'] = (
    merged_df_std['ONG_Avg'] * merged_df_std['Weighted_Std_ONG'] +
    merged_df_std['Coal_Avg'] * merged_df_std['Weighted_Std_coal']
) / (merged_df_std['ONG_Avg'] + merged_df_std['Coal_Avg'])

FF_StdDev_Mean = merged_df_std['FF_dD_stdev'].mean(skipna=True)

# calculate standard errors
merged_df_std['ONG_StdErr']  = np.where(
    merged_df_std['dD_CH4_N_x'] > 0,
    merged_df_std['Weighted_Std_ONG'] / np.sqrt(merged_df_std['dD_CH4_N_x']),
    np.nan
)
merged_df_std['Coal_StdErr'] = np.where(
    merged_df_std['dD_CH4_N_y'] > 0,
    merged_df_std['Weighted_Std_coal'] / np.sqrt(merged_df_std['dD_CH4_N_y']),
    np.nan
)
merged_df_std['FF_StdErr'] = np.where(
    (merged_df_std['dD_CH4_N_x'].fillna(0) + merged_df_std['dD_CH4_N_y'].fillna(0)) > 0,
    merged_df_std['FF_dD_stdev'] / np.sqrt(
        merged_df_std['dD_CH4_N_x'].fillna(0) + merged_df_std['dD_CH4_N_y'].fillna(0)
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



#%% calculate global mean ong and coal dD for each year using edgar emissions

# First sum coal and fossil emissions for each year and then weight each country
Sum_Coal = Country_Coal_emis.iloc[:, 1:].sum()
Sum_ONG = Country_ONG_emis.iloc[:, 1:].sum()
# Then weight each country's dD based on emissions
Country_Coal_weight = Country_Coal_emis.iloc[:,1:]/Sum_Coal
Country_Coal_emis['COUNTRY'] = Country_Coal_emis['COUNTRY'].str.lower() #Convert to lower case
Country_Coal_weight = pd.concat([Country_Coal_emis['COUNTRY'], Country_Coal_weight], axis=1)
Country_ONG_weight = Country_ONG_emis.iloc[:,1:]/Sum_ONG
Country_ONG_emis['COUNTRY'] = Country_ONG_emis['COUNTRY'].str.lower() #Convert to lower case
Country_ONG_weight = pd.concat([Country_ONG_emis['COUNTRY'], Country_ONG_weight], axis=1)

# Now, add in temporal component for US, China, and Canada
# Extract country level mean for ONG and coal
coal_dD_All = merged_df[['COUNTRY', 'dD_CH4_MEAN_y']]
ONG_dD_All = merged_df[['COUNTRY', 'dD_CH4_MEAN_x']]
# Expand dataframe
repeated_values = np.tile(ONG_dD_All['dD_CH4_MEAN_x'].values.reshape(-1, 1), (1, 53))
repeated_df = pd.DataFrame(repeated_values, columns=[f'{1969+i}' for i in range(1, 54)])
ONG_dD_All_new = pd.concat([ONG_dD_All[['COUNTRY']], repeated_df], axis=1)

# Now add temporal components to 3 countries
for i, row in ONG_dD_All_new.iterrows():
    if row['COUNTRY'] == 'china':
        ONG_dD_All_new.iloc[i, 1:] = China_ONG_dD.values
    elif row['COUNTRY'] == 'united states':
        ONG_dD_All_new.iloc[i, 1:] = US_ONG_dD.values
    else:
        ONG_dD_All_new.iloc[i, 1:] = ONG_dD_All.loc[i, 'dD_CH4_MEAN_x']

# Merge with country mean ONG and Coal dD
Country_coal_weight_merge = pd.merge(coal_dD_All, Country_Coal_weight, on='COUNTRY', how='left')
Country_ONG_weight_merge = pd.merge(ONG_dD_All_new, Country_ONG_weight, on='COUNTRY', how='left')
Country_coal_weight_merge.fillna(0, inplace=True)
Country_ONG_weight_merge.fillna(0, inplace=True)

# Calculate weighted average for each year for each country
Country_coal_weight_merge.iloc[:, 2:] = Country_coal_weight_merge.iloc[:, 2:].mul(Country_coal_weight_merge['dD_CH4_MEAN_y'], axis=0)
coal_dD_AnnAvg = Country_coal_weight_merge.iloc[:, 2:].sum()
for i in range(1, 54):
    Country_ONG_weight_merge.iloc[:, i + 53] = Country_ONG_weight_merge.iloc[:, i] * Country_ONG_weight_merge.iloc[:, i + 53]
ONG_dD_AnnAvg = Country_ONG_weight_merge.iloc[:, 54:].sum()


#%% merge data with world geography

# path to the local shapefile
shapefile_path = 'data/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp'
world = gpd.read_file(shapefile_path)
world['SOVEREIGNT'] = world['SOVEREIGNT'].str.lower()
world['SOVEREIGNT'] = world['SOVEREIGNT'].replace({"united states of america": "united states",
    "w. sahara": "western sahara", "bosnia and herz.": "bosnia and herzegovina",
    "central african rep.": "central african republic", "congo": "congo-brazzaville",
    "côte d'ivoire": "cote d'ivoire", "dominican rep.": "dominican republic",
    "dem. rep. congo": "congo-kinshasa", "eq. guinea": "equatorial guinea",
    "falkland is.": "falkland islands", "gambia": "gambia, the",
    "myanmar": "burma", "n. cyprus": "cyprus",
    "palestine": "palestinian territories", "s. sudan": "south sudan",
    "solomon is.": "solomon islands", "somoliland": "somalia",
    "turkey": "turkiye", "russia": "russian federation"})

# merge ff with country geography
merged_FF = world.merge(merged_df_std, how='left', left_on='SOVEREIGNT', right_on='COUNTRY')
merged_FF['FF_dD_mean'].fillna(dD_ONG_glob, inplace=True)

# merge ong with country geography
merged_ONG = world.merge(merged_df_std, how='left', left_on='SOVEREIGNT', right_on='COUNTRY')
merged_ONG['dD_CH4_MEAN_x'].fillna(dD_ONG_glob, inplace=True)

# merge coal with country geography
merged_Coal = world.merge(merged_df_std, how='left', left_on='SOVEREIGNT', right_on='COUNTRY')
merged_Coal['dD_CH4_MEAN_y'].fillna(dD_coal_glob, inplace=True)

# merge standard error and standard deviation
merged_stdev = world.merge(merged_df_std, how='left', left_on='SOVEREIGNT', right_on='COUNTRY')
merged_stdev['Weighted_Std_coal'].fillna(coal_stdev_glob, inplace=True)
merged_stdev['Weighted_Std_ONG'].fillna(ONG_stdev_glob, inplace=True)
merged_stdev['Coal_StdErr'].fillna(Coal_StdErr_Mean, inplace=True)
merged_stdev['ONG_StdErr'].fillna(ONG_StdErr_Mean, inplace=True)
merged_stdev['FF_dD_stdev'].fillna(FF_StdDev_Mean, inplace=True)
merged_stdev['FF_StdErr'].fillna(FF_StdErr_Mean, inplace=True)
merged_stdev['dD_CH4_N_x'].fillna(0, inplace=True)
merged_stdev['dD_CH4_N_y'].fillna(0, inplace=True)


#%% data processing

# shapefile reader
shp_filename = shpreader.natural_earth(resolution='50m', category='cultural', name='admin_0_countries')
shp_reader = shpreader.Reader(shp_filename)

country_names = merged_FF.iloc[:, 3].astype(str) 
country_names = country_names.str.title()
country_names = country_names.replace({
    "Congo-Kinshasa": "Democratic Republic of the Congo",
    "Russia": "Russian Federation",
    "Falkland Islands": "Falkland Islands / Malvinas",
    "Fr. S. Antarctic Lands": "French Southern and Antarctic Lands",
    "Cote D'Ivoire": "Côte d'Ivoire",
    "Congo-Brazzaville": "Republic of the Congo",
    "Eswatini": "Kingdom of eSwatini",
    "Palestinian Territories": "Palestine",
    "Gambia, The": "The Gambia",
    "Laos": "Lao PDR",
    "Burma": "Myanmar",
    "North Korea": "Dem. Rep. Korea",
    "South Korea": "Republic of Korea",
    "Turkiye": "Turkey",
    "Brunei": "Brunei Darussalam",
    "Czechia": "Czech Republic",
    "Bosnia And Herzegovina": "Bosnia and Herzegovina",
    "Trinidad And Tobago": "Trinidad and Tobago",
    "United Republic Of Tanzania": "Tanzania",
    "The Bahamas": "Bahamas",
    "East Timor": "Timor-Leste",
    "Ivory Coast": "Côte d'Ivoire",
    "Democratic Republic Of The Congo": "Democratic Republic of the Congo",
    "Republic Of The Congo": "Republic of the Congo",
    "Republic Of Serbia": "Serbia"
})

data_geometries = [None] * len(country_names)
for i, data_country in enumerate(country_names):
    shp_countries = shp_reader.records()  # need to re-read records for every iteration
    for shp_country in shp_countries:
        if data_country == shp_country.attributes['NAME_LONG']:
            data_geometries[i] = shp_country.geometry
    if data_geometries[i] is None:
        print('ERROR: cannot find ' + data_country + ' in shp_countries')


#%% create 1x1 degree csv files of data      

# define lat and lon grid
lat_grid = np.arange(89.5, -90.5, -1)
lon_grid = np.arange(-179.5, 180.5, 1)

# define grids for each parameter and fill with global averages
data_grid_FF = np.full((len(lat_grid), len(lon_grid)), dD_ONG_glob)
data_grid_Coal = np.full((len(lat_grid), len(lon_grid)), dD_coal_glob)
data_grid_Coal_StdErr = np.full((len(lat_grid), len(lon_grid)), Coal_StdErr_Mean)
data_grid_ONG_StdErr = np.full((len(lat_grid), len(lon_grid)), ONG_StdErr_Mean)
data_grid_FF_StdErr = np.full((len(lat_grid), len(lon_grid)), FF_StdErr_Mean)
data_grid_Coal_StdDev = np.full((len(lat_grid), len(lon_grid)), coal_stdev_glob)
data_grid_ONG_StdDev = np.full((len(lat_grid), len(lon_grid)), ONG_stdev_glob)
data_grid_FF_StdDev = np.full((len(lat_grid), len(lon_grid)), FF_StdDev_Mean)
data_grid_Coal_N = np.zeros((len(lat_grid), len(lon_grid)))
data_grid_ONG_N = np.zeros((len(lat_grid), len(lon_grid)))

# fill grids using country polygons
for i_lat, lat in enumerate(lat_grid):
    for i_lon, lon in enumerate(lon_grid):
        p = shapely.geometry.Point([lon, lat])
        for i, data_geometry in enumerate(data_geometries):
            if data_geometry and p.within(data_geometry):
                data_grid_FF[i_lat, i_lon]       = merged_FF['FF_dD_mean'].iloc[i]
                data_grid_Coal[i_lat, i_lon]     = merged_Coal['dD_CH4_MEAN_y'].iloc[i]
                data_grid_Coal_StdErr[i_lat, i_lon] = merged_stdev['Coal_StdErr'].iloc[i]
                data_grid_ONG_StdErr[i_lat, i_lon]  = merged_stdev['ONG_StdErr'].iloc[i]
                data_grid_FF_StdErr[i_lat, i_lon]   = merged_stdev['FF_StdErr'].iloc[i]
                data_grid_Coal_StdDev[i_lat, i_lon] = merged_stdev['Weighted_Std_coal'].iloc[i]
                data_grid_ONG_StdDev[i_lat, i_lon]  = merged_stdev['Weighted_Std_ONG'].iloc[i]
                data_grid_FF_StdDev[i_lat, i_lon]   = merged_stdev['FF_dD_stdev'].iloc[i]
                data_grid_Coal_N[i_lat, i_lon]      = merged_stdev['dD_CH4_N_y'].iloc[i]
                data_grid_ONG_N[i_lat, i_lon]       = merged_stdev['dD_CH4_N_x'].iloc[i]

# export csvs
np.savetxt('Output/FF_dD_1x1.txt', data_grid_FF, delimiter=',', fmt='%+.3f')
np.savetxt('Output/Coal_dD_1x1.txt', data_grid_Coal, delimiter=',', fmt='%+.3f')
np.savetxt('Output/CoalStdErr_dD_1x1.txt', data_grid_Coal_StdErr, delimiter=',', fmt='%+.3f')
np.savetxt('Output/ONGStdErr_dD_1x1.txt', data_grid_ONG_StdErr, delimiter=',', fmt='%+.3f')
np.savetxt('Output/FFStdErr_dD_1x1.txt', data_grid_FF_StdErr, delimiter=',', fmt='%+.3f')
np.savetxt('Output/CoalStdDev_dD_1x1.txt', data_grid_Coal_StdDev, delimiter=',', fmt='%+.3f')
np.savetxt('Output/ONGStdDev_dD_1x1.txt', data_grid_ONG_StdDev, delimiter=',', fmt='%+.3f')
np.savetxt('Output/FFStdDev_dD_1x1.txt', data_grid_FF_StdDev, delimiter=',', fmt='%+.3f')
np.savetxt('Output/CoalN_dD_1x1.txt', data_grid_Coal_N, delimiter=',', fmt='%+.3f')
np.savetxt('Output/ONGN_dD_1x1.txt', data_grid_ONG_N, delimiter=',', fmt='%+.3f')


#%% Save annual results for ONG and total fossil

# [step 1] make grid points once
grid_points = gpd.GeoDataFrame({
    'geometry': [Point(lon, lat) for lat in lat_grid for lon in lon_grid]
}, crs=world.crs)

# [step 2] spatial join once, to match each grid cell to a country
# we use merged_FF which has all needed geometry + country mapping
grid_with_country = gpd.sjoin(
    grid_points, 
    merged_FF[['SOVEREIGNT', 'geometry']], 
    how='left', 
    predicate='within'
).reset_index(drop=True)

# [step 3] prepare to look up country rows quickly
country_to_row = {country: i for i, country in enumerate(merged_df_std['COUNTRY'])}

# [step 4] loop through each year, update US/China, and fill grid fast
years = np.arange(1970, 2023)
for idx, year in enumerate(years):

    # copy the template dataframe
    year_df = merged_df_std.copy()

    # update ONG_Avg from Country_ONG_Emis
    if str(year) in Country_ONG_emis.columns:
        ong_map = Country_ONG_emis.set_index("COUNTRY")[str(year)].to_dict()
        year_df["ONG_Avg"] = year_df["COUNTRY"].map(ong_map).fillna(year_df["ONG_Avg"])
    else:
        print(f"[warning] year {year} not found in Country_ONG_Emis")

    # update Coal_Avg from Country_Coal_Emis
    if str(year) in Country_Coal_emis.columns:
        coal_map = Country_Coal_emis.set_index("COUNTRY")[str(year)].to_dict()
        year_df["Coal_Avg"] = year_df["COUNTRY"].map(coal_map).fillna(year_df["Coal_Avg"])
    else:
        print(f"[warning] year {year} not found in Country_Coal_Emis")

    # check for missing matches
    missing_ong = year_df.loc[year_df["ONG_Avg"].isna(), "COUNTRY"]
    missing_coal = year_df.loc[year_df["Coal_Avg"].isna(), "COUNTRY"]
    if not missing_ong.empty:
        print(f"[error] ONG emissions not found for: {missing_ong.tolist()} in {year}")
    if not missing_coal.empty:
        print(f"[error] Coal emissions not found for: {missing_coal.tolist()} in {year}")

    # update US/China ONG dD
    year_df.loc[year_df['COUNTRY'] == 'united states', 'dD_CH4_MEAN_x'] = US_ONG_dD.iat[idx]
    year_df.loc[year_df['COUNTRY'] == 'china', 'dD_CH4_MEAN_x'] = China_ONG_dD.iat[idx]

    # recalc FF dD mean
    year_df['FF_dD_mean'] = (
        year_df['ONG_Avg'] * year_df['dD_CH4_MEAN_x'] +
        year_df['Coal_Avg'] * year_df['dD_CH4_MEAN_y']
    ) / (year_df['ONG_Avg'] + year_df['Coal_Avg'])

    year_df['FF_dD_mean'].fillna(dD_FF_glob, inplace=True)
    year_df['dD_CH4_MEAN_x'].fillna(dD_ONG_glob, inplace=True)

    # [step 5] map countries to grids
    ff_values, ong_values = [], []
    for country in grid_with_country['SOVEREIGNT']:
        if pd.isna(country):
            ff_values.append(dD_FF_glob)
            ong_values.append(dD_ONG_glob)
        else:
            row_idx = country_to_row.get(country, None)
            if row_idx is not None:
                ff_values.append(year_df.at[row_idx, 'FF_dD_mean'])
                ong_values.append(year_df.at[row_idx, 'dD_CH4_MEAN_x'])
            else:
                ff_values.append(dD_FF_glob)
                ong_values.append(dD_ONG_glob)

    # [step 6] reshape
    grid_FF = np.array(ff_values).reshape(len(lat_grid), len(lon_grid))
    grid_ONG = np.array(ong_values).reshape(len(lat_grid), len(lon_grid))

    # [step 7] save
    np.savetxt(f"Output/FF_dD_1x1_{year}.txt", grid_FF, delimiter=',', fmt='%+.3f')
    np.savetxt(f"Output/ONG_dD_1x1_{year}.txt", grid_ONG, delimiter=',', fmt='%+.3f')



#%% plot it

lower_left_lat = -90
lower_left_lon = -180
upper_right_lat = 90
upper_right_lon = 180
vmin = -300
vmax = -150

# plot ff
data_gridFF_flip = np.flipud(data_grid_FF)
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
            urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(data_gridFF_flip, cmap='viridis', interpolation='nearest',
         norm=Normalize(vmin=vmin, vmax=vmax))
cbar = plt.colorbar(shrink=.5)
cbar.ax.tick_params(labelsize=8)
cbar.set_label('FF ${\delta}D$-CH$_4$ (‰)', fontsize=8)
m.drawcoastlines(linewidth=.4, color='k')
m.drawcountries(linewidth=.4, color='k')

# plot coal
data_gridCoal_flip = np.flipud(data_grid_Coal)
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
            urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(data_gridCoal_flip, cmap='viridis', interpolation='nearest',
         norm=Normalize(vmin=vmin, vmax=vmax))
cbar = plt.colorbar(shrink=.5)
cbar.ax.tick_params(labelsize=8)
cbar.set_label('Coal ${\delta}D$-CH$_4$ (‰)', fontsize=8)
m.drawcoastlines(linewidth=.4, color='k')
m.drawcountries(linewidth=.4, color='k')


# Now plot maps using exported CSV files
# Load data
ONG_dD_1x1_1970 = pd.read_csv('output/FF_dD_1x1_1970.txt')
ONG_dD_1x1_2021 = pd.read_csv('output/FF_dD_1x1_2021.txt')

# Define vmin and vmax for 13C
vmin13C = -300
vmax13C = -150

# Plot ONG dD for 1970
data_gridONG_flip = np.flipud(ONG_dD_1x1_1970)
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(data_gridONG_flip, cmap='viridis', interpolation='nearest', norm=Normalize(vmin=vmin13C, vmax=vmax13C))
plt.set_cmap('viridis')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('ONG ${\delta}D$-CH$_4$ (‰)', fontsize=8)
m.drawcoastlines(linewidth=.4, color='k')
m.drawcountries(linewidth=.4, color='k')

# Plot ONG dD for 2021
data_gridONG_flip = np.flipud(ONG_dD_1x1_2021)
m = Basemap(llcrnrlat=lower_left_lat, llcrnrlon=lower_left_lon,
              urcrnrlat=upper_right_lat, urcrnrlon=upper_right_lon)
plt.figure(dpi=1000)
m.imshow(data_gridONG_flip, cmap='viridis', interpolation='nearest', norm=Normalize(vmin=vmin13C, vmax=vmax13C))
plt.set_cmap('viridis')
cbar = plt.colorbar(shrink=.5)  # Add a color bar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('2021 ONG ${\delta}D$-CH$_4$ (‰)', fontsize=8)
m.drawcoastlines(linewidth=.4, color='k')
m.drawcountries(linewidth=.4, color='k')


