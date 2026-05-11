#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 11:22:23 2024

@author: ryoung
"""

# This code estimates global and semihemispheric mean dD-CH4 and uncertainties using a 
# compilation of data and NOAA data extension and integration methodology

locals().clear()
#import ccg_filter
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings
import os

# Suppress runtime warnings related to NaN calculations
np.seterr(invalid='ignore', divide='ignore')

# Suppress pandas SettingWithCopyWarning
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)

# Description of uncertainties added: 
    # 1. Network: Each MC iteration removes 2 sites at random and continues with calculations
    # 2. ADDED IN SMOOTHED MBL CALCULATION CODE. Atmospheric: Adds error to monthly means
    # 3. EMBEDDED IN ATMOSPHERIC. Analysis: Adds error to individual measurements
    # 4. Measurement bias: Adds uncertainty to interlaboratory comparisons wrt. MPI

#%% 

# List of input filenames. Data in these files are moved to the MPI scale based on scale comparisons reported in Umezeawa et al., 2018. 
filenames = [
    "alt_01D0_dat.txt", "altMPI_01D0_dat.txt", "asc_01D0_dat.txt", "ato_01D0_dat.txt",
    "azr_01D0_dat.txt", "bal_01D0_dat.txt", "bikMPI_01D0_dat.txt", "brw_01D0_dat.txt",
    "bsc_01D0_dat.txt", "cba_01D0_dat.txt", "cgo_01D0_dat.txt",
    "cvoMPI_01D0_dat.txt", "eom_01D0_dat.txt", "gvnMPI_01D0_dat.txt", "gvnIMAU_01D0_dat.txt",
    "jfjMPI_01D0_dat.txt", "kjnMPI_01D0_dat.txt", "kum_01D0_dat.txt", "lef_01D0_dat.txt",
    "mhd_01D0_dat.txt", "mlo_01D0_dat.txt", "mloIMAU_01D0_dat.txt", "namMPI_01D0_dat.txt",
    "nyaNIPR_01D0_dat.txt", "oxkMPI_01D0_dat.txt", "sisMPI_01D0_dat.txt", "smo_01D0_dat.txt",
    "spo_01D0_dat.txt", "syoNIPR_01D0_dat.txt", "vrsMPI_01D0_dat.txt",
    "zep_01D0_dat.txt",  "zepIMAU_01D0_dat.txt", "zotMPI_01D0_dat.txt"]


# Run a loop where two sites are left out of the DEI calculation at each iteration
# Number of iterations (Increase to 1000 for MC simulation)
iterations = 1000
# Define years
start_year = 2005
end_year = 2024.5
year = np.arange(start_year, (end_year - 0.5), 1)
# Initialize results matrix
# First weekly matrices
WeeklyGlob_matrix = np.zeros((1013, iterations))
WeeklyNH_matrix = np.zeros((1013, iterations))
WeeklySH_matrix = np.zeros((1013, iterations))
WeeklyPN_matrix = np.zeros((1013, iterations))
WeeklyTN_matrix = np.zeros((1013, iterations))
WeeklyPS_matrix = np.zeros((1013, iterations))
WeeklyTS_matrix = np.zeros((1013, iterations))
# Now annual averages 
smoothed_matrix = np.zeros((len(year), iterations))
smoothed_matrixR = np.zeros((len(year), iterations))
smoothedNH_matrix = np.zeros((len(year), iterations))
smoothedNH_matrixR = np.zeros((len(year), iterations))
smoothedSH_matrix = np.zeros((len(year), iterations))
smoothedSH_matrixR = np.zeros((len(year), iterations))
smoothedPN_matrix = np.zeros((len(year), iterations))
smoothedPN_matrixR = np.zeros((len(year), iterations))
smoothedTN_matrix = np.zeros((len(year), iterations))
smoothedTN_matrixR = np.zeros((len(year), iterations))
smoothedTS_matrix = np.zeros((len(year), iterations))
smoothedTS_matrixR = np.zeros((len(year), iterations))
smoothedPS_matrix = np.zeros((len(year), iterations))
smoothedPS_matrixR = np.zeros((len(year), iterations))

# Define intercomparison uncertainties to MPI
IMAU_Unc = 2.2 #New = 1.4; Dasgupta = 2.2
NIPR_Unc = 1.6 #Umezawa = 1.4; Dasgupta = 1.6
INSTAAR_Unc = 1.6 #Umezawa = 1.5; Dasgupta = 1.6


# Load results of MBL_calc_Unc code:
for filename in filenames:
    # Extract the site name by removing the "_dat.txt" suffix
    site_name = filename.replace("_dat.txt", "_smoothedMC")
    # Construct the corresponding smoothedMC file path
    file_path = f'output/{site_name}.txt'
        # Load the matrix from the smoothedMC file
    if os.path.exists(file_path):
        # Load the data and assign it to a variable with the site name
        globals()[site_name] = np.loadtxt(file_path)
        
        
#%% Adjust scaling based on paper analysis (see All_MBL_data_hemplot.xlsx readme file for calculations)

# Define scaling changes:
# These numbers are converting from the scale in the raw data scaling (based on Dasgupta et al., 2025 (preprint)) to the mean Umezawa et al., 2018 scaling. 
# For INSTAAR, the new interocomparison between MPI and INSTAAR reported in Riddell-Young et al., 2025 is used. 
INSTAAR = 1.8 #New = 1.2; Dasgupta = 1.8
IMAU = .5  #Umezawa = 2.33; Dasgupta = 0.5
NIPR = 0  #Umezawa = 1.85; Dasgupta = 0
        
# Load results of MBL_calc_Unc code:
for filename in filenames:
    # Extract the site name by removing the "_dat.txt" suffix
    site_name = filename.replace("_dat.txt", "_smoothedMC")
    # Construct the corresponding smoothedMC file path
    file_path = f'output/{site_name}.txt'
    # Load the matrix from the smoothedMC file
    if os.path.exists(file_path):
        # Load the data and assign it to a variable with the site name
        globals()[site_name] = np.loadtxt(file_path)
        data = globals()[site_name]
        # Modify the dataset based on filename conditions
        if isinstance(data, np.ndarray) and data.shape[1] > 1:
            if "IMAU" in filename:
                data[:, 1:] -= IMAU  # Add 1.33 to the second column for IMAU files
            if "NIPR" in filename:
                data[:, 1:] -= NIPR  # Add 0.85 to the second column for NIPR files
            if all(keyword not in filename for keyword in ["IMAU", "NIPR", "MPI"]):
                data[:, 1:] -= INSTAAR  # Subtract 1 from all columns except the first
            # Store the modified data back in the global variable
            globals()[site_name] = data


#%%
# Run for loop MC error propagation
for k in range(iterations):
    # Load site info
    siteinfo = pd.read_csv('data/siteinfo_all_ch4h2.txt', sep=' \| ', engine='python', header=None)
    siteinfo_filtered = siteinfo[siteinfo[9] != 0].copy() # Removes non-MBL sites
    df_sites_i = siteinfo_filtered[[0, 3]].copy()  # column 0: site name, column 3: latitude
    df_sites_i.columns = ['site', 'latitude']  # rename the columns for clarity
    df_sites = df_sites_i.drop(df_sites_i.sample(2).index).reset_index(drop=True)
    
    # Combine site info with results of curve fitting
    # Initialize a list to store the combined data
    combined_data = []
    # Loop through the sites and match them with the smoothed data
    for site in df_sites['site'].values:
        # Retrieve the corresponding matrix for the site
        site_var_name = site + '_01D0_smoothedMC'
        site_matrix_var = globals().get(site_var_name)
    
        if site_matrix_var is not None:
            # Get latitude from df_sites
            latitude = df_sites.loc[df_sites['site'] == site, 'latitude'].values[0]
    
            # Extract x_interp and smoothed_curve from the matrix
            x_interp = site_matrix_var[:, 0]  # First column is x_interp
            smoothed_curve = site_matrix_var[:, (k+1)]  # (k+1)th column is smoothed_curve for iteration k
    
            # Append all data as a single row in the final DataFrame
            combined_data.append({
                'site': site,
                'latitude': latitude,
                'x_interp': x_interp,
                'smoothed_curve': smoothed_curve})
    # Convert the list of dictionaries into a DataFrame
    final_df = pd.DataFrame(combined_data)
    
    # Add scale interncomparison uncertainty (compared to MPI scale)
    # Find rows where 'IMAU' and 'NIPR' appears in column 1 (assuming column 1 is the 'site' column)
    imau_rows = final_df[final_df['site'].str.contains('IMAU')]
    NIPR_rows = final_df[final_df['site'].str.contains('NIPR')]
    INSTAAR_rows = final_df[final_df['site'].str.len() == 3]
    # Add random values with mean=0 and std=x to column 4 ('smoothed_curve') of those rows
    final_df.loc[final_df['site'].str.contains('IMAU'), 'smoothed_curve'] += np.random.normal(0, IMAU_Unc, size=len(imau_rows))
    final_df.loc[final_df['site'].str.contains('NIPR'), 'smoothed_curve'] += np.random.normal(0, NIPR_Unc, size=len(NIPR_rows))
    final_df.loc[final_df['site'].str.len() == 3, 'smoothed_curve'] += np.random.normal(0, INSTAAR_Unc, size=len(INSTAAR_rows))
    
    # Calculate weekly value for each dataset
    def combine_row_data(row, date, y_column):
        combined_y = []  # Store the combined y values for each interval
        for i in range(len(date) - 1):  # iterate over consecutive pairs in date
            y_vals = []  # Store y values for the current interval
            x_values = row['x_interp']
            y_values = row[y_column]  # 'smoothed_curve' depending on which column to process
            # Find values in x_values between two consecutive numbers in date
            mask = (x_values >= date[i]) & (x_values < date[i+1])
            # Append the corresponding y values
            y_vals.extend(y_values[mask])
            # Store combined y data for the current interval
            combined_y.append(y_vals)
        return combined_y
    # Create the date range for weekly intervals, defining years of analysis
    WeeksInYear = 52
    date = np.arange(start_year, end_year, 1/WeeksInYear)
    # Apply the function to each row of the DataFrame for "smoothed_curve"
    final_df['combined_smoothed'] = final_df.apply(lambda row: combine_row_data(row, date, 'smoothed_curve'), axis=1)
    
    # Separate into four latitudinal blocks
    df_PS = final_df[(final_df['latitude'] >= -90) & (final_df['latitude'] < -30)]
    df_TS = final_df[(final_df['latitude'] >= -30) & (final_df['latitude'] < 0)]
    df_TN = final_df[(final_df['latitude'] >= 0) & (final_df['latitude'] < 30)]
    df_PN = final_df[(final_df['latitude'] >= 30) & (final_df['latitude'] <= 90)]
    
    # Now, we're going to group data by latitudinal block
    date = np.arange(start_year, end_year, 1/WeeksInYear)
    # Function to combine data based on arr intervals and store the site if only one data point exists
    def combine_data_in_blocks(df_block, date):
        combined_smoothed = []
        single_site_weeks = []  # List to store the site for weeks with only one data point
        for i in range(len(date) - 1):  # iterate over consecutive pairs in arr
            smoothed_vals = []
            sites = []  # Store the site names for each week
            for _, row in df_block.iterrows():
                x_interp = row['x_interp']
                smoothed_curve = row['smoothed_curve']
                site = row['site']  # Get the site name
                # Find values in x_interp between two consecutive numbers in date
                mask = (x_interp >= date[i]) & (x_interp < date[i+1])
                # Append the corresponding smoothed_curve values and site
                smoothed_vals.extend(smoothed_curve[mask])
                if mask.any():
                    sites.append(site)
            # Store combined data for the current interval
            combined_smoothed.append(smoothed_vals)
            # If only one data point exists, store the site in the list
            if len(smoothed_vals) == 1:
                single_site_weeks.append(sites[0])  # Store the site name for that week
            else:
                single_site_weeks.append(None)  # Append None if more than one or no data point exists
        return combined_smoothed, single_site_weeks
    
    # Combine data for each latitudinal block
    combined_smoothed_PS, single_site_weeks_PS = combine_data_in_blocks(df_PS, date)
    combined_smoothed_TS, single_site_weeks_TS = combine_data_in_blocks(df_TS, date)
    combined_smoothed_TN, single_site_weeks_TN = combine_data_in_blocks(df_TN, date)
    combined_smoothed_PN, single_site_weeks_PN = combine_data_in_blocks(df_PN, date)
    
    # Calculate average for each week only when there are at least two values
    def average_each_row(combined_smoothed):
        return np.array([np.mean(row) if len(row) >= 1 else np.nan for row in combined_smoothed])
    # Compute the row-wise average for each combined_smoothed variable
    avg_smoothed_PS = average_each_row(combined_smoothed_PS)
    avg_smoothed_TS = average_each_row(combined_smoothed_TS)
    avg_smoothed_TN = average_each_row(combined_smoothed_TN)
    avg_smoothed_PN = average_each_row(combined_smoothed_PN)
    
    # Calculate average for each week only when there is at least one value for comparisons sake
    def average_each_rowB(combined_smoothed):
        return np.array([np.mean(row) if len(row) >= 1 else np.nan for row in combined_smoothed])
    # Compute the row-wise average for each combined_smoothed variable
    avg_smoothed_PS_singlesites = average_each_rowB(combined_smoothed_PS)
    avg_smoothed_TS_singlesites = average_each_rowB(combined_smoothed_TS)
    avg_smoothed_TN_singlesites = average_each_rowB(combined_smoothed_TN)
    avg_smoothed_PN_singlesites = average_each_rowB(combined_smoothed_PN)
    
    # If value is surrounded by NaN on both sides, replace with NaN
    avg_smoothed_PS_sub = np.where((np.isnan(np.roll(avg_smoothed_PS, 1))) & (np.isnan(np.roll(avg_smoothed_PS, -1))), np.nan, avg_smoothed_PS)
    avg_smoothed_TS_sub = np.where((np.isnan(np.roll(avg_smoothed_TS, 1))) & (np.isnan(np.roll(avg_smoothed_TS, -1))), np.nan, avg_smoothed_TS)
    avg_smoothed_TN_sub = np.where((np.isnan(np.roll(avg_smoothed_TN, 1))) & (np.isnan(np.roll(avg_smoothed_TN, -1))), np.nan, avg_smoothed_TN)
    avg_smoothed_PN_sub = np.where((np.isnan(np.roll(avg_smoothed_PN, 1))) & (np.isnan(np.roll(avg_smoothed_PN, -1))), np.nan, avg_smoothed_PN)
    
    # Calculate number of sites for each week in each region
    sites_PS = np.array([sum(1 for item in row if isinstance(item, (int, float))) for row in combined_smoothed_PS])
    sites_TS = np.array([sum(1 for item in row if isinstance(item, (int, float))) for row in combined_smoothed_TS])
    sites_TN = np.array([sum(1 for item in row if isinstance(item, (int, float))) for row in combined_smoothed_TN])
    sites_PN = np.array([sum(1 for item in row if isinstance(item, (int, float))) for row in combined_smoothed_PN])
    # Remove anamolous weeks where two weeks were averaged. 
    for i in range(1, len(sites_PS) - 1):
        if sites_PS[i] == sites_PS[i-1] + 1 == sites_PS[i+1] + 1:
            sites_PS[i] -= 1
        if sites_TS[i] == sites_TS[i-1] + 1 == sites_TS[i+1] + 1:
            sites_TS[i] -= 1
        if sites_TN[i] == sites_TN[i-1] + 1 == sites_TN[i+1] + 1:
            sites_TN[i] -= 1
        if sites_PN[i] == sites_PN[i-1] + 1 == sites_PN[i+1] + 1:
            sites_PN[i] -= 1
    sites_total = sites_PS + sites_TS + sites_TN + sites_PN


    # Fill in data gaps and calculate global average
    
    # Calculate difference between subhemispheric sections and PN
    Smooth_Diff_Smooth_PN_PS = avg_smoothed_PN_sub - avg_smoothed_PS_sub
    Smooth_Diff_Smooth_PN_TS = avg_smoothed_PN_sub - avg_smoothed_TS_sub
    Smooth_Diff_Smooth_PN_TN = avg_smoothed_PN_sub - avg_smoothed_TN_sub
    # Calculate difference between subhemispheric sections and Ps
    Smooth_Diff_Smooth_PS_PN = avg_smoothed_PS_sub - avg_smoothed_PN_sub
    Smooth_Diff_Smooth_PS_TS = avg_smoothed_PS_sub - avg_smoothed_TS_sub
    Smooth_Diff_Smooth_PS_TN = avg_smoothed_PS_sub - avg_smoothed_TN_sub
    
    # Function to calculate weekly averages, ignoring NaN and zeros
    def calculate_weekly_avg_diff(diff_array, num_weeks=WeeksInYear):
        # Create an array to store weekly averages
        weekly_avg_diff = np.zeros(num_weeks)
        # Assume each year's data is grouped into `num_weeks` weeks
        for week in range(num_weeks):
            # Get values for the current week
            week_values = diff_array[week::num_weeks]  # Get every 52nd value, starting at index `week`
            # Mask to ignore NaN and zero values
            mask = ~np.isnan(week_values) & (week_values != 0)
            # Calculate the mean for the current week
            weekly_avg_diff[week] = np.mean(week_values[mask]) if np.any(mask) else np.nan
        return weekly_avg_diff
    
    # Calculate weekly differences for each comparison
    weekly_avg_Smooth_Diff_Smooth_PN_PS = calculate_weekly_avg_diff(Smooth_Diff_Smooth_PN_PS)
    weekly_avg_Smooth_Diff_Smooth_PN_TS = calculate_weekly_avg_diff(Smooth_Diff_Smooth_PN_TS)
    weekly_avg_Smooth_Diff_Smooth_PN_TN = calculate_weekly_avg_diff(Smooth_Diff_Smooth_PN_TN)
    weekly_avg_Smooth_Diff_Smooth_PS_PN = calculate_weekly_avg_diff(Smooth_Diff_Smooth_PS_PN)
    weekly_avg_Smooth_Diff_Smooth_PS_TS = calculate_weekly_avg_diff(Smooth_Diff_Smooth_PS_TS)
    weekly_avg_Smooth_Diff_Smooth_PS_TN = calculate_weekly_avg_diff(Smooth_Diff_Smooth_PS_TN)
    
    # Now, fill in the gaps using the correct weekly differences for the corresponding weeks
    week_indices = np.mod(np.arange(len(avg_smoothed_PS_sub)), WeeksInYear)
    # Calculate PS where there is no data by using mean weekly difference between PN and PS
    avg_smoothed_PS_subs = np.where(np.isnan(avg_smoothed_PS_sub), 
        avg_smoothed_PN_sub - weekly_avg_Smooth_Diff_Smooth_PN_PS[week_indices], avg_smoothed_PS_sub)
    
    def safe_mean(value1, value2):
        mask1_nan = np.isnan(value1)
        mask2_nan = np.isnan(value2)
        # Case when both values are NaN
        if np.all(mask1_nan) and np.all(mask2_nan):
            return np.nan
        # Case when only value1 is NaN
        elif np.all(mask1_nan):
            return value2
        # Case when only value2 is NaN
        elif np.all(mask2_nan):
            return value1
        # Case when neither is NaN, average the two
        else:
            return (value1 + value2) / 2
        
    # Calculate TS where there is no data by using mean weekly difference between PN and TS AND PS and TS
    avg_smoothed_TS_subs = np.where(np.isnan(avg_smoothed_TS_sub), safe_mean(
            avg_smoothed_PN_sub - weekly_avg_Smooth_Diff_Smooth_PN_TS[week_indices],
            avg_smoothed_PS_subs - weekly_avg_Smooth_Diff_Smooth_PS_TS[week_indices]), avg_smoothed_TS_sub)
    # Calculate TN where there is no data by using mean weekly difference between PN and TN AND PS and TN
    avg_smoothed_TN_subs = np.where(np.isnan(avg_smoothed_TN_sub), safe_mean(
            avg_smoothed_PN_sub - weekly_avg_Smooth_Diff_Smooth_PN_TN[week_indices],
            avg_smoothed_PS_subs - weekly_avg_Smooth_Diff_Smooth_PS_TN[week_indices]), avg_smoothed_TN_sub)
    # Calculate PN where there is no data by using mean weekly difference between TN and PN AND TS and PN
    avg_smoothed_PN_subs = np.where(np.isnan(avg_smoothed_PN_sub), safe_mean(
            avg_smoothed_TN_subs + weekly_avg_Smooth_Diff_Smooth_PN_TN[week_indices],
            avg_smoothed_TS_subs + weekly_avg_Smooth_Diff_Smooth_PN_TS[week_indices]), avg_smoothed_PN_sub)
    # Calculate PS where there is no data by using mean weekly difference between TN and PS AND TS and PS
    avg_smoothed_PS_subs2 = np.where(np.isnan(avg_smoothed_PS_subs), safe_mean(
            avg_smoothed_TN_subs + weekly_avg_Smooth_Diff_Smooth_PS_TN[week_indices],
            avg_smoothed_TS_subs + weekly_avg_Smooth_Diff_Smooth_PS_TS[week_indices]), avg_smoothed_PS_subs)
    
    # Calculate global and hemispheric averages
    avg_glob_smoothed = (avg_smoothed_PN_subs + avg_smoothed_PS_subs2 + avg_smoothed_TS_subs + avg_smoothed_TN_subs)/4
    avg_NH_smoothed = (avg_smoothed_PN_subs + avg_smoothed_TN_subs)/2
    avg_SH_smoothed = (avg_smoothed_PS_subs2 + avg_smoothed_TS_subs)/2
    
    # Get rid of anamolous jumps
    for i in range(1, len(avg_glob_smoothed) - 1):
        if abs(avg_glob_smoothed[i - 1] - avg_glob_smoothed[i]) > 0.3 and abs(avg_glob_smoothed[i + 1] - avg_glob_smoothed[i]) > 0.3:
            avg_glob_smoothed[i] = (avg_glob_smoothed[i - 1] + avg_glob_smoothed[i + 1]) / 2
    for i in range(1, len(avg_NH_smoothed) - 1):
        if abs(avg_NH_smoothed[i - 1] - avg_NH_smoothed[i]) > 0.3 and abs(avg_NH_smoothed[i + 1] - avg_NH_smoothed[i]) > 0.3:
            avg_NH_smoothed[i] = (avg_NH_smoothed[i - 1] + avg_NH_smoothed[i + 1]) / 2
    for i in range(1, len(avg_SH_smoothed) - 1):
        if abs(avg_SH_smoothed[i - 1] - avg_SH_smoothed[i]) > 0.3 and abs(avg_SH_smoothed[i + 1] - avg_SH_smoothed[i]) > 0.3:
            avg_SH_smoothed[i] = (avg_SH_smoothed[i - 1] + avg_SH_smoothed[i + 1]) / 2
    
    # Define annual chunks
    full_chunks = len(avg_glob_smoothed) // WeeksInYear  # Calculate how many full 52-element chunks we have
    fullNH_chunks = len(avg_NH_smoothed) // WeeksInYear  # Calculate how many full 52-element chunks we have
    fullSH_chunks = len(avg_SH_smoothed) // WeeksInYear  # Calculate how many full 52-element chunks we have
    fullPN_chunks = len(avg_smoothed_PN_subs) // WeeksInYear  # Calculate how many full 52-element chunks we have
    fullTN_chunks = len(avg_smoothed_TN_subs) // WeeksInYear  # Calculate how many full 52-element chunks we have
    fullTS_chunks = len(avg_smoothed_TS_subs) // WeeksInYear  # Calculate how many full 52-element chunks we have
    fullPS_chunks = len(avg_smoothed_PS_subs2) // WeeksInYear  # Calculate how many full 52-element chunks we have
    # Calculate global annual average using smoothed data
    #First global 
    avg_glob_smoothed_sliced = avg_glob_smoothed[:full_chunks * WeeksInYear]  # Ignore the remainder
    avg_glob_smoothed_reshaped = avg_glob_smoothed_sliced.reshape(full_chunks, WeeksInYear)
    AnnAvg_glob_smoothed = np.nanmean(avg_glob_smoothed_reshaped, axis=1)
    AnnAvg_glob_smoothedR = AnnAvg_glob_smoothed - AnnAvg_glob_smoothed[0]
    #Next NH
    avg_NH_smoothed_sliced = avg_NH_smoothed[:fullNH_chunks * WeeksInYear]  # Ignore the remainder
    avg_NH_smoothed_reshaped = avg_NH_smoothed_sliced.reshape(fullNH_chunks, WeeksInYear)
    AnnAvg_NH_smoothed = np.nanmean(avg_NH_smoothed_reshaped, axis=1)
    AnnAvg_NH_smoothedR = AnnAvg_NH_smoothed - AnnAvg_NH_smoothed[0]
    #Next SH
    avg_SH_smoothed_sliced = avg_SH_smoothed[:fullSH_chunks * WeeksInYear]  # Ignore the remainder
    avg_SH_smoothed_reshaped = avg_SH_smoothed_sliced.reshape(fullSH_chunks, WeeksInYear)
    AnnAvg_SH_smoothed = np.nanmean(avg_SH_smoothed_reshaped, axis=1)
    AnnAvg_SH_smoothedR = AnnAvg_SH_smoothed - AnnAvg_SH_smoothed[0]
    #Next PN
    avg_PN_smoothed_sliced = avg_smoothed_PN_subs[:fullPN_chunks * WeeksInYear]  # Ignore the remainder
    avg_PN_smoothed_reshaped = avg_PN_smoothed_sliced.reshape(fullPN_chunks, WeeksInYear)
    AnnAvg_PN_smoothed = np.nanmean(avg_PN_smoothed_reshaped, axis=1)
    AnnAvg_PN_smoothedR = AnnAvg_PN_smoothed - AnnAvg_PN_smoothed[0]
    #Next TN
    avg_TN_smoothed_sliced = avg_smoothed_TN_subs[:fullTN_chunks * WeeksInYear]  # Ignore the remainder
    avg_TN_smoothed_reshaped = avg_TN_smoothed_sliced.reshape(fullTN_chunks, WeeksInYear)
    AnnAvg_TN_smoothed = np.nanmean(avg_TN_smoothed_reshaped, axis=1)
    AnnAvg_TN_smoothedR = AnnAvg_TN_smoothed - AnnAvg_TN_smoothed[0]
    #Next TS
    avg_TS_smoothed_sliced = avg_smoothed_TS_subs[:fullTS_chunks * WeeksInYear]  # Ignore the remainder
    avg_TS_smoothed_reshaped = avg_TS_smoothed_sliced.reshape(fullTS_chunks, WeeksInYear)
    AnnAvg_TS_smoothed = np.nanmean(avg_TS_smoothed_reshaped, axis=1)
    AnnAvg_TS_smoothedR = AnnAvg_TS_smoothed - AnnAvg_TS_smoothed[0]
    #Next PS
    avg_PS_smoothed_sliced = avg_smoothed_PS_subs2[:fullPS_chunks * WeeksInYear]  # Ignore the remainder
    avg_PS_smoothed_reshaped = avg_PS_smoothed_sliced.reshape(fullPS_chunks, WeeksInYear)
    AnnAvg_PS_smoothed = np.nanmean(avg_PS_smoothed_reshaped, axis=1)
    AnnAvg_PS_smoothedR = AnnAvg_PS_smoothed - AnnAvg_PS_smoothed[0]
    # Aggregate for loop results
    # First weekly means
    WeeklyGlob_matrix[:, k] = avg_glob_smoothed
    WeeklyNH_matrix[:, k] = avg_NH_smoothed
    WeeklySH_matrix[:, k] = avg_SH_smoothed
    WeeklyPN_matrix[:, k] = avg_smoothed_PN_sub
    WeeklyTN_matrix[:, k] = avg_smoothed_TN_subs
    WeeklyPS_matrix[:, k] = avg_smoothed_PS_subs
    WeeklyTS_matrix[:, k] = avg_smoothed_TS_subs
    # Now global averages
    smoothed_matrix[:, k] = AnnAvg_glob_smoothed  # Replace with actual calculation for each iteration
    smoothed_matrixR[:, k] = AnnAvg_glob_smoothedR # Replace with actual calculation for each iteration
    smoothedNH_matrix[:, k] = AnnAvg_NH_smoothed  # Replace with actual calculation for each iteration
    smoothedNH_matrixR[:, k] = AnnAvg_NH_smoothedR # Replace with actual calculation for each iteration
    smoothedSH_matrix[:, k] = AnnAvg_SH_smoothed  # Replace with actual calculation for each iteration
    smoothedSH_matrixR[:, k] = AnnAvg_SH_smoothedR # Replace with actual calculation for each iteration
    smoothedPN_matrix[:, k] = AnnAvg_PN_smoothed  # Replace with actual calculation for each iteration
    smoothedPN_matrixR[:, k] = AnnAvg_PN_smoothedR # Replace with actual calculation for each iteration
    smoothedTN_matrix[:, k] = AnnAvg_TN_smoothed  # Replace with actual calculation for each iteration
    smoothedTN_matrixR[:, k] = AnnAvg_TN_smoothedR # Replace with actual calculation for each iteration
    smoothedTS_matrix[:, k] = AnnAvg_TS_smoothed  # Replace with actual calculation for each iteration
    smoothedTS_matrixR[:, k] = AnnAvg_TS_smoothedR # Replace with actual calculation for each iteration
    smoothedPS_matrix[:, k] = AnnAvg_PS_smoothed  # Replace with actual calculation for each iteration
    smoothedPS_matrixR[:, k] = AnnAvg_PS_smoothedR # Replace with actual calculation for each iteration
    
    print(k)

# Calculate statistics
AnnAvg_glob_smooth_mean = np.nanmean(smoothed_matrix, axis=1)
AnnAvg_glob_smooth_std = np.nanstd(smoothed_matrix, axis=1)
AnnAvg_glob_smooth_meanR = np.nanmean(smoothed_matrixR, axis=1)
AnnAvg_glob_smooth_stdR = np.nanstd(smoothed_matrixR, axis=1)
# Calculate annaul NH and SH statistics
AnnAvg_NH_smooth_mean = np.nanmean(smoothedNH_matrix, axis=1)
AnnAvg_NH_smooth_std = np.nanstd(smoothedNH_matrix, axis=1)
AnnAvg_SH_smooth_mean = np.nanmean(smoothedSH_matrix, axis=1)
AnnAvg_SH_smooth_std = np.nanstd(smoothedSH_matrix, axis=1)
AnnAvg_NH_smooth_meanR = np.nanmean(smoothedNH_matrixR, axis=1)
AnnAvg_NH_smooth_stdR = np.nanstd(smoothedNH_matrixR, axis=1)
AnnAvg_SH_smooth_meanR = np.nanmean(smoothedSH_matrixR, axis=1)
AnnAvg_SH_smooth_stdR = np.nanstd(smoothedSH_matrixR, axis=1)
# Calculate weekly statistics
Weekly_Glob_mean = np.nanmean(WeeklyGlob_matrix, axis=1)
Weekly_NH_mean = np.nanmean(WeeklyNH_matrix, axis=1)
Weekly_SH_mean = np.nanmean(WeeklySH_matrix, axis=1)
Weekly_PN_mean = np.nanmean(WeeklyPN_matrix, axis=1)
Weekly_PS_mean = np.nanmean(WeeklyPS_matrix, axis=1)
Weekly_TN_mean = np.nanmean(WeeklyTN_matrix, axis=1)
Weekly_TS_mean = np.nanmean(WeeklyTS_matrix, axis=1)
# Calculate annaul semi-hem statistics
AnnAvg_PN_smooth_mean = np.nanmean(smoothedPN_matrix, axis=1)
AnnAvg_PN_smooth_std = np.nanstd(smoothedPN_matrix, axis=1)
AnnAvg_TN_smooth_mean = np.nanmean(smoothedTN_matrix, axis=1)
AnnAvg_TN_smooth_std = np.nanstd(smoothedTN_matrix, axis=1)
AnnAvg_PN_smooth_meanR = np.nanmean(smoothedPN_matrixR, axis=1)
AnnAvg_PN_smooth_stdR = np.nanstd(smoothedPN_matrixR, axis=1)
AnnAvg_TN_smooth_meanR = np.nanmean(smoothedTN_matrixR, axis=1)
AnnAvg_TN_smooth_stdR = np.nanstd(smoothedTN_matrixR, axis=1)
AnnAvg_PS_smooth_mean = np.nanmean(smoothedPS_matrix, axis=1)
AnnAvg_PS_smooth_std = np.nanstd(smoothedPS_matrix, axis=1)
AnnAvg_TS_smooth_mean = np.nanmean(smoothedTS_matrix, axis=1)
AnnAvg_TS_smooth_std = np.nanstd(smoothedTS_matrix, axis=1)
AnnAvg_PS_smooth_meanR = np.nanmean(smoothedPS_matrixR, axis=1)
AnnAvg_PS_smooth_stdR = np.nanstd(smoothedPS_matrixR, axis=1)
AnnAvg_TS_smooth_meanR = np.nanmean(smoothedTS_matrixR, axis=1)
AnnAvg_TS_smooth_stdR = np.nanstd(smoothedTS_matrixR, axis=1)

plt.plot(smoothed_matrixR)


#%% Save data

if k >= 100:
    # Save global data
    combined_matrix = np.column_stack((year,AnnAvg_glob_smooth_mean,AnnAvg_glob_smooth_std,AnnAvg_glob_smooth_meanR,AnnAvg_glob_smooth_stdR))
    # Convert to a DataFrame for easy saving
    columns = ['Year', 'AnnAvg_glob_smooth_mean', 'AnnAvg_glob_smooth_std', 'AnnAvg_glob_smooth_meanR', 'AnnAvg_glob_smooth_stdR']
    df_combined = pd.DataFrame(combined_matrix, columns=columns)
    # Save the DataFrame to a CSV file
    output_filename = 'output/GlobMean_dD_dei_DasguptaCal_noBUDS.csv'  # Adding .csv extension
    df_combined.to_csv(output_filename, index=False)

    # Save semi-hemispheric data 
    combined_matrix2 = np.column_stack((date[1:],Weekly_PN_mean,Weekly_TN_mean,Weekly_TS_mean,Weekly_PS_mean))
    # Convert to a DataFrame for easy saving
    columns2 = ['Year', 'PN_smooth_mean', 'TN_smooth_mean', 'TS_smooth_mean', 'PS_smooth_mean']
    df_combined2 = pd.DataFrame(combined_matrix2, columns=columns2)
    # Save the DataFrame to a CSV file
    output_filename2 = 'output/SemiHemMean_dD_dei_DasguptaCal_noBUDS.csv'  # Adding .csv extension
    df_combined2.to_csv(output_filename2, index=False)
    
    # Save hemispheric data 
    combined_matrix3 = np.column_stack((date[1:],Weekly_Glob_mean,Weekly_NH_mean,Weekly_SH_mean))
    # Convert to a DataFrame for easy saving
    columns3 = ['Year', 'Glob_smooth_mean', 'NH_smooth_mean', 'SH_smooth_mean']
    df_combined3 = pd.DataFrame(combined_matrix3, columns=columns3)
    # Save the DataFrame to a CSV file
    output_filename3 = 'output/HemMean_dD_dei_DasguptaCal_noBUDS.csv'  # Adding .csv extension
    df_combined3.to_csv(output_filename3, index=False)
    
    # Save all MC iterations
    combined_matrix4 = np.column_stack((year,smoothed_matrix))
    pd.DataFrame(combined_matrix4).to_excel('output/GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx', index=False, header=False)



#%% Plot results!

# Load global annual dD data 
glob_ann_dD_data = pd.read_excel('data/glob_ann_dD.xlsx').to_numpy()
glob_ann_dD = glob_ann_dD_data[:,1]


# Plot semi-hemispheric averaged sites
plt.figure(1)
plt.gcf().set_facecolor('w')
plt.figure(dpi=500)
plt.gca().tick_params(axis='both', which='both', width=1)
plt.plot(date[1:], avg_smoothed_PN_singlesites, linewidth=1, label='30-90 N')
plt.plot(date[1:], avg_smoothed_PS_singlesites, linewidth=1, label='30-90 S')
plt.plot(date[1:], avg_smoothed_TN_singlesites, linewidth=1, label='0-30 N')
plt.plot(date[1:], avg_smoothed_TS_singlesites, linewidth=1, label='0-30 S')
plt.ylabel('dD (‰)')
plt.ylim(-96, -63)
plt.xlim(2005, 2025)
plt.legend()


# Plot semi-hemispheric averaged sites where times with only 1 site are adjusted
plt.figure(2)
plt.gcf().set_facecolor('w')
plt.figure(dpi=500)
plt.gca().tick_params(axis='both', which='both', width=1)
plt.plot(date[1:], avg_smoothed_PN_sub, linewidth=1, label='30-90 N')
plt.plot(date[1:], avg_smoothed_PS_sub, linewidth=1, label='30-90 S')
plt.plot(date[1:], avg_smoothed_TN_sub, linewidth=1, label='0-30 N')
plt.plot(date[1:], avg_smoothed_TS_sub, linewidth=1, label='0-30 S')
plt.ylabel('dD (‰)')
plt.ylim(-96, -63)
plt.xlim(2005, 2025)
plt.legend()

# Plot filled data gaps 
plt.figure(3)
plt.gcf().set_facecolor('w')
plt.figure(dpi=500)
plt.gca().tick_params(axis='both', which='both', width=1)
plt.plot(date[1:], avg_smoothed_PN_subs, linewidth=1, label='30-90 N')
plt.plot(date[1:], avg_smoothed_PS_subs2, linewidth=1, label='30-90 S')
plt.plot(date[1:], avg_smoothed_TN_subs, linewidth=1, label='0-30 N')
plt.plot(date[1:], avg_smoothed_TS_subs, linewidth=1, label='0-30 S')
plt.ylabel('dD (‰)')
plt.ylim(-96, -63)
plt.xlim(2005, 2025)
plt.legend()


#plot global mean dD
plt.figure(5)
plt.gcf().set_facecolor('w')
plt.figure(dpi=500)
plt.gca().tick_params(axis='both', which='both', width=1)
plt.plot((year+.5),AnnAvg_glob_smooth_mean, linewidth=1, label='Annual Glob Mean', color='blue')
plt.plot((year+.5),(AnnAvg_glob_smooth_mean + 2*AnnAvg_glob_smooth_std), linewidth=1, label='95% confidence', color='blue', linestyle='--')
plt.plot((year+.5),(AnnAvg_glob_smooth_mean - 2*AnnAvg_glob_smooth_std), linewidth=1, color='blue', linestyle='--')
plt.plot(date[1:],Weekly_Glob_mean, linewidth=1, label='Glob Mean')
#plt.plot((glob_ann_dD_data[:,0]+0.5),glob_ann_dD_data[:,1], linewidth=1, label='Glob Mean (old)')
plt.ylabel(r'$\delta D-CH_4$ (‰)')
plt.title(r'Global Mean $\delta D-CH_4$')
plt.grid(axis='x', linestyle='--', linewidth=0.5)  # Add vertical gridlines
plt.xlim(start_year, )
plt.legend()


#plot global mean dD RELATIVE TO 2005
plt.figure(6)
plt.gcf().set_facecolor('w')
plt.figure(dpi=500)
plt.gca().tick_params(axis='both', which='both', width=1)
plt.plot((year+.5),AnnAvg_glob_smooth_meanR, linewidth=1, label='Annual Glob Mean', color='blue')
plt.plot((year+.5),(AnnAvg_glob_smooth_meanR + 2*AnnAvg_glob_smooth_stdR), linewidth=1, label='95% confidence', color='blue', linestyle='--')
plt.plot((year+.5),(AnnAvg_glob_smooth_meanR - 2*AnnAvg_glob_smooth_stdR), linewidth=1, color='blue', linestyle='--')
#plt.plot((glob_ann_dD_data[:,0]+0.5),(glob_ann_dD_data[:,1] - glob_ann_dD_data[10,1]), linewidth=1, label='Glob Mean (old)')
plt.ylabel(r'$\Delta \delta D-CH_4$ (‰)')
plt.title(r'Change in $\delta D-CH_4$')
plt.xlim(start_year, )
plt.legend()


#%% Supplementary figure in main text

# Load site coordinates: 
df_sites_b = siteinfo_filtered[[0, 3, 4]].copy()
df_sites_b.columns = ['site', 'latitude', 'longitude']

from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(12, 10), dpi=500)
gs = gridspec.GridSpec(
    2, 2,
    figure=fig,
    width_ratios=[1, 1],
    height_ratios=[1, 1]
)

# top-left: world map
ax_map = fig.add_subplot(gs[0, 0])
m = Basemap(projection='robin', lon_0=0, resolution='c', ax=ax_map)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='lightgray', lake_color='lightblue')
m.drawmapboundary(fill_color='lightblue')

# adjust map shape here
ax_map.set_aspect(1.5)  # >1 = taller, <1 = wider

# color scheme
color_default = "red"
color_imau    = "blue"
color_nipr    = "green"
color_mpi     = "purple"

# plot sites: first non-default
for idx, row in df_sites_b.iterrows():
    sitename = str(row['site'])
    lat, lon = row['latitude'], row['longitude']
    prefix = sitename[:3].upper()
    suffix = sitename[3:].upper() if len(sitename) > 3 else ""

    if suffix in ["IMAU", "NIPR", "MPI"]:
        if suffix == "IMAU":
            c = color_imau
        elif suffix == "NIPR":
            c = color_nipr
        elif suffix == "MPI":
            c = color_mpi

        x, y = m(lon, lat)
        m.scatter(x, y, c=c, s=50, edgecolor='k', zorder=5)
        ax_map.text(x+300000, y+100000, prefix, fontsize=9,
                    ha='left', va='bottom', fontweight='bold')

# now plot defaults on top
for idx, row in df_sites_b.iterrows():
    sitename = str(row['site'])
    lat, lon = row['latitude'], row['longitude']
    prefix = sitename[:3].upper()
    suffix = sitename[3:].upper() if len(sitename) > 3 else ""

    if suffix not in ["IMAU", "NIPR", "MPI"]:
        x, y = m(lon, lat)
        m.scatter(x, y, c=color_default, s=50, edgecolor='k', zorder=6)  # higher zorder
        ax_map.text(x+300000, y+100000, prefix, fontsize=9,
                    ha='left', va='bottom', fontweight='bold')

ax_map.set_title("Site Locations by Affiliation", fontsize=15, pad=15)

# top-right: weekly number of sites
ax1 = fig.add_subplot(gs[0, 1])
ax1.plot(date[:-1], sites_PS, linewidth=1, label='30-90 S (PS)')
ax1.plot(date[:-1], sites_TS, linewidth=1, label='0-30 S (TS)')
ax1.plot(date[:-1], sites_TN, linewidth=1, label='0-30 N (TN)')
ax1.plot(date[:-1], sites_PN, linewidth=1, label='30-90 N (PN)')
ax1.plot(date[:-1], sites_total, linewidth=2, label='Total # of sites')
ax1.set_ylabel('n', fontsize=15)
ax1.set_xlim(2005, 2024)
ax1.legend(fontsize=10)
ax1.set_title('Weekly Number of Sites', fontsize=15, pad=15)

# bottom-left: hemispheric seasonal cycle
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(date[1:], (Weekly_NH_mean - Weekly_SH_mean), color='blue')
ax2.set_ylabel('NH - SH ${\delta}D$-CH$_4$ (‰)', fontsize=15)
ax2.set_xlim(2005, 2024)
ax2.set_title('Hemispheric Seasonal Cycle', fontsize=15, pad=15)
ax2.set_xlabel('Year', fontsize=15)

# bottom-right: weekly averages
ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(date[1:], Weekly_PS_mean, linewidth=1, label='30-90 S (PS)')
ax3.plot(date[1:], Weekly_TS_mean, linewidth=1, label='0-30 S (TS)')
ax3.plot(date[1:], Weekly_TN_mean, linewidth=1, label='0-30 N (TN)')
ax3.plot(date[1:], Weekly_PN_mean, linewidth=1, label='30-90 N (PN)')
ax3.set_ylabel('${\delta}D$-CH$_4$ (‰)', fontsize=15)
ax3.set_xlim(2005, 2024)
ax3.set_ylim(-97, -66)
ax3.set_title('Weekly Averages of MC Iterations', fontsize=15, pad=15)

# adjust ticks
for ax in [ax1, ax2, ax3]:
    ax.tick_params(axis='both', which='both', width=1, labelsize=15)
    ax.set_xticks(np.arange(2005, 2025, 4))

plt.tight_layout(pad=1.5)
plt.show()



