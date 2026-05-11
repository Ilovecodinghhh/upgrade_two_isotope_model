#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 10:03:55 2024

@author: ryoung
"""

# This code creates an array of smoothed curves for each site used in the global mean calculations by estimating 
# atmospheric/analytical uncertainty related to residual standard deviations from the smoothed curve 


locals().clear()
import ccg_filter
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
    # 2. Atmospheric: Adds error to monthly means
    # 3. No analysis uncertainty needed because this is embedded in atmospheric?

#%% First, calculate MBL values

# List of input filenames
filenames = [
    "alt_01D0_dat.txt", "altMPI_01D0_dat.txt", "asc_01D0_dat.txt", "ato_01D0_dat.txt",
    "azr_01D0_dat.txt", "bal_01D0_dat.txt", "bikMPI_01D0_dat.txt", "brw_01D0_dat.txt",
    "brwIMAU_01D0_dat.txt", "bsc_01D0_dat.txt", "cba_01D0_dat.txt", "cgo_01D0_dat.txt",
    "cgoIMAU_01D0_dat.txt", "cvoMPI_01D0_dat.txt", "eom_01D0_dat.txt", "gvnMPI_01D0_dat.txt", "gvnIMAU_01D0_dat.txt",
    "jfjMPI_01D0_dat.txt", "kjnMPI_01D0_dat.txt", "kum_01D0_dat.txt", "lef_01D0_dat.txt",
    "mhd_01D0_dat.txt", "mlo_01D0_dat.txt", "mloIMAU_01D0_dat.txt", "namMPI_01D0_dat.txt",
    "nyaNIPR_01D0_dat.txt", "oxkMPI_01D0_dat.txt", "sisMPI_01D0_dat.txt", "smo_01D0_dat.txt",
    "smoIMAU_01D0_dat.txt", "spo_01D0_dat.txt", "syoNIPR_01D0_dat.txt", "vrsMPI_01D0_dat.txt",
    "zep_01D0_dat.txt",  "zepIMAU_01D0_dat.txt", "zotMPI_01D0_dat.txt"]

# Dictionaries to store data for the matrices
smoothed_dict = {}
trend_dict = {}
rsd_dict = {} # Residual standard deviation

# Loop through each file
for filename in filenames:
    # Load data
    data = np.loadtxt(f'data/{filename}')
    DecDate = data[:, 0]
    dD = data[:, 1]
    # Create the ccg_filter object
    filt = ccg_filter.ccgFilter(DecDate, dD, shortterm=150, longterm=667, sampleinterval=7, numpolyterms=3, numharmonics=4)
    
    # Get x values, smoothed curve, and trend
    x_interp = filt.xinterp
    smoothed_curve = filt.getSmoothValue(x_interp)
    trend = filt.getTrendValue(x_interp)
    # Extract rsd1 and rsd2
    rsd2 = filt.rsd2
    # Repeat rsd1 and rsd2 to match the length of x_interp, smoothed_curve, and trend
    rsd2_array = np.full_like(x_interp, rsd2)

    # Store the data into the dictionaries
    base_name = filename.replace("_dat.txt", "")
    smoothed_dict[f'{base_name}_x_interp'] = x_interp
    smoothed_dict[f'{base_name}_smoothed_curve'] = smoothed_curve
    trend_dict[f'{base_name}_x_interp'] = x_interp
    trend_dict[f'{base_name}_trend'] = trend
    # Store the rsd1 and rsd2 values
    rsd_dict[f'{base_name}_rsd'] = rsd2

    # Save individual output file with rsd1 and rsd2 included
    output_filename = f'output/{filename.replace("_dat.txt", "_curves_rsd.txt")}'
    np.savetxt(output_filename, np.column_stack((x_interp, smoothed_curve, trend, rsd2_array)), fmt='%.6f', 
               header="x_interp smoothed_curve trend rsd2")


# Now that all data is collected, convert to DataFrames and save
# This will handle varying lengths by padding with NaN
df_smoothed = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in smoothed_dict.items()]))
df_trend = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in trend_dict.items()]))
df_rsd = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in rsd_dict.items()]))


#%% Repeat but save values at sample resolutoin

# Dictionaries to store data for the matrices
smoothed_dict_sampleres = {}
trend_dict_sampleres = {}
rsd_dict_sampleres = {} # Residual standard deviation

# Loop through each file
for filename in filenames:
    # Load data
    data = np.loadtxt(f'data/{filename}')
    DecDate = data[:, 0]
    dD = data[:, 1]
    # Create the ccg_filter object without sampleinterval for exact sample resolution
    filt = ccg_filter.ccgFilter(DecDate, dD, shortterm=150, longterm=667, numpolyterms=3, numharmonics=4)
    
    # Get x values, smoothed curve, and trend at the sample resolution
    x_interp = filt.xp  # Original data points
    smoothed_curve = filt.getSmoothValue(x_interp)
    trend = filt.getTrendValue(x_interp)

    # Store the data into the dictionaries
    base_name = filename.replace("_dat.txt", "")
    smoothed_dict_sampleres[f'{base_name}_x_interp'] = x_interp
    smoothed_dict_sampleres[f'{base_name}_smoothed_curve'] = smoothed_curve
    trend_dict_sampleres[f'{base_name}_x_interp'] = x_interp
    trend_dict_sampleres[f'{base_name}_trend'] = trend

    # Save individual output file
    output_filename = f'MBL_SampleRes/{filename.replace("_dat.txt", "_curves_SampleRes.txt")}'
    np.savetxt(output_filename, np.column_stack((x_interp, smoothed_curve, trend)), fmt='%.6f', 
               header="x_interp smoothed_curve trend")

# Convert data to DataFrames and save
df_smoothed_sampleres = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in smoothed_dict_sampleres.items()]))
df_trend_sampleres = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in trend_dict_sampleres.items()]))

# Save the DataFrames to CSV files
df_smoothed_sampleres.to_csv('MBL_SampleRes/all_smoothed_curves_matrix_sampleres.csv', index=False)
df_trend_sampleres.to_csv('MBL_SampleRes/all_trend_curves_matrix_sampleres.csv', index=False)


#%% Next, repeat code but with MC iterations adding random error based on rsd to the smoothed resolution data.

# Loop through each file
for filename in filenames:
    # Load data
    data = np.loadtxt(f'MBL_SampleRes/{filename.replace("_dat.txt", "_curves_SampleRes.txt")}')
    rds_data = np.loadtxt(f'output/{filename.replace("_dat.txt", "_curves_rsd.txt")}')
    DecDate = data[:, 0]
    dD = data[:, 1]
    rsd = rds_data[3, 3]
    
    # Number of MC iterations
    iterations = 1000
    
    # List to hold smoothed curves for each iteration
    smoothed_curves = []
    
    for k in range(iterations):
        def get_year_month(decimal_date):
            year = int(decimal_date)
            month = int((decimal_date - year) * 12) + 1
            return year, month
        
        atm_error = np.zeros_like(DecDate)
        current_random_value = np.random.normal(0, rsd)
        last_year, last_month = get_year_month(DecDate[0])
        
        for i in range(len(DecDate)):
            year, month = get_year_month(DecDate[i])
            if year != last_year or month != last_month:
                current_random_value = np.random.normal(0, rsd)
                last_year, last_month = year, month
            atm_error[i] = current_random_value
        
        filt = ccg_filter.ccgFilter(DecDate, (dD + atm_error), shortterm=150, longterm=667, sampleinterval=7, numpolyterms=3, numharmonics=4)
        smoothed_curve = filt.getSmoothValue(filt.xinterp)
        
        # Append only smoothed curve for this iteration
        smoothed_curves.append(smoothed_curve)
    
    # Convert the smoothed curves list into a matrix (each column is an iteration)
    smoothed_curves_matrix = np.column_stack(smoothed_curves)
    
    # Concatenate x_interp (DecDate) as the first column, followed by the smoothed curves
    mc_results_matrix = np.column_stack((filt.xinterp, smoothed_curves_matrix))
    
    # Save the matrix to the Output folder
    output_filename = f'Output/{filename.replace("_dat.txt", "_smoothedMC.txt")}'
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    np.savetxt(output_filename, mc_results_matrix, fmt='%.6f', header="DecDate Smoothed_Curve_MC1 Smoothed_Curve_MC2 ... Smoothed_Curve_MC10")
    
    print(f'Saved {output_filename}')

