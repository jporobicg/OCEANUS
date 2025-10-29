# Get_variables_main.py
# Main script for processing variables (temperature, salinity, vertical flux)
# Based on Get_Variables_SS.m

import numpy as np
import os
import glob
import pickle
from read_boxes import read_boxes
from read_faces import read_faces
from box_averages import box_averages
from write_all_variables import write_all_variables

# File paths
BGM_FILE = "Test/SEAP_NESP_final_ll_fixed.bgm"
data_file = "Test/bass2_simple_2017-11.nc"
mesh_file = "Test/bass2_simple_2017-11.nc"

# Parameters
dlev = np.array([0, 25, 50, 100, 250, 400, 700])
varn = ['w', 'salt', 'temp']

print("Reading box geometry...")
nbox, nface, bid, cent, b_area, vert, iface, botz = read_boxes(BGM_FILE)
print("Reading face definitions...")
nulr, nupt1, nupt2 = read_faces(nbox, nface, bid, vert, iface, BGM_FILE)

print(f"Number of boxes: {nbox}")
print(f"Number of faces: {nface}")
print(f"Number of real faces: {len(nulr)}")

print(f"\nProcessing variables: {varn}")



# Process each variable
for year in range(2017, 2025):
    # Initialize arrays to store all variable data
    all_temp_data = None
    all_salt_data = None
    all_w_data = None
    for month in range(1, 13):
        if month < 10:
            month = f"0{month}"
        else:
            month = str(month)
        file_path = f"/datasets/work/nesp-gda-owf-ra/work/data/processed/BASS2_ocean/2017-2024_historical_v2/bass2_simple_{year}-{month}.nc"
        print(f"Processing file: {file_path}")
    
        for v in range(len(varn)):
            avname = varn[v]
            print(f"\nProcessing variable: {avname}")
            
            # Call box_averages for this variable
            box_averages(vert, avname, dlev, file_path, mesh_file, month, year)
        
    ## getting the data for the combined output
    for v in range(len(varn)):
        avname = varn[v]
        print(f"\nProcessing variable: {avname}")
        pattern = f"*_{year}_{avname}_SS_Second_step.npz"
        t_files = glob.glob(pattern)
        if len(t_files) > 0:
            print(f"Found {len(t_files)} processed file(s) for {avname}")
            Var_avg = np.concatenate([np.load(file)['Var_avg'] for file in t_files])
            tims = np.concatenate([np.load(file)['tims'] for file in t_files])
            # Store the data for the combined output
            if avname == 'temp':
                all_temp_data = Var_avg
            elif avname == 'salt':
                all_salt_data = Var_avg
            elif avname == 'w':
                all_w_data = Var_avg
            all_times = tims
    ## Writing the output files
    if all_temp_data is not None and all_salt_data is not None and all_w_data is not None:
        output_file = f"Test/bass2_simple_{year}.nc"
    print(f"\nWriting all variables to: {output_file}")
    write_all_variables(all_times, bid, all_temp_data, all_salt_data, all_w_data, output_file)
    print(f"Variable processing completed for year {year}!")
else:
    print(f"Error: Not all variables were processed successfully for year {year}") 