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

# Initialize arrays to store all variable data
all_temp_data = None
all_salt_data = None
all_w_data = None
all_times = None

# Process each variable
for v in range(len(varn)):
    avname = varn[v]
    print(f"\nProcessing variable: {avname}")
    
    # Call box_averages for this variable
    box_averages(vert, avname, dlev, data_file, mesh_file, 1)
    
    # Look for the processed file
    pattern = f"*{avname}_SS_Second_step.pkl"
    t_files = glob.glob(pattern)
    
    if len(t_files) > 0:
        print(f"Found {len(t_files)} processed file(s) for {avname}")
        
        # Load the processed data
        with open(t_files[0], 'rb') as f_pkl:
            data = pickle.load(f_pkl)
            Av_final = data['Var_avg']
            nctime = data['tims']
        
        print(f"Loaded data shape: {Av_final.shape}")
        print(f"Time steps: {len(nctime)}")
        
        # Store the data for the combined output
        if avname == 'temp':
            all_temp_data = Av_final
        elif avname == 'salt':
            all_salt_data = Av_final
        elif avname == 'w':
            all_w_data = Av_final
        
        # Use the first variable's time array for all
        if all_times is None:
            all_times = nctime
            
    else:
        print(f"No processed files found for {avname}")

# Write all variables to a single NetCDF file
if all_temp_data is not None and all_salt_data is not None and all_w_data is not None:
    output_file = "Test/Variables_combined.nc"
    print(f"\nWriting all variables to: {output_file}")
    write_all_variables(all_times, bid, all_temp_data, all_salt_data, all_w_data, output_file)
    print("Variable processing completed!")
else:
    print("Error: Not all variables were processed successfully") 