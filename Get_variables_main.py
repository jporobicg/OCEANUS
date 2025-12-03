# Get_variables_main.py
# Main script for processing variables (temperature, salinity, vertical flux)
# Based on Get_Variables_SS.m

import numpy as np
import os
import glob
import pickle
import sys
from read_boxes import read_boxes
from read_faces import read_faces
from box_averages import box_averages
from write_all_variables import write_all_variables

# File paths
BGM_FILE = "Test/SEAP_NESP_final_ll_fixed.bgm"
mesh_file = "Test/mesh.nc"

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

# Get single year from CLI and prepare output dir
Year = int(sys.argv[1])
year = Year
output_folder = f"/datasets/work/oa-alantis/work/NESP_hydro/New_version/OCEANUS/temporal/{year}/variables/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder, exist_ok=True)

# Monthly processing: compute and save per-variable NPZs
for m in range(1, 13):
    month_str = f"{m:02d}"
    
    file_path = f"/datasets/work/nesp-gda-owf-ra/work/data/covariate_data/processed/BASS2_ocean/2017-2024_historical_v2/bass2_simple_{year}-{month_str}.nc"
    print(f"Processing file: {file_path}")

    for avname in varn:
        npz_file = os.path.join(output_folder, f"{month_str}_{year}_{avname}_SS_Second_step.npz")
        
        if os.path.exists(npz_file):
            print(f"  -> Skipping {avname} {month_str}-{year} - NPZ exists: {npz_file}")
            continue
        print(f"  -> Processing variable: {avname}")
        box_averages(vert, avname, dlev, file_path, mesh_file, month_str, year, output_folder)

# Year concatenation per variable, then combined write
all_temp_data = None
all_salt_data = None
all_w_data = None
all_times = None

for avname in varn:
    pattern = os.path.join(output_folder, f"*_{year}_{avname}_SS_Second_step.npz")
    t_files = sorted(glob.glob(pattern))
    t_files.sort()
    if len(t_files) == 0:
        print(f"No monthly NPZ files found for variable {avname} in {year}")
        continue
    print(f"Found {len(t_files)} monthly NPZ files for {avname}")
    # Concatenate along time axis (axis=1) for Var_avg (nbox, ntm, nlay)
    var_arrays = []
    time_arrays = []
    for f in t_files:
        data = np.load(f)
        var_arrays.append(data['Var_avg'])
        time_arrays.append(data['tims'])
        data.close()
    Var_avg = np.concatenate(var_arrays, axis=1)
    tims = np.concatenate(time_arrays, axis=0)
    if avname == 'temp':
        all_temp_data = Var_avg
    elif avname == 'salt':
        all_salt_data = Var_avg
    elif avname == 'w':
        all_w_data = Var_avg
    if all_times is None:
        all_times = tims

if all_temp_data is not None and all_salt_data is not None and all_w_data is not None and all_times is not None:
    output_file = os.path.join(output_folder, f"Variables_{year}.nc")
    print(f"\nWriting all variables to: {output_file}")
    write_all_variables(all_times, bid, all_temp_data, all_salt_data, all_w_data, output_file)
    print(f"Variable processing completed for year {year}!")
else:
    print(f"Error: Not all variables were processed successfully for year {year}")
