# Get_physic_main.py
# Main script for OCEANUS transport processing (Python version)

import numpy as np
import os
import glob
import sys
from read_boxes import read_boxes
from read_faces import read_faces
from calculate_transport import calc_transport
from write_nc_transport import write_nc_transport

# Path to the BGM file (update as needed)
BGM_FILE = "Test/SEAP_NESP_final_ll_fixed.bgm"

# Read box geometry from the BGM file
print("Reading box geometry...")
nbox, nface, bid, cent, b_area, vert, iface, botz = read_boxes(BGM_FILE)

# Read face definitions and relationships
print("Reading face definitions...")
nulr, nupt1, nupt2 = read_faces(nbox, nface, bid, vert, iface, BGM_FILE)

# Assign to more readable variable names (as in MATLAB)
lr = nulr      # Neighbor layers
pt1 = nupt1    # Face point 1
pt2 = nupt2    # Face point 2

# Find real faces (not NaN)
irealfaces = np.where(~np.isnan(nupt1[:, 0]))[0]
fcid = irealfaces  # Python is 0-based

# Integration parameters
rimn = 10
dinc = 1

# Depth levels (can be adjusted for your model)
dlev = np.array([0, 20, 50, 100, 250, 700, 4743])

# Print basic results
print(f"Number of boxes: {nbox}")
print(f"Number of faces: {nface}")
print(f"Number of real faces: {len(irealfaces)}")

# Transport processing section
print("\n" + "="*50)
print("TRANSPORT PROCESSING SECTION")
print("="*50)

# Define paths

mesh_file = "Test/mesh.nc"

# Get single year from CLI
Year = int(sys.argv[1])
output_folder = f"/datasets/work/oa-alantis/work/NESP_hydro/New_version/OCEANUS/temporal/{Year}/"
# output folder with the year. create the folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# calculate transport for each month for the provided year
year = Year
for month in range(1, 13):
    if month < 10:
        month = f"0{month}"
    else:
        month = str(month)
    file_path = f"/datasets/work/nesp-gda-owf-ra/work/data/processed/BASS2_ocean/2017-2024_historical_v2/bass2_simple_{year}-{month}.nc"
    ## check if the processed file exists
    processed_file = f"{output_folder}/{month}_{year}_SS_second_Step.npz"
    if os.path.exists(processed_file):
        print(f"Processed file {processed_file} already exists, skipping")
        continue
    else:
        print(f"Processing file: {file_path}")
        calc_transport(vert, pt1, pt2, dlev, dinc, rimn, file_path, mesh_file, month, year, output_folder)


# load and concatenate NPZ files by year
npz_files = glob.glob(f"{output_folder}/*_{year}_SS_second_Step.npz")
npz_files.sort()  # relies on zero-padded month prefix for correct order
T = np.concatenate([np.load(file)['T'] for file in npz_files], axis=1)
tims = np.concatenate([np.load(file)['tims'] for file in npz_files], axis=0)
output_file = f"{output_folder}/bass2_simple_{year}.nc"
write_nc_transport(pt1, pt2, lr, tims, T, fcid, output_file)



    
