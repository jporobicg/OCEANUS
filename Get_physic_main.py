# Get_physic_main.py
# Main script for OCEANUS transport processing (Python version)

import numpy as np
import os
import glob
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
output_folder = "Test/transport"
mesh_file = "Test/mesh.nc"

# Process year by year
for year in range(2017, 2025):  # 2017 to 2024 inclusive
    # calculate transport for each month
    for month in range(1, 13):
        if month < 10:
            month = f"0{month}"
        else:
            month = str(month)
        file_path = f"/datasets/work/nesp-gda-owf-ra/work/data/processed/BASS2_ocean/2017-2024_historical_v2/bass2_simple_{year}-{month}.nc"
        print(f"Processing file: {file_path}")
        calc_transport(vert, pt1, pt2, dlev, dinc, rimn, file_path, mesh_file, month, year)
        
        npz_file = f"{month}{year}_SS_second_Step.npz"
        if os.path.exists(npz_file):
            print(f"  -> Skipping {file_path} - NPZ file exists: {npz_file}")
            continue
    # load and concatenate NPZ files by year
    npz_files = glob.glob(f"{year}_*_SS_second_Step.npz")
    npz_files.sort()
    T = np.concatenate([np.load(file)['T'] for file in npz_files])
    tims = np.concatenate([np.load(file)['tims'] for file in npz_files])
    output_file = f"Test/bass2_simple_{year}.nc"
    write_nc_transport(pt1, pt2, lr, tims, T, fcid, output_file)
    print(f"Loaded {len(npz_files)} NPZ files for year {year}")
    print(f"Transport array shape: {T.shape}")
    print(f"Time steps: {len(tims)}")
    print(f"Output file: {output_file}")


    
