# Get_variables_main.py
# Main script for processing variables (temperature, salinity, vertical flux)
# Based on Get_Variables_SS.m

import numpy as np
import os
import glob
import pickle
import sys
import argparse
import netCDF4 as nc
from read_boxes import read_boxes
from read_faces import read_faces
from box_averages import box_averages
from write_all_variables import write_all_variables

# File paths
BGM_FILE = "Test/SEAP_NESP_final_ll_fixed.bgm"
mesh_file = "Test/mesh.nc"

# Parameters
dlev = np.array([0, 25, 50, 100, 250, 400, 700])

# Variable sets
STANDARD_VARS = ['w', 'salt', 'temp']
EXTRA_VARS = ['DetR_N', 'DetR_N_sed', 'DOR_N', 'DOR_N_sed', 'NH4', 'NO3', 
              'PAR', 'PAR_z', 'PhyL_N', 'PhyS_N', 'porosity', 'ZooL_N', 'ZooS_N']

def extract_variable_metadata(nc_file_path, var_name):
    """
    Extract metadata (long_name, units, FillValue) from source NetCDF file.
    
    Args:
        nc_file_path: Path to source NetCDF file
        var_name: Name of the variable in the NetCDF file
    
    Returns:
        dict: Dictionary with 'long_name', 'units', and 'FillValue_' keys
    """
    try:
        with nc.Dataset(nc_file_path, 'r') as ds:
            if var_name not in ds.variables:
                print(f"Warning: Variable {var_name} not found in {nc_file_path}")
                return {
                    'long_name': var_name,
                    'units': '',
                    'FillValue_': -1e20
                }
            
            var = ds.variables[var_name]
            
            # Extract attributes
            long_name = var.getncattr('long_name') if 'long_name' in var.ncattrs() else var_name
            units = var.getncattr('units') if 'units' in var.ncattrs() else ''
            
            # Try different fill value attribute names
            fill_value = -1e20  # default
            if '_FillValue' in var.ncattrs():
                fill_value = var.getncattr('_FillValue')
            elif 'FillValue_' in var.ncattrs():
                fill_value = var.getncattr('FillValue_')
            elif 'missing_value' in var.ncattrs():
                fill_value = var.getncattr('missing_value')
            
            return {
                'long_name': str(long_name),
                'units': str(units),
                'FillValue_': float(fill_value)
            }
    except Exception as e:
        print(f"Warning: Could not extract metadata for {var_name} from {nc_file_path}: {e}")
        return {
            'long_name': var_name,
            'units': '',
            'FillValue_': -1e20
        }

# Parse command line arguments
parser = argparse.ArgumentParser(description='Process variables for OCEANUS')
parser.add_argument('year', type=int, help='Year to process')
parser.add_argument('--extra', action='store_true', help='Process extra variables instead of standard ones')
args = parser.parse_args()

year = args.year
extra_variables = args.extra

# Select variable set based on flag
if extra_variables:
    varn = EXTRA_VARS
    print("\n=== Processing EXTRA variables ===")
else:
    varn = STANDARD_VARS
    print("\n=== Processing STANDARD variables ===")

print("Reading box geometry...")
nbox, nface, bid, cent, b_area, vert, iface, botz = read_boxes(BGM_FILE)
print("Reading face definitions...")
nulr, nupt1, nupt2 = read_faces(nbox, nface, bid, vert, iface, BGM_FILE)

print(f"Number of boxes: {nbox}")
print(f"Number of faces: {nface}")
print(f"Number of real faces: {len(nulr)}")

print(f"\nProcessing variables: {varn}")

# Prepare output dir
if extra_variables:
    output_folder = f"/datasets/work/oa-alantis/work/NESP_hydro/New_version/OCEANUS/temporal/{year}/variables_extra/"
else:
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
# Collect data and metadata for all variables
variables_info = []
all_times = None

# Get metadata from first available NetCDF file (use first month as reference)
first_month_file = f"/datasets/work/nesp-gda-owf-ra/work/data/covariate_data/processed/BASS2_ocean/2017-2024_historical_v2/bass2_simple_{year}-01.nc"

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
    
    # Extract metadata from source NetCDF file
    metadata = extract_variable_metadata(first_month_file, avname)
    
    # Store variable information
    variables_info.append({
        'name': avname,
        'data': Var_avg,
        'long_name': metadata['long_name'],
        'units': metadata['units'],
        'FillValue_': metadata['FillValue_']
    })
    
    if all_times is None:
        all_times = tims
    elif not np.array_equal(all_times, tims):
        print(f"Warning: Time arrays differ for variable {avname}")

# Write all variables to NetCDF
if len(variables_info) > 0 and all_times is not None:
    if extra_variables:
        output_file = os.path.join(output_folder, f"Variables_extra_{year}.nc")
    else:
        output_file = os.path.join(output_folder, f"Variables_{year}.nc")
    
    print(f"\nWriting {len(variables_info)} variables to: {output_file}")
    write_all_variables(all_times, bid, variables_info, output_file)
    print(f"Variable processing completed for year {year}!")
else:
    print(f"Error: No variables were processed successfully for year {year}")
