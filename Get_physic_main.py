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

test_nc_file = "Test/bass2_simple_2017-11.nc"
mesh_file = "Test/mesh.nc"

if os.path.exists(test_nc_file):
    print(f"Found test NetCDF file: {test_nc_file}")
    data_file = test_nc_file
    mesh_file = test_nc_file  # Using same file for now
    
    # Calculate transport
    print("Calculating transport...")
    T, tims = calc_transport(vert, pt1, pt2, dlev, dinc, rimn, data_file, mesh_file, 1)
    
    # Write transport file
    output_file = "Test/Transport_test.nc"
    print(f"Writing transport file: {output_file}")
    write_nc_transport(pt1, pt2, lr, tims, T, fcid, output_file)
    
    print(f"\nTransport processing completed!")
    print(f"Transport array shape: {T.shape}")
    print(f"Time steps: {len(tims)}")
    print(f"Output file: {output_file}")
    
else:
    print(f"Test NetCDF file not found: {test_nc_file}")
    print("Skipping transport processing - need hydrodynamic data files")
    print("In a real implementation, you would:")
    print("1. Have NetCDF files with velocity data (U, V components)")
    print("2. Have a mesh file with grid coordinates")
    print("3. Process multiple time files in a loop")
    print("4. Combine results into final transport file")
    
