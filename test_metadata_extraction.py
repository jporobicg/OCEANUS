#!/usr/bin/env python3
"""
Test script to verify metadata extraction from NetCDF file
"""

import netCDF4 as nc
import sys
import os

# Define the function directly to avoid argparse conflicts
def extract_variable_metadata(nc_file_path, var_name):
    """
    Extract metadata (long_name, units, FillValue) from source NetCDF file.
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

# Variable sets
STANDARD_VARS = ['w', 'salt', 'temp']
EXTRA_VARS = ['DetR_N', 'DetR_N_sed', 'DOR_N', 'DOR_N_sed', 'NH4', 'NO3', 
              'PAR', 'PAR_z', 'PhyL_N', 'PhyS_N', 'porosity', 'ZooL_N', 'ZooS_N']

# Test file
test_file = "/home/por07g/Documents/Projects/NESP_ParksAustralia/Hydro/bass2_simple_2017-11.nc"

print("="*70)
print("Testing Metadata Extraction")
print("="*70)
print(f"\nTest file: {test_file}")
print(f"File exists: {os.path.exists(test_file)}")

# First, let's see what variables are in the file
print("\n" + "="*70)
print("Variables available in NetCDF file:")
print("="*70)
try:
    with nc.Dataset(test_file, 'r') as ds:
        print(f"\nTotal variables: {len(ds.variables)}")
        print("\nVariable list:")
        for var_name in sorted(ds.variables.keys()):
            var = ds.variables[var_name]
            dims = var.dimensions
            shape = var.shape
            print(f"  - {var_name:20s} : dims={dims}, shape={shape}")
            
            # Show attributes for data variables (not coordinates)
            if len(dims) > 1:  # Likely a data variable
                attrs = var.ncattrs()
                if attrs:
                    print(f"    Attributes: {', '.join(attrs)}")
                    if 'long_name' in attrs:
                        print(f"      long_name: {var.getncattr('long_name')}")
                    if 'units' in attrs:
                        print(f"      units: {var.getncattr('units')}")
                    if '_FillValue' in attrs or 'FillValue_' in attrs or 'missing_value' in attrs:
                        fill_attr = None
                        if '_FillValue' in attrs:
                            fill_attr = '_FillValue'
                        elif 'FillValue_' in attrs:
                            fill_attr = 'FillValue_'
                        elif 'missing_value' in attrs:
                            fill_attr = 'missing_value'
                        if fill_attr:
                            print(f"      {fill_attr}: {var.getncattr(fill_attr)}")
                print()
except Exception as e:
    print(f"Error reading file: {e}")
    sys.exit(1)

# Test metadata extraction for standard variables
print("\n" + "="*70)
print("Testing Standard Variables Metadata Extraction:")
print("="*70)
for var_name in STANDARD_VARS:
    print(f"\nVariable: {var_name}")
    metadata = extract_variable_metadata(test_file, var_name)
    print(f"  long_name: {metadata['long_name']}")
    print(f"  units: {metadata['units']}")
    print(f"  FillValue_: {metadata['FillValue_']}")

# Test metadata extraction for extra variables
print("\n" + "="*70)
print("Testing Extra Variables Metadata Extraction:")
print("="*70)
for var_name in EXTRA_VARS:
    print(f"\nVariable: {var_name}")
    metadata = extract_variable_metadata(test_file, var_name)
    print(f"  long_name: {metadata['long_name']}")
    print(f"  units: {metadata['units']}")
    print(f"  FillValue_: {metadata['FillValue_']}")

print("\n" + "="*70)
print("Test Complete")
print("="*70)
