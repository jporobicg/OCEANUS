#!/usr/bin/env python3
"""
Resample OCEANUS outputs from hourly to 12-hourly resolution
Averages every 12 consecutive time steps and sets correct time variable
"""

import numpy as np
import netCDF4 as nc
import os
import sys
from datetime import datetime, timedelta

def create_12hour_time_array(start_year, ntime_12h):
    """
    Create time array in seconds since 1990-01-01 00:00:00 +10 (AEST)
    
    Args:
        start_year: Starting year (e.g., 2017)
        ntime_12h: Number of 12-hour time steps
    
    Returns:
        time_array: Array of time in seconds since reference date
    """
    # Reference date: 1990-01-01 00:00:00 +10 (AEST)
    # Note: AEST is UTC+10, but we'll work in UTC and adjust
    reference_date = datetime(1990, 1, 1, 0, 0, 0)
    
    # Start date for the year (January 1, 00:00:00)
    start_date = datetime(start_year, 1, 1, 0, 0, 0)
    
    # Calculate seconds from reference to start of year
    seconds_from_ref = (start_date - reference_date).total_seconds()
    
    # Create time array: 0, 43200, 86400, ... (12-hour intervals)
    # dt = 43200 seconds (12 hours)
    dt = 43200  # 12 hours in seconds
    time_array = seconds_from_ref + np.arange(ntime_12h) * dt
    
    return time_array

def resample_transport_file(input_file, output_file, year):
    """
    Resample transport file from hourly to 12-hourly by averaging
    Preserves all variables and dimensions with original names
    
    Args:
        input_file: Path to input transport NetCDF file
        output_file: Path to output transport NetCDF file
        year: Year for time reference
    """
    print(f"Processing transport file: {input_file}")
    
    # Open input file
    with nc.Dataset(input_file, 'r') as ds_in:
        # Get all dimensions and their sizes (preserve names)
        dims_info = {}
        for dim_name, dim_obj in ds_in.dimensions.items():
            if dim_name == 'time':
                dims_info[dim_name] = None  # Unlimited
            else:
                dims_info[dim_name] = len(dim_obj)
        
        # Get time dimension size
        ntime_hourly = len(ds_in.dimensions['time'])
        
        # Calculate number of 12-hour time steps
        ntime_12h = ntime_hourly // 12
        
        if ntime_12h == 0:
            print(f"Warning: Not enough time steps ({ntime_hourly}) for 12-hour averaging")
            return
        
        print(f"  Original time steps: {ntime_hourly}")
        print(f"  New 12-hour time steps: {ntime_12h}")
        
        # Create time array
        time_12h = create_12hour_time_array(year, ntime_12h)
        
        # Create output file
        with nc.Dataset(output_file, 'w', format='NETCDF4') as ds_out:
            # Create all dimensions with original names
            for dim_name, dim_size in dims_info.items():
                ds_out.createDimension(dim_name, dim_size)
            
            # Process all variables
            for var_name, var_in in ds_in.variables.items():
                if var_name == 'time':
                    # Create new time variable
                    var_out = ds_out.createVariable('time', var_in.dtype, var_in.dimensions)
                    var_out[:] = time_12h
                    # Copy all attributes from original
                    for attr in var_in.ncattrs():
                        var_out.setncattr(attr, var_in.getncattr(attr))
                    # Override with correct values
                    var_out.units = 'seconds since 1990-01-01 00:00:00 +10'
                    var_out.dt = 43200.0
                
                elif 'time' in var_in.dimensions:
                    # Resample variables that depend on time
                    var_data = var_in[:]
                    # Get shape and dimensions
                    var_shape = var_data.shape
                    time_idx = var_in.dimensions.index('time')
                    
                    # Reshape to average over time dimension
                    # Calculate new shape
                    new_shape = list(var_shape)
                    new_shape[time_idx] = ntime_12h
                    
                    # Reshape for averaging: [ntime_12h, 12, ...] then average over axis 1
                    reshape_list = list(var_shape)
                    reshape_list[time_idx] = ntime_12h
                    reshape_list.insert(time_idx + 1, 12)
                    
                    # Create indexing to reshape
                    var_reshaped = var_data[:ntime_12h*12].reshape(reshape_list)
                    var_12h = np.nanmean(var_reshaped, axis=time_idx + 1)
                    
                    # Create variable with original dimensions
                    var_out = ds_out.createVariable(var_name, var_in.dtype, var_in.dimensions)
                    var_out[:] = var_12h
                    
                    # Copy all attributes
                    for attr in var_in.ncattrs():
                        if attr != '_FillValue':
                            try:
                                var_out.setncattr(attr, var_in.getncattr(attr))
                            except:
                                pass
                    # Copy fill value if it exists
                    if hasattr(var_in, '_FillValue'):
                        var_out._FillValue = var_in._FillValue
                
                else:
                    # Copy variables that don't depend on time as-is
                    var_out = ds_out.createVariable(var_name, var_in.dtype, var_in.dimensions)
                    var_out[:] = var_in[:]
                    # Copy all attributes
                    for attr in var_in.ncattrs():
                        if attr != '_FillValue':
                            try:
                                var_out.setncattr(attr, var_in.getncattr(attr))
                            except:
                                pass
                    # Copy fill value if it exists
                    if hasattr(var_in, '_FillValue'):
                        var_out._FillValue = var_in._FillValue
            
            # Copy all global attributes
            for attr in ds_in.ncattrs():
                try:
                    ds_out.setncattr(attr, ds_in.getncattr(attr))
                except:
                    pass
            
            # Add processing history
            history = f"Resampled from hourly to 12-hourly resolution by averaging every 12 time steps. Original file: {os.path.basename(input_file)}"
            if 'history' in ds_in.ncattrs():
                history = ds_in.getncattr('history') + '\n' + history
            ds_out.setncattr('history', history)
            ds_out.setncattr('time_resolution', '12-hourly')
            ds_out.setncattr('dt', '43200 seconds')
    
    print(f"  Saved resampled transport file: {output_file}")

def resample_variables_file(input_file, output_file, year):
    """
    Resample variables file from hourly to 12-hourly by averaging
    Preserves all variables and dimensions with original names
    
    Args:
        input_file: Path to input variables NetCDF file
        output_file: Path to output variables NetCDF file
        year: Year for time reference
    """
    print(f"Processing variables file: {input_file}")
    
    # Open input file
    with nc.Dataset(input_file, 'r') as ds_in:
        # Get all dimensions and their sizes (preserve names)
        dims_info = {}
        for dim_name, dim_obj in ds_in.dimensions.items():
            if dim_name == 'time':
                dims_info[dim_name] = None  # Unlimited
            else:
                dims_info[dim_name] = len(dim_obj)
        
        # Get time dimension size
        ntime_hourly = len(ds_in.dimensions['time'])
        
        # Calculate number of 12-hour time steps
        ntime_12h = ntime_hourly // 12
        
        if ntime_12h == 0:
            print(f"Warning: Not enough time steps ({ntime_hourly}) for 12-hour averaging")
            return
        
        print(f"  Original time steps: {ntime_hourly}")
        print(f"  New 12-hour time steps: {ntime_12h}")
        
        # Create time array
        time_12h = create_12hour_time_array(year, ntime_12h)
        
        # Create output file
        with nc.Dataset(output_file, 'w', format='NETCDF4') as ds_out:
            # Create all dimensions with original names
            for dim_name, dim_size in dims_info.items():
                ds_out.createDimension(dim_name, dim_size)
            
            # Process all variables
            for var_name, var_in in ds_in.variables.items():
                if var_name == 'time':
                    # Create new time variable
                    var_out = ds_out.createVariable('time', var_in.dtype, var_in.dimensions)
                    var_out[:] = time_12h
                    # Copy all attributes from original
                    for attr in var_in.ncattrs():
                        var_out.setncattr(attr, var_in.getncattr(attr))
                    # Override with correct values
                    var_out.units = 'seconds since 1990-01-01 00:00:00 +10'
                    var_out.dt = 43200.0
                
                elif 'time' in var_in.dimensions:
                    # Resample variables that depend on time
                    var_data = var_in[:]
                    # Get shape and dimensions
                    var_shape = var_data.shape
                    time_idx = var_in.dimensions.index('time')
                    
                    # Reshape to average over time dimension
                    # Calculate new shape
                    new_shape = list(var_shape)
                    new_shape[time_idx] = ntime_12h
                    
                    # Reshape for averaging: [ntime_12h, 12, ...] then average over axis 1
                    reshape_list = list(var_shape)
                    reshape_list[time_idx] = ntime_12h
                    reshape_list.insert(time_idx + 1, 12)
                    
                    # Create indexing to reshape
                    var_reshaped = var_data[:ntime_12h*12].reshape(reshape_list)
                    var_12h = np.nanmean(var_reshaped, axis=time_idx + 1)
                    
                    # Create variable with original dimensions
                    var_out = ds_out.createVariable(var_name, var_in.dtype, var_in.dimensions)
                    var_out[:] = var_12h
                    
                    # Copy all attributes
                    for attr in var_in.ncattrs():
                        if attr != '_FillValue':
                            try:
                                var_out.setncattr(attr, var_in.getncattr(attr))
                            except:
                                pass
                    # Copy fill value if it exists
                    if hasattr(var_in, '_FillValue'):
                        var_out._FillValue = var_in._FillValue
                
                else:
                    # Copy variables that don't depend on time as-is
                    var_out = ds_out.createVariable(var_name, var_in.dtype, var_in.dimensions)
                    var_out[:] = var_in[:]
                    # Copy all attributes
                    for attr in var_in.ncattrs():
                        if attr != '_FillValue':
                            try:
                                var_out.setncattr(attr, var_in.getncattr(attr))
                            except:
                                pass
                    # Copy fill value if it exists
                    if hasattr(var_in, '_FillValue'):
                        var_out._FillValue = var_in._FillValue
            
            # Copy all global attributes
            for attr in ds_in.ncattrs():
                try:
                    ds_out.setncattr(attr, ds_in.getncattr(attr))
                except:
                    pass
            
            # Add processing history
            history = f"Resampled from hourly to 12-hourly resolution by averaging every 12 time steps. Original file: {os.path.basename(input_file)}"
            if 'history' in ds_in.ncattrs():
                history = ds_in.getncattr('history') + '\n' + history
            ds_out.setncattr('history', history)
            ds_out.setncattr('time_resolution', '12-hourly')
            ds_out.setncattr('dt', '43200 seconds')
    
    print(f"  Saved resampled variables file: {output_file}")

def main():
    """Main function to resample all files"""
    
    # Input and output directories
    input_dir = "/home/por07g/Documents/Projects/NESP_ParksAustralia/Hydro/outputs_final_corrected"
    output_dir = "/home/por07g/Documents/Projects/NESP_ParksAustralia/Hydro/outputs_final_corrected_12h"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Process years 2017-2023
    years = range(2017, 2024)
    
    print("="*60)
    print("Resampling OCEANUS outputs from hourly to 12-hourly")
    print("="*60)
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print()
    
    for year in years:
        print(f"\nProcessing year {year}...")
        
        # Process transport file
        transport_input = os.path.join(input_dir, f"mb_transport_{year}_bass2_simple.nc")
        transport_output = os.path.join(output_dir, f"mb_transport_{year}_bass2_simple.nc")
        
        if os.path.exists(transport_input):
            resample_transport_file(transport_input, transport_output, year)
        else:
            print(f"  Transport file not found: {transport_input}")
        
        # Process variables file
        variables_input = os.path.join(input_dir, f"Variables_{year}.nc")
        variables_output = os.path.join(output_dir, f"Variables_{year}.nc")
        
        if os.path.exists(variables_input):
            resample_variables_file(variables_input, variables_output, year)
        else:
            print(f"  Variables file not found: {variables_input}")
    
    print("\n" + "="*60)
    print("Resampling completed!")
    print(f"Output files saved in: {output_dir}")
    print("="*60)

if __name__ == "__main__":
    main()

