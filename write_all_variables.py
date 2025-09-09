# write_all_variables.py
# Function to write all variables to a single NetCDF file
# Based on write_av_var_new.m

import numpy as np
import netCDF4 as nc
from datetime import datetime

def write_all_variables(tims, bid, temp_data, salt_data, w_data, fnm):
    """
    Write all variables (temperature, salinity, verticalflux) to a single NetCDF file.
    
    Args:
        tims: time array
        bid: box IDs array
        temp_data: temperature data array (nbox, ntime, nlevel)
        salt_data: salinity data array (nbox, ntime, nlevel)
        w_data: vertical flux data array (nbox, ntime, nlevel)
        fnm: output NetCDF file name
    """
    nbx = len(bid)
    ntm = len(tims)
    nlv = temp_data.shape[2]  # Assuming all variables have same shape
    
    # Create NetCDF file
    nc_out = nc.Dataset(fnm, 'w', format='NETCDF4')
    
    # Create dimensions
    time_dim = nc_out.createDimension('time', None)
    level_dim = nc_out.createDimension('level', nlv)
    boxes_dim = nc_out.createDimension('boxes', nbx)
    
    # Create variables
    time_var = nc_out.createVariable('time', 'f8', ('time',))
    boxes_var = nc_out.createVariable('boxes', 'i4', ('boxes',))
    level_var = nc_out.createVariable('level', 'i4', ('level',))
    
    # Create variable data arrays
    var_temp = nc_out.createVariable('temperature', 'f4', ('time', 'boxes', 'level'))
    var_salinity = nc_out.createVariable('salinity', 'f4', ('time', 'boxes', 'level'))
    var_w = nc_out.createVariable('verticalflux', 'f4', ('time', 'boxes', 'level'))
    
    # Set variable attributes
    time_var.setncattr('long_name', 'simulation_time')
    time_var.setncattr('units', 'seconds since 1990-01-01 00:00:00 +10')
    time_var.setncattr('grads_step', 'seconds')
    time_var.setncattr('dt', 86400)  # 86400 seconds = 1 day
    
    boxes_var.setncattr('long_name', 'Box IDs')
    
    level_var.setncattr('long_name', 'layer index; 1=near-surface')
    
    # Temperature attributes
    var_temp.setncattr('long_name', 'sea_water_conservative_temperature')
    var_temp.setncattr('units', 'degrees celsius [c]')
    var_temp.setncattr('FillValue_', -1e20)
    
    # Salinity attributes
    var_salinity.setncattr('long_name', 'sea_water_reference_salinity')
    var_salinity.setncattr('units', 'g kg-1')
    var_salinity.setncattr('FillValue_', -1e20)
    
    # Vertical flux attributes
    var_w.setncattr('long_name', 'ocean vertical velocity')
    var_w.setncattr('units', 'meter second-1')
    var_w.setncattr('field', 'v-velocity, scalar, series')
    
    # Set global attributes
    str_desc = ('Standar Atlantis tempora file - created based on the roms model and the bgm '
                'configuration file. Properties averaged from hydrodynamic model outputs.')
    his = f'Created Using Javier Porobic (more details : https://github.com/jporobicg) codes on {datetime.now().strftime("%Y-%m-%d")}'
    
    nc_out.setncattr('Title', 'Roms files Salish sea ocean output')
    nc_out.setncattr('Comments', str_desc)
    nc_out.setncattr('history', his)
    
    # Leave define mode
    nc_out._enddef()
    
    # Permute data to match NetCDF dimension order (time, boxes, level)
    temp_permuted = np.transpose(temp_data, (1, 0, 2))
    salt_permuted = np.transpose(salt_data, (1, 0, 2))
    w_permuted = np.transpose(w_data, (1, 0, 2))
    
    # Filter extreme values (matching MATLAB code)
    temp_permuted[temp_permuted > 1e15] = 0
    salt_permuted[salt_permuted > 1e15] = 0
    w_permuted[w_permuted > 1e15] = 0
    
    # Write variables
    time_var[:] = tims * 86400  # Convert days to seconds
    level_var[:] = np.arange(1, nlv + 1)
    boxes_var[:] = bid
    var_temp[:] = temp_permuted
    var_salinity[:] = salt_permuted
    var_w[:] = w_permuted
    
    # Close file
    nc_out.close()
    print(f"All variables written to: {fnm}")
