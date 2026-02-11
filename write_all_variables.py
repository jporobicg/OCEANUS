# write_all_variables.py
# Function to write all variables to a single NetCDF file
# Based on write_av_var_new.m

import numpy as np
import netCDF4 as nc
from datetime import datetime

def write_all_variables(tims, bid, variables_info, fnm):
    """
    Write all variables to a single NetCDF file.
    
    Args:
        tims: time array
        bid: box IDs array
        variables_info: list of dictionaries, each containing:
            - 'name': NetCDF variable name (str)
            - 'data': data array (nbox, ntime, nlevel)
            - 'long_name': long descriptive name (str)
            - 'units': units string (str)
            - 'FillValue_': fill value (float)
        fnm: output NetCDF file name
    """
    if len(variables_info) == 0:
        raise ValueError("No variables provided in variables_info")
    
    nbx = len(bid)
    ntm = len(tims)
    nlv = variables_info[0]['data'].shape[2]  # Get nlevel from first variable
    
    # Create NetCDF file
    nc_out = nc.Dataset(fnm, 'w', format='NETCDF4')
    
    # Create dimensions
    time_dim = nc_out.createDimension('time', None)
    level_dim = nc_out.createDimension('level', nlv)
    boxes_dim = nc_out.createDimension('boxes', nbx)
    
    # Create coordinate variables
    time_var = nc_out.createVariable('time', 'f8', ('time',))
    boxes_var = nc_out.createVariable('boxes', 'i4', ('boxes',))
    level_var = nc_out.createVariable('level', 'i4', ('level',))
    
    # Set coordinate variable attributes
    time_var.setncattr('long_name', 'simulation_time')
    time_var.setncattr('units', 'seconds since 1990-01-01 00:00:00 +10')
    time_var.setncattr('grads_step', 'seconds')
    time_var.setncattr('dt', 86400)  # 86400 seconds = 1 day
    
    boxes_var.setncattr('long_name', 'Box IDs')
    
    level_var.setncattr('long_name', 'layer index; 1=near-surface')
    
    # Create data variables dynamically from variables_info
    nc_variables = {}
    for var_info in variables_info:
        var_name = var_info['name']
        nc_var = nc_out.createVariable(var_name, 'f4', ('time', 'boxes', 'level'))
        
        # Set attributes from metadata
        nc_var.setncattr('long_name', var_info['long_name'])
        nc_var.setncattr('units', var_info['units'])
        nc_var.setncattr('FillValue_', var_info['FillValue_'])
        
        # Store reference for later writing
        nc_variables[var_name] = nc_var
    
    # Set global attributes
    str_desc = ('Standar Atlantis tempora file - created based on the roms model and the bgm '
                'configuration file. Properties averaged from hydrodynamic model outputs.')
    his = f'Created Using Javier Porobic (more details : https://github.com/jporobicg) codes on {datetime.now().strftime("%Y-%m-%d")}'
    
    nc_out.setncattr('Title', 'Roms files Salish sea ocean output')
    nc_out.setncattr('Comments', str_desc)
    nc_out.setncattr('history', his)
    
    # Leave define mode
    nc_out._enddef()
    
    # Write coordinate variables
    time_var[:] = tims * 86400  # Convert days to seconds
    level_var[:] = np.arange(1, nlv + 1)
    boxes_var[:] = bid
    
    # Write data variables
    for var_info in variables_info:
        var_name = var_info['name']
        var_data = var_info['data']
        
        # Permute data to match NetCDF dimension order (time, boxes, level)
        var_permuted = np.transpose(var_data, (1, 0, 2))
        
        # Filter extreme values (matching MATLAB code)
        var_permuted[var_permuted > 1e15] = 0
        
        # Write to NetCDF
        nc_variables[var_name][:] = var_permuted
    
    # Close file
    nc_out.close()
    print(f"All variables written to: {fnm}")
