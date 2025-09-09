import numpy as np
import netCDF4 as nc
from datetime import datetime

def write_nc_transport(pt1, pt2, lr, tims, T, fcid, fnm):
    """
    Write transport data to NetCDF file.
    Args:
        pt1: face start points (nface, 2)
        pt2: face end points (nface, 2)
        lr: left/right box IDs for each face (nface, 2)
        tims: time array
        T: transport array (nface, ntime, nlayer)
        fcid: face indices (nface,)
        fnm: output NetCDF file name
    """
    nfc = pt1.shape[0]
    ntm = len(tims)
    nlv = T.shape[2]

    nc_out = nc.Dataset(fnm, 'w', format='NETCDF4')
    nc_out.createDimension('time', None)
    nc_out.createDimension('level', nlv)
    nc_out.createDimension('faces', nfc)

    time_var = nc_out.createVariable('time', 'f8', ('time',))
    level_var = nc_out.createVariable('level', 'i4', ('level',))
    faces_var = nc_out.createVariable('faces', 'i4', ('faces',))
    pt1x_var = nc_out.createVariable('pt1_x', 'f4', ('faces',))
    pt1y_var = nc_out.createVariable('pt1_y', 'f4', ('faces',))
    pt2x_var = nc_out.createVariable('pt2_x', 'f4', ('faces',))
    pt2y_var = nc_out.createVariable('pt2_y', 'f4', ('faces',))
    dest_var = nc_out.createVariable('dest_boxid', 'i4', ('faces',))
    sour_var = nc_out.createVariable('source_boxid', 'i4', ('faces',))
    transp_var = nc_out.createVariable('transport', 'f4', ('time', 'faces', 'level'))

    # Add variable attributes (matching MATLAB code)
    # Time variable attributes
    time_var.setncattr('long_name', 'simulation_time')
    time_var.setncattr('units', 'seconds since 1990-01-01 00:00:00 +10')
    time_var.setncattr('grads_step', 'seconds')
    time_var.setncattr('dt', 86400)  # 86400 seconds = 1 day
    
    # Level variable attributes
    level_var.setncattr('long_name', 'layer index; 1=near-surface')
    
    # Faces variable attributes
    faces_var.setncattr('long_name', 'face IDs')
    
    # Point 1 x-coordinate attributes
    pt1x_var.setncattr('long_name', 'x coordinate of point 1 of face')
    pt1x_var.setncattr('units', 'degree_east')
    
    # Point 1 y-coordinate attributes
    pt1y_var.setncattr('long_name', 'y coordinate of point 1 of face')
    pt1y_var.setncattr('units', 'degree_north')
    
    # Point 2 x-coordinate attributes
    pt2x_var.setncattr('long_name', 'x coordinate of point 2 of face')
    pt2x_var.setncattr('units', 'degree_east')
    
    # Point 2 y-coordinate attributes
    pt2y_var.setncattr('long_name', 'y coordinate of point 2 of face')
    pt2y_var.setncattr('units', 'degree_north')
    
    # Destination box attributes
    dest_var.setncattr('long_name', 'ID of destination box')
    dest_var.setncattr('units', 'IDs')
    
    # Source box attributes
    sour_var.setncattr('long_name', 'ID of source box')
    sour_var.setncattr('units', 'IDs')
    
    # Transport variable attributes
    transp_var.setncattr('long_name', 'flux across face')
    transp_var.setncattr('units', 'meter second-1')
    transp_var.setncattr('comment', 'm/s is to left, viewing from pt1 to pt2')
    transp_var.setncattr('FillValue_', -1e20)

    # Add global attributes (matching MATLAB code)
    str_desc = ('Transports across faces of box model derived from hydrodynamic'
                ' model current field. Positive numbers are water gain in dest'
                ' box, loss in source box.')
    his = f'Created Using Javier Porobic (more details : https://github.com/jporobicg) codes on {datetime.now().strftime("%Y-%m-%d")}'
    
    nc_out.setncattr('Title', 'Transport from NEMO salish sea model')
    nc_out.setncattr('Comments', str_desc)
    nc_out.setncattr('Institution', 'CSIRO')
    nc_out.setncattr('Conventions', 'Standar')
    nc_out.setncattr('history', his)

    # Leave define mode (matching MATLAB code)
    nc_out._enddef()

    time_var[:] = tims * 86400  # Convert days to seconds
    level_var[:] = np.arange(1, nlv+1)
    faces_var[:] = fcid
    pt1x_var[:] = pt1[:, 0]
    pt1y_var[:] = pt1[:, 1]
    pt2x_var[:] = pt2[:, 0]
    pt2y_var[:] = pt2[:, 1]
    dest_var[:] = lr[:, 0]
    sour_var[:] = lr[:, 1]
    T_permuted = np.transpose(T, (1, 0, 2))
    
    # Data validation: filter out extreme values (matching MATLAB code)
    T_permuted[T_permuted > 1e15] = 0
    
    transp_var[:] = T_permuted
    nc_out.close() 
    print(f"Transport file written: {fnm}") 