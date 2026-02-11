# box_averages.py
# Function to calculate average properties for modeling polygons
# Based on box_av_SS.m

import numpy as np
import netCDF4 as nc
import os
import pickle
from scipy.interpolate import griddata
from matplotlib.path import Path

def latcor(lat):
    """Calculate latitude correction factor."""
    return np.cos(np.radians(lat))

def polyarea(x, y):
    """Calculate polygon area using shoelace formula."""
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

def box_averages(vert, varn, dlev, fnm, fll, month_str, year, output_folder=None):
    """
    Calculate average properties for modeling polygons.
    
    Args:
        vert: list of box vertices
        varn: variable name to process
        dlev: depth levels array
        fnm: NetCDF data file path
        fll: NetCDF mesh file path
        month_str: two-digit month string (e.g., "01", "12")
        year: year number
    """
    # Open NetCDF files
    nc_data = nc.Dataset(fnm, 'r')
    nc_mesh = nc.Dataset(fll, 'r')
    
    # Calculate depth intervals
    dint = np.diff(dlev)
    nlay = len(dint)
    
    # Get basic dimensions
    tims = nc_data.variables['time'][:]
    lon = nc_mesh.variables['longitude'][:].transpose(1, 0)
    lat = nc_mesh.variables['latitude'][:].transpose(1, 0)
    
    # Get depth coordinate - try to get from data file, fallback to mesh
    if 'zc' in nc_data.variables:
        zc = nc_data.variables['zc'][:]  # Depth coordinate from data file
    elif 'zc' in nc_mesh.variables:
        zc = nc_mesh.variables['zc'][:]
    else:
        raise ValueError("Depth coordinate 'zc' not found in data or mesh file")
    
    ntm = len(tims)
    nbox = len(vert)
    
    # Check variable dimensions to handle different depth structures
    if varn not in nc_data.variables:
        raise ValueError(f"Variable {varn} not found in {fnm}")
    
    var_dims = nc_data.variables[varn].dimensions
    var_shape = nc_data.variables[varn].shape
    
    # Determine if variable has depth dimension and get its size
    has_depth = False
    var_depth_size = 1
    depth_dim_idx = None
    
    # Common depth dimension names
    depth_dim_names = ['k', 'depth', 'z', 'level', 'zc']
    
    for i, dim_name in enumerate(var_dims):
        if dim_name in depth_dim_names or 'depth' in dim_name.lower() or 'level' in dim_name.lower():
            has_depth = True
            depth_dim_idx = i
            var_depth_size = var_shape[i]
            break
    
    # If variable has depth, try to get its depth coordinate
    var_zc = zc  # Default to using zc from file
    if has_depth and depth_dim_idx is not None:
        # Try to find depth coordinate variable that matches variable depth dimension
        depth_coord_name = None
        for coord_name in ['zc', 'depth', 'z', 'level']:
            if coord_name in nc_data.variables:
                coord_data = nc_data.variables[coord_name][:]
                # Check if this coordinate matches the variable's depth dimension
                if len(coord_data) == var_depth_size:
                    var_zc = coord_data
                    depth_coord_name = coord_name
                    break
        
        # If no matching depth coordinate found, use subset of zc or create one
        if depth_coord_name is None:
            if len(zc) >= var_depth_size:
                # Use first var_depth_size levels from zc
                var_zc = zc[:var_depth_size]
                print(f"Info: Using first {var_depth_size} levels from zc for variable {varn}")
            else:
                # Variable has more depth levels than zc - use zc and pad if needed
                var_zc = zc
                print(f"Warning: Variable {varn} has {var_depth_size} depth levels but zc has only {len(zc)}")
    
    # Ensure output folder
    if output_folder is None:
        output_folder = "."
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)

    # month_str is provided by caller (two-digit, e.g., 01, 12)

    # First step file (grid preparation) shared across all months/variables/years in this folder
    file1 = os.path.join(output_folder, "SS_First_Step.pkl")
    
    if not os.path.exists(file1):
        print(f"Creating first step file: {file1}")
        
        # Initialize arrays
        barea = np.zeros(nbox)
        ngrd = 0
        idx = {}
        flo = []
        fla = []
        ix = {}
        iy = {}
        nav = np.zeros(nbox, dtype=int)
        
        # Prepare integration grid for all boxes
        for ibx in range(nbox):
            x = vert[ibx][:, 0]
            y = vert[ibx][:, 1]
            
            # Create polygon path for point-in-polygon test
            poly_path = Path(np.column_stack([x, y]))
            
            # Find grid points inside polygon
            lon_flat = lon.ravel()
            lat_flat = lat.ravel()
            points = np.column_stack([lon_flat, lat_flat])
            ig = poly_path.contains_points(points)
            ig = np.where(ig)[0]  # Get indices of points inside polygon
            
            nav[ibx] = len(ig)
            ii = list(range(ngrd, ngrd + nav[ibx]))
            ngrd = ngrd + nav[ibx]
            idx[ibx] = ii
            
            # Get coordinates of points inside polygon
            flo.extend(lon_flat[ig])
            fla.extend(lat_flat[ig])
            
            # Get grid indices
            lon_indices, lat_indices = np.unravel_index(ig, lon.shape)
            ix[ibx] = lon_indices
            iy[ibx] = lat_indices
            
            # Calculate polygon area (important for production)
            if varn == 'w':
                # Convert lat, lon to x, y for area calculation
                xx = (x - np.mean(x)) * latcor(np.mean(y)) * 111119
                yy = y * 111119
                barea[ibx] = polyarea(xx, yy) / 1000000
            else:
                barea[ibx] = 1.0
        
        # Save first step data
        first_step_data = {
            'zc': zc,
            'idx': idx,
            'flo': np.array(flo),
            'fla': np.array(fla),
            'nav': nav,
            'iy': iy,
            'ix': ix,
            'barea': barea
        }
        with open(file1, 'wb') as f_pkl:
            pickle.dump(first_step_data, f_pkl)
            
    else:
        print(f"Loading existing first step file: {file1}")
        with open(file1, 'rb') as f_pkl:
            first_step_data = pickle.load(f_pkl)
            
        zc = first_step_data['zc']
        idx = first_step_data['idx']
        flo = first_step_data['flo']
        fla = first_step_data['fla']
        nav = first_step_data['nav']
        iy = first_step_data['iy']
        ix = first_step_data['ix']
        barea = first_step_data['barea']
    
    # Initialize variable average array
    Var_avg = np.full((nbox, ntm, nlay), np.nan)
    
    # Second step file (variable processing) per variable
    file2 = os.path.join(output_folder, f"{month_str}_{year}_{varn}_SS_Second_step.npz")
    
    if not os.path.exists(file2):
        print(f"Creating second step file: {file2}")
        
        # Process each time step
        for id in range(ntm):
              
            # Read variable data for this time step
            # Handle different variable dimensions (2D or 3D)
            var_full = nc_data.variables[varn]
            var_full_shape = var_full.shape
            
            if has_depth and depth_dim_idx is not None:
                # 3D variable: need to extract and rearrange to (depth, lat, lon)
                # Get all dimensions for this time step
                if len(var_full_shape) == 4:  # (time, depth, lat, lon) or (time, lat, lon, depth)
                    if depth_dim_idx == 1:  # depth is second dimension: (time, depth, lat, lon)
                        varData = var_full[id, :, :, :]  # (depth, lat, lon)
                    elif depth_dim_idx == 3:  # depth is last: (time, lat, lon, depth)
                        varData = var_full[id, :, :, :].transpose(2, 0, 1)  # Rearrange to (depth, lat, lon)
                    else:
                        # Try to read and figure out dimensions
                        varData = var_full[id, :, :, :]  # Default assumption
                        # Check if we need to transpose
                        if varData.shape[0] != var_depth_size:
                            varData = varData.transpose(2, 0, 1)  # Try transposing
                elif len(var_full_shape) == 3:  # (time, lat, lon) - no depth, but has_depth was True
                    # This shouldn't happen, but handle it
                    varData_2d = var_full[id, :, :]  # (lat, lon)
                    varData = varData_2d[np.newaxis, :, :]  # (1, lat, lon)
                    var_zc = np.array([0])  # Surface level
                else:
                    raise ValueError(f"Unexpected variable shape for {varn}: {var_full_shape}")
            else:
                # 2D variable: (time, lat, lon) - expand to have depth dimension
                if len(var_full_shape) == 3:
                    varData_2d = var_full[id, :, :]  # (lat, lon)
                    varData = varData_2d[np.newaxis, :, :]  # (1, lat, lon)
                    var_zc = np.array([0])  # Surface level
                else:
                    raise ValueError(f"Unexpected variable shape for {varn}: {var_full_shape}")
            
            # Verify varData depth dimension matches var_zc
            if varData.shape[0] != len(var_zc):
                print(f"Warning: varData depth dimension ({varData.shape[0]}) doesn't match var_zc length ({len(var_zc)})")
                # Adjust var_zc to match varData
                if len(var_zc) > varData.shape[0]:
                    var_zc = var_zc[:varData.shape[0]]
                elif len(zc) >= varData.shape[0]:
                    var_zc = zc[:varData.shape[0]]
                else:
                    # Create a simple depth array
                    var_zc = np.linspace(zc[0], zc[-1], varData.shape[0])
            
            # Handle missing values
            varData[varData >= 1.0e16] = np.nan
            
            # Find valid grid points
            lon_flat = lon.ravel()
            lat_flat = lat.ravel()
            imodxy = ~(np.isnan(lon_flat) | (lon_flat == 0))
            xlon = lon_flat[imodxy]
            ylat = lat_flat[imodxy]
            
            # Process by layers
            interpolatedVarData = np.full((nlay, len(flo)), np.nan)
            
            for layer in range(nlay):
                # Average variable data for this depth layer
                # Use var_zc (variable's depth coordinate) instead of zc
                layer_mask = (var_zc <= -dlev[layer]) & (var_zc >= -dlev[layer+1])
                if np.any(layer_mask):
                    # Ensure layer_mask matches varData depth dimension
                    if len(layer_mask) == varData.shape[0]:
                        layer_var_data = np.nanmean(varData[layer_mask, :, :], axis=0)
                    else:
                        # If dimensions don't match, use all available depth levels
                        print(f"Warning: Depth dimension mismatch for {varn}, using all available levels")
                        layer_var_data = np.nanmean(varData, axis=0)
                    
                    # Remove NaNs and interpolate
                    layer_var_data_flat = layer_var_data.ravel()
                    valid_data = layer_var_data_flat[imodxy]
                    valid_mask = ~np.isnan(valid_data)
                    
                    if np.sum(valid_mask) > 0:
                        valid_xlon = xlon[valid_mask]
                        valid_ylat = ylat[valid_mask]
                        valid_values = valid_data[valid_mask]
                        
                        # Interpolate to face integration points
                        interpolated = griddata((valid_xlon, valid_ylat), valid_values, 
                                               (flo, fla), method='linear', fill_value=np.nan)
                        interpolatedVarData[layer, :] = interpolated
            
            # Calculate averages for each box
            for jj in range(nbox):
                # 3D variables
                tmp = np.zeros((nlay, nav[jj]))
                
                for lvs in range(nlay):
                    for ij in range(nav[jj]):
                        # Average by box by grid cell by layer
                        tmp[lvs, ij] = barea[jj] * np.nanmean(interpolatedVarData[lvs, idx[jj][ij]])
                    
                    Var_avg[jj, id, lvs] = np.nanmean(tmp[lvs, :])
        
        # Save second step data
        np.savez(file2, Var_avg=Var_avg, tims=tims)            
    
    # Close NetCDF files
    nc_data.close()
    nc_mesh.close()
 
 