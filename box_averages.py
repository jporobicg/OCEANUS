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

def box_averages(vert, varn, dlev, fnm, fll, month, year, output_folder=None):
    """
    Calculate average properties for modeling polygons.
    
    Args:
        vert: list of box vertices
        varn: variable name to process
        dlev: depth levels array
        fnm: NetCDF data file path
        fll: NetCDF mesh file path
        month: month number
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
    zc = nc_data.variables['zc'][:]  # Depth coordinate
    
    ntm = len(tims)
    nbox = len(vert)
    
    # Ensure output folder
    if output_folder is None:
        output_folder = "."
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)

    # First step file (grid preparation)
    file1 = os.path.join(output_folder, f"{month}_{year}_{varn}_SS_First_Step.pkl")
    
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
    
    # Second step file (variable processing)
    file2 = os.path.join(output_folder, f"{month}_{year}_{varn}_SS_Second_step.npz")
    
    if not os.path.exists(file2):
        print(f"Creating second step file: {file2}")
        
        # Process each time step
        for id in range(ntm):
            print(f"Processing time step {id+1}/{ntm}")
            
            # Read variable data for this time step
            # Variables have shape (time, depth, lat, lon)
            varData = nc_data.variables[varn][id, :, :, :]  # (depth, lat, lon)
            
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
                layer_mask = (zc <= -dlev[layer]) & (zc >= -dlev[layer+1])
                if np.any(layer_mask):
                    layer_var_data = np.nanmean(varData[layer_mask, :, :], axis=0)
                    
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
 
 