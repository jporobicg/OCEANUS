import numpy as np
import netCDF4 as nc
import os
from scipy.interpolate import griddata
import pickle


def cart2pol(x, y):
    """Convert Cartesian coordinates to polar coordinates."""
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return theta, r


def calc_transport(vert, pt1, pt2, dlev, dinc, rimn, fnm, fll, month, year):
    """
    Calculate transport through faces from hydrodynamic data.
    Complete Python version of transport_SS.m
    Generates intermediate files but does not return values.

    Args:
        vert: box vertices
        pt1: face start points (nface, 2)
        pt2: face end points (nface, 2)
        dlev: depth levels (array)
        dinc: integration step length
        rimn: minimum integration steps per face
        fnm: hydrodynamic data file (NetCDF)
        fll: mesh file (NetCDF)
        f: file index (for output naming)
    """
    # Open NetCDF files
    nc_data = nc.Dataset(fnm, 'r')
    nc_mesh = nc.Dataset(fll, 'r')

    # Get basic dimensions
    dint = np.diff(dlev)
    nlay = len(dint)
    tims = nc_data.variables['time'][:]
    ntm = len(tims)
    rimx = 400
    # Get grid coordinates
    lon = nc_mesh.variables['longitude'][:].transpose(1, 0)
    lat = nc_mesh.variables['latitude'][:].transpose(1, 0)
    zc = nc_data.variables['zc'][:]
    # print(zc)

    nfc = pt1.shape[0]
    ndps = len(zc)

    # Initialize transport array
    T = np.full((nfc, ntm, nlay), np.nan)

    # Check if first step file exists
    file1 = 'SS_first_Step.pkl'
    if not os.path.exists(file1):
        # Prepare the integration grid for each face
        ngrd = 0
        idx = {}
        flo = []
        fla = []
        dirn = []
        vint = {}
        abovebot = np.ones((ndps, 0))

        for ifc in range(nfc):  # Number of faces
            x = [pt1[ifc, 0], pt2[ifc, 0]]
            y = [pt1[ifc, 1], pt2[ifc, 1]]
            yd = y[1] - y[0]  # Distance between points
            xd = x[1] - x[0]

            if abs(xd) > 0.00001:
                lcor = np.cos(np.radians(np.mean(y)))
                xd = lcor * xd
            else:
                lcor = 1

            rdist = 111191 * np.sqrt(yd**2 + xd**2)  # distance in kilometres
            dirn.append(np.arctan2(yd, xd))  # cartesian to polar coordinates

            rinc = int(np.ceil(rdist / dinc))  # To cope with Large domain
            if rinc < rimn:
                rinc = rimn
            if rinc > rimx:
                rinc = rimx
                yi = yd / rinc
                Y = np.linspace(y[0] + yi/2, y[1] - yi/2, rinc)
                xi = xd / rinc
                X = np.linspace(x[0] + xi/2, x[1] - xi/2, rinc)
                ninc = len(Y)
            if ninc == 0:
                ninc = len(X)
                Y = np.full_like(X, np.mean(y))
            else:
                if len(X) == 0:
                    X = np.full_like(Y, np.mean(x))
                elif ninc > 2:
                    # Adjustment for lat correction variation along face
                    lenX = X[-1] - X[0]
                    xi = xd / rinc
                    xn = xi / np.cos(np.radians(Y[0]))
                    X[0] = x[0] + xn / 2
                    # compare length of X and Y

                    for gg in range(1, ninc):
                        #print(f"indexing {gg}")
                        X[gg] = X[gg-1] + xi / np.cos(np.radians(Y[gg]))
                        # Correcting the last X
                    lenXnew = X[-1] - X[0]
                    tmp = X - X[0]
                    X = X[0] + tmp * lenX / lenXnew

            # Saving result in a list
            ii = list(range(ngrd, ngrd + ninc))
            ngrd = ngrd + ninc
            idx[ifc] = ii
            flo.extend(X)
            fla.extend(Y)

            # In this data all levels are above the bottom
            abovebot = np.column_stack([abovebot, np.ones((ndps, ninc))])

            # Volume interval is distance (changed to m*10^6 so transport in Sv)
            # by depth intervals. (The distance need to be changed from km to m)
            vint[ifc] = dint * rdist * 1000

        # Save first step data
        first_step_data = {
            'vint': vint,
            'abovebot': abovebot,
            'idx': idx,
            'flo': np.array(flo),
            'fla': np.array(fla),
            'ngrd': ngrd,
            'dirn': np.array(dirn),
            'lon': lon,
            'lat': lat,
            'zc': zc,
            'nfc': nfc
        }
        with open(file1, 'wb') as f_pkl:
            pickle.dump(first_step_data, f_pkl)

        nc_data.close()
        nc_mesh.close()

    else:
        with open(file1, 'rb') as f_pkl:
            first_step_data = pickle.load(f_pkl)

        vint = first_step_data['vint']
        abovebot = first_step_data['abovebot']
        idx = first_step_data['idx']
        flo = first_step_data['flo']
        fla = first_step_data['fla']
        ngrd = first_step_data['ngrd']
        dirn = first_step_data['dirn']
        lon = first_step_data['lon']
        lat = first_step_data['lat']
        zc = first_step_data['zc']
        nfc = first_step_data['nfc']
        
        # Reopen NetCDF files for second step processing
        nc_data = nc.Dataset(fnm, 'r')
        nc_mesh = nc.Dataset(fll, 'r')
        # Creation of the Final File

    file2 = f"{month}_{year}_SS_second_Step.npz"
    if not os.path.exists(file2):
        # Read velocity data for each time step
        for id in range(ntm):
            # Read velocity data for this time step
            u = nc_data.variables['u'][id, :, :, :]
            v = nc_data.variables['v'][id, :, :, :]

            # Handle missing values
            u[u >= 1e20] = np.nan
            v[v >= 1e20] = np.nan

            U = np.zeros((nlay, ngrd))
            V = np.zeros((nlay, ngrd))

            # Interpolation by layers

            for layer in range(nlay):
                # Average velocity data for this depth layer
                layer_mask = (zc <= -dlev[layer]) & (zc >= -dlev[layer+1])
                if np.any(layer_mask):
                    layer_u_data = np.nanmean(u[layer_mask, :, :], axis=0)
                    layer_v_data = np.nanmean(v[layer_mask, :, :], axis=0)
                    # Flatten arrays for masking and interpolation
                    layer_u_data_flat = layer_u_data.ravel()
                    layer_v_data_flat = layer_v_data.ravel()
                    lon_flat = lon.ravel()
                    lat_flat = lat.ravel()
                    valid_mask = ~(np.isnan(layer_u_data_flat) | np.isnan(layer_v_data_flat))
                    valid_lon = lon_flat[valid_mask]
                    valid_lat = lat_flat[valid_mask]
                    valid_u = layer_u_data_flat[valid_mask]
                    valid_v = layer_v_data_flat[valid_mask]

                    if len(valid_u) > 0:
                        # Interpolate to face integration points
                        U[layer, :] = griddata((valid_lon, valid_lat), valid_u, (flo, fla), method='linear', fill_value=np.nan)
                        V[layer, :] = griddata((valid_lon, valid_lat), valid_v, (flo, fla), method='linear', fill_value=np.nan)
                    else:
                        # If no valid data, fill with NaN
                        U[layer, :] = np.nan
                        V[layer, :] = np.nan
                        # Calculate transport for each face
            dir2, spd = cart2pol(U, V)
            for jj in range(nfc):
                ii = idx[jj]
                # Calculate velocity component along face
                tt = spd[:, ii] * np.sin(dir2[:, ii] - dirn[jj])

                # For each layer
                for lvs in range(nlay):
                    value = np.nanmean(tt[lvs, :])
                    T[jj, id, lvs] = vint[jj][lvs] * value

        # Save second step data as NPZ file
        np.savez_compressed(file2, T=T, tims=tims)

    else:
        # File already exists, just load to verify
        data = np.load(file2)
        T = data['T']
        tims = data['tims']
        data.close()

    nc_data.close()
    nc_mesh.close()
    
    # Function no longer returns values - only generates files
