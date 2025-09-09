#!/usr/bin/env python3
"""
Mass Balance Correction for OCEANUS Transport
=============================================

This script corrects transport fluxes to ensure mass balance (net flux = 0) for each box.
It only adjusts faces that currently have no flux values (NaN or zero).

Author: OCEANUS Development Team
Date: 2024
"""

import numpy as np
import netCDF4 as nc
import os
from datetime import datetime


def load_transport_data(transport_file):
    """
    Load transport data from NetCDF file.
    
    Parameters:
    -----------
    transport_file : str
        Path to the transport NetCDF file
        
    Returns:
    --------
    dict
        Dictionary containing transport data and metadata
    """
    print(f"Loading transport data from: {transport_file}")
    
    with nc.Dataset(transport_file, 'r') as ds:
        # Read transport data
        transport = ds.variables['transport'][:]  # (time, faces, level)
        
        # Read face connectivity
        dest = ds.variables['dest_boxid'][:]  # destination box for each face
        sour = ds.variables['source_boxid'][:]  # source box for each face
        
        # Read face coordinates
        pt1x = ds.variables['pt1_x'][:]
        pt1y = ds.variables['pt1_y'][:]
        pt2x = ds.variables['pt2_x'][:]
        pt2y = ds.variables['pt2_y'][:]
        
        # Read time and level data
        time = ds.variables['time'][:]
        level = ds.variables['level'][:]
        
        print(f"Transport data shape: {transport.shape}")
        print(f"Number of faces: {len(dest)}")
        print(f"Number of time steps: {len(time)}")
        print(f"Number of levels: {len(level)}")
        
        return {
            'transport': transport,
            'dest': dest,
            'sour': sour,
            'pt1x': pt1x,
            'pt1y': pt1y,
            'pt2x': pt2x,
            'pt2y': pt2y,
            'time': time,
            'level': level
        }


def calculate_box_net_flux(transport_data, box_id, time_step, level):
    """
    Calculate net flux for a specific box at a given time and level.
    
    Parameters:
    -----------
    transport_data : dict
        Transport data dictionary
    box_id : int
        Box ID to calculate net flux for
    time_step : int
        Time step index
    level : int
        Level index
        
    Returns:
    --------
    float
        Net flux for the box (positive = net inflow, negative = net outflow)
    """
    dest = transport_data['dest']
    sour = transport_data['sour']
    transport = transport_data['transport']
    
    net_flux = 0.0
    
    # Check all faces
    for face_id in range(len(dest)):
        flux = transport[time_step, face_id, level]
        
        # Skip NaN values
        if np.isnan(flux):
            continue
            
        # If box is destination, flux is positive (inflow)
        if dest[face_id] == box_id:
            net_flux += flux
        # If box is source, flux is negative (outflow)
        elif sour[face_id] == box_id:
            net_flux -= flux
    
    return net_flux


def find_empty_faces(transport_data, time_step, level):
    """
    Find faces that have no flux values (NaN or zero).
    
    Parameters:
    -----------
    transport_data : dict
        Transport data dictionary
    time_step : int
        Time step index
    level : int
        Level index
        
    Returns:
    --------
    list
        List of face IDs with no flux values
    """
    transport = transport_data['transport']
    empty_faces = []
    
    for face_id in range(transport.shape[1]):
        flux = transport[time_step, face_id, level]
        if np.isnan(flux) or flux == 0.0:
            empty_faces.append(face_id)
    
    return empty_faces


def find_connected_boxes(transport_data, box_id):
    """
    Find boxes connected to a given box through faces.
    
    Parameters:
    -----------
    transport_data : dict
        Transport data dictionary
    box_id : int
        Box ID to find connections for
        
    Returns:
    --------
    list
        List of connected box IDs
    """
    dest = transport_data['dest']
    sour = transport_data['sour']
    connected_boxes = []
    
    for face_id in range(len(dest)):
        if dest[face_id] == box_id:
            connected_boxes.append(sour[face_id])
        elif sour[face_id] == box_id:
            connected_boxes.append(dest[face_id])
    
    return list(set(connected_boxes))  # Remove duplicates


def apply_mass_balance_correction(transport_data, tolerance=0.01):
    """
    Apply mass balance correction to transport data.
    
    Parameters:
    -----------
    transport_data : dict
        Transport data dictionary
    tolerance : float
        Tolerance for considering a box balanced (default: 1%)
        
    Returns:
    --------
    dict
        Corrected transport data
    """
    print("Applying mass balance correction...")
    
    # Create a copy of transport data
    corrected_transport = transport_data['transport'].copy()
    dest = transport_data['dest']
    sour = transport_data['sour']
    
    n_time_steps = corrected_transport.shape[0]
    n_faces = corrected_transport.shape[1]
    n_levels = corrected_transport.shape[2]
    
    # Get unique box IDs
    all_boxes = set(dest) | set(sour)
    all_boxes.discard(-9)  # Remove invalid box ID
    all_boxes = sorted(list(all_boxes))
    
    print(f"Processing {len(all_boxes)} boxes across {n_time_steps} time steps and {n_levels} levels")
    
    corrections_applied = 0
    
    # Process each time step and level
    for time_step in range(n_time_steps):
        for level in range(n_levels):
            print(f"Processing time step {time_step+1}/{n_time_steps}, level {level+1}/{n_levels}")
            
            # Find empty faces for this time/level
            empty_faces = find_empty_faces(transport_data, time_step, level)
            
            if len(empty_faces) == 0:
                continue  # No empty faces to work with
            
            # Calculate net flux for each box
            box_imbalances = {}
            for box_id in all_boxes:
                net_flux = calculate_box_net_flux(transport_data, box_id, time_step, level)
                if abs(net_flux) > tolerance:
                    box_imbalances[box_id] = net_flux
            
            # Apply corrections using empty faces
            for box_id, imbalance in box_imbalances.items():
                if len(empty_faces) == 0:
                    break  # No more empty faces to use
                
                # Find empty faces connected to this box
                connected_empty_faces = []
                for face_id in empty_faces:
                    if dest[face_id] == box_id or sour[face_id] == box_id:
                        connected_empty_faces.append(face_id)
                
                if len(connected_empty_faces) == 0:
                    continue  # No empty faces connected to this box
                
                # Apply correction to the first available empty face
                face_id = connected_empty_faces[0]
                
                # Determine correction direction
                if dest[face_id] == box_id:
                    # Box is destination, positive flux means inflow
                    correction = -imbalance
                else:
                    # Box is source, negative flux means outflow
                    correction = imbalance
                
                # Apply correction
                corrected_transport[time_step, face_id, level] = correction
                
                # Remove this face from empty faces list
                empty_faces.remove(face_id)
                
                corrections_applied += 1
                
                print(f"  Applied correction: Box {box_id}, Face {face_id}, Correction: {correction:.6f} Sv")
    
    print(f"Total corrections applied: {corrections_applied}")
    
    # Create corrected transport data
    corrected_data = transport_data.copy()
    corrected_data['transport'] = corrected_transport
    
    return corrected_data


def validate_mass_balance(transport_data, tolerance=0.01):
    """
    Validate mass balance for all boxes.
    
    Parameters:
    -----------
    transport_data : dict
        Transport data dictionary
    tolerance : float
        Tolerance for considering a box balanced
        
    Returns:
    --------
    dict
        Validation results
    """
    print("Validating mass balance...")
    
    dest = transport_data['dest']
    sour = transport_data['sour']
    transport = transport_data['transport']
    
    # Get unique box IDs
    all_boxes = set(dest) | set(sour)
    all_boxes.discard(-9)  # Remove invalid box ID
    all_boxes = sorted(list(all_boxes))
    
    n_time_steps = transport.shape[0]
    n_levels = transport.shape[2]
    
    validation_results = {
        'total_checks': 0,
        'balanced_boxes': 0,
        'unbalanced_boxes': 0,
        'max_imbalance': 0.0,
        'avg_imbalance': 0.0,
        'imbalances': []
    }
    
    total_imbalance = 0.0
    
    for time_step in range(n_time_steps):
        for level in range(n_levels):
            for box_id in all_boxes:
                net_flux = calculate_box_net_flux(transport_data, box_id, time_step, level)
                
                validation_results['total_checks'] += 1
                validation_results['imbalances'].append(abs(net_flux))
                
                if abs(net_flux) <= tolerance:
                    validation_results['balanced_boxes'] += 1
                else:
                    validation_results['unbalanced_boxes'] += 1
                    total_imbalance += abs(net_flux)
                
                validation_results['max_imbalance'] = max(validation_results['max_imbalance'], abs(net_flux))
    
    validation_results['avg_imbalance'] = total_imbalance / validation_results['total_checks']
    
    print(f"Validation results:")
    print(f"  Total checks: {validation_results['total_checks']}")
    print(f"  Balanced boxes: {validation_results['balanced_boxes']}")
    print(f"  Unbalanced boxes: {validation_results['unbalanced_boxes']}")
    print(f"  Balance percentage: {(validation_results['balanced_boxes']/validation_results['total_checks'])*100:.1f}%")
    print(f"  Maximum imbalance: {validation_results['max_imbalance']:.6f} Sv")
    print(f"  Average imbalance: {validation_results['avg_imbalance']:.6f} Sv")
    
    return validation_results


def write_corrected_transport(transport_data, output_file):
    """
    Write corrected transport data to NetCDF file.
    
    Parameters:
    -----------
    transport_data : dict
        Corrected transport data
    output_file : str
        Output file path
    """
    print(f"Writing corrected transport data to: {output_file}")
    
    with nc.Dataset(output_file, 'w', format='NETCDF4') as ds:
        # Create dimensions
        time_dim = ds.createDimension('time', None)  # UNLIMITED dimension
        faces_dim = ds.createDimension('faces', len(transport_data['dest']))
        level_dim = ds.createDimension('level', len(transport_data['level']))
        
        # Create variables
        time_var = ds.createVariable('time', 'f8', ('time',))
        level_var = ds.createVariable('level', 'i4', ('level',))
        faces_var = ds.createVariable('faces', 'i4', ('faces',))
        
        pt1x_var = ds.createVariable('pt1_x', 'f4', ('faces',))
        pt1y_var = ds.createVariable('pt1_y', 'f4', ('faces',))
        pt2x_var = ds.createVariable('pt2_x', 'f4', ('faces',))
        pt2y_var = ds.createVariable('pt2_y', 'f4', ('faces',))
        
        dest_var = ds.createVariable('dest_boxid', 'i4', ('faces',))
        sour_var = ds.createVariable('source_boxid', 'i4', ('faces',))
        
        transp_var = ds.createVariable('transport', 'f4', ('time', 'faces', 'level'))
        
        # Set variable attributes
        time_var.setncattr('long_name', 'simulation_time')
        time_var.setncattr('units', 'seconds since 1990-01-01 00:00:00 +10')
        time_var.setncattr('grads_step', 'seconds')
        time_var.setncattr('dt', 86400)
        
        level_var.setncattr('long_name', 'layer index; 1=near-surface')
        
        faces_var.setncattr('long_name', 'Face IDs')
        
        pt1x_var.setncattr('long_name', 'x coordinate of point 1 of face')
        pt1x_var.setncattr('units', 'degree_east')
        pt1y_var.setncattr('long_name', 'y coordinate of point 1 of face')
        pt1y_var.setncattr('units', 'degree_north')
        pt2x_var.setncattr('long_name', 'x coordinate of point 2 of face')
        pt2x_var.setncattr('units', 'degree_east')
        pt2y_var.setncattr('long_name', 'y coordinate of point 2 of face')
        pt2y_var.setncattr('units', 'degree_north')
        
        dest_var.setncattr('long_name', 'ID of destination box')
        dest_var.setncattr('units', 'IDs')
        sour_var.setncattr('long_name', 'ID of source box')
        sour_var.setncattr('units', 'IDs')
        
        transp_var.setncattr('long_name', 'flux across face')
        transp_var.setncattr('units', 'meter second-1')
        transp_var.setncattr('comment', 'm/s is to left, viewing from pt1 to pt2')
        transp_var.setncattr('FillValue_', -1e20)
        
        # Set global attributes
        str_desc = ('Transports across faces of box model derived from hydrodynamic'
                    ' model current field. Mass balance corrected.')
        his = f'Created Using Javier Porobic (more details : https://github.com/jporobicg) codes on {datetime.now().strftime("%Y-%m-%d")}'

        ds.setncattr('Title', 'Transport Bass Strait Australian Model BASS2 (Mass Balance Corrected)')
        ds.setncattr('Comments', str_desc)
        ds.setncattr('Institution', 'CSIRO')
        ds.setncattr('Conventions', 'Standar')
        ds.setncattr('history', his)
        
        # Leave define mode
        ds._enddef()
        
        # Write data
        time_var[:] = transport_data['time']
        level_var[:] = transport_data['level']
        faces_var[:] = np.arange(len(transport_data['dest']))
        
        pt1x_var[:] = transport_data['pt1x']
        pt1y_var[:] = transport_data['pt1y']
        pt2x_var[:] = transport_data['pt2x']
        pt2y_var[:] = transport_data['pt2y']
        
        dest_var[:] = transport_data['dest']
        sour_var[:] = transport_data['sour']
        
        transp_var[:] = transport_data['transport']
    
    print(f"Corrected transport file written: {output_file}")


def main():
    """
    Main function to run mass balance correction.
    """
    print("OCEANUS Mass Balance Correction")
    print("=" * 35)
    
    # File paths
    input_file = "Test/Transport_test.nc"
    output_file = "Test/Transport_test_mass_balanced.nc"
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Transport file not found: {input_file}")
        return
    
    # Load transport data
    transport_data = load_transport_data(input_file)
    
    # Validate initial mass balance
    print("\n" + "="*50)
    print("INITIAL MASS BALANCE VALIDATION")
    print("="*50)
    initial_validation = validate_mass_balance(transport_data)
    
    # Apply mass balance correction
    print("\n" + "="*50)
    print("APPLYING MASS BALANCE CORRECTION")
    print("="*50)
    corrected_data = apply_mass_balance_correction(transport_data, tolerance=0.01)
    
    # Validate final mass balance
    print("\n" + "="*50)
    print("FINAL MASS BALANCE VALIDATION")
    print("="*50)
    final_validation = validate_mass_balance(corrected_data)
    
    # Write corrected transport file
    print("\n" + "="*50)
    print("WRITING CORRECTED TRANSPORT FILE")
    print("="*50)
    write_corrected_transport(corrected_data, output_file)
    
    # Summary
    print("\n" + "="*50)
    print("MASS BALANCE CORRECTION SUMMARY")
    print("="*50)
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print(f"Initial balance: {initial_validation['balanced_boxes']/initial_validation['total_checks']*100:.1f}%")
    print(f"Final balance: {final_validation['balanced_boxes']/final_validation['total_checks']*100:.1f}%")
    print(f"Improvement: {final_validation['balanced_boxes']/final_validation['total_checks']*100 - initial_validation['balanced_boxes']/initial_validation['total_checks']*100:.1f}%")
    
    print("\nMass balance correction completed!")


if __name__ == "__main__":
    main()
