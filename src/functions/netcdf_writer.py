"""
NetCDF writer for the Salish Sea Atlantis Model.

This module handles writing transport and variable data to NetCDF files
in the format required by the Atlantis model.
"""

import numpy as np
import xarray as xr
from typing import List, Dict, Optional, Union
import logging
from datetime import datetime
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def write_transport_file(output_file: str,
                        pt1: np.ndarray,
                        pt2: np.ndarray,
                        layer_connections: np.ndarray,
                        timestamps: np.ndarray,
                        transport: np.ndarray,
                        face_ids: np.ndarray) -> None:
    """
    Write transport data to a NetCDF file.
    
    Args:
        output_file: Path to output NetCDF file
        pt1: Face point 1 coordinates
        pt2: Face point 2 coordinates
        layer_connections: Array of layer connections
        timestamps: Array of timestamps
        transport: Transport data array
        face_ids: Array of face IDs
    """
    # Create dataset
    ds = xr.Dataset(
        data_vars={
            'transport': (('time', 'face', 'depth'), transport),
            'face_pt1': (('face', 'coord'), pt1),
            'face_pt2': (('face', 'coord'), pt2),
            'layer_connections': (('face', 'connection'), layer_connections),
            'face_id': ('face', face_ids)
        },
        coords={
            'time': timestamps,
            'face': np.arange(len(face_ids)),
            'depth': np.arange(transport.shape[2]),
            'coord': ['lon', 'lat'],
            'connection': ['box1', 'box2']
        }
    )
    
    # Add metadata
    ds.transport.attrs.update({
        'units': 'm3/s',
        'long_name': 'Volume transport across faces',
        'description': 'Positive transport is from box1 to box2'
    })
    
    ds.face_pt1.attrs.update({
        'units': 'degrees',
        'long_name': 'Face point 1 coordinates'
    })
    
    ds.face_pt2.attrs.update({
        'units': 'degrees',
        'long_name': 'Face point 2 coordinates'
    })
    
    ds.layer_connections.attrs.update({
        'long_name': 'Connected box indices',
        'description': 'Indices of boxes connected by each face'
    })
    
    ds.attrs.update({
        'title': 'Salish Sea Atlantis Model Transport Data',
        'creation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'contact': 'CSIRO Marine Research'
    })
    
    # Write to file
    try:
        ds.to_netcdf(output_file)
        logger.info(f"Successfully wrote transport data to {output_file}")
    except Exception as e:
        logger.error(f"Error writing transport data to {output_file}: {e}")
        raise

def write_variable_file(output_file: str,
                       timestamps: np.ndarray,
                       box_ids: np.ndarray,
                       variables: Dict[str, np.ndarray]) -> None:
    """
    Write variable data to a NetCDF file.
    
    Args:
        output_file: Path to output NetCDF file
        timestamps: Array of timestamps
        box_ids: Array of box IDs
        variables: Dictionary of variable arrays
    """
    # Create dataset with coordinates
    ds = xr.Dataset(
        coords={
            'time': timestamps,
            'box': box_ids,
            'depth': np.arange(next(iter(variables.values())).shape[2])
        }
    )
    
    # Add each variable to the dataset
    for var_name, var_data in variables.items():
        ds[var_name] = (('time', 'box', 'depth'), var_data)
        
        # Add variable metadata
        ds[var_name].attrs.update({
            'units': get_variable_units(var_name),
            'long_name': get_variable_long_name(var_name),
            'description': get_variable_description(var_name)
        })
    
    # Add global metadata
    ds.attrs.update({
        'title': 'Salish Sea Atlantis Model Variable Data',
        'creation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'contact': 'CSIRO Marine Research'
    })
    
    # Write to file
    try:
        ds.to_netcdf(output_file)
        logger.info(f"Successfully wrote variable data to {output_file}")
    except Exception as e:
        logger.error(f"Error writing variable data to {output_file}: {e}")
        raise

def get_variable_units(var_name: str) -> str:
    """
    Get the appropriate units for a variable.
    
    Args:
        var_name: Name of the variable
        
    Returns:
        String containing the units
    """
    units_dict = {
        'temperature': 'degC',
        'salinity': 'PSU',
        'wVelocity': 'm/s',
        'microzooplankton': 'mg N/m3',
        'mesozooplankton': 'mg N/m3',
        'diatoms': 'mg N/m3',
        'dissolved_organic_nitrogen': 'mg N/m3',
        'ammonium': 'mg N/m3',
        'nitrate': 'mg N/m3',
        'silicon': 'mg Si/m3',
        'biogenic_silicon': 'mg Si/m3'
    }
    
    return units_dict.get(var_name, 'dimensionless')

def get_variable_long_name(var_name: str) -> str:
    """
    Get the long name for a variable.
    
    Args:
        var_name: Name of the variable
        
    Returns:
        String containing the long name
    """
    long_names = {
        'temperature': 'Water Temperature',
        'salinity': 'Salinity',
        'wVelocity': 'Vertical Velocity',
        'microzooplankton': 'Microzooplankton Concentration',
        'mesozooplankton': 'Mesozooplankton Concentration',
        'diatoms': 'Diatom Concentration',
        'dissolved_organic_nitrogen': 'Dissolved Organic Nitrogen',
        'ammonium': 'Ammonium Concentration',
        'nitrate': 'Nitrate Concentration',
        'silicon': 'Silicon Concentration',
        'biogenic_silicon': 'Biogenic Silicon Concentration'
    }
    
    return long_names.get(var_name, var_name.replace('_', ' ').title())

def get_variable_description(var_name: str) -> str:
    """
    Get the description for a variable.
    
    Args:
        var_name: Name of the variable
        
    Returns:
        String containing the description
    """
    descriptions = {
        'temperature': 'Water temperature averaged within each box',
        'salinity': 'Salinity averaged within each box',
        'wVelocity': 'Vertical velocity component averaged within each box',
        'microzooplankton': 'Concentration of microzooplankton biomass',
        'mesozooplankton': 'Concentration of mesozooplankton biomass',
        'diatoms': 'Concentration of diatom biomass',
        'dissolved_organic_nitrogen': 'Concentration of dissolved organic nitrogen',
        'ammonium': 'Concentration of ammonium',
        'nitrate': 'Concentration of nitrate',
        'silicon': 'Concentration of dissolved silicon',
        'biogenic_silicon': 'Concentration of biogenic silicon'
    }
    
    return descriptions.get(var_name, f'Average {var_name.replace("_", " ")} within each box') 