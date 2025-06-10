"""
Coordinate utilities for the Salish Sea Atlantis Model.

This module provides functions for coordinate transformations and corrections,
particularly for latitude-dependent calculations.
"""

import numpy as np
from typing import Union, Tuple

def lat_correction(lat: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate the latitude correction factor for distance calculations.
    
    This function computes the correction factor needed to account for the 
    Earth's sphericity when calculating distances at different latitudes.
    
    Args:
        lat: Latitude in degrees (scalar or array)
        
    Returns:
        Correction factor(s) for the given latitude(s)
    """
    # Convert latitude to radians
    lat_rad = np.deg2rad(lat)
    
    # Calculate correction factor
    return np.cos(lat_rad)

def distance_between_points(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Calculate the great circle distance between two points on Earth.
    
    Args:
        lat1: Latitude of first point in degrees
        lon1: Longitude of first point in degrees
        lat2: Latitude of second point in degrees
        lon2: Longitude of second point in degrees
        
    Returns:
        Distance in kilometers
    """
    # Earth's radius in kilometers
    R = 6371.0
    
    # Convert coordinates to radians
    lat1, lon1, lat2, lon2 = map(np.deg2rad, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    return R * c

def interpolate_at_latitude(lat: float, values: np.ndarray, lats: np.ndarray) -> float:
    """
    Interpolate a value at a specific latitude using surrounding data points.
    
    Args:
        lat: Target latitude for interpolation
        values: Array of values corresponding to latitudes
        lats: Array of latitudes where values are known
        
    Returns:
        Interpolated value at the target latitude
    """
    return np.interp(lat, lats, values)

def coordinate_to_index(coord: float, coord_array: np.ndarray) -> Tuple[int, float]:
    """
    Find the index and fractional position of a coordinate in an array.
    
    Args:
        coord: Target coordinate value
        coord_array: Array of coordinate values
        
    Returns:
        Tuple containing:
        - Index of the lower bound
        - Fractional position between lower and upper bound (0-1)
    """
    if coord <= coord_array[0]:
        return 0, 0.0
    if coord >= coord_array[-1]:
        return len(coord_array) - 2, 1.0
    
    idx = np.searchsorted(coord_array, coord) - 1
    frac = (coord - coord_array[idx]) / (coord_array[idx + 1] - coord_array[idx])
    
    return idx, frac 