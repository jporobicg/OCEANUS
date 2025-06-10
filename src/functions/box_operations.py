"""
Box operations for the Salish Sea Atlantis Model.

This module handles calculations and operations performed within model boxes,
including averaging and interpolation of variables.
"""

import numpy as np
import xarray as xr
from typing import Tuple, List, Dict, Optional
import logging
from pathlib import Path

from .coordinate_utils import lat_correction, coordinate_to_index

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class BoxProcessor:
    def __init__(self,
                 vertices: List[np.ndarray],
                 depth_levels: np.ndarray):
        """
        Initialize the box processor.
        
        Args:
            vertices: List of vertex arrays defining box boundaries
            depth_levels: Array of depth levels for vertical discretization
        """
        self.vertices = vertices
        self.depth_levels = depth_levels
        self.n_boxes = len(vertices)
        self.n_depths = len(depth_levels) - 1
        
    def calculate_box_averages(self,
                             variable_data: xr.Dataset,
                             variable_name: str,
                             mesh_file: str,
                             timestep: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate spatial averages of a variable within each box.
        
        Args:
            variable_data: Dataset containing the variable
            variable_name: Name of the variable to process
            mesh_file: Path to the mesh file
            timestep: Current timestep to process
            
        Returns:
            Tuple containing:
            - Array of box averages
            - Timestamps
        """
        # Load mesh data
        mesh = xr.open_dataset(mesh_file)
        
        # Extract variable for current timestep
        var = variable_data[variable_name].isel(time=timestep)
        
        # Initialize output array
        box_averages = np.zeros((self.n_boxes, self.n_depths))
        
        # Process each box
        for box_idx in range(self.n_boxes):
            box_averages[box_idx, :] = self._calculate_box_average(
                box_idx, var, mesh
            )
            
        return box_averages, variable_data.time.values[timestep]
    
    def _calculate_box_average(self,
                             box_idx: int,
                             variable: xr.DataArray,
                             mesh: xr.Dataset) -> np.ndarray:
        """
        Calculate average values within a single box.
        
        Args:
            box_idx: Index of the box
            variable: Variable data array
            mesh: Mesh dataset
            
        Returns:
            Array of averaged values for each depth level
        """
        # Get box vertices
        verts = self.vertices[box_idx]
        
        # Find bounding box
        lon_min, lon_max = np.min(verts[:, 0]), np.max(verts[:, 0])
        lat_min, lat_max = np.min(verts[:, 1]), np.max(verts[:, 1])
        
        # Extract relevant portion of the variable field
        var_subset = variable.sel(
            longitude=slice(lon_min, lon_max),
            latitude=slice(lat_min, lat_max)
        )
        
        # Create mask for points inside the box
        mask = self._create_box_mask(var_subset.longitude, var_subset.latitude, verts)
        
        # Initialize output array
        averages = np.zeros(self.n_depths)
        
        # Calculate averages for each depth level
        for d in range(self.n_depths):
            depth_mid = (self.depth_levels[d] + self.depth_levels[d + 1]) / 2
            
            # Get variable values at this depth
            var_depth = var_subset.sel(depth=depth_mid, method='nearest')
            
            # Apply mask and calculate average
            valid_points = var_depth.values[mask]
            if len(valid_points) > 0:
                averages[d] = np.nanmean(valid_points)
            else:
                averages[d] = np.nan
                
        return averages
    
    def _create_box_mask(self,
                        lons: xr.DataArray,
                        lats: xr.DataArray,
                        vertices: np.ndarray) -> np.ndarray:
        """
        Create a boolean mask for points inside a polygon.
        
        Args:
            lons: Longitude coordinates
            lats: Latitude coordinates
            vertices: Vertices of the polygon
            
        Returns:
            Boolean mask array
        """
        # Create coordinate meshgrid
        lon_grid, lat_grid = np.meshgrid(lons, lats)
        
        # Initialize mask
        mask = np.zeros_like(lon_grid, dtype=bool)
        
        # Use ray casting algorithm to determine points inside polygon
        for i in range(len(lon_grid)):
            for j in range(len(lon_grid[0])):
                if self._point_in_polygon(lon_grid[i, j], lat_grid[i, j], vertices):
                    mask[i, j] = True
                    
        return mask
    
    @staticmethod
    def _point_in_polygon(x: float, y: float, vertices: np.ndarray) -> bool:
        """
        Determine if a point is inside a polygon using ray casting algorithm.
        
        Args:
            x: X-coordinate of the point
            y: Y-coordinate of the point
            vertices: Array of polygon vertices
            
        Returns:
            True if point is inside polygon, False otherwise
        """
        n = len(vertices)
        inside = False
        
        p1x, p1y = vertices[0]
        for i in range(n + 1):
            p2x, p2y = vertices[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
            
        return inside

def process_box_averages(vertices: List[np.ndarray],
                        depth_levels: np.ndarray,
                        variable_file: str,
                        variable_name: str,
                        mesh_file: str,
                        timestep: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Process box averages for a variable at a single timestep.
    
    Args:
        vertices: List of vertex arrays defining box boundaries
        depth_levels: Array of depth levels
        variable_file: Path to variable data file
        variable_name: Name of the variable to process
        mesh_file: Path to mesh file
        timestep: Timestep to process
        
    Returns:
        Tuple containing:
        - Array of box averages
        - Timestamps
    """
    # Create box processor
    processor = BoxProcessor(
        vertices=vertices,
        depth_levels=depth_levels
    )
    
    # Load variable data
    try:
        variable_data = xr.open_dataset(variable_file)
    except Exception as e:
        logger.error(f"Error loading variable data from {variable_file}: {e}")
        raise
        
    # Calculate averages
    try:
        averages, timestamp = processor.calculate_box_averages(
            variable_data=variable_data,
            variable_name=variable_name,
            mesh_file=mesh_file,
            timestep=timestep
        )
    except Exception as e:
        logger.error(f"Error calculating box averages for timestep {timestep}: {e}")
        raise
    finally:
        variable_data.close()
        
    return averages, timestamp 