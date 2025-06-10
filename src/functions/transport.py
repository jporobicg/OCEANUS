"""
Transport calculations for the Salish Sea Atlantis Model.

This module handles the calculation of transport between model boxes,
including vertical and horizontal fluxes.
"""

import numpy as np
import xarray as xr
from typing import Tuple, List, Dict, Optional
from pathlib import Path
import logging

from .coordinate_utils import lat_correction, distance_between_points
from .box_reader import read_boxes, read_faces

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TransportCalculator:
    def __init__(self, 
                 vert: List[np.ndarray],
                 pt1: np.ndarray,
                 pt2: np.ndarray,
                 depth_levels: np.ndarray,
                 distance_increment: float,
                 rim_points: int):
        """
        Initialize the transport calculator.
        
        Args:
            vert: List of vertex arrays for each box
            pt1: Face point 1 coordinates
            pt2: Face point 2 coordinates
            depth_levels: Array of depth levels
            distance_increment: Increment for distance calculations
            rim_points: Number of rim points for integration
        """
        self.vert = vert
        self.pt1 = pt1
        self.pt2 = pt2
        self.depth_levels = depth_levels
        self.dinc = distance_increment
        self.rimn = rim_points
        
    def calculate_transport(self, 
                          velocity_data: xr.Dataset,
                          mesh_file: str,
                          timestep: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate transport between boxes for a given timestep.
        
        Args:
            velocity_data: Dataset containing velocity components
            mesh_file: Path to the mesh file
            timestep: Current timestep to process
            
        Returns:
            Tuple containing:
            - Transport matrix
            - Timestamps
        """
        # Load mesh data
        mesh = xr.open_dataset(mesh_file)
        
        # Extract velocity components for the current timestep
        u = velocity_data['u'].isel(time=timestep)
        v = velocity_data['v'].isel(time=timestep)
        w = velocity_data.get('w', None)  # Vertical velocity if available
        
        # Initialize transport arrays
        n_faces = len(self.pt1)
        n_depths = len(self.depth_levels) - 1
        transport = np.zeros((n_faces, n_depths))
        
        # Calculate transport for each face
        for face_idx in range(n_faces):
            if np.any(np.isnan(self.pt1[face_idx])) or np.any(np.isnan(self.pt2[face_idx])):
                continue
                
            transport[face_idx, :] = self._calculate_face_transport(
                face_idx, u, v, w, mesh
            )
        
        return transport, velocity_data.time.values[timestep]
    
    def _calculate_face_transport(self,
                                face_idx: int,
                                u: xr.DataArray,
                                v: xr.DataArray,
                                w: Optional[xr.DataArray],
                                mesh: xr.Dataset) -> np.ndarray:
        """
        Calculate transport through a single face.
        
        Args:
            face_idx: Index of the face
            u: Eastward velocity component
            v: Northward velocity component
            w: Vertical velocity component (optional)
            mesh: Mesh dataset
            
        Returns:
            Array of transport values for each depth level
        """
        # Get face points
        p1 = self.pt1[face_idx]
        p2 = self.pt2[face_idx]
        
        # Calculate face properties
        face_length = distance_between_points(p1[1], p1[0], p2[1], p2[0])
        face_angle = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])
        
        # Create integration points along the face
        n_points = max(2, int(face_length / self.dinc))
        x_points = np.linspace(p1[0], p2[0], n_points)
        y_points = np.linspace(p1[1], p2[1], n_points)
        
        # Initialize transport array for this face
        transport = np.zeros(len(self.depth_levels) - 1)
        
        # Calculate transport at each depth level
        for d in range(len(self.depth_levels) - 1):
            depth_mid = (self.depth_levels[d] + self.depth_levels[d + 1]) / 2
            
            # Interpolate velocities to face points
            u_face = self._interpolate_to_points(u.sel(depth=depth_mid), x_points, y_points)
            v_face = self._interpolate_to_points(v.sel(depth=depth_mid), x_points, y_points)
            
            # Calculate normal velocity component
            normal_vel = (u_face * np.cos(face_angle + np.pi/2) + 
                        v_face * np.sin(face_angle + np.pi/2))
            
            # Calculate transport
            depth_layer = self.depth_levels[d + 1] - self.depth_levels[d]
            transport[d] = np.mean(normal_vel) * face_length * depth_layer
            
            # Add vertical transport if available
            if w is not None:
                w_face = self._interpolate_to_points(w.sel(depth=depth_mid), x_points, y_points)
                transport[d] += np.mean(w_face) * face_length * self.dinc
        
        return transport
    
    def _interpolate_to_points(self,
                             field: xr.DataArray,
                             x_points: np.ndarray,
                             y_points: np.ndarray) -> np.ndarray:
        """
        Interpolate a field to a set of points.
        
        Args:
            field: Field to interpolate
            x_points: X coordinates of target points
            y_points: Y coordinates of target points
            
        Returns:
            Interpolated values at target points
        """
        return field.interp(
            longitude=xr.DataArray(x_points),
            latitude=xr.DataArray(y_points),
            method='linear'
        ).values

def process_transport(vert: List[np.ndarray],
                     pt1: np.ndarray,
                     pt2: np.ndarray,
                     depth_levels: np.ndarray,
                     velocity_file: str,
                     mesh_file: str,
                     timestep: int,
                     distance_increment: float = 1.0,
                     rim_points: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    Process transport calculations for a single timestep.
    
    Args:
        vert: List of vertex arrays for each box
        pt1: Face point 1 coordinates
        pt2: Face point 2 coordinates
        depth_levels: Array of depth levels
        velocity_file: Path to velocity data file
        mesh_file: Path to mesh file
        timestep: Timestep to process
        distance_increment: Increment for distance calculations
        rim_points: Number of rim points for integration
        
    Returns:
        Tuple containing:
        - Transport matrix
        - Timestamps
    """
    # Create transport calculator
    calculator = TransportCalculator(
        vert=vert,
        pt1=pt1,
        pt2=pt2,
        depth_levels=depth_levels,
        distance_increment=distance_increment,
        rim_points=rim_points
    )
    
    # Load velocity data
    try:
        velocity_data = xr.open_dataset(velocity_file)
    except Exception as e:
        logger.error(f"Error loading velocity data from {velocity_file}: {e}")
        raise
        
    # Calculate transport
    try:
        transport, timestamp = calculator.calculate_transport(
            velocity_data=velocity_data,
            mesh_file=mesh_file,
            timestep=timestep
        )
    except Exception as e:
        logger.error(f"Error calculating transport for timestep {timestep}: {e}")
        raise
    finally:
        velocity_data.close()
        
    return transport, timestamp 