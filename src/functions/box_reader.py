"""
Box geometry reader for the Salish Sea Atlantis Model.

This module provides functionality to read and process box geometry model (BGM) files,
extracting box definitions, face relationships, and geometric properties.
"""

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Optional
import re

def read_boxes(bgm_file: str) -> Tuple[int, int, np.ndarray, np.ndarray, np.ndarray, List[np.ndarray], np.ndarray, np.ndarray]:
    """
    Read box definitions from a BGM file.
    
    Args:
        bgm_file: Path to the BGM file
        
    Returns:
        Tuple containing:
        - nbox: Number of boxes
        - nface: Number of faces
        - bid: Box IDs
        - cent: Centroids
        - b_area: Box areas
        - vert: List of vertex arrays for each box
        - iface: Face IDs
        - botz: Bottom depths
    """
    with open(bgm_file, 'r') as f:
        lines = f.readlines()
    
    # Extract projection information (if needed)
    proj_line = next(line for line in lines if line.startswith('projection'))
    
    # Find number of boxes and faces
    nbox = int(next(line for line in lines if line.startswith('nbox')).split()[1])
    nface = int(next(line for line in lines if line.startswith('nface')).split()[1])
    
    # Initialize arrays
    bid = np.zeros(nbox, dtype=int)
    cent = np.zeros((nbox, 2))
    b_area = np.zeros(nbox)
    vert = []
    botz = np.zeros(nbox)
    iface = np.zeros((nface, 2), dtype=int)
    
    box_pattern = re.compile(r'box(\d+)')
    face_pattern = re.compile(r'face(\d+)')
    
    current_box = -1
    vertices = []
    
    for line in lines:
        if box_pattern.match(line):
            if current_box >= 0 and vertices:
                vert.append(np.array(vertices))
                vertices = []
            current_box = int(box_pattern.match(line).group(1))
            bid[current_box] = current_box
        elif line.strip().startswith('area'):
            b_area[current_box] = float(line.split()[1])
        elif line.strip().startswith('centroid'):
            parts = line.split()
            cent[current_box] = [float(parts[1]), float(parts[2])]
        elif line.strip().startswith('botz'):
            botz[current_box] = float(line.split()[1])
        elif line.strip() and not any(p in line for p in ['projection', 'nbox', 'nface', 'box', 'face']):
            try:
                x, y = map(float, line.strip().split()[:2])
                vertices.append([x, y])
            except ValueError:
                continue
    
    # Add the last box's vertices
    if vertices:
        vert.append(np.array(vertices))
    
    return nbox, nface, bid, cent, b_area, vert, iface, botz

def read_faces(nbox: int, nface: int, bid: np.ndarray, vert: List[np.ndarray], 
               iface: np.ndarray, bgm_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read face definitions and compute face relationships.
    
    Args:
        nbox: Number of boxes
        nface: Number of faces
        bid: Box IDs
        vert: List of vertex arrays for each box
        iface: Face IDs
        bgm_file: Path to the BGM file
        
    Returns:
        Tuple containing:
        - nulr: Neighboring layers
        - nupt1: Face point 1
        - nupt2: Face point 2
    """
    # Initialize arrays
    nulr = np.full((nface, 2), np.nan)  # Neighboring layers
    nupt1 = np.full((nface, 2), np.nan)  # Face point 1
    nupt2 = np.full((nface, 2), np.nan)  # Face point 2
    
    with open(bgm_file, 'r') as f:
        lines = f.readlines()
    
    face_pattern = re.compile(r'face(\d+)')
    current_face = -1
    
    for i, line in enumerate(lines):
        if face_pattern.match(line):
            current_face = int(face_pattern.match(line).group(1))
            # Get the next line containing the face points
            points_line = lines[i + 1].strip()
            try:
                x1, y1, x2, y2 = map(float, points_line.split()[:4])
                nupt1[current_face] = [x1, y1]
                nupt2[current_face] = [x2, y2]
            except ValueError:
                continue
            
            # Find connecting boxes for this face
            for b1 in range(nbox):
                for b2 in range(b1 + 1, nbox):
                    if check_face_connection(vert[b1], vert[b2], nupt1[current_face], nupt2[current_face]):
                        nulr[current_face] = [b1, b2]
                        break
                if not np.isnan(nulr[current_face, 0]):
                    break
    
    return nulr, nupt1, nupt2

def check_face_connection(box1_verts: np.ndarray, box2_verts: np.ndarray, 
                         pt1: np.ndarray, pt2: np.ndarray, tolerance: float = 1e-10) -> bool:
    """
    Check if two boxes share a face defined by two points.
    
    Args:
        box1_verts: Vertices of first box
        box2_verts: Vertices of second box
        pt1: First point of the face
        pt2: Second point of the face
        tolerance: Numerical tolerance for floating point comparison
        
    Returns:
        bool: True if boxes share the face, False otherwise
    """
    def point_on_edge(point: np.ndarray, box_verts: np.ndarray) -> bool:
        for i in range(len(box_verts)):
            v1 = box_verts[i]
            v2 = box_verts[(i + 1) % len(box_verts)]
            
            # Check if point lies on the edge
            d = np.abs(np.cross(v2 - v1, point - v1)) / np.linalg.norm(v2 - v1)
            if d < tolerance:
                # Check if point lies between vertices
                t = np.dot(point - v1, v2 - v1) / np.dot(v2 - v1, v2 - v1)
                if 0 <= t <= 1:
                    return True
        return False
    
    # Check if both points lie on edges of both boxes
    return (point_on_edge(pt1, box1_verts) and point_on_edge(pt2, box1_verts) and
            point_on_edge(pt1, box2_verts) and point_on_edge(pt2, box2_verts)) 