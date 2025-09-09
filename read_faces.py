# read_faces.py
# Function to read face definitions from a BGM file

import numpy as np

def read_faces(nbox, nface, bid, verts, iface, fname):
    """
    Read face definitions from BGM file - simplified version matching MATLAB read_faces2.m
    
    Returns:
        nulr: array of left/right box IDs for each face
        nupt1: array of first point coordinates for each face  
        nupt2: array of second point coordinates for each face
    """
    # Initialize arrays like MATLAB
    nulr = np.full((nface, 2), -9)  # left/right box IDs
    nupt1 = np.full((nface, 2), np.nan)  # first point
    nupt2 = np.full((nface, 2), np.nan)  # second point
    
    current_face = None
    
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Find face number
            if '# Face' in line:
                face_num = int(line.split('Face number')[1].strip())
                current_face = face_num  # MATLAB uses 1-based indexing
            
            # Extract first point coordinates
            elif '.p1' in line:
                if current_face is not None:
                    coords = line.split('.p1')[1].strip().split()
                    nupt1[current_face] = [float(coords[0]), float(coords[1])]
            
            # Extract second point coordinates  
            elif '.p2' in line:
                if current_face is not None:
                    coords = line.split('.p2')[1].strip().split()
                    nupt2[current_face] = [float(coords[0]), float(coords[1])]
            
            # Extract left/right box IDs
            elif '.lr' in line:
                if current_face is not None:
                    lr_ids = line.split('.lr')[1].strip().split()
                    nulr[current_face] = [int(lr_ids[0]), int(lr_ids[1])]
    
    return nulr, nupt1, nupt2 