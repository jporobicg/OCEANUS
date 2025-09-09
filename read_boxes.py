# read_boxes.py
# Function to read box geometry from a BGM file

import numpy as np
import re

def read_boxes(bgm_file):
    """
    Reads the BGM file and extracts box geometry information.
    Returns:
        nbox: number of boxes
        nface: number of faces
        bid: numpy array of box IDs
        cent: numpy array of centroids (nbox, 2)
        b_area: numpy array of box areas
        vert: list of numpy arrays of vertices for each box
        iface: list of face IDs for each box
        botz: numpy array of bottom depths
    """
    nbox = None
    nface = None
    # First pass: get nbox and nface
    with open(bgm_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('nbox'):
                nbox = int(line.split()[1])
            elif line.startswith('nface'):
                nface = int(line.split()[1])
            if nbox is not None and nface is not None:
                break
    # Prepare storage
    bid = np.arange(nbox, dtype=int)
    cent = np.zeros((nbox, 2))
    b_area = np.zeros(nbox)
    vert = [[] for _ in range(nbox)]
    iface = [[] for _ in range(nbox)]
    botz = np.zeros(nbox)
    # Second pass: parse box info
    current_box = None
    with open(bgm_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Detect start of a new box
            m = re.match(r'box(\d+)\.label', line)
            if m:
                current_box = int(m.group(1))
                continue
            if current_box is not None:
                # Centroid (inside)
                if line.startswith(f'box{current_box}.inside'):
                    parts = line.split()
                    cent[current_box, 0] = float(parts[1])
                    cent[current_box, 1] = float(parts[2])
                # Area
                elif line.startswith(f'box{current_box}.area'):
                    b_area[current_box] = float(line.split()[1])
                # Bottom depth
                elif line.startswith(f'box{current_box}.botz'):
                    botz[current_box] = float(line.split()[1])
                # Vertices
                elif line.startswith(f'box{current_box}.vert ') and not line.startswith(f'box{current_box}.vertmix'):
                    parts = line.split()
                    if len(parts) >= 3:
                        x, y = float(parts[1]), float(parts[2])
                        vert[current_box].append([x, y])
                    else:
                        print(f"Warning: Malformed vert line for box {current_box}: {line}")
                # Face IDs
                elif line.startswith(f'box{current_box}.iface'):
                    # Face IDs are all on one line after the keyword
                    parts = line.split()
                    if len(parts) > 1:
                        iface[current_box] = [int(fid) for fid in parts[1:]]
                    else:
                        # If no face IDs on this line, it might be empty or on next line
                        print(f"Info: Empty iface for box {current_box}")
    # Convert vertices to numpy arrays
    vert = [np.array(v) if len(v) > 0 else np.zeros((0, 2)) for v in vert]
    return nbox, nface, bid, cent, b_area, vert, iface, botz 