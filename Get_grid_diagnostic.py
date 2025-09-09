#!/usr/bin/env python3
"""
Grid Diagnostic Analysis for OCEANUS
====================================

This script analyzes the coverage of BGM boxes and faces against the hydrodynamic model domain.
It identifies boxes and faces that fall outside the valid hydrodynamic domain and generates
diagnostic reports and visualizations.

Author: OCEANUS Development Team
Date: 2024
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import os
from datetime import datetime
from read_boxes import read_boxes
from read_faces import read_faces


def get_hydrodynamic_domain(mesh_file):
    """
    Extract the hydrodynamic model domain from mesh.nc file.
    
    Parameters:
    -----------
    mesh_file : str
        Path to the mesh.nc file
        
    Returns:
    --------
    dict
        Dictionary containing domain information:
        - lon_min, lon_max: longitude bounds
        - lat_min, lat_max: latitude bounds
        - lon_grid, lat_grid: full grid coordinates
        - domain_mask: boolean mask of valid ocean points
    """
    print(f"Reading hydrodynamic domain from: {mesh_file}")
    
    with nc.Dataset(mesh_file, 'r') as ds:
        # Read longitude and latitude
        lon = ds.variables['longitude'][:]
        lat = ds.variables['latitude'][:]
        
        # Get domain bounds
        lon_min, lon_max = np.nanmin(lon), np.nanmax(lon)
        lat_min, lat_max = np.nanmin(lat), np.nanmax(lat)
        
        # Create domain mask (exclude NaN values)
        domain_mask = ~(np.isnan(lon) | np.isnan(lat))
        
        print(f"Hydrodynamic domain bounds:")
        print(f"  Longitude: {lon_min:.4f} to {lon_max:.4f}")
        print(f"  Latitude:  {lat_min:.4f} to {lat_max:.4f}")
        print(f"  Grid size: {lon.shape}")
        
        return {
            'lon_min': lon_min,
            'lon_max': lon_max,
            'lat_min': lat_min,
            'lat_max': lat_max,
            'lon_grid': lon,
            'lat_grid': lat,
            'domain_mask': domain_mask
        }


def check_box_coverage(box_centroids, domain_info):
    """
    Check which BGM boxes fall within the hydrodynamic domain.
    
    Parameters:
    -----------
    box_centroids : np.ndarray
        Array of box centroids (n_boxes, 2) with [lon, lat]
    domain_info : dict
        Domain information from get_hydrodynamic_domain()
        
    Returns:
    --------
    dict
        Coverage analysis results for boxes
    """
    print("Analyzing box coverage...")
    
    n_boxes = len(box_centroids)
    inside_domain = []
    outside_domain = []
    
    for i, centroid in enumerate(box_centroids):
        lon, lat = centroid[0], centroid[1]
        
        # Check if centroid is within domain bounds
        in_bounds = (domain_info['lon_min'] <= lon <= domain_info['lon_max'] and
                    domain_info['lat_min'] <= lat <= domain_info['lat_max'])
        
        if in_bounds:
            inside_domain.append(i)
        else:
            outside_domain.append(i)
    
    coverage_percentage = (len(inside_domain) / n_boxes) * 100
    
    print(f"Box coverage analysis:")
    print(f"  Total boxes: {n_boxes}")
    print(f"  Inside domain: {len(inside_domain)} ({coverage_percentage:.1f}%)")
    print(f"  Outside domain: {len(outside_domain)} ({100-coverage_percentage:.1f}%)")
    
    return {
        'total_boxes': n_boxes,
        'inside_domain': inside_domain,
        'outside_domain': outside_domain,
        'coverage_percentage': coverage_percentage
    }


def check_face_coverage(face_points, domain_info):
    """
    Check which BGM faces fall within the hydrodynamic domain.
    
    Parameters:
    -----------
    face_points : tuple
        Tuple of (pt1, pt2) arrays with face endpoints
    domain_info : dict
        Domain information from get_hydrodynamic_domain()
        
    Returns:
    --------
    dict
        Coverage analysis results for faces
    """
    print("Analyzing face coverage...")
    
    pt1, pt2 = face_points
    n_faces = len(pt1)
    inside_domain = []
    outside_domain = []
    partially_inside = []
    
    for i in range(n_faces):
        # Check if face endpoints are valid (not NaN)
        if np.isnan(pt1[i, 0]) or np.isnan(pt2[i, 0]):
            outside_domain.append(i)
            continue
            
        lon1, lat1 = pt1[i, 0], pt1[i, 1]
        lon2, lat2 = pt2[i, 0], pt2[i, 1]
        
        # Check if both endpoints are within domain bounds
        p1_in_bounds = (domain_info['lon_min'] <= lon1 <= domain_info['lon_max'] and
                       domain_info['lat_min'] <= lat1 <= domain_info['lat_max'])
        p2_in_bounds = (domain_info['lon_min'] <= lon2 <= domain_info['lon_max'] and
                       domain_info['lat_min'] <= lat2 <= domain_info['lat_max'])
        
        if p1_in_bounds and p2_in_bounds:
            inside_domain.append(i)
        elif p1_in_bounds or p2_in_bounds:
            partially_inside.append(i)
        else:
            outside_domain.append(i)
    
    total_valid_faces = len(inside_domain) + len(partially_inside) + len(outside_domain)
    coverage_percentage = (len(inside_domain) / total_valid_faces) * 100 if total_valid_faces > 0 else 0
    
    print(f"Face coverage analysis:")
    print(f"  Total faces: {n_faces}")
    print(f"  Inside domain: {len(inside_domain)} ({coverage_percentage:.1f}%)")
    print(f"  Partially inside: {len(partially_inside)}")
    print(f"  Outside domain: {len(outside_domain)} ({100-coverage_percentage:.1f}%)")
    
    return {
        'total_faces': n_faces,
        'inside_domain': inside_domain,
        'partially_inside': partially_inside,
        'outside_domain': outside_domain,
        'coverage_percentage': coverage_percentage
    }


def plot_box_coverage(box_centroids, box_analysis, domain_info, box_vertices, output_file):
    """
    Create visualization of box coverage against hydrodynamic domain.
    
    Parameters:
    -----------
    box_centroids : np.ndarray
        Array of box centroids
    box_analysis : dict
        Results from check_box_coverage()
    domain_info : dict
        Domain information
    box_vertices : list
        List of box vertex arrays
    output_file : str
        Output file path for the plot
    """
    print(f"Creating box coverage plot: {output_file}")
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Calculate BGM domain bounds
    all_lons = []
    all_lats = []
    for vertices in box_vertices:
        if len(vertices) > 0:
            all_lons.extend(vertices[:, 0])
            all_lats.extend(vertices[:, 1])
    
    bgm_lon_min, bgm_lon_max = np.min(all_lons), np.max(all_lons)
    bgm_lat_min, bgm_lat_max = np.min(all_lats), np.max(all_lats)
    
    # Plot hydrodynamic domain grid
    lon_grid = domain_info['lon_grid']
    lat_grid = domain_info['lat_grid']
    domain_mask = domain_info['domain_mask']
    
    # Plot valid ocean points
    valid_lon = lon_grid[domain_mask]
    valid_lat = lat_grid[domain_mask]
    ax.scatter(valid_lon, valid_lat, c='lightblue', s=1, alpha=0.3, label='Hydrodynamic Grid')
    
    # Plot BGM box polygons
    for i, vertices in enumerate(box_vertices):
        if len(vertices) > 0:
            # Determine color based on whether box is inside domain
            if i in box_analysis['inside_domain']:
                color = 'lightgreen'
                alpha = 0.3
            else:
                color = 'lightcoral'
                alpha = 0.3
            
            # Plot polygon
            polygon = plt.Polygon(vertices, facecolor=color, edgecolor='gray', 
                                linewidth=0.5, alpha=alpha)
            ax.add_patch(polygon)
    
    # Plot box centroids
    inside_boxes = box_centroids[box_analysis['inside_domain']]
    outside_boxes = box_centroids[box_analysis['outside_domain']]
    
    if len(inside_boxes) > 0:
        ax.scatter(inside_boxes[:, 0], inside_boxes[:, 1], 
                  c='green', s=50, alpha=0.7, label='Boxes Inside Domain')
        
        # Add box numbers for inside boxes
        for i, box_id in enumerate(box_analysis['inside_domain']):
            ax.annotate(f'B{box_id}', (inside_boxes[i, 0], inside_boxes[i, 1]), 
                       fontsize=8, ha='center', va='bottom')
    
    if len(outside_boxes) > 0:
        ax.scatter(outside_boxes[:, 0], outside_boxes[:, 1], 
                  c='red', s=50, alpha=0.7, label='Boxes Outside Domain')
        
        # Add box numbers for outside boxes
        for i, box_id in enumerate(box_analysis['outside_domain']):
            ax.annotate(f'B{box_id}', (outside_boxes[i, 0], outside_boxes[i, 1]), 
                       fontsize=8, ha='center', va='bottom', color='red')
    
    # Set plot properties
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('BGM Box Coverage vs Hydrodynamic Domain')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Set axis limits to BGM bounds with some padding
    padding_lon = (bgm_lon_max - bgm_lon_min) * 0.05
    padding_lat = (bgm_lat_max - bgm_lat_min) * 0.05
    ax.set_xlim(bgm_lon_min - padding_lon, bgm_lon_max + padding_lon)
    ax.set_ylim(bgm_lat_min - padding_lat, bgm_lat_max + padding_lat)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Box coverage plot saved: {output_file}")


def plot_face_coverage(face_points, face_analysis, domain_info, box_vertices, output_file):
    """
    Create visualization of face coverage against hydrodynamic domain.
    
    Parameters:
    -----------
    face_points : tuple
        Tuple of (pt1, pt2) arrays
    face_analysis : dict
        Results from check_face_coverage()
    domain_info : dict
        Domain information
    box_vertices : list
        List of box vertex arrays (for domain bounds calculation)
    output_file : str
        Output file path for the plot
    """
    print(f"Creating face coverage plot: {output_file}")
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Calculate BGM domain bounds
    all_lons = []
    all_lats = []
    for vertices in box_vertices:
        if len(vertices) > 0:
            all_lons.extend(vertices[:, 0])
            all_lats.extend(vertices[:, 1])
    
    bgm_lon_min, bgm_lon_max = np.min(all_lons), np.max(all_lons)
    bgm_lat_min, bgm_lat_max = np.min(all_lats), np.max(all_lats)
    
    # Plot hydrodynamic domain grid
    lon_grid = domain_info['lon_grid']
    lat_grid = domain_info['lat_grid']
    domain_mask = domain_info['domain_mask']
    
    # Plot valid ocean points
    valid_lon = lon_grid[domain_mask]
    valid_lat = lat_grid[domain_mask]
    ax.scatter(valid_lon, valid_lat, c='lightblue', s=1, alpha=0.3, label='Hydrodynamic Grid')
    
    pt1, pt2 = face_points
    
    # Plot faces by category
    colors = {'inside': 'green', 'partial': 'orange', 'outside': 'red'}
    labels = {'inside': 'Faces Inside Domain', 'partial': 'Faces Partially Inside', 'outside': 'Faces Outside Domain'}
    
    for category, face_ids in [('inside', face_analysis['inside_domain']),
                              ('partial', face_analysis['partially_inside']),
                              ('outside', face_analysis['outside_domain'])]:
        
        if len(face_ids) > 0:
            for face_id in face_ids:
                if not (np.isnan(pt1[face_id, 0]) or np.isnan(pt2[face_id, 0])):
                    x = [pt1[face_id, 0], pt2[face_id, 0]]
                    y = [pt1[face_id, 1], pt2[face_id, 1]]
                    ax.plot(x, y, c=colors[category], alpha=0.7, linewidth=1)
                    
                    # Add face number at midpoint
                    mid_x, mid_y = np.mean(x), np.mean(y)
                    ax.annotate(f'F{face_id}', (mid_x, mid_y), 
                               fontsize=6, ha='center', va='bottom', 
                               color=colors[category] if category != 'inside' else 'black')
    
    # Add legend
    for category, color in colors.items():
        key = f'{category}_domain' if category != 'partial' else 'partially_inside'
        if len(face_analysis[key]) > 0:
            ax.plot([], [], c=color, label=labels[category])
    
    # Set plot properties
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('BGM Face Coverage vs Hydrodynamic Domain')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Set axis limits to BGM bounds with some padding
    padding_lon = (bgm_lon_max - bgm_lon_min) * 0.05
    padding_lat = (bgm_lat_max - bgm_lat_min) * 0.05
    ax.set_xlim(bgm_lon_min - padding_lon, bgm_lon_max + padding_lon)
    ax.set_ylim(bgm_lat_min - padding_lat, bgm_lat_max + padding_lat)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Face coverage plot saved: {output_file}")


def write_diagnostic_report(box_analysis, face_analysis, domain_info, output_file):
    """
    Write detailed diagnostic report to .diag file.
    
    Parameters:
    -----------
    box_analysis : dict
        Results from check_box_coverage()
    face_analysis : dict
        Results from check_face_coverage()
    domain_info : dict
        Domain information
    output_file : str
        Output file path for the .diag file
    """
    print(f"Writing diagnostic report: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("OCEANUS Grid Diagnostic Report\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Domain Information
        f.write("HYDRODYNAMIC DOMAIN INFORMATION\n")
        f.write("-" * 35 + "\n")
        f.write(f"Longitude bounds: {domain_info['lon_min']:.4f} to {domain_info['lon_max']:.4f}\n")
        f.write(f"Latitude bounds:  {domain_info['lat_min']:.4f} to {domain_info['lat_max']:.4f}\n")
        f.write(f"Grid dimensions: {domain_info['lon_grid'].shape}\n\n")
        
        # Box Analysis
        f.write("BOX COVERAGE ANALYSIS\n")
        f.write("-" * 25 + "\n")
        f.write(f"Total boxes: {box_analysis['total_boxes']}\n")
        f.write(f"Inside domain: {len(box_analysis['inside_domain'])} ({box_analysis['coverage_percentage']:.1f}%)\n")
        f.write(f"Outside domain: {len(box_analysis['outside_domain'])} ({100-box_analysis['coverage_percentage']:.1f}%)\n\n")
        
        if box_analysis['outside_domain']:
            f.write("Boxes outside domain:\n")
            for box_id in box_analysis['outside_domain']:
                f.write(f"  Box {box_id}\n")
            f.write("\n")
        
        # Face Analysis
        f.write("FACE COVERAGE ANALYSIS\n")
        f.write("-" * 26 + "\n")
        f.write(f"Total faces: {face_analysis['total_faces']}\n")
        f.write(f"Inside domain: {len(face_analysis['inside_domain'])} ({face_analysis['coverage_percentage']:.1f}%)\n")
        f.write(f"Partially inside: {len(face_analysis['partially_inside'])}\n")
        f.write(f"Outside domain: {len(face_analysis['outside_domain'])}\n\n")
        
        if face_analysis['outside_domain']:
            f.write("Faces outside domain:\n")
            for face_id in face_analysis['outside_domain']:
                f.write(f"  Face {face_id}\n")
            f.write("\n")
        
        if face_analysis['partially_inside']:
            f.write("Faces partially inside domain:\n")
            for face_id in face_analysis['partially_inside']:
                f.write(f"  Face {face_id}\n")
            f.write("\n")
        
        # Recommendations
        f.write("RECOMMENDATIONS\n")
        f.write("-" * 15 + "\n")
        
        if box_analysis['outside_domain']:
            f.write("1. Consider removing or adjusting boxes outside the hydrodynamic domain\n")
            f.write("2. These boxes may cause issues in transport calculations\n")
        
        if face_analysis['outside_domain']:
            f.write("3. Faces outside domain will have zero or invalid transport values\n")
            f.write("4. Consider implementing interpolation for these faces\n")
        
        if face_analysis['partially_inside']:
            f.write("5. Partially inside faces may need special handling in transport calculations\n")
        
        f.write("\n6. Run mass balance correction after addressing domain coverage issues\n")
        
        print(f"Diagnostic report saved: {output_file}")


def main():
    """
    Main function to run the grid diagnostic analysis.
    """
    print("OCEANUS Grid Diagnostic Analysis")
    print("=" * 40)
    
    # File paths
    bgm_file = "Test/SEAP_NESP_final_ll_fixed.bgm"
    mesh_file = "Test/mesh.nc"
    
    # Check if files exist
    if not os.path.exists(bgm_file):
        print(f"Error: BGM file not found: {bgm_file}")
        return
    
    if not os.path.exists(mesh_file):
        print(f"Error: Mesh file not found: {mesh_file}")
        return
    
    # Create Diagnostics directory if it doesn't exist
    os.makedirs("Diagnostics", exist_ok=True)
    
    # Read hydrodynamic domain
    domain_info = get_hydrodynamic_domain(mesh_file)
    
    # Read BGM geometry
    print(f"\nReading BGM geometry from: {bgm_file}")
    nbox, nface, bid, cent, b_area, vert, iface, botz = read_boxes(bgm_file)
    nulr, nupt1, nupt2 = read_faces(nbox, nface, bid, vert, iface, bgm_file)
    
    # Analyze coverage
    box_analysis = check_box_coverage(cent, domain_info)
    face_analysis = check_face_coverage((nupt1, nupt2), domain_info)
    
    # Generate visualizations
    plot_box_coverage(cent, box_analysis, domain_info, vert, "Diagnostics/grid_coverage_boxes.png")
    plot_face_coverage((nupt1, nupt2), face_analysis, domain_info, vert, "Diagnostics/grid_coverage_faces.png")
    
    # Write diagnostic report
    write_diagnostic_report(box_analysis, face_analysis, domain_info, "Diagnostics/domain_analysis.diag")
    
    print("\nGrid diagnostic analysis completed!")
    print("Output files:")
    print("  - Diagnostics/grid_coverage_boxes.png")
    print("  - Diagnostics/grid_coverage_faces.png")
    print("  - Diagnostics/domain_analysis.diag")


if __name__ == "__main__":
    main()
