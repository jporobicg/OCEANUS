#!/usr/bin/env python3
"""
Grid Comparison Script for OCEANUS
Compares BGM grid (boxes/faces) with hydrodynamic model grid
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import sys
import geopandas as gpd
from shapely.geometry import Point, LineString

# Add parent directory to path to import our functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from read_boxes import read_boxes
from read_faces import read_faces

def load_hydrodynamic_grid(nc_file):
    """Load hydrodynamic grid from NetCDF file"""
    print(f"Loading hydrodynamic grid from: {nc_file}")
    
    with xr.open_dataset(nc_file) as ds:
        # Get grid coordinates
        if 'lon' in ds.variables and 'lat' in ds.variables:
            lon = ds.lon.values
            lat = ds.lat.values
        elif 'longitude' in ds.variables and 'latitude' in ds.variables:
            lon = ds.longitude.values
            lat = ds.latitude.values
        else:
            # Try to find 2D coordinate variables
            coords = list(ds.coords.keys())
            print(f"Available coordinates: {coords}")
            
            # Look for common 2D coordinate names
            lon_names = ['lon', 'longitude', 'lon_rho', 'lon_u', 'lon_v']
            lat_names = ['lat', 'latitude', 'lat_rho', 'lat_u', 'lat_v']
            
            lon_var = None
            lat_var = None
            
            for name in lon_names:
                if name in ds.variables:
                    lon_var = name
                    break
                    
            for name in lat_names:
                if name in ds.variables:
                    lat_var = name
                    break
            
            if lon_var and lat_var:
                lon = ds[lon_var].values
                lat = ds[lat_var].values
            else:
                raise ValueError("Could not find longitude/latitude variables")
        
        # Get depth information if available
        depth_info = {}
        if 'depth' in ds.variables:
            depth_info['depth'] = ds.depth.values
        elif 'z' in ds.variables:
            depth_info['z'] = ds.z.values
        elif 'zc' in ds.variables:
            depth_info['zc'] = ds.zc.values
            
        print(f"Hydrodynamic grid shape: lon {lon.shape}, lat {lat.shape}")
        print(f"Lon range: {np.nanmin(lon):.4f} to {np.nanmax(lon):.4f}")
        print(f"Lat range: {np.nanmin(lat):.4f} to {np.nanmax(lat):.4f}")
        
        return lon, lat, depth_info

def load_shapefile_grid(shp_file):
    """Load grid from shapefile"""
    print(f"Loading shapefile grid from: {shp_file}")
    
    # Read the shapefile
    gdf = gpd.read_file(shp_file)
    
    # Extract centroids
    centroids = np.array([geom.centroid.coords[0] for geom in gdf.geometry])
    
    # Calculate areas (in square degrees)
    areas = np.array([geom.area for geom in gdf.geometry])
    
    # Extract boundaries for face lines
    face_pt1 = []
    face_pt2 = []
    
    for geom in gdf.geometry:
        if geom.geom_type == 'Polygon':
            # Get the exterior boundary
            coords = list(geom.exterior.coords)
            # Create face lines from consecutive points
            for i in range(len(coords) - 1):
                face_pt1.append(coords[i])
                face_pt2.append(coords[i + 1])
        elif geom.geom_type == 'MultiPolygon':
            for poly in geom.geoms:
                coords = list(poly.exterior.coords)
                for i in range(len(coords) - 1):
                    face_pt1.append(coords[i])
                    face_pt2.append(coords[i + 1])
    
    face_pt1 = np.array(face_pt1)
    face_pt2 = np.array(face_pt2)
    
    # For shapefile, we don't have depth information, so create dummy values
    depths = np.full(len(centroids), -50.0)  # Default depth
    
    print(f"Shapefile grid: {len(centroids)} polygons, {len(face_pt1)} face segments")
    print(f"Centroids range: lon {np.min(centroids[:, 0]):.4f} to {np.max(centroids[:, 0]):.4f}")
    print(f"Centroids range: lat {np.min(centroids[:, 1]):.4f} to {np.max(centroids[:, 1]):.4f}")
    
    return {
        'centroids': centroids,
        'areas': areas,
        'depths': depths,
        'face_pt1': face_pt1,
        'face_pt2': face_pt2,
        'nbox': len(centroids),
        'nface': len(face_pt1)
    }

def load_bgm_grid(bgm_file):
    """Load BGM grid (boxes and faces)"""
    print(f"Loading BGM grid from: {bgm_file}")
    
    # Read boxes first to get nbox
    nbox, nface, bid, centroids, areas, verts, iface, depths = read_boxes(bgm_file)
    
    # Create data for read_faces
    # Read faces
    nulr, nupt1, nupt2 = read_faces(nbox, nface, bid, verts, iface, bgm_file)
    
    # Convert to the format we need for plotting
    # Filter out invalid faces (where coordinates are NaN)
    valid_faces = ~np.isnan(nupt1[:, 0])
    face_pt1 = nupt1[valid_faces]
    face_pt2 = nupt2[valid_faces]
    
    print(f"BGM grid: {nbox} boxes, {np.sum(valid_faces)} valid faces")
    print(f"Centroids range: lon {np.min(centroids[:, 0]):.4f} to {np.max(centroids[:, 0]):.4f}")
    print(f"Centroids range: lat {np.min(centroids[:, 1]):.4f} to {np.max(centroids[:, 1]):.4f}")
    
    return {
        'centroids': centroids,
        'areas': areas,
        'depths': depths,
        'face_pt1': face_pt1,
        'face_pt2': face_pt2,
        'nulr': nulr[valid_faces],
        'nbox': nbox,
        'nface': nface
    }

def plot_grid_comparison(hydro_lon, hydro_lat, bgm_data, output_file=None):
    """Create a comparison plot of both grids"""
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: Hydrodynamic grid
    ax1.scatter(hydro_lon.flatten(), hydro_lat.flatten(), 
                c='blue', alpha=0.6, s=1, label='Hydrodynamic grid')
    ax1.set_title('Hydrodynamic Model Grid')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: BGM grid (boxes and faces)
    # Plot centroids
    ax2.scatter(bgm_data['centroids'][:, 0], bgm_data['centroids'][:, 1], 
                c='red', s=20, alpha=0.8, label='Box centroids')
    
    # Plot face lines
    for i in range(len(bgm_data['face_pt1'])):
        ax2.plot([bgm_data['face_pt1'][i, 0], bgm_data['face_pt2'][i, 0]],
                 [bgm_data['face_pt1'][i, 1], bgm_data['face_pt2'][i, 1]], 
                 'k-', alpha=0.3, linewidth=0.5)
    
    ax2.set_title('BGM Grid (Boxes and Faces)')
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Plot 3: Overlay comparison
    # Hydrodynamic grid (smaller, more transparent)
    ax3.scatter(hydro_lon.flatten(), hydro_lat.flatten(), 
                c='blue', alpha=0.3, s=0.5, label='Hydrodynamic grid')
    
    # BGM centroids
    ax3.scatter(bgm_data['centroids'][:, 0], bgm_data['centroids'][:, 1], 
                c='red', s=30, alpha=0.8, label='Box centroids')
    
    # BGM faces
    for i in range(len(bgm_data['face_pt1'])):
        ax3.plot([bgm_data['face_pt1'][i, 0], bgm_data['face_pt2'][i, 0]],
                 [bgm_data['face_pt1'][i, 1], bgm_data['face_pt2'][i, 1]], 
                 'k-', alpha=0.5, linewidth=1)
    
    ax3.set_title('Grid Overlay Comparison')
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    
    plt.show()

def analyze_spatial_coverage(hydro_lon, hydro_lat, bgm_data):
    """Analyze spatial coverage and overlap"""
    
    print("\n" + "="*50)
    print("SPATIAL COVERAGE ANALYSIS")
    print("="*50)
    
    # Hydrodynamic grid bounds
    hydro_lon_min, hydro_lon_max = np.nanmin(hydro_lon), np.nanmax(hydro_lon)
    hydro_lat_min, hydro_lat_max = np.nanmin(hydro_lat), np.nanmax(hydro_lat)
    
    # BGM grid bounds
    bgm_lon_min = np.min(bgm_data['centroids'][:, 0])
    bgm_lon_max = np.max(bgm_data['centroids'][:, 0])
    bgm_lat_min = np.min(bgm_data['centroids'][:, 1])
    bgm_lat_max = np.max(bgm_data['centroids'][:, 1])
    
    print(f"Hydrodynamic grid bounds:")
    print(f"  Longitude: {hydro_lon_min:.4f} to {hydro_lon_max:.4f}")
    print(f"  Latitude:  {hydro_lat_min:.4f} to {hydro_lat_max:.4f}")
    
    print(f"\nBGM grid bounds:")
    print(f"  Longitude: {bgm_lon_min:.4f} to {bgm_lon_max:.4f}")
    print(f"  Latitude:  {bgm_lat_min:.4f} to {bgm_lat_max:.4f}")
    
    # Check overlap
    lon_overlap = (bgm_lon_min <= hydro_lon_max) and (bgm_lon_max >= hydro_lon_min)
    lat_overlap = (bgm_lat_min <= hydro_lat_max) and (bgm_lat_max >= hydro_lat_min)
    
    print(f"\nSpatial overlap analysis:")
    print(f"  Longitude overlap: {'YES' if lon_overlap else 'NO'}")
    print(f"  Latitude overlap:  {'YES' if lat_overlap else 'NO'}")
    print(f"  Full overlap:      {'YES' if (lon_overlap and lat_overlap) else 'NO'}")
    
    if lon_overlap and lat_overlap:
        # Calculate overlap area
        overlap_lon_min = max(hydro_lon_min, bgm_lon_min)
        overlap_lon_max = min(hydro_lon_max, bgm_lon_max)
        overlap_lat_min = max(hydro_lat_min, bgm_lat_min)
        overlap_lat_max = min(hydro_lat_max, bgm_lat_max)
        
        print(f"\nOverlap region:")
        print(f"  Longitude: {overlap_lon_min:.4f} to {overlap_lon_max:.4f}")
        print(f"  Latitude:  {overlap_lat_min:.4f} to {overlap_lat_max:.4f}")
    
    # Count boxes within hydrodynamic domain
    within_domain = ((bgm_data['centroids'][:, 0] >= hydro_lon_min) & 
                    (bgm_data['centroids'][:, 0] <= hydro_lon_max) &
                    (bgm_data['centroids'][:, 1] >= hydro_lat_min) & 
                    (bgm_data['centroids'][:, 1] <= hydro_lat_max))
    
    print(f"\nBox coverage:")
    print(f"  Total boxes: {len(bgm_data['centroids'])}")
    print(f"  Boxes within hydrodynamic domain: {np.sum(within_domain)}")
    print(f"  Coverage percentage: {100 * np.sum(within_domain) / len(bgm_data['centroids']):.1f}%")

def main():
    """Main function to run grid comparison"""
    
    # File paths
    nc_file = "bass2_simple_2017-11.nc"
    shp_file = "/home/por07g/Documents/Projects/NESP_ParksAustralia/Hydro/shp_file/Clean_version/SEAP_NESP_final_ll.shp"
    output_plot = "grid_comparison_shapefile.png"
    
    print("OCEANUS Grid Comparison Tool (Shapefile)")
    print("="*40)
    
    # Check if files exist
    if not os.path.exists(nc_file):
        print(f"Error: Hydrodynamic file not found: {nc_file}")
        return
    
    if not os.path.exists(shp_file):
        print(f"Error: Shapefile not found: {shp_file}")
        return
    
    try:
        # Load grids
        hydro_lon, hydro_lat, depth_info = load_hydrodynamic_grid(nc_file)
        grid_data = load_shapefile_grid(shp_file)
        
        # Analyze spatial coverage
        analyze_spatial_coverage(hydro_lon, hydro_lat, grid_data)
        
        # Create comparison plot
        plot_grid_comparison(hydro_lon, hydro_lat, grid_data, output_plot)
        
        print(f"\nGrid comparison completed!")
        print(f"Plot saved to: {output_plot}")
        
    except Exception as e:
        print(f"Error during grid comparison: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 