#!/usr/bin/env python3
"""
Surface Temperature and Salinity Visualization with Atlantis Model Polygons
=======================================================================

This script creates side-by-side plots of surface temperature and salinity
from the original hydrodynamic data with Atlantis model box polygons overlaid.

Author: OCEANUS Development Team
Date: 2024
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from read_boxes import read_boxes

def load_hydrodynamic_data(data_file):
    """
    Load hydrodynamic data from NetCDF file.
    
    Parameters:
    -----------
    data_file : str
        Path to the NetCDF data file
        
    Returns:
    --------
    dict
        Dictionary containing data arrays and coordinates
    """
    print(f"Loading hydrodynamic data from: {data_file}")
    
    with nc.Dataset(data_file, 'r') as ds:
        # Read coordinates
        lon = ds.variables['longitude'][:]
        lat = ds.variables['latitude'][:]
        zc = ds.variables['zc'][:]  # depth levels
        
        # Read surface data (first time step, surface level)
        temp = ds.variables['temp'][0, 58, :, :]  # [time, level, lat, lon]
        salt = ds.variables['salt'][0, 58, :, :]  # [time, level, lat, lon]
        
        # Handle missing values
        temp[temp >= 1e20] = np.nan
        salt[salt >= 1e20] = np.nan
        
        print(f"Data shape: {temp.shape}")
        print(f"Longitude range: {lon.min():.3f} to {lon.max():.3f}")
        print(f"Latitude range: {lat.min():.3f} to {lat.max():.3f}")
        print(f"Temperature range: {np.nanmin(temp):.2f} to {np.nanmax(temp):.2f} °C")
        print(f"Salinity range: {np.nanmin(salt):.2f} to {np.nanmax(salt):.2f} PSU")
        
        return {
            'lon': lon,
            'lat': lat,
            'temp': temp,
            'salt': salt,
            'zc': zc
        }

def load_atlantis_boxes(bgm_file):
    """
    Load Atlantis model box geometry.
    
    Parameters:
    -----------
    bgm_file : str
        Path to the BGM file
        
    Returns:
    --------
    dict
        Dictionary containing box geometry
    """
    print(f"Loading Atlantis model boxes from: {bgm_file}")
    
    nbox, nface, bid, cent, b_area, vert, iface, botz = read_boxes(bgm_file)
    
    print(f"Number of boxes: {nbox}")
    print(f"Number of faces: {nface}")
    
    return {
        'nbox': nbox,
        'nface': nface,
        'bid': bid,
        'cent': cent,
        'b_area': b_area,
        'vert': vert,
        'iface': iface,
        'botz': botz
    }

def create_surface_plots(data, boxes, output_file='surface_temperature_salinity.png'):
    """
    Create side-by-side plots of surface temperature and salinity with Atlantis polygons.
    
    Parameters:
    -----------
    data : dict
        Hydrodynamic data dictionary
    boxes : dict
        Atlantis box geometry dictionary
    output_file : str
        Output filename for the plot
    """
    print("Creating surface plots...")
    
    # Set up the figure with two subplots (vertically arranged)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 16), 
                                   subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Define color maps
    temp_cmap = plt.cm.RdYlBu_r
    salt_cmap = plt.cm.viridis
    
    # Plot 1: Surface Temperature
    print("Plotting surface temperature...")
    im1 = ax1.pcolormesh(data['lon'], data['lat'], data['temp'], 
                        cmap=temp_cmap, transform=ccrs.PlateCarree(),
                        shading='auto', alpha=0.8)
    
    # Add coastlines and features
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.3)
    ax1.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
    ax1.add_feature(cfeature.LAND, color='lightgray', alpha=0.3)
    
    # Add Atlantis model polygons
    for i in range(boxes['nbox']):
        if len(boxes['vert'][i]) > 0:
            # Create polygon from vertices
            polygon = patches.Polygon(boxes['vert'][i], 
                                    closed=True, 
                                    linewidth=0.8, 
                                    edgecolor='black', 
                                    facecolor='none',
                                    alpha=0.8)
            ax1.add_patch(polygon)
            
            # Add box ID labels at centroids
            if not np.isnan(boxes['cent'][i, 0]) and not np.isnan(boxes['cent'][i, 1]):
                ax1.text(boxes['cent'][i, 0], boxes['cent'][i, 1], str(i), 
                        fontsize=6, ha='center', va='center', 
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    
    ax1.set_title('Surface Temperature (°C)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Longitude', fontsize=12)
    ax1.set_ylabel('Latitude', fontsize=12)
    ax1.gridlines(draw_labels=True, alpha=0.5)
    
    # Add colorbar for temperature (on the right side)
    cbar1 = plt.colorbar(im1, ax=ax1, orientation='vertical', pad=0.05, shrink=0.8)
    cbar1.set_label('Temperature (°C)', fontsize=10)
    
    # Plot 2: Surface Salinity
    print("Plotting surface salinity...")
    im2 = ax2.pcolormesh(data['lon'], data['lat'], data['salt'], 
                        cmap=salt_cmap, transform=ccrs.PlateCarree(),
                        shading='auto', alpha=0.8)
    
    # Add coastlines and features
    ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.3)
    ax2.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
    ax2.add_feature(cfeature.LAND, color='lightgray', alpha=0.3)
    
    # Add Atlantis model polygons
    for i in range(boxes['nbox']):
        if len(boxes['vert'][i]) > 0:
            # Create polygon from vertices
            polygon = patches.Polygon(boxes['vert'][i], 
                                    closed=True, 
                                    linewidth=0.8, 
                                    edgecolor='black', 
                                    facecolor='none',
                                    alpha=0.8)
            ax2.add_patch(polygon)
            
            # Add box ID labels at centroids
            if not np.isnan(boxes['cent'][i, 0]) and not np.isnan(boxes['cent'][i, 1]):
                ax2.text(boxes['cent'][i, 0], boxes['cent'][i, 1], str(i), 
                        fontsize=6, ha='center', va='center', 
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    
    ax2.set_title('Surface Salinity (PSU)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Longitude', fontsize=12)
    ax2.set_ylabel('Latitude', fontsize=12)
    ax2.gridlines(draw_labels=True, alpha=0.5)
    
    # Add colorbar for salinity (on the right side)
    cbar2 = plt.colorbar(im2, ax=ax2, orientation='vertical', pad=0.05, shrink=0.8)
    cbar2.set_label('Salinity (PSU)', fontsize=10)
    
    # Set the extent to focus on the data region
    lon_min, lon_max = data['lon'].min(), data['lon'].max()
    lat_min, lat_max = data['lat'].min(), data['lat'].max()
    
    # Add some padding
    lon_pad = (lon_max - lon_min) * 0.05
    lat_pad = (lat_max - lat_min) * 0.05
    
    ax1.set_extent([lon_min - lon_pad, lon_max + lon_pad, 
                    lat_min - lat_pad, lat_max + lat_pad], ccrs.PlateCarree())
    ax2.set_extent([lon_min - lon_pad, lon_max + lon_pad, 
                    lat_min - lat_pad, lat_max + lat_pad], ccrs.PlateCarree())
    
    # Add overall title
    fig.suptitle('OCEANUS: Surface Temperature and Salinity with Atlantis Model Boxes', 
                 fontsize=16, fontweight='bold', y=0.95)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.subplots_adjust(top=0.95, hspace=0.3)
    
    print(f"Saving plot to: {output_file}")
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print("Plot creation completed!")

def main():
    """
    Main function to create surface temperature and salinity plots.
    """
    print("OCEANUS Surface Temperature and Salinity Visualization")
    print("=" * 55)
    
    # File paths
    data_file = "Test/bass2_simple_2017-11.nc"
    bgm_file = "Test/SEAP_NESP_final_ll_fixed.bgm"
    output_file = "surface_temperature_salinity_atlantis.png"
    
    try:
        # Load hydrodynamic data
        data = load_hydrodynamic_data(data_file)
        
        # Load Atlantis model boxes
        boxes = load_atlantis_boxes(bgm_file)
        
        # Create plots
        create_surface_plots(data, boxes, output_file)
        
        print(f"\nVisualization completed successfully!")
        print(f"Output file: {output_file}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
