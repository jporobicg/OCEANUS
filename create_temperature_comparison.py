#!/usr/bin/env python3
"""
Create temperature comparison plot: Original BASS-2 vs Transformed Atlantis
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import matplotlib.patches as patches
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec

def load_bass2_temperature():
    """Load original BASS-2 temperature data"""
    print("Loading BASS-2 temperature data...")
    
    # Load BASS-2 data
    bass2_file = 'Test/bass2_simple_2017-11.nc'
    ds = Dataset(bass2_file, 'r')
    
    # Extract temperature at surface (level 58) for time step 0
    temp = ds.variables['temp'][0, 58, :, :]  # [time, level, lat, lon]
    
    # Get coordinates
    lons = ds.variables['longitude'][:]
    lats = ds.variables['latitude'][:]
    
    ds.close()
    
    return temp, lons, lats

def load_atlantis_temperature():
    """Load transformed Atlantis temperature data"""
    print("Loading Atlantis temperature data...")
    
    # Load Atlantis variables data
    atlantis_file = 'Test/Variables_combined.nc'
    ds = Dataset(atlantis_file, 'r')
    
    # Extract temperature for time step 0, surface layer (index 0)
    # Shape is [time, boxes, layers], so we need [0, :, 0]
    temp = ds.variables['temperature'][0, :, 0]  # [time, box, layer]
    
    ds.close()
    
    return temp

def load_atlantis_geometry():
    """Load Atlantis model geometry from BGM file"""
    print("Loading Atlantis geometry from BGM file...")
    
    # Load BGM file to get box coordinates
    bgm_file = 'Test/SEAP_NESP_final_ll_fixed.bgm'
    
    boxes = []
    with open(bgm_file, 'r') as f:
        lines = f.readlines()
    
    # Parse BGM file to extract box coordinates
    current_box = []
    box_count = 0
    
    for line in lines:
        line = line.strip()
        
        # Check if this is a vertex line for a box
        if line.startswith('box') and '.vert ' in line:
            # Extract coordinates from vertex line
            parts = line.split()
            if len(parts) >= 3:
                try:
                    lon = float(parts[1])
                    lat = float(parts[2])
                    current_box.append([lon, lat])
                except ValueError:
                    continue
        
        # Check if this is the start of a new box (but not a vertex line)
        elif line.startswith('box') and '.vert ' not in line and '.label' not in line:
            # Save previous box if it exists
            if current_box:
                boxes.append(current_box)
                box_count += 1
            current_box = []
    
    # Add the last box if exists
    if current_box:
        boxes.append(current_box)
        box_count += 1
    
    print(f"Loaded {box_count} boxes from BGM file")
    return boxes

def create_temperature_comparison():
    """Create side-by-side temperature comparison plot"""
    print("Creating temperature comparison plot...")
    
    # Load data
    bass2_temp, bass2_lons, bass2_lats = load_bass2_temperature()
    atlantis_temp = load_atlantis_temperature()
    boxes = load_atlantis_geometry()
    
    # Create figure with subplots - vertical layout
    fig = plt.figure(figsize=(12, 16))
    
    # Create grid layout with shared colorbar - vertical arrangement
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 0.05], width_ratios=[1, 1], hspace=0.3, wspace=0.1)
    
    # Define projection
    proj = ccrs.PlateCarree()
    
    # Find common temperature range for shared colorbar
    temp_min = min(np.nanmin(bass2_temp), np.nanmin(atlantis_temp))
    temp_max = max(np.nanmax(bass2_temp), np.nanmax(atlantis_temp))
    
    # Create normalization for shared colorbar
    norm = Normalize(vmin=temp_min, vmax=temp_max)
    
    # Plot 1: Original BASS-2 temperature
    ax1 = fig.add_subplot(gs[0, :], projection=proj)
    
    # Plot BASS-2 temperature
    im1 = ax1.pcolormesh(bass2_lons, bass2_lats, bass2_temp, 
                        transform=proj, cmap='RdYlBu_r', norm=norm, shading='nearest')
    
    # Add coastlines and features
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.3)
    ax1.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
    ax1.add_feature(cfeature.LAND, color='lightgray', alpha=0.3)
    
    # Set extent - expanded to show all polygons (based on BGM bounds: 140.08-152.36°E, -41.89 to -36.99°N)
    ax1.set_extent([139.8, 152.6, -42.2, -36.7], crs=proj)
    
    # Add title
    ax1.set_title('Original BASS-2 Temperature\nTime Step 0, Surface Layer', 
                 fontsize=12, fontweight='bold')
    
    # Add gridlines
    ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, 
                 alpha=0.5, linestyle='--')
    
    # Plot 2: Transformed Atlantis temperature
    ax2 = fig.add_subplot(gs[1, :], projection=proj)
    
    # Plot Atlantis temperature as polygons using BGM geometry
    plotted_count = 0
    outside_domain_count = 0
    
    for i, box_temp in enumerate(atlantis_temp):
        if i < len(boxes) and len(boxes[i]) > 0:
            if not np.isnan(box_temp):
                # Polygon with temperature data - use color based on temperature
                polygon = patches.Polygon(boxes[i], closed=True, 
                                        facecolor=plt.cm.RdYlBu_r(norm(box_temp)), 
                                        edgecolor='black', linewidth=0.5, alpha=0.8)
                plotted_count += 1
            else:
                # Polygon outside domain - use gray color with same line style
                polygon = patches.Polygon(boxes[i], closed=True, 
                                        facecolor='lightgray', 
                                        edgecolor='black', linewidth=0.5, alpha=0.6)
                outside_domain_count += 1
            
            ax2.add_patch(polygon)
    
    print(f"Plotted {plotted_count} polygons with temperature data")
    print(f"Plotted {outside_domain_count} polygons outside domain (gray)")
    print(f"Total polygons: {plotted_count + outside_domain_count} out of {len(atlantis_temp)}")
    
    # Add coastlines and features
    ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.3)
    ax2.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
    ax2.add_feature(cfeature.LAND, color='lightgray', alpha=0.3)
    
    # Set extent - expanded to show all polygons (based on BGM bounds: 140.08-152.36°E, -41.89 to -36.99°N)
    ax2.set_extent([139.8, 152.6, -42.2, -36.7], crs=proj)
    
    # Add title
    ax2.set_title('Transformed Atlantis Temperature\nTime Step 0, Surface Layer', 
                 fontsize=12, fontweight='bold')
    
    # Add gridlines
    ax2.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, 
                 alpha=0.5, linestyle='--')
    
    # Add shared colorbar
    cbar_ax = fig.add_subplot(gs[2, :])
    cbar = plt.colorbar(im1, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Temperature (°C)', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # Add overall title
    fig.suptitle('Temperature Comparison: BASS-2 vs Atlantis Transformation', 
                fontsize=16, fontweight='bold', y=0.95)
    
    # Add statistics text
    stats_text = f'BASS-2 Range: {temp_min:.2f} - {temp_max:.2f}°C\n' \
                f'Atlantis Range: {np.nanmin(atlantis_temp):.2f} - {np.nanmax(atlantis_temp):.2f}°C\n' \
                f'Atlantis Boxes: {len(atlantis_temp)} | BGM Polygons: {len(boxes)}'
    
    fig.text(0.5, 0.02, stats_text, ha='center', va='bottom', 
             fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, bottom=0.1)
    
    # Save figure
    output_file = 'NESP-Atlantis-Pilot/figures/temperature_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Temperature comparison plot saved: {output_file}")
    
    plt.show()

if __name__ == "__main__":
    create_temperature_comparison()
