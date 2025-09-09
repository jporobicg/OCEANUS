#!/usr/bin/env python3
"""
Transport Data Visualization Script
Plots transport data from OCEANUS output to assess model performance
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd
from datetime import datetime, timedelta

def load_transport_data(filename):
    """Load transport data from NetCDF file"""
    print(f"Loading transport data from: {filename}")
    
    with nc.Dataset(filename, 'r') as ds:
        # Load basic dimensions
        time = ds.variables['time'][:]
        level = ds.variables['level'][:]
        faces = ds.variables['faces'][:]
        
        # Load transport data
        transport = ds.variables['transport'][:]  # (level, faces, time)
        
        # Load face coordinates
        pt1_x = ds.variables['pt1_x'][:]
        pt1_y = ds.variables['pt1_y'][:]
        pt2_x = ds.variables['pt2_x'][:]
        pt2_y = ds.variables['pt2_y'][:]
        
        # Load box relationships
        dest_boxid = ds.variables['dest_boxid'][:]
        source_boxid = ds.variables['source_boxid'][:]
        
        print(f"Data loaded:")
        print(f"  Time steps: {len(time)}")
        print(f"  Levels: {len(level)}")
        print(f"  Faces: {len(faces)}")
        print(f"  Transport shape: {transport.shape}")
        
        return {
            'time': time,
            'level': level,
            'faces': faces,
            'transport': transport,
            'pt1_x': pt1_x,
            'pt1_y': pt1_y,
            'pt2_x': pt2_x,
            'pt2_y': pt2_y,
            'dest_boxid': dest_boxid,
            'source_boxid': source_boxid
        }

def plot_transport_time_series(data, output_prefix="transport"):
    """Plot transport time series for different levels"""
    print("Creating time series plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, level_idx in enumerate(range(data['transport'].shape[2])):  # Changed from shape[0] to shape[2]
        ax = axes[i]
        
        # Calculate mean transport across all faces for this level
        mean_transport = np.nanmean(data['transport'][:, :, level_idx], axis=1)  # Changed indexing
        std_transport = np.nanstd(data['transport'][:, :, level_idx], axis=1)    # Changed indexing
        
        # Convert time to days (assuming seconds since epoch)
        time_days = data['time'] / (24 * 3600)
        
        ax.plot(time_days, mean_transport, 'b-', linewidth=2, label='Mean')
        ax.fill_between(time_days, 
                       mean_transport - std_transport, 
                       mean_transport + std_transport, 
                       alpha=0.3, color='blue', label='±1σ')
        
        ax.set_title(f'Level {data["level"][level_idx]} (Depth Layer {level_idx+1})')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Transport (m³/s)')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add statistics
        stats_text = f'Mean: {np.nanmean(mean_transport):.2f}\nStd: {np.nanstd(mean_transport):.2f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_time_series.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_transport_spatial_distribution(data, output_prefix="transport"):
    """Plot spatial distribution of transport"""
    print("Creating spatial distribution plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, level_idx in enumerate(range(data['transport'].shape[2])):  # Changed from shape[0] to shape[2]
        ax = axes[i]
        
        # Calculate mean transport across time for this level
        mean_transport = np.nanmean(data['transport'][:, :, level_idx], axis=0)  # Changed indexing
        
        # Create scatter plot of face centers with transport magnitude
        face_centers_x = (data['pt1_x'] + data['pt2_x']) / 2
        face_centers_y = (data['pt1_y'] + data['pt2_y']) / 2
        
        # Color by transport magnitude
        scatter = ax.scatter(face_centers_x, face_centers_y, 
                           c=mean_transport, cmap='RdBu_r', 
                           s=20, alpha=0.7, edgecolors='black', linewidth=0.5)
        
        ax.set_title(f'Level {data["level"][level_idx]} - Mean Transport')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Transport (m³/s)')
        
        # Add statistics
        stats_text = f'Mean: {np.nanmean(mean_transport):.2f}\nMax: {np.nanmax(mean_transport):.2f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_spatial_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_transport_statistics(data, output_prefix="transport"):
    """Plot statistical analysis of transport data"""
    print("Creating statistical analysis plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Transport distribution by level
    ax1 = axes[0, 0]
    for i in range(data['transport'].shape[2]):  # Changed from shape[0] to shape[2]
        level_data = data['transport'][:, :, i].flatten()  # Changed indexing
        level_data = level_data[~np.isnan(level_data)]
        ax1.hist(level_data, bins=50, alpha=0.6, label=f'Level {data["level"][i]}')
    
    ax1.set_xlabel('Transport (m³/s)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Transport Distribution by Level')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Box plot of transport by level
    ax2 = axes[0, 1]
    level_data_list = []
    level_labels = []
    for i in range(data['transport'].shape[2]):  # Changed from shape[0] to shape[2]
        level_data = data['transport'][:, :, i].flatten()  # Changed indexing
        level_data = level_data[~np.isnan(level_data)]
        level_data_list.append(level_data)
        level_labels.append(f'Level {data["level"][i]}')
    
    ax2.boxplot(level_data_list, labels=level_labels)
    ax2.set_ylabel('Transport (m³/s)')
    ax2.set_title('Transport Distribution by Level (Box Plot)')
    ax2.grid(True, alpha=0.3)
    
    # 3. Time series of total transport
    ax3 = axes[1, 0]
    total_transport = np.nansum(data['transport'], axis=(1, 2))  # Changed from (0, 1) to (1, 2)
    time_days = data['time'] / (24 * 3600)
    
    ax3.plot(time_days, total_transport, 'g-', linewidth=2)
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Total Transport (m³/s)')
    ax3.set_title('Total Transport Over Time')
    ax3.grid(True, alpha=0.3)
    
    # 4. Transport magnitude vs face index
    ax4 = axes[1, 1]
    mean_transport_by_face = np.nanmean(data['transport'], axis=(0, 2))  # Changed from (0, 2) to (0, 2) - this one is correct
    face_indices = np.arange(len(mean_transport_by_face))
    
    ax4.scatter(face_indices, mean_transport_by_face, alpha=0.6, s=10)
    ax4.set_xlabel('Face Index')
    ax4.set_ylabel('Mean Transport (m³/s)')
    ax4.set_title('Mean Transport by Face')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_statistics.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_transport_heatmap(data, output_prefix="transport"):
    """Plot transport heatmap (faces vs time)"""
    print("Creating transport heatmap...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, level_idx in enumerate(range(data['transport'].shape[2])):  # Changed from shape[0] to shape[2]
        ax = axes[i]
        
        # Get transport data for this level
        transport_level = data['transport'][:, :, level_idx]  # Changed indexing
        
        # Create heatmap
        im = ax.imshow(transport_level, aspect='auto', cmap='RdBu_r', 
                      extent=[0, len(data['faces']), 0, len(data['time'])])
        
        ax.set_title(f'Level {data["level"][level_idx]} - Transport Heatmap')
        ax.set_xlabel('Face Index')
        ax.set_ylabel('Time Step')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Transport (m³/s)')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()

def generate_summary_report(data, output_prefix="transport"):
    """Generate a summary report of transport data"""
    print("Generating summary report...")
    
    report = []
    report.append("OCEANUS Transport Data Summary Report")
    report.append("=" * 50)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    # Basic statistics
    report.append("Basic Statistics:")
    report.append(f"  Total time steps: {len(data['time'])}")
    report.append(f"  Total levels: {len(data['level'])}")
    report.append(f"  Total faces: {len(data['faces'])}")
    report.append(f"  Time range: {data['time'][0]:.1f} to {data['time'][-1]:.1f} seconds")
    report.append("")
    
    # Transport statistics by level
    report.append("Transport Statistics by Level:")
    for i in range(data['transport'].shape[0]):
        level_data = data['transport'][i, :, :].flatten()
        level_data = level_data[~np.isnan(level_data)]
        
        report.append(f"  Level {data['level'][i]}:")
        report.append(f"    Mean: {np.nanmean(level_data):.4f} m³/s")
        report.append(f"    Std:  {np.nanstd(level_data):.4f} m³/s")
        report.append(f"    Min:  {np.nanmin(level_data):.4f} m³/s")
        report.append(f"    Max:  {np.nanmax(level_data):.4f} m³/s")
        report.append(f"    Non-zero: {np.sum(level_data != 0)}/{len(level_data)} ({100*np.sum(level_data != 0)/len(level_data):.1f}%)")
        report.append("")
    
    # Overall statistics
    total_transport = np.nansum(data['transport'], axis=(0, 1))
    report.append("Overall Statistics:")
    report.append(f"  Total transport mean: {np.nanmean(total_transport):.4f} m³/s")
    report.append(f"  Total transport std:  {np.nanstd(total_transport):.4f} m³/s")
    report.append(f"  Total transport min:  {np.nanmin(total_transport):.4f} m³/s")
    report.append(f"  Total transport max:  {np.nanmax(total_transport):.4f} m³/s")
    
    # Write report to file
    with open(f'{output_prefix}_summary_report.txt', 'w') as f:
        f.write('\n'.join(report))
    
    print('\n'.join(report))
    print(f"\nSummary report saved to: {output_prefix}_summary_report.txt")

def main():
    """Main function to run all plotting routines"""
    print("OCEANUS Transport Data Visualization")
    print("=" * 40)
    
    # Load data
    try:
        data = load_transport_data('Transport_test.nc')
    except FileNotFoundError:
        print("Error: Transport_test.nc not found in current directory")
        print("Please ensure the transport file is in the Test directory")
        return
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    # Set up plotting style
    plt.style.use('default')
    
    # Generate all plots
    try:
        plot_transport_time_series(data)
        plot_transport_spatial_distribution(data)
        plot_transport_statistics(data)
        plot_transport_heatmap(data)
        generate_summary_report(data)
        
        print("\nAll plots generated successfully!")
        print("Files created:")
        print("  - transport_time_series.png")
        print("  - transport_spatial_distribution.png")
        print("  - transport_statistics.png")
        print("  - transport_heatmap.png")
        print("  - transport_summary_report.txt")
        
    except Exception as e:
        print(f"Error generating plots: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 