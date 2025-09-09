#!/usr/bin/env python3
"""
Transport Correction Comparison Plot
===================================

This script compares the original and mass-balanced transport outputs,
focusing on faces that were corrected during the mass balance process.

Author: OCEANUS Development Team
Date: 2024
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import os
from datetime import datetime


def load_transport_data(transport_file):
    """
    Load transport data from NetCDF file.
    
    Parameters:
    -----------
    transport_file : str
        Path to the transport NetCDF file
        
    Returns:
    --------
    dict
        Dictionary containing transport data
    """
    print(f"Loading transport data from: {transport_file}")
    
    with nc.Dataset(transport_file, 'r') as ds:
        transport = ds.variables['transport'][:]  # (time, faces, level)
        # Try different possible variable names
        try:
            dest = ds.variables['dest_boxid'][:]
            sour = ds.variables['source_boxid'][:]
        except KeyError:
            try:
                dest = ds.variables['dest'][:]
                sour = ds.variables['sour'][:]
            except KeyError:
                print(f"Error: Could not find dest/sour variables in {transport_file}")
                raise
        time = ds.variables['time'][:]
        level = ds.variables['level'][:]
        
        return {
            'transport': transport,
            'dest': dest,
            'sour': sour,
            'time': time,
            'level': level
        }


def find_corrected_faces(original_data, corrected_data, tolerance=1e-6):
    """
    Find faces that were corrected during mass balance process.
    
    Parameters:
    -----------
    original_data : dict
        Original transport data
    corrected_data : dict
        Corrected transport data
    tolerance : float
        Tolerance for considering values different
        
    Returns:
    --------
    list
        List of face IDs that were corrected
    """
    print("Finding corrected faces...")
    
    original_transport = original_data['transport']
    corrected_transport = corrected_data['transport']
    
    corrected_faces = []
    
    for face_id in range(original_transport.shape[1]):
        original_face = original_transport[:, face_id, :]
        corrected_face = corrected_transport[:, face_id, :]
        
        # Check if any values changed significantly
        if np.any(np.abs(original_face - corrected_face) > tolerance):
            corrected_faces.append(face_id)
    
    print(f"Found {len(corrected_faces)} corrected faces")
    return corrected_faces


def plot_face_comparison(original_data, corrected_data, corrected_faces, output_dir):
    """
    Create comparison plots for corrected faces.
    
    Parameters:
    -----------
    original_data : dict
        Original transport data
    corrected_data : dict
        Corrected transport data
    corrected_faces : list
        List of corrected face IDs
    output_dir : str
        Output directory for plots
    """
    print(f"Creating comparison plots for {len(corrected_faces)} corrected faces...")
    
    original_transport = original_data['transport']
    corrected_transport = corrected_data['transport']
    time = original_data['time']
    level = original_data['level']
    dest = original_data['dest']
    sour = original_data['sour']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot 1: Time series comparison for each corrected face
    n_faces = len(corrected_faces)
    n_levels = len(level)
    
    # Create subplots for each corrected face
    fig, axes = plt.subplots(n_faces, n_levels, figsize=(4*n_levels, 3*n_faces))
    if n_faces == 1:
        axes = axes.reshape(1, -1)
    
    for i, face_id in enumerate(corrected_faces):
        for j in range(n_levels):
            ax = axes[i, j]
            
            # Get time series for this face and level
            original_ts = original_transport[:, face_id, j]
            corrected_ts = corrected_transport[:, face_id, j]
            
            # Plot time series
            ax.plot(time, original_ts, 'b-', alpha=0.7, label='Original', linewidth=1)
            ax.plot(time, corrected_ts, 'r-', alpha=0.7, label='Corrected', linewidth=1)
            
            # Calculate statistics
            original_mean = np.nanmean(original_ts)
            corrected_mean = np.nanmean(corrected_ts)
            original_std = np.nanstd(original_ts)
            corrected_std = np.nanstd(corrected_ts)
            
            # Add statistics text
            stats_text = f'Face {face_id}\nBox {sour[face_id]}→{dest[face_id]}\nLevel {level[j]}\n'
            stats_text += f'Original: {original_mean:.2e}±{original_std:.2e}\n'
            stats_text += f'Corrected: {corrected_mean:.2e}±{corrected_std:.2e}'
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=8, 
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            ax.set_xlabel('Time')
            ax.set_ylabel('Transport (m/s)')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            
            # Rotate x-axis labels for better readability
            ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/corrected_faces_timeseries.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Time series comparison saved: {output_dir}/corrected_faces_timeseries.png")
    
    # Plot 2: Scatter plot of original vs corrected values
    fig, ax = plt.subplots(figsize=(10, 8))
    
    all_original = []
    all_corrected = []
    
    for face_id in corrected_faces:
        original_face = original_transport[:, face_id, :].flatten()
        corrected_face = corrected_transport[:, face_id, :].flatten()
        
        # Remove NaN values
        valid_mask = ~(np.isnan(original_face) | np.isnan(corrected_face))
        if np.any(valid_mask):
            all_original.extend(original_face[valid_mask])
            all_corrected.extend(corrected_face[valid_mask])
    
    all_original = np.array(all_original)
    all_corrected = np.array(all_corrected)
    
    # Create scatter plot
    ax.scatter(all_original, all_corrected, alpha=0.6, s=20)
    
    # Add 1:1 line
    min_val = min(np.nanmin(all_original), np.nanmin(all_corrected))
    max_val = max(np.nanmax(all_original), np.nanmax(all_corrected))
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='1:1 line')
    
    ax.set_xlabel('Original Transport (m/s)')
    ax.set_ylabel('Corrected Transport (m/s)')
    ax.set_title('Original vs Corrected Transport Values\n(Corrected Faces Only)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    correlation = np.corrcoef(all_original, all_corrected)[0, 1]
    rmse = np.sqrt(np.mean((all_original - all_corrected)**2))
    
    stats_text = f'Correlation: {correlation:.3f}\nRMSE: {rmse:.2e}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
           verticalalignment='top', fontsize=12,
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/original_vs_corrected_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Scatter plot saved: {output_dir}/original_vs_corrected_scatter.png")
    
    # Plot 3: Correction magnitude histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    
    corrections = all_corrected - all_original
    
    # Create histogram of correction magnitudes
    ax.hist(np.abs(corrections), bins=50, alpha=0.7, edgecolor='black')
    ax.set_xlabel('Absolute Correction Magnitude (m/s)')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Correction Magnitudes\n(Corrected Faces Only)')
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    mean_correction = np.mean(np.abs(corrections))
    max_correction = np.max(np.abs(corrections))
    
    stats_text = f'Mean correction: {mean_correction:.2e}\nMax correction: {max_correction:.2e}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
           verticalalignment='top', fontsize=12,
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/correction_magnitude_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Correction histogram saved: {output_dir}/correction_magnitude_histogram.png")


def create_summary_report(original_data, corrected_data, corrected_faces, output_dir):
    """
    Create a summary report of the corrections.
    
    Parameters:
    -----------
    original_data : dict
        Original transport data
    corrected_data : dict
        Corrected transport data
    corrected_faces : list
        List of corrected face IDs
    output_dir : str
        Output directory for report
    """
    print("Creating summary report...")
    
    original_transport = original_data['transport']
    corrected_transport = corrected_data['transport']
    dest = original_data['dest']
    sour = original_data['sour']
    level = original_data['level']
    
    report_file = f'{output_dir}/correction_summary.txt'
    
    with open(report_file, 'w') as f:
        f.write("Transport Correction Summary Report\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"Total corrected faces: {len(corrected_faces)}\n")
        f.write(f"Total faces in dataset: {original_transport.shape[1]}\n")
        f.write(f"Correction percentage: {(len(corrected_faces)/original_transport.shape[1])*100:.1f}%\n\n")
        
        f.write("CORRECTED FACES DETAILS\n")
        f.write("-" * 30 + "\n")
        
        for face_id in corrected_faces:
            f.write(f"\nFace {face_id} (Box {sour[face_id]} → Box {dest[face_id]}):\n")
            
            original_face = original_transport[:, face_id, :]
            corrected_face = corrected_transport[:, face_id, :]
            
            for level_idx in range(len(level)):
                original_level = original_face[:, level_idx]
                corrected_level = corrected_face[:, level_idx]
                
                # Calculate statistics
                original_mean = np.nanmean(original_level)
                corrected_mean = np.nanmean(corrected_level)
                correction = corrected_mean - original_mean
                
                f.write(f"  Level {level[level_idx]}: ")
                f.write(f"Original={original_mean:.2e}, ")
                f.write(f"Corrected={corrected_mean:.2e}, ")
                f.write(f"Correction={correction:.2e}\n")
        
        # Overall statistics
        f.write(f"\nOVERALL STATISTICS\n")
        f.write("-" * 20 + "\n")
        
        all_original = original_transport[:, corrected_faces, :].flatten()
        all_corrected = corrected_transport[:, corrected_faces, :].flatten()
        
        # Remove NaN values
        valid_mask = ~(np.isnan(all_original) | np.isnan(all_corrected))
        all_original = all_original[valid_mask]
        all_corrected = all_corrected[valid_mask]
        
        corrections = all_corrected - all_original
        
        f.write(f"Total corrected values: {len(corrections)}\n")
        f.write(f"Mean correction magnitude: {np.mean(np.abs(corrections)):.2e}\n")
        f.write(f"Max correction magnitude: {np.max(np.abs(corrections)):.2e}\n")
        f.write(f"Standard deviation of corrections: {np.std(corrections):.2e}\n")
        f.write(f"Correlation (original vs corrected): {np.corrcoef(all_original, all_corrected)[0, 1]:.3f}\n")
        
        print(f"Summary report saved: {report_file}")


def main():
    """
    Main function to run the transport comparison analysis.
    """
    print("Transport Correction Comparison Analysis")
    print("=" * 45)
    
    # File paths
    original_file = "Transport_test.nc"
    corrected_file = "Transport_test_mass_balanced.nc"
    output_dir = "Diagnostics/transport_comparison"
    
    # Check if files exist
    if not os.path.exists(original_file):
        print(f"Error: Original transport file not found: {original_file}")
        return
    
    if not os.path.exists(corrected_file):
        print(f"Error: Corrected transport file not found: {corrected_file}")
        return
    
    # Load transport data
    print("Loading transport data...")
    original_data = load_transport_data(original_file)
    corrected_data = load_transport_data(corrected_file)
    
    # Find corrected faces
    corrected_faces = find_corrected_faces(original_data, corrected_data)
    
    if len(corrected_faces) == 0:
        print("No corrected faces found. No comparison plots will be created.")
        return
    
    # Create comparison plots
    plot_face_comparison(original_data, corrected_data, corrected_faces, output_dir)
    
    # Create summary report
    create_summary_report(original_data, corrected_data, corrected_faces, output_dir)
    
    print("\nTransport comparison analysis completed!")
    print(f"Output files saved in: {output_dir}")


if __name__ == "__main__":
    main()
