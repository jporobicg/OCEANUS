#!/usr/bin/env python3
"""
Transport Correction Summary
===========================

This script creates a summary visualization showing the key differences
between the original and mass-balanced transport files.
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import os


def create_correction_summary():
    """
    Create a summary visualization of the transport corrections.
    """
    print("Creating transport correction summary...")
    
    # Load both files
    with nc.Dataset('Transport_test.nc', 'r') as ds1:
        orig = ds1.variables['transport'][:]
        dest = ds1.variables['dest_boxid'][:]
        sour = ds1.variables['source_boxid'][:]
        
    with nc.Dataset('Transport_test_mass_balanced.nc', 'r') as ds2:
        corr = ds2.variables['transport'][:]
    
    # Create output directory
    os.makedirs("Diagnostics/correction_summary", exist_ok=True)
    
    # Create summary figure
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: NaN reduction
    ax1 = axes[0, 0]
    
    orig_nan_count = np.sum(np.isnan(orig), axis=(0, 2))  # NaN count per face
    corr_nan_count = np.sum(np.isnan(corr), axis=(0, 2))
    
    faces_with_nan = np.where(orig_nan_count > 0)[0]
    
    if len(faces_with_nan) > 0:
        x_pos = np.arange(len(faces_with_nan))
        width = 0.35
        
        ax1.bar(x_pos - width/2, orig_nan_count[faces_with_nan], width, 
               label='Original', alpha=0.7, color='red')
        ax1.bar(x_pos + width/2, corr_nan_count[faces_with_nan], width, 
               label='Corrected', alpha=0.7, color='green')
        
        ax1.set_xlabel('Face ID')
        ax1.set_ylabel('Number of NaN Values')
        ax1.set_title('NaN Values Before and After Correction')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels([f'F{fid}' for fid in faces_with_nan], rotation=45)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    
    # Plot 2: Correction magnitude by face
    ax2 = axes[0, 1]
    
    # Find faces with corrections
    corrected_faces = []
    correction_magnitudes = []
    
    for face_id in range(orig.shape[1]):
        orig_face = orig[:, face_id, :]
        corr_face = corr[:, face_id, :]
        
        # Find where NaN was corrected
        nan_corrected = np.isnan(orig_face) & ~np.isnan(corr_face)
        
        if np.any(nan_corrected):
            corrected_faces.append(face_id)
            correction_vals = corr_face[nan_corrected]
            correction_magnitudes.append(np.mean(np.abs(correction_vals)))
    
    if corrected_faces:
        ax2.bar(range(len(corrected_faces)), correction_magnitudes, alpha=0.7, color='blue')
        ax2.set_xlabel('Face ID')
        ax2.set_ylabel('Mean Correction Magnitude (m/s)')
        ax2.set_title('Mean Correction Magnitude by Face')
        ax2.set_xticks(range(len(corrected_faces)))
        ax2.set_xticklabels([f'F{fid}' for fid in corrected_faces], rotation=45)
        ax2.grid(True, alpha=0.3)
        ax2.set_yscale('log')
    
    # Plot 3: Time series of total NaN count
    ax3 = axes[1, 0]
    
    orig_nan_time = np.sum(np.isnan(orig), axis=(1, 2))  # NaN count per time step
    corr_nan_time = np.sum(np.isnan(corr), axis=(1, 2))
    
    time_steps = np.arange(len(orig_nan_time))
    
    ax3.plot(time_steps, orig_nan_time, 'r-', label='Original', linewidth=2)
    ax3.plot(time_steps, corr_nan_time, 'g-', label='Corrected', linewidth=2)
    ax3.fill_between(time_steps, orig_nan_time, corr_nan_time, alpha=0.3, color='blue')
    
    ax3.set_xlabel('Time Step')
    ax3.set_ylabel('Total NaN Count')
    ax3.set_title('NaN Count Over Time')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Level-wise NaN reduction
    ax4 = axes[1, 1]
    
    orig_nan_level = np.sum(np.isnan(orig), axis=(0, 1))  # NaN count per level
    corr_nan_level = np.sum(np.isnan(corr), axis=(0, 1))
    
    levels = np.arange(len(orig_nan_level))
    x_pos = np.arange(len(levels))
    width = 0.35
    
    ax4.bar(x_pos - width/2, orig_nan_level, width, label='Original', alpha=0.7, color='red')
    ax4.bar(x_pos + width/2, corr_nan_level, width, label='Corrected', alpha=0.7, color='green')
    
    ax4.set_xlabel('Level')
    ax4.set_ylabel('Number of NaN Values')
    ax4.set_title('NaN Values by Level')
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels([f'L{l}' for l in levels])
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Diagnostics/correction_summary/transport_correction_summary.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create summary statistics
    summary_file = 'Diagnostics/correction_summary/correction_summary.txt'
    
    with open(summary_file, 'w') as f:
        f.write("Transport Correction Summary\n")
        f.write("=" * 35 + "\n\n")
        
        f.write("OVERALL STATISTICS\n")
        f.write("-" * 20 + "\n")
        f.write(f"Original NaN count: {np.sum(np.isnan(orig))}\n")
        f.write(f"Corrected NaN count: {np.sum(np.isnan(corr))}\n")
        f.write(f"NaN reduction: {np.sum(np.isnan(orig)) - np.sum(np.isnan(corr))}\n")
        f.write(f"Reduction percentage: {(np.sum(np.isnan(orig)) - np.sum(np.isnan(corr))) / np.sum(np.isnan(orig)) * 100:.1f}%\n\n")
        
        f.write("CORRECTED FACES\n")
        f.write("-" * 15 + "\n")
        if corrected_faces:
            for i, face_id in enumerate(corrected_faces):
                f.write(f"Face {face_id}: {sour[face_id]} → {dest[face_id]}, "
                       f"Mean correction: {correction_magnitudes[i]:.2e} m/s\n")
        else:
            f.write("No faces were corrected.\n")
        
        f.write(f"\nTIME SERIES ANALYSIS\n")
        f.write("-" * 22 + "\n")
        f.write(f"Original NaN range: {np.min(orig_nan_time)} - {np.max(orig_nan_time)}\n")
        f.write(f"Corrected NaN range: {np.min(corr_nan_time)} - {np.max(corr_nan_time)}\n")
        f.write(f"Average NaN reduction per time step: {np.mean(orig_nan_time - corr_nan_time):.1f}\n")
        
        f.write(f"\nLEVEL ANALYSIS\n")
        f.write("-" * 15 + "\n")
        for level_idx in range(len(levels)):
            reduction = orig_nan_level[level_idx] - corr_nan_level[level_idx]
            f.write(f"Level {level_idx}: {orig_nan_level[level_idx]} → {corr_nan_level[level_idx]} "
                   f"(reduction: {reduction})\n")
    
    print(f"Summary created: Diagnostics/correction_summary/")
    print(f"Summary file: {summary_file}")


if __name__ == "__main__":
    create_correction_summary()
