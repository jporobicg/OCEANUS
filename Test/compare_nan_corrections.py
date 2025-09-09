#!/usr/bin/env python3
"""
NaN Correction Analysis
======================

This script analyzes the NaN corrections made during the mass balance process.
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import os


def analyze_nan_corrections():
    """
    Analyze the NaN corrections made during mass balance process.
    """
    print("Analyzing NaN corrections...")
    
    # Load both files
    with nc.Dataset('Transport_test.nc', 'r') as ds1:
        orig = ds1.variables['transport'][:]
        dest = ds1.variables['dest_boxid'][:]
        sour = ds1.variables['source_boxid'][:]
        
    with nc.Dataset('Transport_test_mass_balanced.nc', 'r') as ds2:
        corr = ds2.variables['transport'][:]
    
    print(f"Data shape: {orig.shape}")
    
    # Find faces with NaN corrections
    corrected_faces = []
    correction_details = []
    
    for face_id in range(orig.shape[1]):
        orig_face = orig[:, face_id, :]
        corr_face = corr[:, face_id, :]
        
        # Find where NaN was corrected
        nan_corrected = np.isnan(orig_face) & ~np.isnan(corr_face)
        
        if np.any(nan_corrected):
            corrected_faces.append(face_id)
            
            # Get correction details
            for time_idx in range(orig.shape[0]):
                for level_idx in range(orig.shape[2]):
                    if nan_corrected[time_idx, level_idx]:
                        correction_value = corr_face[time_idx, level_idx]
                        correction_details.append({
                            'face_id': face_id,
                            'time_idx': time_idx,
                            'level_idx': level_idx,
                            'correction_value': correction_value,
                            'box_connection': f"{sour[face_id]}→{dest[face_id]}"
                        })
    
    print(f"Found {len(corrected_faces)} faces with NaN corrections")
    print(f"Total NaN corrections: {len(correction_details)}")
    
    if len(correction_details) == 0:
        print("No NaN corrections found.")
        return
    
    # Create output directory
    os.makedirs("Diagnostics/nan_corrections", exist_ok=True)
    
    # Plot 1: Correction values histogram
    correction_values = [d['correction_value'] for d in correction_details]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(correction_values, bins=30, alpha=0.7, edgecolor='black')
    ax.set_xlabel('Correction Value (m/s)')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of NaN Correction Values')
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    mean_corr = np.mean(np.abs(correction_values))
    max_corr = np.max(np.abs(correction_values))
    
    stats_text = f'Mean correction: {mean_corr:.2e}\nMax correction: {max_corr:.2e}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
           verticalalignment='top', fontsize=12,
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('Diagnostics/nan_corrections/correction_values_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Corrections by face
    face_corrections = {}
    for detail in correction_details:
        face_id = detail['face_id']
        if face_id not in face_corrections:
            face_corrections[face_id] = []
        face_corrections[face_id].append(detail['correction_value'])
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    face_ids = list(face_corrections.keys())
    correction_counts = [len(face_corrections[fid]) for fid in face_ids]
    
    bars = ax.bar(range(len(face_ids)), correction_counts)
    ax.set_xlabel('Face ID')
    ax.set_ylabel('Number of NaN Corrections')
    ax.set_title('Number of NaN Corrections by Face')
    ax.set_xticks(range(len(face_ids)))
    ax.set_xticklabels([f'F{fid}' for fid in face_ids], rotation=45)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Diagnostics/nan_corrections/corrections_by_face.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create detailed report
    report_file = 'Diagnostics/nan_corrections/nan_correction_report.txt'
    
    with open(report_file, 'w') as f:
        f.write("NaN Correction Analysis Report\n")
        f.write("=" * 40 + "\n\n")
        
        f.write(f"Total faces with corrections: {len(corrected_faces)}\n")
        f.write(f"Total NaN corrections: {len(correction_details)}\n")
        f.write(f"Original NaN count: {np.sum(np.isnan(orig))}\n")
        f.write(f"Final NaN count: {np.sum(np.isnan(corr))}\n")
        f.write(f"NaN reduction: {np.sum(np.isnan(orig)) - np.sum(np.isnan(corr))}\n\n")
        
        f.write("CORRECTED FACES DETAILS\n")
        f.write("-" * 30 + "\n")
        
        for face_id in corrected_faces:
            face_details = [d for d in correction_details if d['face_id'] == face_id]
            f.write(f"\nFace {face_id} ({face_details[0]['box_connection']}):\n")
            f.write(f"  Corrections: {len(face_details)}\n")
            
            correction_vals = [d['correction_value'] for d in face_details]
            f.write(f"  Mean correction: {np.mean(correction_vals):.2e}\n")
            f.write(f"  Max correction: {np.max(np.abs(correction_vals)):.2e}\n")
            
            # Show time/level distribution
            time_indices = [d['time_idx'] for d in face_details]
            level_indices = [d['level_idx'] for d in face_details]
            f.write(f"  Time range: {min(time_indices)} to {max(time_indices)}\n")
            f.write(f"  Levels: {sorted(set(level_indices))}\n")
        
        f.write(f"\nOVERALL STATISTICS\n")
        f.write("-" * 20 + "\n")
        f.write(f"Mean correction magnitude: {np.mean(np.abs(correction_values)):.2e}\n")
        f.write(f"Max correction magnitude: {np.max(np.abs(correction_values)):.2e}\n")
        f.write(f"Standard deviation: {np.std(correction_values):.2e}\n")
        
        # Time and level distribution
        time_indices = [d['time_idx'] for d in correction_details]
        level_indices = [d['level_idx'] for d in correction_details]
        
        f.write(f"\nTime distribution:\n")
        for t in sorted(set(time_indices)):
            count = time_indices.count(t)
            f.write(f"  Time {t}: {count} corrections\n")
        
        f.write(f"\nLevel distribution:\n")
        for l in sorted(set(level_indices)):
            count = level_indices.count(l)
            f.write(f"  Level {l}: {count} corrections\n")
    
    print(f"Analysis complete. Results saved in Diagnostics/nan_corrections/")
    print(f"Report: {report_file}")


if __name__ == "__main__":
    analyze_nan_corrections()
