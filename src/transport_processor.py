"""
Main transport processing script for the Salish Sea Atlantis Model.

This script processes transport data between model boxes, calculating fluxes
and writing the results to NetCDF files.
"""

import os
import glob
import numpy as np
import xarray as xr
from pathlib import Path
import logging
from typing import List, Tuple
import argparse

from functions.box_reader import read_boxes, read_faces
from functions.transport import process_transport
from functions.netcdf_writer import write_transport_file

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def process_transport_data(bgm_file: str,
                         velocity_dir: str,
                         mesh_file: str,
                         output_file: str,
                         depth_levels: List[float] = None) -> None:
    """
    Process transport data for the Salish Sea model.
    
    Args:
        bgm_file: Path to the BGM file
        velocity_dir: Directory containing velocity data files
        mesh_file: Path to the mesh file
        output_file: Path for output NetCDF file
        depth_levels: List of depth levels (optional)
    """
    # Default depth levels if not provided
    if depth_levels is None:
        depth_levels = [0, 25, 50, 100, 250, 400, 700]
    depth_levels = np.array(depth_levels)
    
    logger.info("Reading box geometry...")
    # Read box geometry
    nbox, nface, bid, cent, b_area, vert, iface, botz = read_boxes(bgm_file)
    
    # Read face relationships
    logger.info("Reading face relationships...")
    nulr, nupt1, nupt2 = read_faces(nbox, nface, bid, vert, iface, bgm_file)
    
    # Get list of velocity files
    velocity_files = sorted(glob.glob(os.path.join(velocity_dir, 'MERGED_*.nc')))
    if not velocity_files:
        raise FileNotFoundError(f"No velocity files found in {velocity_dir}")
    
    logger.info(f"Found {len(velocity_files)} velocity files to process")
    
    # Process each velocity file
    all_transport = []
    all_timestamps = []
    
    for i, velocity_file in enumerate(velocity_files, 1):
        logger.info(f"Processing file {i}/{len(velocity_files)}: {os.path.basename(velocity_file)}")
        
        # Open velocity dataset to get number of timesteps
        with xr.open_dataset(velocity_file) as ds:
            n_timesteps = len(ds.time)
        
        # Process each timestep
        for timestep in range(n_timesteps):
            try:
                transport, timestamp = process_transport(
                    vert=vert,
                    pt1=nupt1,
                    pt2=nupt2,
                    depth_levels=depth_levels,
                    velocity_file=velocity_file,
                    mesh_file=mesh_file,
                    timestep=timestep
                )
                
                all_transport.append(transport)
                all_timestamps.append(timestamp)
                
            except Exception as e:
                logger.error(f"Error processing timestep {timestep} of {os.path.basename(velocity_file)}: {e}")
                continue
    
    if not all_transport:
        raise RuntimeError("No transport data was successfully processed")
    
    # Combine all timesteps
    transport_array = np.stack(all_transport)
    timestamp_array = np.array(all_timestamps)
    
    # Get face IDs (0-based indices of real faces)
    real_faces = np.where(~np.isnan(nupt1[:, 0]))[0]
    
    # Write output file
    logger.info(f"Writing transport data to {output_file}")
    write_transport_file(
        output_file=output_file,
        pt1=nupt1,
        pt2=nupt2,
        layer_connections=nulr,
        timestamps=timestamp_array,
        transport=transport_array,
        face_ids=real_faces
    )
    
    logger.info("Transport processing completed successfully")

def main():
    """Main entry point with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Process transport data for the Salish Sea Atlantis Model"
    )
    
    parser.add_argument(
        '--bgm-file',
        type=str,
        default='20190812_Salish_Sea_ll_fixed.bgm',
        help='Path to the BGM file'
    )
    
    parser.add_argument(
        '--velocity-dir',
        type=str,
        required=True,
        help='Directory containing velocity data files'
    )
    
    parser.add_argument(
        '--mesh-file',
        type=str,
        required=True,
        help='Path to the mesh file'
    )
    
    parser.add_argument(
        '--output-file',
        type=str,
        default='SS_Transport.nc',
        help='Path for output NetCDF file'
    )
    
    parser.add_argument(
        '--depth-levels',
        type=float,
        nargs='+',
        help='List of depth levels (optional)'
    )
    
    args = parser.parse_args()
    
    try:
        process_transport_data(
            bgm_file=args.bgm_file,
            velocity_dir=args.velocity_dir,
            mesh_file=args.mesh_file,
            output_file=args.output_file,
            depth_levels=args.depth_levels
        )
    except Exception as e:
        logger.error(f"Error during transport processing: {e}")
        raise

if __name__ == '__main__':
    main() 