"""
Main variables processing script for the Salish Sea Atlantis Model.

This script processes oceanographic and biological variables, calculating
box averages and writing the results to NetCDF files.
"""

import os
import glob
import numpy as np
import xarray as xr
from pathlib import Path
import logging
from typing import List, Dict, Optional
import argparse

from functions.box_reader import read_boxes, read_faces
from functions.box_operations import process_box_averages
from functions.netcdf_writer import write_variable_file

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def process_variables(bgm_file: str,
                     variable_files: Dict[str, str],
                     mesh_file: str,
                     output_file: str,
                     depth_levels: List[float] = None) -> None:
    """
    Process variables for the Salish Sea model.
    
    Args:
        bgm_file: Path to the BGM file
        variable_files: Dictionary mapping variable names to file paths
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
    
    # Initialize storage for all variables
    all_variables = {}
    all_timestamps = []
    
    # Process each variable
    for var_name, var_file in variable_files.items():
        logger.info(f"Processing variable: {var_name}")
        
        try:
            # Open variable dataset to get number of timesteps
            with xr.open_dataset(var_file) as ds:
                n_timesteps = len(ds.time)
            
            # Initialize storage for this variable
            var_data = []
            var_timestamps = []
            
            # Process each timestep
            for timestep in range(n_timesteps):
                try:
                    averages, timestamp = process_box_averages(
                        vertices=vert,
                        depth_levels=depth_levels,
                        variable_file=var_file,
                        variable_name=var_name,
                        mesh_file=mesh_file,
                        timestep=timestep
                    )
                    
                    var_data.append(averages)
                    var_timestamps.append(timestamp)
                    
                except Exception as e:
                    logger.error(f"Error processing timestep {timestep} for {var_name}: {e}")
                    continue
            
            if var_data:
                # Stack all timesteps for this variable
                all_variables[var_name] = np.stack(var_data)
                
                # Store timestamps (only need to do this once)
                if not all_timestamps:
                    all_timestamps = np.array(var_timestamps)
                    
            else:
                logger.warning(f"No data was successfully processed for {var_name}")
                
        except Exception as e:
            logger.error(f"Error processing variable {var_name}: {e}")
            continue
    
    if not all_variables:
        raise RuntimeError("No variables were successfully processed")
    
    # Write output file
    logger.info(f"Writing variable data to {output_file}")
    write_variable_file(
        output_file=output_file,
        timestamps=all_timestamps,
        box_ids=bid,
        variables=all_variables
    )
    
    logger.info("Variable processing completed successfully")

def main():
    """Main entry point with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Process variables for the Salish Sea Atlantis Model"
    )
    
    parser.add_argument(
        '--bgm-file',
        type=str,
        default='20190812_Salish_Sea_ll_fixed.bgm',
        help='Path to the BGM file'
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
        default='SS_Variables.nc',
        help='Path for output NetCDF file'
    )
    
    parser.add_argument(
        '--depth-levels',
        type=float,
        nargs='+',
        help='List of depth levels (optional)'
    )
    
    parser.add_argument(
        '--variable-files',
        type=str,
        nargs='+',
        required=True,
        help='Space-separated list of variable files in format: name:path'
    )
    
    args = parser.parse_args()
    
    # Parse variable files argument
    variable_files = {}
    for var_file in args.variable_files:
        try:
            name, path = var_file.split(':')
            variable_files[name] = path
        except ValueError:
            parser.error(f"Invalid variable file format: {var_file}. Use name:path format.")
    
    try:
        process_variables(
            bgm_file=args.bgm_file,
            variable_files=variable_files,
            mesh_file=args.mesh_file,
            output_file=args.output_file,
            depth_levels=args.depth_levels
        )
    except Exception as e:
        logger.error(f"Error during variable processing: {e}")
        raise

if __name__ == '__main__':
    main() 