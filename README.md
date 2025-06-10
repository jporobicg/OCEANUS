# OCEANUS

**O**ceanographic **C**onverter & **E**xtractor for **A**tla**N**tis **U**sing **S**imulations

```
    ╔════════════════════════════════════════════╗
    ║               O C E A N U S               ║
    ║                                            ║
    ║    Oceanographic Converter & Extractor     ║
    ║     for Atlantis Using Simulations         ║
    ╚════════════════════════════════════════════╝
```

A powerful Python tool for transforming oceanographic data into Atlantis ecosystem model inputs. OCEANUS streamlines the preparation of hydrodynamic data, ensuring seamless integration with Atlantis simulations.

## Overview
OCEANUS masters the flow of oceanographic data into Atlantis by:
- Converting complex hydrodynamic data into Atlantis-compatible formats
- Processing transport calculations between model boxes
- Handling variable averaging within boxes
- Managing vertical layer distributions
- Generating precise NetCDF outputs

## Features
- Transport data processing from NetCDF files
- Variable data processing and averaging
- Box model geometry handling
- NetCDF output generation

## Installation
1. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Linux/Mac
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Project Structure
- `src/transport_processor.py`: Main transport processing script
- `src/variables_processor.py`: Variables processing script
- `src/functions/`: Helper functions for data processing

## Model Workflow

```
                                    [Input Files]
                                         │
                    ┌───────────────────┴───────────────────┐
                    │                                       │
            [Transport Processing]               [Variables Processing]
                    │                                       │
                    ▼                                       ▼
        ┌─────────────────────┐               ┌─────────────────────┐
        │ transport_processor │               │ variables_processor  │
        └─────────────────────┘               └─────────────────────┘
                    │                                       │
            ┌───────┴───────┐                       ┌──────┴──────┐
            ▼               ▼                       ▼             ▼
    [Read Geometry]  [Process Transport]    [Read Geometry] [Process Variables]
         │                   │                    │               │
         ▼                   ▼                    ▼               ▼
   ┌──────────┐    ┌────────────────┐    ┌──────────┐    ┌─────────────┐
   │box_reader│    │   transport    │    │box_reader│    │box_operations│
   └──────────┘    └────────────────┘    └──────────┘    └─────────────┘
         │                   │                    │               │
         └───────┐    ┌─────┘                    └───────┐    ┌─┘
                 ▼    ▼                                   ▼    ▼
            [Transport NetCDF]                      [Variables NetCDF]
                 │                                         │
                 └─────────────────┐       ┌──────────────┘
                                   ▼       ▼
                               [Output Files]

```

## Function Dependencies

### Transport Processing Chain:
```
transport_processor.py
├── box_reader.py
│   └── coordinate_utils.py
├── transport.py
│   ├── coordinate_utils.py
│   └── box_reader.py
└── netcdf_writer.py

```

### Variables Processing Chain:
```
variables_processor.py
├── box_reader.py
│   └── coordinate_utils.py
├── box_operations.py
│   └── coordinate_utils.py
└── netcdf_writer.py
```

## Core Functions Description

### Geometry Processing
- `read_boxes()`: Reads box definitions from BGM file
- `read_faces()`: Processes face relationships between boxes
- `check_face_connection()`: Validates face connections between boxes

### Transport Calculations
- `process_transport()`: Main transport processing function
- `calculate_transport()`: Computes transport between boxes
- `calculate_face_transport()`: Handles individual face calculations

### Variable Processing
- `process_variables()`: Main variable processing function
- `calculate_box_averages()`: Computes spatial averages
- `process_box_averages()`: Handles box-level calculations

### Coordinate Operations
- `lat_correction()`: Applies latitude-dependent corrections
- `distance_between_points()`: Calculates distances
- `interpolate_at_latitude()`: Handles latitude interpolation

### Output Generation
- `write_transport_file()`: Creates transport NetCDF files
- `write_variable_file()`: Creates variable NetCDF files

## Usage
The project contains two main processing scripts:

1. Transport Processing:
```bash
python src/transport_processor.py
```

2. Variables Processing:
```bash
python src/variables_processor.py
```

## Data Requirements
- Input NetCDF files containing transport data
- BGM (Box Geometry Model) files
- Mesh configuration file
- Variable data files 