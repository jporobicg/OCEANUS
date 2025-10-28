# NESP-Atlantis Pilot Test Documentation

This folder contains the complete documentation and results from the pilot test of the OCEANUS package for transforming BASS-2 hydrodynamic model output to Atlantis model structure.

## Contents

- **`NESP-Atlantis_Pilot_Report.md`** - Complete technical report (Markdown)
- **`NESP-Atlantis_Pilot_Report.tex`** - LaTeX source document
- **`figures/`** - All generated visualizations and diagnostic plots
- **`Makefile`** - LaTeX compilation makefile
- **`compile_report.sh`** - Compilation script

## Generated Figures

1. **`workflow_diagram.png`** - Data transformation workflow
2. **`grid_comparison.png`** - Structured vs unstructured grid comparison
3. **`summary_table.png`** - Pilot run statistics summary
4. **`correction_statistics.png`** - Mass balance correction diagnostics
5. **`surface_temperature_salinity_atlantis.png`** - Surface fields with polygon overlay
6. **`temperature_comparison.png`** - Direct comparison: BASS-2 vs Atlantis temperature

## Quick Summary

- **Domain:** Bass Strait region (140.4°E - 151.9°E, -41.4°N - -36.8°N)
- **Boxes:** 147 Atlantis model boxes
- **Corrections:** 1,440 mass balance corrections applied
- **Coverage:** 95.2% domain coverage achieved
- **Status:** Pilot test completed successfully

## Usage

### Viewing the Report

**Markdown version:** `NESP-Atlantis_Pilot_Report.md`

**LaTeX version:** Compile to PDF using one of these methods:

#### Method 1: Using the compilation script
```bash
./compile_report.sh
```

#### Method 2: Using make
```bash
make all
```

#### Method 3: Direct LaTeX compilation
```bash
pdflatex NESP-Atlantis_Pilot_Report.tex
pdflatex NESP-Atlantis_Pilot_Report.tex  # Run twice for references
```

### Requirements

- LaTeX distribution (TeX Live, MiKTeX, etc.)
- `pdflatex` command available

All figures are referenced within the report and can be viewed independently in the `figures/` subdirectory.
