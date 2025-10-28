#!/bin/bash
# Compilation script for NESP-Atlantis Pilot Report
# ================================================

echo "NESP-Atlantis Pilot Report Compilation"
echo "======================================"

# Check if pdflatex is available
if ! command -v pdflatex &> /dev/null; then
    echo "Error: pdflatex is not installed or not in PATH"
    echo "Please install a LaTeX distribution (e.g., TeX Live, MiKTeX)"
    exit 1
fi

# Check if the main tex file exists
if [ ! -f "NESP-Atlantis_Pilot_Report.tex" ]; then
    echo "Error: NESP-Atlantis_Pilot_Report.tex not found"
    exit 1
fi

# Check if figures directory exists
if [ ! -d "figures" ]; then
    echo "Error: figures directory not found"
    exit 1
fi

echo "Compiling LaTeX document..."
echo "Running pdflatex (first pass)..."
pdflatex NESP-Atlantis_Pilot_Report.tex

echo "Running pdflatex (second pass for references)..."
pdflatex NESP-Atlantis_Pilot_Report.tex

# Check if compilation was successful
if [ -f "NESP-Atlantis_Pilot_Report.pdf" ]; then
    echo ""
    echo "✓ Compilation successful!"
    echo "✓ PDF generated: NESP-Atlantis_Pilot_Report.pdf"
    
    # Get file size
    SIZE=$(du -h NESP-Atlantis_Pilot_Report.pdf | cut -f1)
    echo "✓ File size: $SIZE"
    
    # Clean up auxiliary files
    echo "Cleaning up auxiliary files..."
    rm -f *.aux *.log *.out *.toc *.fdb_latexmk *.fls *.synctex.gz
    
    echo ""
    echo "Report ready for viewing!"
    echo "To view: xdg-open NESP-Atlantis_Pilot_Report.pdf (Linux)"
    echo "         open NESP-Atlantis_Pilot_Report.pdf (macOS)"
    
else
    echo ""
    echo "✗ Compilation failed!"
    echo "Check the .log file for error details"
    exit 1
fi



