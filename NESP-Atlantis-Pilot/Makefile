# Makefile for NESP-Atlantis Pilot Report
# ========================================

# LaTeX compiler
LATEX = pdflatex
BIBTEX = bibtex

# Main document
MAIN = NESP-Atlantis_Pilot_Report
TEX_FILE = $(MAIN).tex
PDF_FILE = $(MAIN).pdf

# Default target
all: $(PDF_FILE)

# Compile LaTeX document
$(PDF_FILE): $(TEX_FILE)
	$(LATEX) $(TEX_FILE)
	$(LATEX) $(TEX_FILE)  # Run twice for proper references

# Clean up auxiliary files
clean:
	rm -f *.aux *.log *.out *.toc *.fdb_latexmk *.fls *.synctex.gz

# Clean everything including PDF
distclean: clean
	rm -f $(PDF_FILE)

# View the PDF (Linux)
view: $(PDF_FILE)
	xdg-open $(PDF_FILE)

# View the PDF (macOS)
view-mac: $(PDF_FILE)
	open $(PDF_FILE)

# Help
help:
	@echo "Available targets:"
	@echo "  all        - Compile the LaTeX document to PDF"
	@echo "  clean      - Remove auxiliary LaTeX files"
	@echo "  distclean  - Remove all generated files including PDF"
	@echo "  view       - Open the PDF (Linux)"
	@echo "  view-mac   - Open the PDF (macOS)"
	@echo "  help       - Show this help message"

.PHONY: all clean distclean view view-mac help
