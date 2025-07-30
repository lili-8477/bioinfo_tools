# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics tools repository containing multiple computational tools for genomics and single-cell analysis:

1. **CellChat** - R package for cell-cell communication analysis from single-cell and spatially resolved transcriptomics data
2. **PFS_scanner** - Python tool for finding target sequences in RNA with pattern matching
3. **Genomics utilities** - Collection of shell scripts and Python tools for ChIP-seq, peak calling, and gene annotation processing

## Core Components

### CellChat R Package (`CellChat/`)
- **Main Class**: CellChat object-oriented system for storing and analyzing cell communication data
- **Key R files**:
  - `R/CellChat_class.R` - Core S4 class definition
  - `R/modeling.R` - Communication probability computation
  - `R/analysis.R` - Network analysis and pattern identification
  - `R/visualization.R` - Plotting and visualization functions
  - `R/database.R` - CellChatDB interaction database management
- **Dependencies**: Uses NMF (>= 0.23.0), circlize (>= 0.4.12), ComplexHeatmap for core functionality
- **Data**: Contains curated ligand-receptor databases for human, mouse, and zebrafish

### Python Tools
- **PFS_scanner** (`PFS_scanner/PFS_Scanner.py`): RNA sequence pattern matching with customizable patterns (N=[ACGU], V=[ACG])
- **Gene extraction utilities**: Scripts for processing GTF files and extracting gene positions using pybedtools

### Shell Scripts for Genomics Pipelines
- **MACS2 pipeline** (`macs2_pipeline.sh`): ChIP-seq peak calling with automatic BAM renaming and bigWig conversion
- **Peak calling from BigWig** (`callPeaksFromBigwig.sh`): Direct peak calling from signal tracks
- **BigWig utilities**: Conversion tools using UCSC genome browser utilities

## Common Development Commands

### CellChat R Package
```bash
# Install from source
R CMD INSTALL CellChat/

# Build package documentation
cd CellChat && R -e "roxygen2::roxygenise()"

# Run R package checks
cd CellChat && R CMD check .
```

### Python Tools
```bash
# Run PFS scanner
python PFS_scanner/PFS_Scanner.py

# Gene extraction example
python extract_all_genes.py input.gtf --output_bed output.bed
```

### Genomics Pipeline Usage
```bash
# MACS2 peak calling pipeline
./macs2_pipeline.sh

# Peak calling from BigWig files
./callPeaksFromBigwig.sh -c 5.0 -l 500 sample.bigwig
```

## Key Architecture Patterns

### CellChat Framework
- **Object-oriented design**: Central CellChat S4 class stores expression data, metadata, and communication networks
- **Modular analysis pipeline**: Sequential workflow from data preprocessing → communication inference → network analysis → visualization
- **Spatial analysis support**: Handles both single-cell RNA-seq and spatially resolved transcriptomics with distance-based communication modeling
- **Database integration**: Uses CellChatDB for curated ligand-receptor interactions with support for custom databases

### Data Processing Flow
1. **Input**: Expression matrices, cell metadata, spatial coordinates (optional)
2. **Preprocessing**: Data normalization, gene filtering, cell type identification
3. **Communication inference**: Probability calculation based on ligand-receptor expression and spatial proximity
4. **Network analysis**: Pattern recognition, centrality measures, pathway enrichment
5. **Visualization**: Multiple plot types including chord diagrams, heatmaps, spatial plots

### Genome Data Handling
- **Reference files**: Chromosome sizes for hg19, hg38, mm10, mm39 genomes
- **GTF processing**: Standard transcript/gene extraction with chromosome filtering
- **BigWig workflows**: Signal processing and peak identification pipelines

## File Naming Conventions
- **R functions**: camelCase (e.g., `computeCommunProb`, `netVisual_aggregate`)
- **Python scripts**: snake_case with descriptive names
- **Shell scripts**: descriptive names with `.sh` extension
- **Data files**: Genome build prefixes (hg38, mm10) for reference files