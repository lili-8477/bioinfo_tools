# Bioinfo Tools

A collection of computational tools for bioinformatics analysis, designed to support daily work in computational biology. This repository contains utilities for ChIP-seq analysis, genomics data processing, CRISPR analysis, and more.

## Repository Structure

```
bioinfo_tools/
├── bin/                    # Binary executables and external tools
├── scripts/               # Analysis scripts organized by domain
│   ├── chipseq/          # ChIP-seq analysis tools
│   ├── genomics/         # General genomics utilities
│   ├── RNAseq/           # RNA-seq analysis tools
│   └── crispr/           # CRISPR/RNA analysis tools
├── reference/             # Reference data and annotations
│   ├── genome_sizes/     # Chromosome size files for various genomes
│   └── annotations/      # GTF/GFF annotation files
├── config/               # Configuration files
└── docs/                 # Documentation and guides
```

## Tools Overview

### ChIP-seq Analysis (`scripts/chipseq/`)
- **macs2_pipeline.sh** - Automated ChIP-seq peak calling pipeline with MACS2
- **callPeaksFromBigwig.sh** - Direct peak calling from BigWig signal files
- **merge_bigwigs.sh** - Utility for merging multiple BigWig files

### Genomics Utilities (`scripts/genomics/`)
- **extract_all_genes.py** - Extract gene positions from GTF files to BED format
- **extract_genes_position.py** - Gene position extraction with filtering
- **extract_promoters.py** - Extract promoter regions from gene annotations
- **bed_split_proximal_distal.py** - Split BED files into proximal/distal regions

### RNAseq Analysis (`scripts/RNAseq/`)
- **homolog_converter.py** - Convert gene symbols between mouse and human using GSEAPY Biomart

### CRISPR/RNA Analysis (`scripts/crispr/`)
- **PFS_Scanner** - Python tool for finding target sequences in RNA with customizable pattern matching

### Binary Tools (`bin/`)
- UCSC Genome Browser utilities (bedGraphToBigWig, bigWigToBedGraph, etc.)
- Nextflow workflow manager

### Reference Data (`reference/`)
- **genome_sizes/** - Chromosome size files for hg19, hg38, mm10, mm39
- **annotations/** - GTF annotation files and intron databases

## Quick Start

### ChIP-seq Peak Calling
```bash
# Basic MACS2 pipeline
./scripts/chipseq/macs2_pipeline.sh

# Peak calling from BigWig files
./scripts/chipseq/callPeaksFromBigwig.sh -c 5.0 -l 500 sample.bigwig
```

### Gene Annotation Processing
```bash
# Extract all genes from GTF to BED format
python scripts/genomics/extract_all_genes.py annotation.gtf --output_bed genes.bed

# Extract promoter regions
python scripts/genomics/extract_promoters.py annotation.gtf
```

### Gene Homolog Conversion
```bash
# Convert mouse genes to human
./scripts/RNAseq/homolog_converter.py -g Actb Kif11 Tpx2 -d m2h

# Convert human genes to mouse from file
./scripts/RNAseq/homolog_converter.py -f human_genes.txt -d h2m -o mouse_orthologs.tsv
```

### CRISPR Analysis
```bash
# Run PFS scanner for target sequence identification
cd scripts/crispr/PFS_scanner
python PFS_Scanner.py
```

## Requirements

### Python Dependencies
- pybedtools
- pandas
- gseapy (for homolog conversion)
- argparse
- re (standard library)

### R Dependencies
- Not applicable (CellChat removed)

### External Tools
- MACS2 for peak calling
- UCSC Genome Browser utilities (included in bin/)
- Nextflow (optional, for workflow management)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/lili-8477/bioinfo_tools.git
cd bioinfo_tools
```

2. Make scripts executable:
```bash
chmod +x scripts/chipseq/*.sh
chmod +x bin/*
```

3. Add bin directory to your PATH (optional):
```bash
export PATH=$PATH:$(pwd)/bin
```

## Usage Examples

### ChIP-seq Analysis Workflow
1. Prepare your BAM files and rename table
2. Run the MACS2 pipeline: `./scripts/chipseq/macs2_pipeline.sh`
3. Convert results to BigWig format (handled automatically)

### Gene Annotation Processing
1. Download GTF annotations for your genome of interest
2. Extract gene positions: `python scripts/genomics/extract_all_genes.py genome.gtf`
3. Process specific regions as needed

## Contributing

This is a personal toolkit repository. Feel free to fork and adapt for your own use.

## License

Tools are provided as-is for research use. Please check individual tool licenses where applicable.

## Support

For issues or questions, please open an issue on GitHub.