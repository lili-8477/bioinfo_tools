# Genomics Utilities

Python scripts for processing genomic annotations and coordinate data.

## Tools

### extract_all_genes.py
Extract all gene positions from GTF files and convert to BED format.

**Features:**
- Filters out non-standard chromosomes
- Removes duplicate genes
- Provides processing statistics
- Creates standard 6-column BED output

**Usage:**
```bash
python extract_all_genes.py input.gtf --output_bed genes.bed
```

### extract_genes_position.py
Enhanced gene position extraction with additional filtering options.

**Usage:**
```bash
python extract_genes_position.py input.gtf
```

### extract_promoters.py
Extract promoter regions from gene annotations.

**Features:**
- Configurable promoter window sizes
- Strand-aware promoter definition
- TSS-based coordinate calculation

**Usage:**
```bash
python extract_promoters.py input.gtf
```

### bed_split_proximal_distal.py
Split BED files into proximal and distal regions based on genomic features.

**Usage:**
```bash
python bed_split_proximal_distal.py input.bed
```

## Dependencies

- **pybedtools** - Required for all scripts
- **argparse** - Command line parsing
- **collections** - Data structures

## Input Formats

- **GTF files** - Standard GTF/GFF format with gene_name and gene_id attributes
- **BED files** - Standard BED format (tab-separated)

## Output Formats

- **BED files** - 6-column format: chr, start, end, name, score, strand
- Processing logs and statistics printed to stdout