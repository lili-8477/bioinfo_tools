# ChIP-seq Analysis Tools

This directory contains scripts for ChIP-seq data analysis and peak calling.

## Tools

### macs2_pipeline.sh
Automated ChIP-seq peak calling pipeline using MACS2.

**Features:**
- Automatic BAM file renaming based on sample mapping
- Interactive ChIP target and peak mode selection
- Supports both narrow and broad peak calling
- Automatic bedGraph to BigWig conversion
- Batch processing of multiple samples

**Usage:**
```bash
./macs2_pipeline.sh
```

**Requirements:**
- rename_table.txt file with old->new filename mappings
- MACS2 installed
- bedGraphToBigWig tool (included in ../bin/)

### callPeaksFromBigwig.sh
Direct peak calling from BigWig signal files using MACS2 bdgpeakcall.

**Features:**
- Converts BigWig to bedGraph automatically
- Customizable signal cutoff, minimum length, and gap parameters
- Batch processing support
- Summary statistics output

**Usage:**
```bash
./callPeaksFromBigwig.sh [-c cutoff] [-l min_length] [-g max_gap] file1.bigwig [file2.bigwig ...]
```

### merge_bigwigs.sh
Utility for merging multiple BigWig files.

**Usage:**
```bash
./merge_bigwigs.sh
```

## Configuration

Default parameters can be modified in each script:
- Genome size files: Located in ../../reference/genome_sizes/
- Output directories: Typically macs2/ or peaks_output/
- Signal thresholds: Adjustable via command line options