# Reference Annotations

This directory contains genome annotation files for various species.

## Files

### refgene/
- **hg38.refGene.gtf.gz** - Human genome (hg38) gene annotations (included)
- **mm10.refGene.gtf** - Mouse genome (mm10) gene annotations (download separately - 145MB)

### Other Files
- **mm10_introns.tab.gz** - Mouse intron annotations

## Large File Downloads

Due to GitHub file size limits, some large annotation files need to be downloaded separately:

### mm10.refGene.gtf (145MB)
Download from UCSC Genome Browser or other genome annotation sources:
```bash
# Example download (adjust URL as needed)
wget -O reference/annotations/refgene/mm10.refGene.gtf [URL_TO_MM10_GTF]
```

## File Formats

- **GTF files** - Gene Transfer Format with gene and transcript annotations
- **TAB files** - Tab-separated annotation files
- **GZ files** - Gzip compressed annotation files

## Usage

These annotation files are used by the genomics utilities in `scripts/genomics/` for:
- Gene position extraction
- Promoter region identification
- Intron/exon analysis
- Coordinate conversion