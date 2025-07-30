#!/usr/bin/env python3

import argparse
import pybedtools
from collections import defaultdict

def extract_all_gene_positions(gtf_file, output_bed='all_genes.bed'):
    """
    Extract all gene positions from a GTF file and write to a BED file.
    Removes entries with unknown or non-standard chromosomes.
    
    Parameters:
        gtf_file (str): Path to input GTF file.
        output_bed (str): Output BED file path.
    """
    gtf = pybedtools.BedTool(gtf_file)
    gene_regions = {}

    # Define standard chromosomes (you can adjust for your use case)
    standard_chroms = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}

    genes_processed = 0
    non_standard_chroms = 0

    for feature in gtf:
        if feature[2] != 'transcript':
            continue

        chrom = feature.chrom
        if chrom not in standard_chroms:
            non_standard_chroms += 1
            continue

        # Parse attributes from GTF
        try:
            attrs = dict(item.strip().replace('"', '').split(' ') for item in feature[8].strip(';').split('; '))
        except Exception as e:
            print(f"Error parsing attributes: {e}")
            continue

        gene_name = attrs.get('gene_name') or attrs.get('gene_id')
        if not gene_name:
            continue

        genes_processed += 1
        strand = feature.strand
        start = int(feature.start)
        end = int(feature.end)

        # Avoid duplicates by gene
        if gene_name not in gene_regions:
            gene_regions[gene_name] = pybedtools.create_interval_from_list([
                chrom,
                str(start),
                str(end),
                gene_name,
                '0',
                strand
            ])

    print(f"\nProcessing Information:")
    print(f"Total transcripts processed: {genes_processed}")
    print(f"Total unique genes found: {len(gene_regions)}")
    print(f"Total non-standard chromosomes skipped: {non_standard_chroms}")

    # Create BED file
    gene_bed = pybedtools.BedTool((str(iv) for iv in gene_regions.values()))
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract all gene positions from a GTF file and write to a BED file.")
    parser.add_argument("gtf_file", help="Path to the input GTF file.")
    parser.add_argument("--output_bed", default="all_genes.bed", help="Output BED file name (default: all_genes.bed)")
    args = parser.parse_args()

    extract_all_gene_positions(args.gtf_file, args.output_bed)