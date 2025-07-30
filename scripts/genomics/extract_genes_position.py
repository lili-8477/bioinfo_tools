#!/usr/bin/env python3

import argparse
import pybedtools
from collections import defaultdict

def extract_gene_positions_from_gtf(gtf_file, gene_list_file, output_bed='genes.bed'):
    """
    Extract gene positions from a GTF file for genes listed in a given text file and write to a BED file.
    Removes entries with unknown or non-standard chromosomes.
    """
    gtf = pybedtools.BedTool(gtf_file)
    gene_regions = {}

    # Define standard chromosomes (you can adjust for your use case)
    standard_chroms = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}

    # Read gene list with explicit encoding and strip any BOM or whitespace
    try:
        with open(gene_list_file, 'r', encoding='utf-8-sig') as f:
            gene_list = [line.strip() for line in f if line.strip()]  # Changed to list
            gene_set = set(gene_list)  # Create a set for faster lookups
    except UnicodeDecodeError:
        # If utf-8-sig fails, try regular utf-8
        with open(gene_list_file, 'r', encoding='utf-8') as f:
            gene_list = [line.strip() for line in f if line.strip()]  # Changed to list
            gene_set = set(gene_list)  # Create a set for faster lookups
    
    print(f"Number of genes in list: {len(gene_list)}")
    print("First few genes in list:", gene_list[:5])  # Now using list directly
    print("Last few genes in list:", gene_list[-5:])  # Now using list directly

    genes_found = 0
    genes_not_found = 0
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

        if gene_name in gene_set:  # Using set for faster lookups
            genes_found += 1
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
        else:
            genes_not_found += 1

    print(f"\nDebugging Information:")
    print(f"Total genes found in GTF: {genes_found}")
    print(f"Total genes not found in GTF: {genes_not_found}")
    print(f"Total non-standard chromosomes skipped: {non_standard_chroms}")
    print(f"Total unique genes in output: {len(gene_regions)}")

    if len(gene_regions) == 0:
        print("\nWARNING: No genes were found! Possible issues:")
        print("1. Gene names in your list might not match the GTF file format")
        print("2. The chromosome format might be different (e.g., 'chr1' vs '1')")
        
        # Print a sample of the GTF file
        print("\nSample of GTF file:")
        for feature in gtf[:5]:
            print(feature)

    gene_bed = pybedtools.BedTool((str(iv) for iv in gene_regions.values()))
    gene_bed.saveas(output_bed)
    print(f"\nGene BED saved to: {output_bed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract gene positions from a GTF file for genes listed in a given text file and write to a BED file.")
    parser.add_argument("gtf_file", help="Path to the input GTF file.")
    parser.add_argument("gene_list_file", help="Path to a text file containing gene names (one per line).")
    parser.add_argument("--output_bed", default="genes.bed", help="Output BED file name (default: genes.bed)")
    args = parser.parse_args()

    extract_gene_positions_from_gtf(args.gtf_file, args.gene_list_file, args.output_bed)