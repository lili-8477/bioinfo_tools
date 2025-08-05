#!/usr/bin/env python3

import argparse
from pybedtools import BedTool
import pandas as pd
import os

def extract_genes_from_gtf(gtf_file):
    """Extract gene features from GTF"""
    gtf = BedTool(gtf_file)
    genes = gtf.filter(lambda x: x[2] == 'transcript').saveas()
    return genes

def find_nearest_genes(bed_file, gtf_file, output_tsv, gene_list_file):
    regions = BedTool(bed_file)
    genes = extract_genes_from_gtf(gtf_file)

    # Determine number of columns in BED file
    bed_columns = 0
    with open(bed_file, 'r') as f:
        first_line = f.readline().strip()
        if first_line:
            bed_columns = len(first_line.split('\t'))
    
    print(f"Detected {bed_columns} columns in BED file")

    # Find nearest gene
    nearest = regions.sort().closest(genes.sort(), d=True)

    results = []
    gene_set = set()

    for feature in nearest:
        chrom, start, end = feature[0], feature[1], feature[2]
        
        # After closest operation: bed_columns + gtf fields (9)
        # GTF attributes column is at index: bed_columns + 8 (0-indexed)
        gtf_attributes_index = bed_columns + 8
        
        if len(feature) > gtf_attributes_index:
            gene_info = feature[gtf_attributes_index]
        else:
            gene_info = ""
        
        distance = feature[-1]
        gene_name = "NA"

        if gene_info:
            for item in gene_info.split(';'):
                item = item.strip()
                if item.startswith("gene_name"):
                    gene_name = item.split('"')[1]
                    break

        results.append([chrom, start, end, gene_name, distance])
        if gene_name != "NA":
            gene_set.add(gene_name)

    # Write full TSV
    df = pd.DataFrame(results, columns=["chrom", "start", "end", "nearest_gene", "distance"])
    df.to_csv(output_tsv, sep='\t', index=False)

    # Write unique gene list
    with open(gene_list_file, 'w') as f:
        for gene in sorted(gene_set):
            f.write(f"{gene}\n")

def main():
    parser = argparse.ArgumentParser(description="Find nearest genes using pybedtools and GTF")
    parser.add_argument("--bed", required=True, help="Input BED file")
    parser.add_argument("--gtf", required=True, help="RefSeq GTF file")
    parser.add_argument("--out", required=True, help="Output TSV file for nearest genes")
    parser.add_argument("--gene-list", required=True, help="Output TXT file for unique gene list")
    args = parser.parse_args()

    if not os.path.exists(args.bed):
        raise FileNotFoundError(f"BED file not found: {args.bed}")
    if not os.path.exists(args.gtf):
        raise FileNotFoundError(f"GTF file not found: {args.gtf}")

    find_nearest_genes(args.bed, args.gtf, args.out, args.gene_list)

if __name__ == "__main__":
    main()
