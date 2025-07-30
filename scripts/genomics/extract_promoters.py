#!/usr/bin/env python3
import argparse
import pybedtools
from collections import defaultdict

def extract_promoters_from_gtf(gtf_file, upstream=2000, downstream=2000, output_bed='promoters.bed'):
    """
    Extract promoter regions from a GTF file and write to a BED file.
    """
    gtf = pybedtools.BedTool(gtf_file)
    promoter_regions = {}

    for feature in gtf:
        if feature[2] != 'transcript':
            continue

        # Parse attributes from GTF
        try:
            attrs = dict(item.strip().replace('"', '').split(' ') for item in feature[8].strip(';').split('; '))
        except Exception as e:
            continue  # skip malformed lines

        gene_name = attrs.get('gene_name') or attrs.get('gene_id')
        if not gene_name:
            continue

        chrom = feature.chrom
        strand = feature.strand
        start = int(feature.start)
        end = int(feature.end)

        tss = start if strand == '+' else end

        if strand == '+':
            promoter_start = max(0, tss - upstream)
            promoter_end = tss + downstream
        else:
            promoter_start = max(0, tss - downstream)
            promoter_end = tss + upstream

        # Avoid duplicates by gene
        if gene_name not in promoter_regions:
            promoter_regions[gene_name] = pybedtools.create_interval_from_list([
                chrom,
                str(promoter_start),
                str(promoter_end),
                gene_name,
                '0',
                strand
            ])

    # Convert to BedTool and save
    promoter_bed = pybedtools.BedTool(list(promoter_regions.values()))
    promoter_bed.saveas(output_bed)
    print(f"Promoter BED saved to: {output_bed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract promoter regions from a GTF file and write to a BED file.")
    parser.add_argument("gtf_file", help="Path to the input GTF file.")
    parser.add_argument("--upstream", type=int, default=2000, help="Bases upstream of TSS (default: 2000)")
    parser.add_argument("--downstream", type=int, default=2000, help="Bases downstream of TSS (default: 2000)")
    parser.add_argument("--output_bed", default="promoters.bed", help="Output BED file name (default: promoters.bed)")
    args = parser.parse_args()

    extract_promoters_from_gtf(args.gtf_file, args.upstream, args.downstream, args.output_bed) 