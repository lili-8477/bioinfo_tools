#!/usr/bin/env python3

import argparse
import os
import subprocess
import tempfile
import gzip


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split a BED file into proximal and distal regions based on RefSeq TSSs."
    )
    parser.add_argument("-i", "--input", required=True, help="Input BED file to split.")
    parser.add_argument("-r", "--refseq", required=True, help="RefSeq GTF or BED file.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output files.")
    parser.add_argument("-p", "--promoter_window", type=int, default=2000, help="Promoter window size (default: ±2000 bp).")
    parser.add_argument("--gtf", action="store_true", help="Specify if refseq input is in GTF format (default: BED).")
    return parser.parse_args()


def smart_open(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def extract_promoters_from_gtf(gtf_file, promoter_window, output_file):
    gene_coords = {}  # gene_name -> (chrom, start, end, strand)

    def parse_attributes(attr_string):
        attrs = {}
        for attr in attr_string.strip().split(";"):
            if not attr.strip():
                continue
            key_value = attr.strip().split(" ", 1)
            if len(key_value) == 2:
                key, value = key_value
                attrs[key] = value.strip('"')
        return attrs

    with open(output_file, "w") as out:
        with smart_open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 9 or cols[2] != "transcript":
                    continue
                chrom = cols[0]
                start = int(cols[3]) - 1  # GTF to BED
                end = int(cols[4])
                strand = cols[6]
                attrs = parse_attributes(cols[8])
                gene_name = attrs.get("gene_name", ".")

                # Store gene position (not promoter window)
                if gene_name != ".":
                    gene_coords[gene_name] = (chrom, start, end, strand)

                # Write promoter window for use in intersection
                tss0 = start if strand == "+" else end
                prom_start = max(0, tss0 - promoter_window)
                prom_end = tss0 + promoter_window
                out.write(f"{chrom}\t{prom_start}\t{prom_end}\t{gene_name}\t.\t{strand}\n")

    return gene_coords

def split_regions(input_bed, promoter_bed, output_prefix, is_gtf, gene_coords=None):
    proximal_bed = f"{output_prefix}_proximal.bed"
    distal_bed = f"{output_prefix}_distal.bed"

    subprocess.run([
        "bedtools", "intersect",
        "-a", input_bed,
        "-b", promoter_bed,
        "-wa", "-u"
    ], stdout=open(proximal_bed, "w"), check=True)

    subprocess.run([
        "bedtools", "intersect",
        "-a", input_bed,
        "-b", promoter_bed,
        "-v"
    ], stdout=open(distal_bed, "w"), check=True)

    print(f"Proximal BED saved to: {proximal_bed}")
    print(f"Distal BED saved to:   {distal_bed}")

    if is_gtf and gene_coords:
        intersected_bed = f"{output_prefix}_intersected_genes.bed"
        intersected_genes_txt = f"{output_prefix}_intersected_genes.txt"

        seen = set()
        with open(intersected_bed, "w") as out_bed, open(intersected_genes_txt, "w") as out_txt:
            result = subprocess.run([
                "bedtools", "intersect",
                "-a", promoter_bed,
                "-b", input_bed,
                "-wa"
            ], capture_output=True, text=True, check=True)

            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                fields = line.strip().split("\t")
                gene_name = fields[3]
                if gene_name != "." and gene_name not in seen and gene_name in gene_coords:
                    seen.add(gene_name)
                    chrom, start, end, strand = gene_coords[gene_name]
                    out_bed.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}\n")
                    out_txt.write(gene_name + "\n")

        print(f"Intersected gene BED saved to: {intersected_bed}")
        print(f"Gene name list saved to:       {intersected_genes_txt}")


def main():
    args = parse_args()

    if args.gtf:
        promoter_bed = f"{os.path.basename(args.refseq).split('.')[0]}_promoter_{args.promoter_window}.bed"
        gene_coords = extract_promoters_from_gtf(args.refseq, args.promoter_window, promoter_bed)
        print(f"Promoter BED saved: {promoter_bed}")
    else:
        promoter_bed = args.refseq
        gene_coords = {}

    # Perform splitting and (in GTF mode) gene intersection
    split_regions(args.input, promoter_bed, args.output_prefix, is_gtf=args.gtf, gene_coords=gene_coords)

    


if __name__ == "__main__":
    main()

