#!/bin/bash
set -euo pipefail

# === Usage ===
# ./merge_bigwigs.sh output_prefix chrom_sizes.bed sample1.bw sample2.bw ...

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <output_prefix> <chrom_sizes_file> <bigwig1> <bigwig2> [bigwig3 ...]"
  exit 1
fi

# === Inputs ===
OUT_PREFIX="$1"
CHROM_SIZES="$2"
shift 2
BW_FILES=("$@")

# === Output files ===
MERGED_BG="${OUT_PREFIX}.bedGraph"
SORTED_BG="${OUT_PREFIX}.sorted.bedGraph"
MERGED_BW="${OUT_PREFIX}.bw"
# You can later run bedGraphToBigWig manually if needed

# === Check required tools ===
command -v bigWigMerge >/dev/null 2>&1 || { echo >&2 "bigWigMerge not found"; exit 1; }
command -v sort >/dev/null 2>&1 || { echo >&2 "GNU sort not found"; exit 1; }

# === Merge BigWig files into raw bedGraph ===
echo "Merging ${#BW_FILES[@]} BigWig files into ${MERGED_BG}..."
bigWigMerge "${BW_FILES[@]}" "$MERGED_BG"

# === Sort BedGraph for downstream tools (bedGraphToBigWig) ===
echo "Sorting BedGraph to ${SORTED_BG}..."
LC_COLLATE=C sort -k1,1 -k2,2n --stable "$MERGED_BG" > "$SORTED_BG" 
# Step 2: Convert to BigWig (requires chrom.sizes file)
bedGraphToBigWig "$SORTED_BG" /Users/li/bioinfo_tools/hg19.chrom.sizes "$MERGED_BW"

echo "Done. Sorted bedGraph saved to: $SORTED_BG"
