#!/bin/bash
set -euo pipefail

# === Default settings ===
CHROM_SIZES_DEFAULT="/Users/li/bioinfo_tools/hg19.chrom.sizes"

# === Usage ===
# ./merge_bigwigs.sh output_prefix [chrom_sizes.bed] sample1.bw sample2.bw ...

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <output_prefix> [chrom_sizes_file] <bigwig1> <bigwig2> [bigwig3 ...]"
  echo "Note: If chrom_sizes_file is omitted, default is: $CHROM_SIZES_DEFAULT"
  exit 1
fi

# === Determine if chrom_sizes_file is provided ===
OUT_PREFIX="$1"
shift

# Check if the second argument ends with ".bw"
if [[ "$1" == *.bw ]]; then
  CHROM_SIZES="$CHROM_SIZES_DEFAULT"
else
  CHROM_SIZES="$1"
  shift
fi

# === Remaining arguments are BigWig files ===
BW_FILES=("$@")

# === Output files ===
MERGED_BG="${OUT_PREFIX}.bedGraph"
SORTED_BG="${OUT_PREFIX}.sorted.bedGraph"
MERGED_BW="${OUT_PREFIX}.bw"

# === Check required tools ===
for tool in bigWigMerge bedGraphToBigWig sort; do
  command -v $tool >/dev/null 2>&1 || { echo >&2 "$tool not found in PATH"; exit 1; }
done

# === Merge BigWig files ===
echo "Merging ${#BW_FILES[@]} BigWig files into ${MERGED_BG}..."
bigWigMerge "${BW_FILES[@]}" "$MERGED_BG"

# === Sort bedGraph ===
echo "Sorting BedGraph to ${SORTED_BG}..."
LC_COLLATE=C sort -k1,1 -k2,2n --stable "$MERGED_BG" > "$SORTED_BG"

# === Convert to BigWig ===
echo "Converting to BigWig: ${MERGED_BW} (chrom sizes: $CHROM_SIZES)"
bedGraphToBigWig "$SORTED_BG" "$CHROM_SIZES" "$MERGED_BW"

echo "Done. Outputs:"
echo "  BedGraph:   $MERGED_BG"
echo "  Sorted BG:  $SORTED_BG"
echo "  Merged BW:  $MERGED_BW"
