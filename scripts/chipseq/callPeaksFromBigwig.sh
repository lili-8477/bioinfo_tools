#!/bin/bash

# Generalized script to call peaks from BigWig files using MACS2
# Converts BigWig to bedGraph and then calls peaks with macs2 bdgpeakcall

set -e  # Exit on any error

# Default parameters
CUTOFF=3.0
LENGTH=300
GAP=100

# Usage message
usage() {
    echo "Usage: $0 [-c cutoff] [-l min_length] [-g max_gap] bigwig1 [bigwig2 ...]"
    echo ""
    echo "Options:"
    echo "  -c CUTOFF        Signal cutoff for peak calling (default: $CUTOFF)"
    echo "  -l LENGTH        Minimum peak length (default: $LENGTH)"
    echo "  -g GAP           Maximum gap to merge nearby signals (default: $GAP)"
    echo ""
    echo "Example:"
    echo "  $0 -c 5.0 -l 500 -g 200 sample1.bigwig sample2.bigwig"
    exit 1
}

# Parse options
while getopts ":c:l:g:" opt; do
    case $opt in
        c) CUTOFF="$OPTARG" ;;
        l) LENGTH="$OPTARG" ;;
        g) GAP="$OPTARG" ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

# Check for at least one BigWig file
if [ "$#" -lt 1 ]; then
    usage
fi

echo "Starting generalized peak calling from BigWig files..."

# Check dependencies
command -v bigWigToBedGraph >/dev/null 2>&1 || { echo "bigWigToBedGraph is required but not installed. Please install UCSC tools."; exit 1; }
command -v macs2 >/dev/null 2>&1 || { echo "macs2 is required but not installed. Please install MACS2."; exit 1; }

# Create output directories
mkdir -p peaks_output bedgraph_files

# Process each BigWig file
for BW in "$@"; do
    if [ ! -f "$BW" ]; then
        echo "File not found: $BW"
        continue
    fi

    BASE=$(basename "$BW" .bigwig)
    BEDGRAPH="bedgraph_files/${BASE}.bedgraph"
    PEAKS="peaks_output/${BASE}_peaks.bed"
    
    echo "Processing $BW..."
    if [ -f "$BEDGRAPH" ]; then
        echo "  Skipping bedGraph conversion: $BEDGRAPH already exists"
    else
        echo "  Converting to bedGraph: $BEDGRAPH"
        bigWigToBedGraph "$BW" "$BEDGRAPH"
    fi

    echo "  Calling peaks with MACS2: $PEAKS"
    macs2 bdgpeakcall -i "$BEDGRAPH" -c "$CUTOFF" -l "$LENGTH" -g "$GAP" -o "$PEAKS"
done

# Summary
echo ""
echo "Peak calling completed!"
echo "Results saved in 'peaks_output/' directory."

echo ""
echo "Peak count summary:"
for FILE in peaks_output/*.bed; do
    COUNT=$(wc -l < "$FILE")
    echo "$(basename "$FILE"): $COUNT peaks"
done

echo ""
echo "Parameters used:"
echo "- Cutoff (-c): $CUTOFF"
echo "- Minimum length (-l): $LENGTH"
echo "- Maximum gap (-g): $GAP"
