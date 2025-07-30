#!/bin/bash

#SBATCH --job-name="tunnel"
#SBATCH --time=4:00:00 # walltime
#SBATCH --cpus-per-task=12 # number of cores
#SBATCH --mem-per-cpu=12G # memory per CPU core
#SBATCH --partition=jonesk-gpu-np # partition, abbreviated by -p
#SBATCH --account=jonesk-gpu-np

module use $HOME/MyModules
module load miniconda3/latest
source activate bioinfo_env

set -euo pipefail

# === USER CONFIGURATION ===
RENAME_TABLE="rename21004.txt"
GENOME_SIZE_FILE="/uufs/chpc.utah.edu/common/home/u6025146/software/genome_ref_data/hg38.chrom.sizes"
MACS2_OUTDIR="macs2"
BEDGRAPH2BIGWIG="/uufs/chpc.utah.edu/common/home/u6025146/software/genome_ref_data/bedGraphToBigWig"

# === Create output directory ===
mkdir -p "$MACS2_OUTDIR"

#set chip and peak mode!!!!
CHIP=${CHIP:-GLTSCR1}
PEAKMODE=${PEAKMODE:-narrow}
if [[ -z "$CHIP" || -z "$PEAKMODE" ]]; then
  echo "❌ CHIP and PEAKMODE must be set."
  echo "Usage: sbatch --export=CHIP=GLTSCR1,PEAKMODE=narrow your_script.sh"
  exit 1
fi

# === Step 1: Rename BAM files ===
echo "[Step 1] Renaming BAM files..."

# If your rename table might contain Windows CRLF line endings, uncomment the next line:
# dos2unix -q "$RENAME_TABLE" 2>/dev/null || true

awk 'NF>=2 {gsub("\r",""); print $1, $2}' "$RENAME_TABLE" | while read -r old new; do
    src="${old}.bam"
    dst="${new}.bam"

    if [[ -e "$src" ]]; then
        echo "Renaming $src -> $dst"
        mv -v -- "$src" "$dst"

        # Rename index if present (handles either foo.bam.bai or foo.bai)
        for idx in "${old}.bam.bai" "${old}.bai"; do
            if [[ -e "$idx" ]]; then
                newidx="${dst}.bai"
                echo "Renaming $idx -> $newidx"
                mv -v -- "$idx" "$newidx"
                break
            fi
        done
    else
        echo "⚠️  Skipping ${src} (not found)"
    fi
done

echo "[Step 3] Running MACS2 in $PEAKMODE mode for chip: $CHIP..."

# === Improved Step 3: Run MACS2 ===
for TREATMENT_BAM in *_${CHIP}_REP*.bam; do
  if [[ ! -f "$TREATMENT_BAM" ]]; then
    echo "⚠️ No BAM files found for CHIP=$CHIP"
    continue
  fi

  BASE=$(echo "$TREATMENT_BAM" | sed "s/_${CHIP}_REP[0-9]*\.bam//")
  REP=$(echo "$TREATMENT_BAM" | grep -oP "REP[0-9]+")

  CONTROL_BAM="${BASE}_INPUT_${REP}.bam"

  if [[ -f "$CONTROL_BAM" ]]; then
    echo "✅ Running MACS2 for $TREATMENT_BAM vs $CONTROL_BAM"

    OUTNAME="${BASE}_${CHIP}_${REP}"
    GENOME="${GENOME:-hs}"        
    QVAL="${QVAL:-0.001}"          
    BROAD_QVAL="${BROAD_QVAL:-0.01}"

    if [[ "$PEAKMODE" == "narrow" ]]; then
      macs2 callpeak \
        -t "$TREATMENT_BAM" \
        -c "$CONTROL_BAM" \
        -g "$GENOME" \
        -n "$OUTNAME" \
        --outdir "$MACS2_OUTDIR" \
        -B --SPMR \
        --call-summits \
        -q "$QVAL" \
        --keep-dup auto \
        2> "${MACS2_OUTDIR}/${OUTNAME}.log"

    elif [[ "$PEAKMODE" == "broad" ]]; then
      macs2 callpeak \
        -t "$TREATMENT_BAM" \
        -c "$CONTROL_BAM" \
        -g "$GENOME" \
        -n "$OUTNAME" \
        --outdir "$MACS2_OUTDIR" \
        -B --SPMR \
        --broad \
        --broad-cutoff "$BROAD_QVAL" \
        --keep-dup auto \
        2> "${MACS2_OUTDIR}/${OUTNAME}.log"

    else
      echo "❌ Unknown PEAKMODE: $PEAKMODE (expected 'narrow' or 'broad')" >&2
      exit 1
    fi

  else
    echo "⚠️ Skipping ${TREATMENT_BAM}: missing control ${CONTROL_BAM}"
  fi
done



# === Step 4: Convert bedGraph to bigWig ===
echo "[Step 4] Converting bedGraph to bigWig..."

cd "$MACS2_OUTDIR"
for bdg in *_treat_pileup.bdg; do
    [[ -f "$bdg" ]] || continue
    prefix="${bdg%_treat_pileup.bdg}"
    echo "Converting $bdg -> ${prefix}.bw"
    "$BEDGRAPH2BIGWIG" "$bdg" "$GENOME_SIZE_FILE" "${prefix}.bw"
done
cd ..

echo "✅ Pipeline finished!"

