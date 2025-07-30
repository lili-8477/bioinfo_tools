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
RENAME_TABLE="rename20828.txt"
GENOME_SIZE_FILE="~/software/genome_ref_data/hg38.chrom.sizes"
MACS2_OUTDIR="macs2"
BEDGRAPH2BIGWIG="~/software/genome_ref_data/bedGraphToBigWig"

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

echo "[Step 2] Running MACS2 in $PEAKMODE mode for chip: $CHIP..."

# === Step 3: Run MACS2 ===
# Get all unique replicate groups
cut -f2 "$RENAME_TABLE" | awk -v chip="$CHIP" '
{
  split($1, parts, "_");
  rep = parts[length(parts)];
  base = parts[1];
  for (i = 2; i < length(parts) - 1; i++) {
    base = base "_" parts[i];
  }
  if (index($1, chip) > 0) {
    print base, rep;
  }
}' | sort -u | while read -r BASE REP; do
  TREATMENT_BAM=$(ls ${BASE}_${CHIP}_${REP}.bam 2>/dev/null || true)
  CONTROL_BAM=$(ls ${BASE}_INPUT_${REP}.bam 2>/dev/null || true)

  if [[ -n "$TREATMENT_BAM" && -n "$CONTROL_BAM" ]]; then
    echo "✅ Running MACS2 for $TREATMENT_BAM vs $CONTROL_BAM"
    # your MACS2 call here
    # Derive a readable name and allow overrides for genome and thresholds
    OUTNAME="${BASE}_${CHIP}_${REP}"
    GENOME="${GENOME:-hs}"        # use "hs" for human, "mm" for mouse, or an effective size (e.g., 2.7e9)
    QVAL="${QVAL:-0.001}"          # narrow peak q-value cutoff
    BROAD_QVAL="${BROAD_QVAL:-0.01}"  # broad peak q-value cutoff

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
    echo "⚠️ Skipping ${BASE}_${CHIP}_${REP}: missing treatment or control"
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

