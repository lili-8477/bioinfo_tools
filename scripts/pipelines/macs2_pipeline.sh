#!/bin/bash
#SBATCH --job-name="tunnel"
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=jonesk-gpu-np
#SBATCH --account=jonesk-gpu-np

module use $HOME/MyModules
module load miniconda3/latest
source activate bioinfo_env

set -euo pipefail

# --- USER SETTINGS (minimal) ---
MAPPING_FILE="${MAPPING_FILE:-list.txt}"   # 3 cols: ID  sample_name  genome
OUTROOT="${OUTROOT:-macs2}"
PEAKMODE="${PEAKMODE:-narrow}"                    # narrow | broad
QVAL="${QVAL:-0.001}"                             # for narrow
BROAD_QVAL="${BROAD_QVAL:-0.01}"                  # for broad
KEEPDUP="${KEEPDUP:-auto}"

# bedGraphToBigWig + chrom.sizes (set both if you want mm10 bigWigs too)
BEDGRAPH2BIGWIG="${BEDGRAPH2BIGWIG:-/uufs/chpc.utah.edu/common/home/u6025146/software/genome_ref_data/bedGraphToBigWig}"
CHROMSIZES_HG38="${CHROMSIZES_HG38:-/uufs/chpc.utah.edu/common/home/u6025146/software/genome_ref_data/hg38.chrom.sizes}"
CHROMSIZES_MM10="${CHROMSIZES_MM10:-/uufs/chpc.utah.edu/common/home/u6025146/software/genome_ref_data/mm10.chrom.sizes}"  # e.g. /path/to/mm10.chrom.sizes

mkdir -p "$OUTROOT"
command -v macs2 >/dev/null || { echo "macs2 not found"; exit 1; }

# --- Read mapping once ---
mapfile -t MAPLINES < <(grep -v '^[[:space:]]*$' "$MAPPING_FILE" | grep -v '^[#]')

# --- Step 1: Rename BAMs & indexes ---
echo "[1/3] Renaming BAMs"
for L in "${MAPLINES[@]}"; do
  IFS=$'\t' read -r ID NAME GENOME <<<"$L"
  [[ -z "${ID:-}" || -z "${NAME:-}" ]] && continue
  
  # Rename BAM file
  if [[ -e "${ID}.bam" ]]; then
    echo "  ${ID}.bam -> ${NAME}.bam"
    mv -- "${ID}.bam" "${NAME}.bam" || { echo "ERROR: Failed to rename ${ID}.bam"; exit 1; }
  fi
  
  # Rename index file (prefer .bam.bai, fallback to .bai)
  if [[ -e "${ID}.bam.bai" ]]; then
    echo "  ${ID}.bam.bai -> ${NAME}.bam.bai"
    mv -- "${ID}.bam.bai" "${NAME}.bam.bai" || { echo "ERROR: Failed to rename ${ID}.bam.bai"; exit 1; }
  elif [[ -e "${ID}.bai" ]]; then
    echo "  ${ID}.bai -> ${NAME}.bam.bai"
    mv -- "${ID}.bai" "${NAME}.bam.bai" || { echo "ERROR: Failed to rename ${ID}.bai"; exit 1; }
  fi
  
  # Verify renamed files exist
  [[ -f "${NAME}.bam" ]] || { echo "ERROR: ${NAME}.bam not found after rename"; exit 1; }
  [[ -f "${NAME}.bam.bai" ]] || echo "  WARNING: No index found for ${NAME}.bam"
done

# --- Find input per genome (sample_name contains 'input') ---
declare -A INPUT_FOR_GENOME
declare -A GENOMES_FOUND
shopt -s nocasematch
for L in "${MAPLINES[@]}"; do
  IFS=$'\t' read -r _ NAME GENOME <<<"$L"
  GENOMES_FOUND["$GENOME"]=1
  [[ "$NAME" =~ input ]] && INPUT_FOR_GENOME["$GENOME"]="$NAME"
done
shopt -u nocasematch

# Validate that each genome has an input control
echo "Validating input controls..."
for GENOME in "${!GENOMES_FOUND[@]}"; do
  if [[ -z "${INPUT_FOR_GENOME[$GENOME]:-}" ]]; then
    echo "ERROR: No input control found for genome $GENOME"
    echo "       Add a sample with 'input' in the name for genome $GENOME"
    exit 1
  else
    echo "  $GENOME: using ${INPUT_FOR_GENOME[$GENOME]} as control"
  fi
done

# --- Step 2: MACS2 peak calling ---
echo "[2/3] MACS2 ($PEAKMODE) ..."
for L in "${MAPLINES[@]}"; do
  IFS=$'\t' read -r _ NAME GENOME <<<"$L"
  # skip inputs
  shopt -s nocasematch; [[ "$NAME" =~ input ]] && { shopt -u nocasematch; continue; }; shopt -u nocasematch

  TREAT="${NAME}.bam"
  CTRL="${INPUT_FOR_GENOME[$GENOME]:-}.bam"
  # Verify treatment and control files exist
  if [[ ! -f "$TREAT" ]]; then
    echo "  ERROR: Treatment file $TREAT not found"
    exit 1
  fi
  if [[ ! -f "$CTRL" ]]; then
    echo "  ERROR: Control file $CTRL not found"
    exit 1
  fi
  
  echo "  Processing: $TREAT vs $CTRL (genome: $GENOME)"

  # genome flag for macs2
  case "$GENOME" in
    hg38|GRCh38) GFLAG="hs" ;;
    mm10|GRCm38) GFLAG="mm" ;;
    *)           GFLAG="" ; echo "  warn: unknown genome $GENOME (no -g)";;
  esac

  OUTDIR="${OUTROOT}/${GENOME}"
  mkdir -p "$OUTDIR"
  PREFIX="$NAME"

  if [[ "$PEAKMODE" == "narrow" ]]; then
    macs2 callpeak -t "$TREAT" -c "$CTRL" ${GFLAG:+-g "$GFLAG"} \
      -n "$PREFIX" --outdir "$OUTDIR" -B --SPMR --call-summits \
      -q "$QVAL" --keep-dup "$KEEPDUP" 2> "${OUTDIR}/${PREFIX}.log"
  else
    macs2 callpeak -t "$TREAT" -c "$CTRL" ${GFLAG:+-g "$GFLAG"} \
      -n "$PREFIX" --outdir "$OUTDIR" -B --SPMR --broad --broad-cutoff "$BROAD_QVAL" \
      --keep-dup "$KEEPDUP" 2> "${OUTDIR}/${PREFIX}.log"
  fi
  if [[ $? -eq 0 ]]; then
    echo "  ✅ Success: ${OUTDIR}/${PREFIX}_peaks.*"
  else
    echo "  ❌ MACS2 failed for $NAME"
    exit 1
  fi
done

# --- Step 3: bedGraph -> bigWig (per genome if chrom.sizes available) ---
echo "[3/3] bedGraph → bigWig"
if [[ -x "$BEDGRAPH2BIGWIG" ]]; then
  # hg38
  if [[ -n "$CHROMSIZES_HG38" && -f "$CHROMSIZES_HG38" && -d "${OUTROOT}/hg38" ]]; then
    for B in "${OUTROOT}/hg38/"*_treat_pileup.bdg; do
      [[ -e "$B" ]] || continue
      "$BEDGRAPH2BIGWIG" "$B" "$CHROMSIZES_HG38" "${B%_treat_pileup.bdg}.bw"
    done
  fi
  # mm10
  if [[ -n "$CHROMSIZES_MM10" && -f "$CHROMSIZES_MM10" && -d "${OUTROOT}/mm10" ]]; then
    for B in "${OUTROOT}/mm10/"*_treat_pileup.bdg; do
      [[ -e "$B" ]] || continue
      "$BEDGRAPH2BIGWIG" "$B" "$CHROMSIZES_MM10" "${B%_treat_pileup.bdg}.bw"
    done
  fi
else
  echo "  bedGraphToBigWig not found; skipping."
fi

echo "✅ Done."
