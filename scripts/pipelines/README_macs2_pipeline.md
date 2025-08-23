# MACS2 ChIP-seq Peak Calling Pipeline

A comprehensive SLURM-compatible pipeline for automated MACS2 peak calling with BAM file renaming and BigWig generation.

## Features

- **Automated BAM renaming**: Converts ID-based filenames to sample names
- **Dual genome support**: Handles hg38 and mm10 samples simultaneously
- **Input validation**: Ensures each genome has appropriate control samples
- **Peak calling modes**: Supports both narrow and broad peak calling
- **BigWig generation**: Automatic conversion from bedGraph to BigWig
- **Error handling**: Robust error checking with informative messages
- **SLURM integration**: Ready for cluster submission

## Requirements

### Software Dependencies
- MACS2 (tested with version 2.2.7.1)
- bedGraphToBigWig (UCSC tools)
- Bash 4.0+

### Environment Setup
- Conda environment with MACS2 installed
- Access to chromosome size files (hg38.chrom.sizes, mm10.chrom.sizes)
- SLURM cluster (optional, can run locally)

## Input Format

### Sample Mapping File (`list.txt`)
Tab-separated file with 3 columns:
```
ID          sample_name    genome
26238X1     SUCCSFPQ      hg38
26238X2     SUCCNONO      hg38
26238X3     SUCCinput     hg38
26238X4     3777SFPQ      mm10
26238X5     3777NONO      mm10
26238X6     3777input     mm10
```

**Requirements:**
- Input/control samples must contain "input" in the sample_name (case-insensitive)
- Each genome must have exactly one input sample
- BAM files should be named as `{ID}.bam` with corresponding `.bam.bai` or `.bai` index files

### Directory Structure
```
your_data_directory/
├── list.txt
├── 26238X1.bam
├── 26238X1.bam.bai
├── 26238X2.bam
├── 26238X2.bam.bai
├── ...
└── 26238X6.bam.bai
```

## Usage

### Basic SLURM Submission
```bash
cd /path/to/your/data
sbatch /path/to/macs2_pipeline.sh
```

### Local Execution
```bash
cd /path/to/your/data
# Load required modules
module use $HOME/MyModules
module load miniconda3/latest
source activate bioinfo_env

# Run pipeline
bash /path/to/macs2_pipeline.sh
```

### Custom Parameters
```bash
# Narrow peaks with stricter q-value
QVAL=0.0001 sbatch macs2_pipeline.sh

# Broad peak calling
PEAKMODE=broad BROAD_QVAL=0.05 sbatch macs2_pipeline.sh

# Custom output directory
OUTROOT=my_peaks sbatch macs2_pipeline.sh

# Different mapping file
MAPPING_FILE=samples.tsv sbatch macs2_pipeline.sh

# Skip BigWig generation
BEDGRAPH2BIGWIG="" sbatch macs2_pipeline.sh
```

## Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `MAPPING_FILE` | `list.txt` | Sample mapping file |
| `OUTROOT` | `macs2` | Output root directory |
| `PEAKMODE` | `narrow` | Peak calling mode (`narrow` or `broad`) |
| `QVAL` | `0.001` | Q-value cutoff for narrow peaks |
| `BROAD_QVAL` | `0.01` | Q-value cutoff for broad peaks |
| `KEEPDUP` | `auto` | MACS2 duplicate handling |
| `BEDGRAPH2BIGWIG` | `/uufs/chpc.utah.edu/...` | Path to bedGraphToBigWig |
| `CHROMSIZES_HG38` | `/uufs/chpc.utah.edu/...` | hg38 chromosome sizes |
| `CHROMSIZES_MM10` | `/uufs/chpc.utah.edu/...` | mm10 chromosome sizes |

### SLURM Parameters
```bash
#SBATCH --job-name="macs2_peaks"
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=your-partition
#SBATCH --account=your-account
```

## Output Structure

```
macs2/
├── hg38/
│   ├── SUCCSFPQ_peaks.narrowPeak
│   ├── SUCCSFPQ_peaks.xls
│   ├── SUCCSFPQ_summits.bed
│   ├── SUCCSFPQ_treat_pileup.bdg
│   ├── SUCCSFPQ.bw
│   ├── SUCCNONO_peaks.narrowPeak
│   ├── SUCCNONO.bw
│   └── ...
└── mm10/
    ├── 3777SFPQ_peaks.narrowPeak
    ├── 3777SFPQ.bw
    ├── 3777NONO_peaks.narrowPeak
    ├── 3777NONO.bw
    └── ...
```

### Output Files Per Sample
- `*_peaks.narrowPeak` or `*_peaks.broadPeak`: Peak locations
- `*_peaks.xls`: Detailed peak information
- `*_summits.bed`: Peak summits (narrow mode only)
- `*_treat_pileup.bdg`: Coverage bedGraph
- `*.bw`: BigWig coverage track
- `*.log`: MACS2 log file

## Pipeline Steps

1. **BAM Renaming**: Renames BAM files from IDs to sample names
2. **Input Validation**: Verifies input controls exist for each genome
3. **Peak Calling**: Runs MACS2 with appropriate parameters
4. **BigWig Generation**: Converts bedGraph files to BigWig format

## Error Handling

The pipeline includes comprehensive error checking:
- Validates input file existence
- Ensures each genome has input controls
- Checks MACS2 execution success
- Verifies file operations complete successfully

## Troubleshooting

### Common Issues

**"No input control found for genome X"**
- Ensure your input samples contain "input" in the sample_name
- Verify genome names match exactly between treatment and input samples

**"Treatment file X.bam not found"**
- Check that BAM files exist and follow the naming convention
- Verify the mapping file uses the correct IDs

**"MACS2 failed for sample X"**
- Check the individual log files in the output directory
- Verify BAM files are properly indexed
- Ensure sufficient memory allocation

### Performance Tuning
- Adjust `--cpus-per-task` based on available resources
- Increase `--mem-per-cpu` for large BAM files
- Use `KEEPDUP=1` to retain all reads for shallow sequencing

## Example Workflow

```bash
# 1. Prepare your data
cd /scratch/my_chipseq_data
ls *.bam *.bai list.txt

# 2. Submit job
sbatch /path/to/macs2_pipeline.sh

# 3. Monitor progress
squeue -u $USER
tail -f slurm-JOBID.out

# 4. Check results
ls macs2/*/
```

## License

This pipeline is part of the bioinfo_tools repository.