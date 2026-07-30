#!/bin/bash
#SBATCH --time=0-00:30:00
#SBATCH --mem=4GB
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/filter_bam_%a.out
#SBATCH --error=logs/filter_bam_%a.err

set -euo pipefail

# Usage: sbatch --array=1-N filter_bam.sh lookup.txt output_dir include_regions.bed [--single_end] [--keep_both_reads]
#
# By default, assumes paired-end data and keeps R1 reads only (0x40)
# Use --keep_both_reads to retain both R1 and R2 reads
# Use --single_end flag for single-end data processing

eval $(spack load --sh samtools)

# Default: paired-end mode, R1 only
SINGLE_END=false
KEEP_BOTH_READS=false

# Parse required arguments
if [[ $# -lt 3 ]]; then
    echo "Error: Insufficient arguments" >&2
    echo "Usage: $0 lookup.txt output_dir include_regions.bed [--single_end]" >&2
    exit 1
fi

LOOKUP=$1
OUTPUT_DIR=$2
INCLUDE_REGIONS=$3
shift 3

# Parse optional flags
while [[ $# -gt 0 ]]; do
    case $1 in
        --single_end)
            SINGLE_END=true
            shift
            ;;
        --keep_both_reads)
            KEEP_BOTH_READS=true
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Valid options: --single_end, --keep_both_reads" >&2
            exit 1
            ;;
    esac
done

# ==========================================
# Validate inputs
# ==========================================
if [[ ! -f "$INCLUDE_REGIONS" ]]; then
    echo "Error: Include regions file not found: $INCLUDE_REGIONS" >&2
    exit 1
fi

if ! read -r bam_file < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$LOOKUP"); then
    echo "Error: Unable to read BAM file from lookup at line ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

if [[ ! -f "$bam_file" ]]; then
    echo "Error: BAM file not found: $bam_file" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR" logs

BASENAME="$(basename "$bam_file" .bam)"
OUTPUT_BAM="$OUTPUT_DIR/${BASENAME}.filtered.bam"

# ==========================================
# Configure filtering based on data type
# ==========================================
if [[ "$SINGLE_END" == "true" ]]; then
    # Single-end: basic quality filtering only
    # -F 0x904: exclude unmapped (0x4), secondary (0x100), supplementary (0x800)
    FILTER_FLAGS=0x904
    REQUIRED_FLAGS=""
    DATA_TYPE="single-end"
    echo "=========================================="
    echo "Processing single-end data"
    echo "BAM file: $bam_file"
    echo "=========================================="
else
    # Paired-end: require proper pairing
    # -F 0x090C: exclude unmapped (0x4), mate unmapped (0x8), secondary (0x100), supplementary (0x800)
    # -f 0x0002: require properly paired
    # -f 0x0040: require R1 only (default; omitted with --keep_both_reads)
    FILTER_FLAGS=0x090C
    if [[ "$KEEP_BOTH_READS" == "true" ]]; then
        REQUIRED_FLAGS=0x0002
        DATA_TYPE="paired-end (both reads)"
        echo "=========================================="
        echo "Processing paired-end data"
        echo "BAM file: $bam_file"
        echo "Mode: Keeping both R1 and R2 reads"
        echo "=========================================="
    else
        # 0x0042 = properly paired (0x2) + read 1 (0x40)
        REQUIRED_FLAGS=0x0042
        DATA_TYPE="paired-end (R1 only)"
        echo "=========================================="
        echo "Processing paired-end data"
        echo "BAM file: $bam_file"
        echo "Mode: R1 only (use --keep_both_reads to retain R2)"
        echo "=========================================="
    fi
fi

# ==========================================
# Filter, sort, and index
# ==========================================
echo "Filtering $DATA_TYPE reads..."

if [[ -n "$REQUIRED_FLAGS" ]]; then
    # Paired-end mode
    samtools view -h -b \
        -F "$FILTER_FLAGS" \
        -f "$REQUIRED_FLAGS" \
        -q 10 \
        -L "$INCLUDE_REGIONS" \
        "$bam_file" | \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "$OUTPUT_BAM" -
else
    # Single-end mode
    samtools view -h -b \
        -F "$FILTER_FLAGS" \
        -q 10 \
        -L "$INCLUDE_REGIONS" \
        "$bam_file" | \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "$OUTPUT_BAM" -
fi

echo "Indexing: $OUTPUT_BAM"
samtools index -@ ${SLURM_CPUS_PER_TASK} "$OUTPUT_BAM"

# ==========================================
# Generate statistics
# ==========================================
TOTAL=$(samtools view -c "$bam_file")
FILTERED=$(samtools view -c "$OUTPUT_BAM")
PCT=$(awk "BEGIN {printf \"%.2f\", ($FILTERED/$TOTAL)*100}")

echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Data type:        $DATA_TYPE"
echo "Input:            $bam_file"
echo "Output:           $OUTPUT_BAM"
echo "Include regions:  $INCLUDE_REGIONS"
echo ""
echo "Total reads:      $TOTAL"
echo "Filtered reads:   $FILTERED"
echo "Retained:         ${PCT}%"
echo "=========================================="
echo "✓ Complete!"
