#!/bin/bash
#SBATCH --job-name=chec_filter_bam
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH -o logs/filter_bam_%a.log
#SBATCH -e logs/filter_bam_%a.log
#SBATCH --container=docker://quay.io/biocontainers/samtools:1.21--h50ea8bc_0

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Independent of the dmel spike-in branch (01a/01b) - only needs the
# S. cerevisiae alignment from 01_align.sh. Runs whether or not
# --align_dmel is used; this is a separate quantification pathway.
#
# Two calling conventions, same pattern as 01b_dmel_normalized_coverage.sh:
#
# 1. LOOKUP MODE (array job, per-sample):
#      sbatch --array=1-N 01c_filter_bam.sh <lookup_file> <include_regions.bed> [bam_type] [--keep_both_reads]
#    bam_type: "nuclear" (default) or "full" - same convention as 02/03/01b.
#
# 2. DIRECT MODE (single job, no sample sheet - e.g. a manually-combined
#    control BAM, same situation maketagdir_control.sh handles):
#      sbatch 01c_filter_bam.sh --bam <bam_file> --include-regions <bed_file> \
#          --output-name <name> [--keep_both_reads]
#
# Filters to properly-paired (or, with --keep_both_reads, still requires
# properly paired but keeps both mates), MAPQ>=10 reads overlapping the
# given include_regions.bed. Default keeps R1 only, since R1's 5' end is
# what 01d_genomecov_5p.sh treats as the cut site - keeping R2 as well would
# double-count/misrepresent cut sites unless you specifically want both
# mates for some other purpose.
#
# Output filename is derived from the INPUT BAM's own basename (not
# reconstructed from regulator/replicate), so a bam_type-specific suffix
# like "_nuclear" is preserved rather than lost, e.g.:
#   RME1_A_nuclear.bam -> RME1_A_nuclear_filtered.bam
#   RME1_A.bam         -> RME1_A_filtered.bam
# 01d_genomecov_5p.sh needs to know the same bam_type to reconstruct this
# exact filename in lookup mode - see that script's usage.
OUTPUT_DIR="results"
LOG_DIR="logs"
KEEP_BOTH_READS=false

# ============================================================================
# MODE DETECTION AND ARGUMENT PARSING
# ============================================================================
if [[ "${1:-}" == --bam || "${1:-}" == --include-regions || "${1:-}" == --output-name ]]; then
    MODE="direct"
else
    MODE="lookup"
fi

if [[ "${MODE}" == "direct" ]]; then
    DIRECT_BAM=""
    DIRECT_INCLUDE_REGIONS=""
    DIRECT_OUTPUT_NAME=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --bam)               DIRECT_BAM="$2"; shift 2 ;;
            --include-regions)    DIRECT_INCLUDE_REGIONS="$2"; shift 2 ;;
            --output-name)        DIRECT_OUTPUT_NAME="$2"; shift 2 ;;
            --keep_both_reads)    KEEP_BOTH_READS=true; shift ;;
            *)
                echo "ERROR: Unrecognized argument in direct mode: $1"
                exit 1
                ;;
        esac
    done

    if [[ -z "${DIRECT_BAM}" || -z "${DIRECT_INCLUDE_REGIONS}" || -z "${DIRECT_OUTPUT_NAME}" ]]; then
        echo "ERROR: Direct mode requires --bam, --include-regions, and --output-name"
        echo "Usage: 01c_filter_bam.sh --bam <bam_file> --include-regions <bed_file> --output-name <name> [--keep_both_reads]"
        exit 1
    fi

    INPUT_BAM="${DIRECT_BAM}"
    INCLUDE_REGIONS="${DIRECT_INCLUDE_REGIONS}"
    SAMPLE_LABEL="${DIRECT_OUTPUT_NAME}"
    OUT_DIR="${OUTPUT_DIR}/filtered_bams/${DIRECT_OUTPUT_NAME}"
    INPUT_BASENAME="$(basename "${INPUT_BAM}" .bam)"
    OUTPUT_BAM="${OUT_DIR}/${INPUT_BASENAME}_filtered.bam"

else
    LOOKUP_FILE="$1"
    INCLUDE_REGIONS="$2"
    shift 2
    BAM_TYPE="nuclear"

    # Optional bam_type positional, then flags
    if [[ $# -gt 0 && "$1" != --* ]]; then
        BAM_TYPE="$1"
        shift
    fi
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --keep_both_reads) KEEP_BOTH_READS=true; shift ;;
            *)
                echo "ERROR: Unrecognized argument: $1"
                exit 1
                ;;
        esac
    done

    case "${BAM_TYPE}" in
        nuclear) BAM_SUFFIX="_nuclear.bam" ;;
        full)    BAM_SUFFIX=".bam" ;;
        *)
            echo "ERROR: Invalid bam_type '${BAM_TYPE}' - must be 'nuclear' or 'full'"
            exit 1
            ;;
    esac

    LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")
    if [[ -z "$LINE" ]]; then
        echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
        exit 1
    fi
    read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

    INPUT_BAM="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}${BAM_SUFFIX}"
    SAMPLE_LABEL="${REGULATOR}_${REPLICATE}"
    OUT_DIR="${OUTPUT_DIR}/filtered_bams/${REGULATOR}/${REPLICATE}"
    INPUT_BASENAME="$(basename "${INPUT_BAM}" .bam)"
    OUTPUT_BAM="${OUT_DIR}/${INPUT_BASENAME}_filtered.bam"
fi

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM not found: ${INPUT_BAM}"
    exit 1
fi

if [[ ! -f "${INCLUDE_REGIONS}" ]]; then
    echo "ERROR: Include regions file not found: ${INCLUDE_REGIONS}"
    exit 1
fi

echo "Filtering ${SAMPLE_LABEL}"
echo "  Input BAM: ${INPUT_BAM}"
echo "  Include regions: ${INCLUDE_REGIONS}"
echo "  Keep both reads: ${KEEP_BOTH_READS}"
echo "  Output: ${OUTPUT_BAM}"

# ============================================================================
# FILTER, SORT, INDEX
# ============================================================================
# -F 0x090C: exclude unmapped (0x4), mate unmapped (0x8), secondary (0x100), supplementary (0x800)
# -f 0x0002: require properly paired
# -f 0x0040: require R1 only (default; omitted with --keep_both_reads)
FILTER_FLAGS=0x090C
if [[ "${KEEP_BOTH_READS}" == "true" ]]; then
    REQUIRED_FLAGS=0x0002
else
    REQUIRED_FLAGS=0x0042   # properly paired (0x2) + read1 (0x40)
fi

samtools view -h -b \
    -F "${FILTER_FLAGS}" \
    -f "${REQUIRED_FLAGS}" \
    -q 10 \
    -L "${INCLUDE_REGIONS}" \
    "${INPUT_BAM}" | \
samtools sort -@ "${SLURM_CPUS_PER_TASK:-1}" -o "${OUTPUT_BAM}" -

samtools index -@ "${SLURM_CPUS_PER_TASK:-1}" "${OUTPUT_BAM}"

# ============================================================================
# SUMMARY
# ============================================================================
TOTAL=$(samtools view -c "${INPUT_BAM}")
FILTERED=$(samtools view -c "${OUTPUT_BAM}")
PCT=$(awk "BEGIN {printf \"%.2f\", ($FILTERED/$TOTAL)*100}")

echo ""
echo "Total reads:    ${TOTAL}"
echo "Filtered reads: ${FILTERED}"
echo "Retained:       ${PCT}%"
echo ""
echo "✓ Filter complete: ${OUTPUT_BAM}"
