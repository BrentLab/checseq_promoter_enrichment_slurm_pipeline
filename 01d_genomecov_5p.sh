#!/bin/bash
#SBATCH --job-name=chec_genomecov_5p
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=00:30:00
#SBATCH -o logs/genomecov_5p_%a.log
#SBATCH -e logs/genomecov_5p_%a.log
#SBATCH --container=docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Runs on the output of 01c_filter_bam.sh. Computes per-base, per-strand 5'
# read-end ("cut site") coverage - the standard ChEC-seq/MNase cut-site
# quantification, distinct from 01b_dmel_normalized_coverage.sh's whole-
# fragment coverage (see that script's header for the -pc vs -5 discussion).
#
# Two calling conventions, same pattern as the other 01x scripts:
#
# 1. LOOKUP MODE (array job, per-sample):
#      sbatch --array=1-N 01d_genomecov_5p.sh <lookup_file> [bam_type]
#    bam_type: "nuclear" (default) or "full" - must match whatever bam_type
#    01c_filter_bam.sh was run with for this sample, since the filtered
#    BAM's filename is derived from the INPUT bam's own basename (e.g.
#    RME1_A_nuclear_filtered.bam vs RME1_A_filtered.bam) - this script needs
#    to know which one to look for.
#
# 2. DIRECT MODE (single job, no sample sheet - e.g. a manually-combined
#    control BAM, already filtered by 01c_filter_bam.sh in direct mode):
#      sbatch 01d_genomecov_5p.sh --bam <filtered_bam_file> --output-name <name>
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# MODE DETECTION AND ARGUMENT PARSING
# ============================================================================
if [[ "${1:-}" == --bam || "${1:-}" == --output-name ]]; then
    MODE="direct"
else
    MODE="lookup"
fi

if [[ "${MODE}" == "direct" ]]; then
    DIRECT_BAM=""
    DIRECT_OUTPUT_NAME=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --bam)          DIRECT_BAM="$2"; shift 2 ;;
            --output-name)  DIRECT_OUTPUT_NAME="$2"; shift 2 ;;
            *)
                echo "ERROR: Unrecognized argument in direct mode: $1"
                exit 1
                ;;
        esac
    done

    if [[ -z "${DIRECT_BAM}" || -z "${DIRECT_OUTPUT_NAME}" ]]; then
        echo "ERROR: Direct mode requires --bam and --output-name"
        echo "Usage: 01d_genomecov_5p.sh --bam <filtered_bam_file> --output-name <name>"
        exit 1
    fi

    INPUT_BAM="${DIRECT_BAM}"
    SAMPLE_LABEL="${DIRECT_OUTPUT_NAME}"
    COV_DIR="${OUTPUT_DIR}/genomecov_5p/${DIRECT_OUTPUT_NAME}"
    OUTPUT_BED="${COV_DIR}/${DIRECT_OUTPUT_NAME}_r1_5p.bed"

else
    LOOKUP_FILE="$1"
    BAM_TYPE="${2:-nuclear}"

    case "${BAM_TYPE}" in
        nuclear) BAM_SUFFIX="_nuclear" ;;
        full)    BAM_SUFFIX="" ;;
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

    # Matches the filename 01c_filter_bam.sh produces: derived from the
    # input BAM's own basename (${REGULATOR}_${REPLICATE}${BAM_SUFFIX}),
    # not reconstructed independently.
    INPUT_BAM="${OUTPUT_DIR}/filtered_bams/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}${BAM_SUFFIX}_filtered.bam"
    SAMPLE_LABEL="${REGULATOR}_${REPLICATE}"
    COV_DIR="${OUTPUT_DIR}/genomecov_5p/${REGULATOR}/${REPLICATE}"
    OUTPUT_BED="${COV_DIR}/${REGULATOR}_${REPLICATE}_r1_5p.bed"
fi

mkdir -p "${COV_DIR}" "${LOG_DIR}"

if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM not found: ${INPUT_BAM}"
    echo "Make sure 01c_filter_bam.sh completed successfully for this sample"
    exit 1
fi

echo "Computing 5' cut-site coverage for ${SAMPLE_LABEL}"
echo "  Input BAM: ${INPUT_BAM}"
echo "  Output: ${OUTPUT_BED}"

TMP_DIR=$(mktemp -d -t genomecov_XXXXXX)
trap "rm -rf ${TMP_DIR}" EXIT

echo "Generating plus strand coverage..."
bedtools genomecov -ibam "${INPUT_BAM}" -5 -dz -strand + > "${TMP_DIR}/plus.cov"

echo "Generating minus strand coverage..."
bedtools genomecov -ibam "${INPUT_BAM}" -5 -dz -strand - > "${TMP_DIR}/minus.cov"

awk 'OFS="\t" {print $1, $2, $2+1, ".", $3, "+"}' "${TMP_DIR}/plus.cov" > "${TMP_DIR}/plus.bed"
awk 'OFS="\t" {print $1, $2, $2+1, ".", $3, "-"}' "${TMP_DIR}/minus.cov" > "${TMP_DIR}/minus.bed"

cat "${TMP_DIR}/plus.bed" "${TMP_DIR}/minus.bed" > "${TMP_DIR}/combined_unsorted.bed"
bedtools sort -i "${TMP_DIR}/combined_unsorted.bed" > "${OUTPUT_BED}"

PLUS_COUNT=$(wc -l < "${TMP_DIR}/plus.bed")
MINUS_COUNT=$(wc -l < "${TMP_DIR}/minus.bed")
TOTAL_COUNT=$(wc -l < "${OUTPUT_BED}")

echo ""
echo "Plus strand positions:  ${PLUS_COUNT}"
echo "Minus strand positions: ${MINUS_COUNT}"
echo "Total positions:        ${TOTAL_COUNT}"
echo ""
echo "✓ 5' cut-site coverage complete: ${OUTPUT_BED}"
