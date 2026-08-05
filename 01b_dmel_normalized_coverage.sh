#!/bin/bash
#SBATCH --job-name=chec_dmel_cov
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH -o logs/dmel_cov_%a.log
#SBATCH -e logs/dmel_cov_%a.log
#SBATCH --container=docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Two calling conventions:
#
# 1. LOOKUP MODE (array job, per-sample - the normal pipeline path):
#      sbatch --array=1-N 01b_dmel_normalized_coverage.sh <lookup_file> [bam_type]
#    bam_type: "nuclear" (default) or "full"
#      nuclear -> {regulator}_{replicate}_nuclear.bam (chrM filtered out)
#      full    -> {regulator}_{replicate}.bam         (all chromosomes, incl. chrM)
#    Derives everything (input BAM, dmel counts file, output paths) from the
#    lookup row at SLURM_ARRAY_TASK_ID, same as the other per-sample scripts.
#
# 2. DIRECT MODE (single job, no sample sheet - e.g. a manually-combined
#    replicate BAM, same situation maketagdir_control.sh handles for tag
#    directories):
#      sbatch 01b_dmel_normalized_coverage.sh --bam <bam_file> \
#          --dmel-counts <count> --output-name <name>
#    Writes to results/coverage/<name>/<name>_dmel_norm.bedgraph. Use this
#    when there's no regulator/replicate row to derive paths from - e.g. a
#    combined_freemnase.bam built by merging replicate BAMs by hand, the
#    same way maketagdir_control.sh's CONTROL_BAM is built.
#    <count> is a plain integer - the total dmel-mapped reads to normalize
#    against (e.g. summed across the replicates you combined into the BAM).
#    If you have 01a_map_to_dmel.sh's *_dmel_counts.txt file(s) instead of
#    just the number, extract it yourself, e.g.:
#      awk -F'\t' '{sum+=$NF} END{print sum}' rep1_dmel_counts.txt rep2_dmel_counts.txt
#
# Computes per-base genome coverage from the S. cerevisiae alignment,
# normalized to the sample's D. melanogaster spike-in read count:
#
#     normalized_coverage(bp) = raw_coverage(bp) / dmel_mapped_reads * 10000
#
# This follows the coverage normalization described in Mahendrawada et al.
# 2024 (and Donczew & Hahn 2018, elifesciences.org/articles/50109, whose
# methods this paper's own methods cite and reuse): "the number of reads
# that mapped at that position divided by the number of all D. melanogaster
# reads mapped in the sample and multiplied by 10000".
#
# Neither methods reference specifies exactly what "coverage at a position"
# means for paired-end data - whole-fragment depth, per-read footprint depth,
# or 5'-end/cut-site depth are all consistent with the wording. Given that
# ambiguity, this uses `bedtools genomecov -pc` (paired-end fragment
# coverage: counts each base between R1's start and R2's end as covered by
# that fragment, rather than only the two read footprints with a gap between
# them). This is the standard, well-defined choice for paired-end coverage
# tracks generally, and doesn't require guessing a fragment-length estimate
# the way single-end extension-based methods (HOMER/MACS) do, since the true
# fragment span is directly observed from both mates. NOTE: this is a
# distinct, different convention from the -5 (5'-end/cut-site) coverage used
# in promoter_enrichment/genomecov.sh elsewhere in this repo - that script
# and this one deliberately answer different questions (cut-site density at
# specific loci vs. genome-wide fragment-level occupancy).
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# MODE DETECTION AND ARGUMENT PARSING
# ============================================================================
if [[ "${1:-}" == --bam || "${1:-}" == --dmel-counts || "${1:-}" == --output-name ]]; then
    MODE="direct"
else
    MODE="lookup"
fi

if [[ "${MODE}" == "direct" ]]; then
    DIRECT_BAM=""
    DIRECT_DMEL_COUNT=""
    DIRECT_OUTPUT_NAME=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --bam)
                DIRECT_BAM="$2"; shift 2 ;;
            --dmel-counts)
                DIRECT_DMEL_COUNT="$2"; shift 2 ;;
            --output-name)
                DIRECT_OUTPUT_NAME="$2"; shift 2 ;;
            *)
                echo "ERROR: Unrecognized argument in direct mode: $1"
                exit 1
                ;;
        esac
    done

    if [[ -z "${DIRECT_BAM}" || -z "${DIRECT_DMEL_COUNT}" || -z "${DIRECT_OUTPUT_NAME}" ]]; then
        echo "ERROR: Direct mode requires --bam, --dmel-counts, and --output-name"
        echo "Usage: 01b_dmel_normalized_coverage.sh --bam <bam_file> --dmel-counts <count> --output-name <name>"
        exit 1
    fi

    if ! [[ "${DIRECT_DMEL_COUNT}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: --dmel-counts must be a plain non-negative integer (got: '${DIRECT_DMEL_COUNT}')"
        exit 1
    fi

    INPUT_BAM="${DIRECT_BAM}"
    DMEL_COUNT="${DIRECT_DMEL_COUNT}"
    SAMPLE_LABEL="${DIRECT_OUTPUT_NAME}"

    COV_DIR="${OUTPUT_DIR}/coverage/${DIRECT_OUTPUT_NAME}"
    RAW_BEDGRAPH="${COV_DIR}/${DIRECT_OUTPUT_NAME}_raw.bedgraph"
    NORM_BEDGRAPH="${COV_DIR}/${DIRECT_OUTPUT_NAME}_dmel_norm.bedgraph"

    echo "Direct mode: Computing dmel-normalized coverage for ${SAMPLE_LABEL}"
    echo "  Input BAM: ${INPUT_BAM}"
    echo "  Output: ${NORM_BEDGRAPH}"

else
    LOOKUP_FILE="$1"
    BAM_TYPE="${2:-nuclear}"

    case "${BAM_TYPE}" in
        nuclear) BAM_SUFFIX="_nuclear.bam" ;;
        full)    BAM_SUFFIX=".bam" ;;
        *)
            echo "ERROR: Invalid bam_type '${BAM_TYPE}' - must be 'nuclear' or 'full'"
            exit 1
            ;;
    esac

    # Get the task line from lookup file (skip header)
    # Since SLURM_ARRAY_TASK_ID starts at 1, add 1 to skip the header line
    LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

    if [[ -z "$LINE" ]]; then
        echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
        exit 1
    fi

    # Parse TSV: regulator_symbol, replicate, fastq_1, fastq_2
    read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

    SAMPLE_DIR="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}"
    INPUT_BAM="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}${BAM_SUFFIX}"
    DMEL_COUNTS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel_counts.txt"
    SAMPLE_LABEL="${REGULATOR}_${REPLICATE}"

    COV_DIR="${OUTPUT_DIR}/coverage/${REGULATOR}/${REPLICATE}"
    RAW_BEDGRAPH="${COV_DIR}/${REGULATOR}_${REPLICATE}_raw.bedgraph"
    NORM_BEDGRAPH="${COV_DIR}/${REGULATOR}_${REPLICATE}_dmel_norm.bedgraph"

    echo "Task ${SLURM_ARRAY_TASK_ID}: Computing dmel-normalized coverage for ${SAMPLE_LABEL}"
    echo "  BAM type: ${BAM_TYPE}"
    echo "  Input BAM: ${INPUT_BAM}"
    echo "  Dmel counts file: ${DMEL_COUNTS_FILE}"
    echo "  Output: ${NORM_BEDGRAPH}"
fi

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${COV_DIR}" "${LOG_DIR}"

# ============================================================================
# VALIDATE INPUT
# ============================================================================
if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM not found: ${INPUT_BAM}"
    exit 1
fi

if [[ "${MODE}" == "lookup" ]]; then
    if [[ ! -f "${DMEL_COUNTS_FILE}" ]]; then
        echo "ERROR: Dmel counts file not found: ${DMEL_COUNTS_FILE}"
        exit 1
    fi

    # Read count from the last tab-separated field - works for both a bare
    # single-number file and the standard regulator/replicate/count format
    DMEL_COUNT=$(awk -F'\t' '{print $NF}' "${DMEL_COUNTS_FILE}")

    if ! [[ "${DMEL_COUNT}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Could not parse a numeric dmel count from ${DMEL_COUNTS_FILE} (got: '${DMEL_COUNT}')"
        exit 1
    fi
fi

# DMEL_COUNT is set either above (lookup mode) or directly from --dmel-counts
# during argument parsing (direct mode) - validated as numeric in both cases.
if [[ "${DMEL_COUNT}" -eq 0 ]]; then
    echo "ERROR: Dmel mapped read count is 0 for ${SAMPLE_LABEL} - cannot normalize"
    echo "  (no dmel spike-in reads recovered; check the upstream dmel alignment for this sample)"
    exit 1
fi

echo "  Dmel mapped reads (normalization denominator): ${DMEL_COUNT}"

# ============================================================================
# GENOME-WIDE COVERAGE (paired-end fragment depth, all positions incl. zero)
# ============================================================================
# -bga: report all positions, including zero-coverage regions
# -pc:  paired-end fragment coverage - count each base from R1's start
#       through R2's end as covered (the whole sequenced fragment), not just
#       the two read footprints with an uncovered gap between them
echo "Computing raw genome coverage (paired-end fragment depth)..."
bedtools genomecov -ibam "${INPUT_BAM}" -bga -pc > "${RAW_BEDGRAPH}"

# ============================================================================
# NORMALIZE: coverage / dmel_count * 10000
# ============================================================================
echo "Normalizing to dmel spike-in (x 10000)..."
awk -v dmel="${DMEL_COUNT}" 'BEGIN{OFS="\t"} {
    printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, ($4 / dmel) * 10000
}' "${RAW_BEDGRAPH}" > "${NORM_BEDGRAPH}"

# Raw bedgraph is an intermediate; keep it only for QC/debugging purposes.
# Comment out the next line if you'd rather keep it around.
rm -f "${RAW_BEDGRAPH}"

# ============================================================================
# VALIDATE OUTPUT
# ============================================================================
if [[ ! -s "${NORM_BEDGRAPH}" ]]; then
    echo "ERROR: Normalized bedgraph is empty or missing: ${NORM_BEDGRAPH}"
    exit 1
fi

NUM_INTERVALS=$(wc -l < "${NORM_BEDGRAPH}")

echo ""
echo "✓ Dmel-normalized coverage complete for ${SAMPLE_LABEL}"
echo "  Intervals: ${NUM_INTERVALS}"
echo "  Output: ${NORM_BEDGRAPH}"
