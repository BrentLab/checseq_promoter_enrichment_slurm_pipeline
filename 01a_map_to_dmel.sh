#!/bin/bash
#SBATCH --job-name=chec_map_dmel
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o logs/map_dmel_%a.log
#SBATCH -e logs/map_dmel_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/bowtie2_samtools_gzip:0c7fd2d5085a1394

set -uo pipefail

# ============================================================================
# CONFIGURATION - Fill in these paths before running
# ============================================================================
# Two calling conventions:
#
# 1. LOOKUP MODE (array job, per-sample - the normal pipeline path):
#      sbatch --array=1-N 01a_map_to_dmel.sh <lookup_file>
#    Derives the unmapped-reads FASTQs (produced by 01_align.sh) and all
#    output paths from the lookup row at SLURM_ARRAY_TASK_ID.
#
# 2. DIRECT MODE (single job, no sample sheet - e.g. manually-combined
#    replicate FASTQs, same situation 01b_dmel_normalized_coverage.sh and
#    maketagdir_control.sh handle for their own inputs):
#      sbatch 01a_map_to_dmel.sh --fastq-r1 <r1.fastq.gz> \
#          --fastq-r2 <r2.fastq.gz> --output-name <name>
#    Writes to results/bams/<name>/<name>_dmel.bam (+ stats/flagstats/
#    coverage/idxstats/counts alongside it, all named <name>_dmel_*). Use
#    this when there's no regulator/replicate row to derive paths from -
#    e.g. unmapped R1/R2 FASTQs you've concatenated across replicates by
#    hand.
DMEL_BOWTIE_INDEX="/ref/mblab/data/dmelanogaster/bowtie2_index/dmel-all-chromosome-r6.65"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# BOWTIE PARAMETERS (same as S. cerevisiae)
# ============================================================================
BOWTIE_PARAMS="-I 10 -X 700 --local --very-sensitive-local --no-unal --no-mixed --no-discordant"

# ============================================================================
# MODE DETECTION AND ARGUMENT PARSING
# ============================================================================
if [[ "${1:-}" == --fastq-r1 || "${1:-}" == --fastq-r2 || "${1:-}" == --output-name ]]; then
    MODE="direct"
else
    MODE="lookup"
fi

if [[ "${MODE}" == "direct" ]]; then
    DIRECT_R1=""
    DIRECT_R2=""
    DIRECT_OUTPUT_NAME=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --fastq-r1)
                DIRECT_R1="$2"; shift 2 ;;
            --fastq-r2)
                DIRECT_R2="$2"; shift 2 ;;
            --output-name)
                DIRECT_OUTPUT_NAME="$2"; shift 2 ;;
            *)
                echo "ERROR: Unrecognized argument in direct mode: $1"
                exit 1
                ;;
        esac
    done

    if [[ -z "${DIRECT_R1}" || -z "${DIRECT_R2}" || -z "${DIRECT_OUTPUT_NAME}" ]]; then
        echo "ERROR: Direct mode requires --fastq-r1, --fastq-r2, and --output-name"
        echo "Usage: 01a_map_to_dmel.sh --fastq-r1 <r1.fastq.gz> --fastq-r2 <r2.fastq.gz> --output-name <name>"
        exit 1
    fi

    UNMAPPED_R1="${DIRECT_R1}"
    UNMAPPED_R2="${DIRECT_R2}"
    SAMPLE_LABEL="${DIRECT_OUTPUT_NAME}"

    SAMPLE_DIR="${OUTPUT_DIR}/bams/${DIRECT_OUTPUT_NAME}"
    DMEL_BAM="${SAMPLE_DIR}/${DIRECT_OUTPUT_NAME}_dmel.bam"
    DMEL_COUNTS_FILE="${SAMPLE_DIR}/${DIRECT_OUTPUT_NAME}_dmel_counts.txt"

    echo "Direct mode: Mapping unmapped reads to D. melanogaster for ${SAMPLE_LABEL}"
    echo "  Unmapped R1: ${UNMAPPED_R1}"
    echo "  Unmapped R2: ${UNMAPPED_R2}"

else
    LOOKUP_FILE="${1:-samples.tsv}"

    # Get the task line from lookup file (skip header)
    LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

    if [[ -z "$LINE" ]]; then
        echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
        exit 1
    fi

    # Parse TSV: regulator_symbol, replicate, fastq_1, fastq_2
    read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

    # All dmel outputs live alongside the other BAMs in bams/{reg}/{rep}/
    SAMPLE_DIR="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}"
    UNMAPPED_R1="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_unmapped_R1.fastq.gz"
    UNMAPPED_R2="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_unmapped_R2.fastq.gz"
    DMEL_BAM="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel.bam"
    DMEL_COUNTS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel_counts.txt"
    SAMPLE_LABEL="${REGULATOR}_${REPLICATE}"

    echo "Task ${SLURM_ARRAY_TASK_ID}: Mapping unmapped reads to D. melanogaster for ${SAMPLE_LABEL}"
    echo "  Unmapped R1: ${UNMAPPED_R1}"
    echo "  Unmapped R2: ${UNMAPPED_R2}"
fi

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${SAMPLE_DIR}" "${LOG_DIR}"

# ============================================================================
# VALIDATE INPUT
# ============================================================================
if [[ ! -f "${UNMAPPED_R1}" ]] || [[ ! -f "${UNMAPPED_R2}" ]]; then
    echo "WARNING: Unmapped FASTQ files not found"
    [[ ! -f "${UNMAPPED_R1}" ]] && echo "  Missing: ${UNMAPPED_R1}"
    [[ ! -f "${UNMAPPED_R2}" ]] && echo "  Missing: ${UNMAPPED_R2}"
    echo "  Skipping D. melanogaster mapping for ${SAMPLE_LABEL}"
    echo "0" > "${DMEL_COUNTS_FILE}"
    exit 0
fi

# Count unmapped reads
UNMAPPED_COUNT=$(gunzip -c "${UNMAPPED_R1}" 2>/dev/null | grep -c "^@" || echo "0")

if [[ "$UNMAPPED_COUNT" -eq 0 ]]; then
    echo "No unmapped reads found, skipping D. melanogaster alignment"
    echo "0" > "${DMEL_COUNTS_FILE}"
    exit 0
fi

echo "  Unmapped read pairs: ${UNMAPPED_COUNT}"

# ============================================================================
# ALIGN TO D. MELANOGASTER
# ============================================================================
echo "Aligning unmapped reads to D. melanogaster..."
bowtie2 \
    -p 8 \
    -q \
    -x "${DMEL_BOWTIE_INDEX}" \
    ${BOWTIE_PARAMS} \
    -1 <(gunzip -c "${UNMAPPED_R1}") \
    -2 <(gunzip -c "${UNMAPPED_R2}") \
    2> "${LOG_DIR}/${SAMPLE_LABEL}_dmel_bowtie.log" | \
    samtools view -b -h -F 4 - | \
    samtools sort -@ 4 -o "${DMEL_BAM}" -

samtools index "${DMEL_BAM}"

# ============================================================================
# GENERATE STATISTICS
# ============================================================================
DMEL_STATS_FILE="${SAMPLE_DIR}/${SAMPLE_LABEL}_dmel_samtools_stats.txt"
DMEL_FLAGSTATS_FILE="${SAMPLE_DIR}/${SAMPLE_LABEL}_dmel_samtools_flagstats.txt"
DMEL_COVERAGE_FILE="${SAMPLE_DIR}/${SAMPLE_LABEL}_dmel_coverage.txt"
DMEL_IDXSTATS_FILE="${SAMPLE_DIR}/${SAMPLE_LABEL}_dmel_idxstats.txt"

echo "Generating D. melanogaster alignment statistics..."
samtools stats "${DMEL_BAM}" > "${DMEL_STATS_FILE}"
samtools flagstats "${DMEL_BAM}" > "${DMEL_FLAGSTATS_FILE}"
samtools coverage "${DMEL_BAM}" > "${DMEL_COVERAGE_FILE}"
samtools idxstats "${DMEL_BAM}" > "${DMEL_IDXSTATS_FILE}"

# ============================================================================
# COUNT D. MELANOGASTER READS
# ============================================================================
DMEL_MAPPED_COUNT=$(samtools view -c "${DMEL_BAM}")

echo ""
echo "D. melanogaster mapping results for ${SAMPLE_LABEL}:"
echo "  Unmapped S. cerevisiae reads: ${UNMAPPED_COUNT}"
echo "  Mapped to D. melanogaster: ${DMEL_MAPPED_COUNT}"
echo ""

# Write count to file for use in downstream spike-in normalization.
# Kept as a 3-column regulator/replicate/count row in lookup mode (matching
# the original format 01b_dmel_normalized_coverage.sh expects); in direct
# mode there's no regulator/replicate, so SAMPLE_LABEL is written into both
# of those fields - 01b's dmel-count parser reads the LAST tab-separated
# field regardless, so this is compatible either way.
if [[ "${MODE}" == "direct" ]]; then
    echo -e "${SAMPLE_LABEL}\t${SAMPLE_LABEL}\t${DMEL_MAPPED_COUNT}" > "${DMEL_COUNTS_FILE}"
else
    echo -e "${REGULATOR}\t${REPLICATE}\t${DMEL_MAPPED_COUNT}" > "${DMEL_COUNTS_FILE}"
fi

echo "✓ D. melanogaster mapping complete"
echo "  BAM: ${DMEL_BAM}"
echo "  Stats: ${DMEL_STATS_FILE}"
echo "  Flagstats: ${DMEL_FLAGSTATS_FILE}"
echo "  Coverage: ${DMEL_COVERAGE_FILE}"
echo "  Idxstats: ${DMEL_IDXSTATS_FILE}"
echo "  Counts file: ${DMEL_COUNTS_FILE}"
