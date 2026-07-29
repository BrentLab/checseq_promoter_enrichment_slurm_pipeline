#!/bin/bash
#SBATCH --job-name=chec_tagdir_control
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH -o logs/tagdir_control.log
#SBATCH -e logs/tagdir_control.log
#SBATCH --container=oras://community.wave.seqera.io/library/homer_samtools:0e83b23821fcb7e6

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# This is a STANDALONE script, run manually before the automated pipeline
# (submit_pipeline.sh). It is not part of the automated dependency chain,
# since the control tag directory is typically built once and reused across
# pipeline runs rather than regenerated every time.
#
# Usage: maketagdir_control.sh <control_bam> <control_stats_file>
#   <control_stats_file> is the samtools_stats.txt produced by 01_align.sh for
#   the control sample. Fragment length is derived at runtime from its
#   "insert size average" line, since HOMER's own autocorrelation-based
#   estimate is unreliable for a background/input control (see QC discussion).
#
# Run this after 01_align.sh has completed for the control sample(s), and
# before running submit_pipeline.sh (which expects results/tag_dirs/control_MNase
# to already exist -- see 03_findpeaks.sh, which errors clearly if it's missing).
CONTROL_BAM="${1:-/path/to/free_MNase_control.bam}"
CONTROL_STATS_FILE="${2:?ERROR: must provide path to control samtools_stats.txt as \$2}"
GENOME_FASTA="/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529_chr_normalized.fa"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}/tag_dirs" "${LOG_DIR}"

CONTROL_TAG_DIR="${OUTPUT_DIR}/tag_dirs/control_MNase"

echo "Creating HOMER tag directory for control MNase dataset..."
echo "  Input BAM: ${CONTROL_BAM}"
echo "  Stats file: ${CONTROL_STATS_FILE}"
echo "  Output tag directory: ${CONTROL_TAG_DIR}"

# Validate control BAM exists and is indexed
if [[ ! -f "${CONTROL_BAM}" ]]; then
    echo "ERROR: Control BAM file not found: ${CONTROL_BAM}"
    exit 1
fi

if [[ ! -f "${CONTROL_BAM}.bai" ]]; then
    echo "Indexing control BAM..."
    samtools index "${CONTROL_BAM}"
fi

# ============================================================================
# DETERMINE FRAGMENT LENGTH FROM SAMTOOLS STATS
# ============================================================================
if [[ ! -f "${CONTROL_STATS_FILE}" ]]; then
    echo "ERROR: Control stats file not found: ${CONTROL_STATS_FILE}"
    echo "  This file is produced by 01_align.sh - make sure alignment completed first."
    exit 1
fi

FRAGMENT_LENGTH=$(grep "^SN.*insert size average" "${CONTROL_STATS_FILE}" | awk '{print int($NF + 0.5)}')

if [[ -z "${FRAGMENT_LENGTH}" ]] || [[ "${FRAGMENT_LENGTH}" -lt 1 ]]; then
    echo "ERROR: Could not determine fragment length from ${CONTROL_STATS_FILE}"
    exit 1
fi

echo "  Fragment length (from control insert size average): ${FRAGMENT_LENGTH} bp"

# ============================================================================
# MAKE TAG DIRECTORY
# ============================================================================
# Convert BAM to SAM (makeTagDirectory expects SAM format)
CONTROL_SAM="${OUTPUT_DIR}/tag_dirs/control_MNase_temp.sam"
samtools view -h "${CONTROL_BAM}" > "${CONTROL_SAM}"

# -format sam       : explicit format (avoids relying on autodetection)
# -unique -mapq 10   : HOMER's default filtering, made explicit for clarity
# -fragLength        : override HOMER's autocorrelation estimate with the
#                      directly measured insert size average, since HOMER
#                      could not produce a reliable peak-size estimate on
#                      its own for this background control
# -genome -checkGC   : run sequence bias (GC%) diagnostics - important to
#                      check on a background control, since a GC-shifted
#                      library would look like fake "enrichment" at GC-rich
#                      regions (e.g. TSS/CpG-like regions) downstream
makeTagDirectory "${CONTROL_TAG_DIR}" "${CONTROL_SAM}" \
    -format sam -unique -mapq 10 \
    -fragLength "${FRAGMENT_LENGTH}" \
    -genome "${GENOME_FASTA}" -checkGC \
    2>&1 | tee "${LOG_DIR}/tagdir_control_makeTagDirectory.log"

# Clean up temporary SAM
rm -f "${CONTROL_SAM}"

# ============================================================================
# VALIDATE OUTPUT
# ============================================================================
if [[ ! -d "${CONTROL_TAG_DIR}" ]]; then
    echo "ERROR: Failed to create tag directory"
    exit 1
fi

# Report summary from tag directory
echo ""
echo "Tag directory created successfully:"
ls -lh "${CONTROL_TAG_DIR}" | head -n 5

if [[ -f "${CONTROL_TAG_DIR}/tagInfo.txt" ]]; then
    echo ""
    echo "Tag directory statistics:"
    grep "^genome" "${CONTROL_TAG_DIR}/tagInfo.txt" || true
    grep "fragmentLengthEstimate=" "${CONTROL_TAG_DIR}/tagInfo.txt" || true
fi

echo ""
echo "✓ Control tag directory ready: ${CONTROL_TAG_DIR}"
