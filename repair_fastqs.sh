#!/bin/bash
#SBATCH --job-name=repair_fastq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o logs/repair_%a.log
#SBATCH -e logs/repair_%a.log
#SBATCH --container=docker://quay.io/staphb/bbtools:latest

set -uo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Usage: repair_fastq_pairs.sh <lookup_file>
#   TSV with header: regulator_symbol, replicate, fastq_1, fastq_2
#
# NOTE: repair.sh only fixes genuinely-paired-but-desynced FASTQs (same
# sequencing run, mismatched read counts). If fastq_1/fastq_2 in a row come
# from two DIFFERENT SRX/SRR accessions, this will not produce a meaningful
# fix - read names won't overlap between unrelated runs, so output will be
# tiny or empty. That pattern indicates a lookup-table pairing bug upstream,
# not something repair.sh can solve. Check the read counts reported below.
LOOKUP_FILE="$1"
OUTPUT_DIR="repaired_fastq"
LOG_DIR="logs"

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

# Get the task line from lookup file (skip header)
# Since SLURM_ARRAY_TASK_ID starts at 1, add 1 to skip the header line
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

if [[ -z "$LINE" ]]; then
    echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    exit 1
fi

# Parse TSV: regulator_symbol, replicate, fastq_1, fastq_2
read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: Repairing pair for ${REGULATOR}_${REPLICATE}"
echo "  R1: ${FASTQ_R1}"
echo "  R2: ${FASTQ_R2}"

# ============================================================================
# VALIDATE INPUT
# ============================================================================
if [[ ! -f "$FASTQ_R1" ]]; then
    echo "ERROR: FASTQ R1 not found: $FASTQ_R1"
    exit 1
fi

if [[ ! -f "$FASTQ_R2" ]]; then
    echo "ERROR: FASTQ R2 not found: $FASTQ_R2"
    exit 1
fi

# Output filenames keep the original basenames, written into repaired_fastq/
R1_BASENAME="$(basename "${FASTQ_R1}")"
R2_BASENAME="$(basename "${FASTQ_R2}")"
OUT_R1="${OUTPUT_DIR}/${R1_BASENAME}"
OUT_R2="${OUTPUT_DIR}/${R2_BASENAME}"
OUT_SINGLETONS="${OUTPUT_DIR}/${REGULATOR}_${REPLICATE}_singletons.fastq.gz"

# ============================================================================
# REPORT INPUT READ COUNTS (before repair, for comparison)
# ============================================================================
IN_R1_COUNT=$(( $(zcat "${FASTQ_R1}" | wc -l) / 4 ))
IN_R2_COUNT=$(( $(zcat "${FASTQ_R2}" | wc -l) / 4 ))
echo "  Input read counts: R1=${IN_R1_COUNT}  R2=${IN_R2_COUNT}"

# ============================================================================
# RUN repair.sh
# ============================================================================
echo "Running repair.sh..."
repair.sh \
    in1="${FASTQ_R1}" \
    in2="${FASTQ_R2}" \
    out1="${OUT_R1}" \
    out2="${OUT_R2}" \
    outs="${OUT_SINGLETONS}" \
    overwrite=true

REPAIR_EXIT=$?
if [[ $REPAIR_EXIT -ne 0 ]]; then
    echo "ERROR: repair.sh failed with exit code $REPAIR_EXIT"
    exit 1
fi

# ============================================================================
# VALIDATE AND REPORT OUTPUT
# ============================================================================
if [[ ! -f "${OUT_R1}" ]] || [[ ! -f "${OUT_R2}" ]]; then
    echo "ERROR: repair.sh did not produce expected output files"
    exit 1
fi

OUT_R1_COUNT=$(( $(zcat "${OUT_R1}" | wc -l) / 4 ))
OUT_R2_COUNT=$(( $(zcat "${OUT_R2}" | wc -l) / 4 ))
SINGLETON_COUNT=0
if [[ -f "${OUT_SINGLETONS}" ]]; then
    SINGLETON_COUNT=$(( $(zcat "${OUT_SINGLETONS}" | wc -l) / 4 ))
fi

echo ""
echo "Repair summary for ${REGULATOR}_${REPLICATE}:"
echo "  Output read counts: R1=${OUT_R1_COUNT}  R2=${OUT_R2_COUNT}  singletons=${SINGLETON_COUNT}"

# Flag the case that indicates a mispaired (not just desynced) input: if the
# repaired output retains only a tiny fraction of the smaller input file's
# reads, R1/R2 likely came from unrelated sequencing runs, not a fixable desync.
if [[ "${IN_R1_COUNT}" -gt 0 ]]; then
    RETAINED_PCT=$(awk "BEGIN {printf \"%.1f\", (${OUT_R1_COUNT} / ${IN_R1_COUNT}) * 100}")
    echo "  Retained ${RETAINED_PCT}% of R1's original reads"
    if (( $(awk "BEGIN {print (${RETAINED_PCT} < 10)}") )); then
        echo "  WARNING: <10% of reads retained - R1/R2 may be from unrelated"
        echo "           sequencing runs (mispaired in the lookup table), not a"
        echo "           genuine desync. Check the source accessions for this row."
    fi
fi

echo ""
echo "✓ Repair complete"
echo "  R1: ${OUT_R1}"
echo "  R2: ${OUT_R2}"
echo "  Singletons: ${OUT_SINGLETONS}"
