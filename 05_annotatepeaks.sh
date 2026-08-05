#!/bin/bash
#SBATCH --job-name=chec_annotatePeaks
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o logs/annotatePeaks_%a.log
#SBATCH -e logs/annotatePeaks_%a.log
#SBATCH --container=docker://quay.io/biocontainers/homer:5.1--pl5321hc52dbad_1

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Usage: 05_annotatepeaks.sh <lookup_file>
LOOKUP_FILE="$1"
GENOME_FASTA="/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529_chr_normalized.fa"
GTF_FILE="/ref/mblab/data/yeast_data/reprocess_mahendrawada/mahendrawada_slurm_pipeline/sacCer3.ensGene.gtf"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================
if [[ ! -f "${GENOME_FASTA}" ]]; then
    echo "ERROR: Genome FASTA not found: ${GENOME_FASTA}"
    exit 1
fi

if [[ ! -f "${GTF_FILE}" ]]; then
    echo "ERROR: GTF file not found: ${GTF_FILE}"
    exit 1
fi

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}/peaks" "${LOG_DIR}"

# Get the task line from lookup file (skip header)
# Since SLURM_ARRAY_TASK_ID starts at 1, add 1 to skip the header line
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

if [[ -z "$LINE" ]]; then
    echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    exit 1
fi

# Parse TSV: regulator_symbol, replicate, fastq_1, fastq_2
read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

# HOMER native peak file from 03_findpeaks.sh (annotatePeaks.pl accepts this
# directly - no need to go through the pos2bed.sh BED conversion first)
PEAK_FILE="${OUTPUT_DIR}/peaks/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_peaks.txt"

ANNOT_DIR="${OUTPUT_DIR}/peaks/${REGULATOR}/${REPLICATE}"
OUTPUT_FILE="${ANNOT_DIR}/${REGULATOR}_${REPLICATE}_annotated_peaks.txt"
ERR_FILE="${ANNOT_DIR}/${REGULATOR}_${REPLICATE}_annotatePeaks.err"

echo "Task ${SLURM_ARRAY_TASK_ID}: Annotating peaks for ${REGULATOR}_${REPLICATE}"
echo "  Peak file: ${PEAK_FILE}"
echo "  Genome FASTA: ${GENOME_FASTA}"
echo "  GTF: ${GTF_FILE}"
echo "  Output: ${OUTPUT_FILE}"

# ============================================================================
# VALIDATE INPUT
# ============================================================================
if [[ ! -f "${PEAK_FILE}" ]]; then
    echo "ERROR: Peak file not found: ${PEAK_FILE}"
    echo "Make sure 03_findpeaks.sh completed successfully"
    exit 1
fi

mkdir -p "${ANNOT_DIR}"

# ============================================================================
# RUN annotatePeaks.pl
# ============================================================================
echo "Running annotatePeaks.pl..."
annotatePeaks.pl \
    "${PEAK_FILE}" \
    "${GENOME_FASTA}" \
    -gtf "${GTF_FILE}" \
    > "${OUTPUT_FILE}" 2> "${ERR_FILE}"

# ============================================================================
# VALIDATE OUTPUT
# ============================================================================
if [[ ! -s "${OUTPUT_FILE}" ]]; then
    echo "ERROR: annotatePeaks.pl produced an empty or missing output file"
    echo "Check stderr log: ${ERR_FILE}"
    cat "${ERR_FILE}"
    exit 1
fi

PEAK_COUNT=$(tail -n +2 "${OUTPUT_FILE}" | wc -l)

echo ""
echo "✓ Annotation complete for ${REGULATOR}_${REPLICATE}"
echo "  Annotated peaks: ${PEAK_COUNT}"
echo "  Output: ${OUTPUT_FILE}"
