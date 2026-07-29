#!/bin/bash
#SBATCH --job-name=chec_pos2bed
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/pos2bed_%a.log
#SBATCH --error=logs/pos2bed_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/homer_samtools:0e83b23821fcb7e6

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
LOOKUP_FILE="$1"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}/peaks" "${LOG_DIR}"

# Create sample-specific peak directory
PEAK_DIR="${OUTPUT_DIR}/peaks"
mkdir -p "${PEAK_DIR}"

# Get the task line from lookup file (skip header)
# Since SLURM_ARRAY_TASK_ID starts at 1, add 1 to skip the header line
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

if [[ -z "$LINE" ]]; then
    echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    exit 1
fi

# Parse TSV: regulator_symbol, replicate, fastq_1, fastq_2
read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

INPUT_PEAK_FILE="${OUTPUT_DIR}/peaks/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_peaks.txt"
OUTPUT_BED_FILE="${OUTPUT_DIR}/peaks/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_peaks.bed"
PEAK_SUMMIT_FILE="${OUTPUT_DIR}/peaks/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_peaks_summits.txt"

echo "Task ${SLURM_ARRAY_TASK_ID}: Converting peaks to BED for ${REGULATOR}_${REPLICATE}"
echo "  Input: ${INPUT_PEAK_FILE}"
echo "  Output BED: ${OUTPUT_BED_FILE}"

# ============================================================================
# VALIDATE INPUT
# ============================================================================
if [[ ! -f "${INPUT_PEAK_FILE}" ]]; then
    echo "ERROR: Peak file not found: ${INPUT_PEAK_FILE}"
    echo "Make sure 03_findpeaks.sh completed successfully"
    exit 1
fi

# ============================================================================
# CONVERT TO BED FORMAT
# ============================================================================
# pos2bed.pl converts HOMER peak format to standard BED format
echo "Running pos2bed.pl..."
pos2bed.pl "${INPUT_PEAK_FILE}" > "${OUTPUT_BED_FILE}"

# ============================================================================
# EXTRACT PEAK SUMMITS
# ============================================================================
# Calculate peak summit as mid-range between peak borders
# HOMER peak file format: PeakID chr start end strand ...
# BED format: chr start end name score strand
echo "Extracting peak summits..."
awk 'NR > 1 && !/^#/ {
    # Skip header lines and comments
    chr = $2
    start = $3
    end = $4
    strand = $5
    peakid = $1
    # Summit is midpoint between start and end
    summit = int((start + end) / 2)
    # Output: chr summit summit name score strand
    print chr "\t" summit "\t" (summit + 1) "\t" peakid "\t0\t" strand
}' "${INPUT_PEAK_FILE}" > "${PEAK_SUMMIT_FILE}"

# ============================================================================
# VALIDATE OUTPUT
# ============================================================================
if [[ ! -f "${OUTPUT_BED_FILE}" ]]; then
    echo "ERROR: BED conversion failed"
    exit 1
fi

NUM_PEAKS=$(grep -v "^#" "${OUTPUT_BED_FILE}" | wc -l)
NUM_SUMMITS=$(grep -v "^#" "${PEAK_SUMMIT_FILE}" | wc -l)

echo ""
echo "Conversion summary for ${REGULATOR}_${REPLICATE}:"
echo "  Peaks in BED: ${NUM_PEAKS}"
echo "  Peak summits: ${NUM_SUMMITS}"
echo ""
echo "First few peaks:"
head -n 5 "${OUTPUT_BED_FILE}"

echo ""
echo "✓ BED conversion complete"
echo "  Peaks: ${OUTPUT_BED_FILE}"
echo "  Summits: ${PEAK_SUMMIT_FILE}"
