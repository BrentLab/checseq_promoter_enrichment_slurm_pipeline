#!/bin/bash
#SBATCH --job-name=chec_findpeaks
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH -o logs/findpeaks_%a.log
#SBATCH -e logs/findpeaks_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/homer_samtools:0e83b23821fcb7e6

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Usage: 03_findpeaks.sh <lookup_file> [bam_type]
#   bam_type: "nuclear" (default) or "full" - must match whatever BAM_TYPE
#   was used to build the sample's tag directory in 02_maketagdir_samples.sh,
#   since the genome size below is the statistical background denominator
#   for HOMER's Poisson model and needs to reflect the actual sequence space
#   the tag directory was built over.
LOOKUP_FILE="$1"
BAM_TYPE="${2:-nuclear}"
OUTPUT_DIR="results"
LOG_DIR="logs"
CONTROL_TAG_DIR="${OUTPUT_DIR}/tag_dirs/control_MNase"

# sacCer3 genome size, chosen based on bam_type:
#   nuclear -> 12071326 (sum of chrI-chrXVI only, excludes chrM/85779bp)
#              230218+813184+316620+1531933+576874+270161+1090940+562643
#              +439888+745751+666816+1078177+924431+784333+1091291+948066
#   full    -> 12157105 (chrI-chrXVI + chrM, i.e. 12071326 + 85779)
case "${BAM_TYPE}" in
    nuclear) GENOME_SIZE=12071326 ;;
    full)    GENOME_SIZE=12157105 ;;
    *)
        echo "ERROR: Invalid bam_type '${BAM_TYPE}' - must be 'nuclear' or 'full'"
        exit 1
        ;;
esac

# HOMER findPeaks parameters (from methods)
# -C 0: disable clonal filtering (appropriate for MNase-treated data)
# -L 6: 6-fold enrichment over local background (vs. default 4-fold)
# -F 10: 10-fold enrichment over control (vs. default 4-fold)
# -gsize: explicit genome size (nuclear- or full-genome, per bam_type above),
#         rather than each sample re-estimating its own slightly different
#         value from tag coverage
FINDPEAKS_PARAMS="-o auto -C 0 -L 6 -F 10 -gsize ${GENOME_SIZE}"

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

SAMPLE_TAG_DIR="${OUTPUT_DIR}/tag_dirs/${REGULATOR}_${REPLICATE}"
OUTPUT_PEAK_FILE="${OUTPUT_DIR}/peaks/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_peaks.txt"

# Create sample-specific peak directory (must come after OUTPUT_PEAK_FILE is set)
PEAK_DIR="$(dirname "${OUTPUT_PEAK_FILE}")"
mkdir -p "${PEAK_DIR}"

echo "Task ${SLURM_ARRAY_TASK_ID}: Calling peaks for ${REGULATOR}_${REPLICATE}"
echo "  BAM type: ${BAM_TYPE} (genome size: ${GENOME_SIZE})"
echo "  Sample tag directory: ${SAMPLE_TAG_DIR}"
echo "  Control tag directory: ${CONTROL_TAG_DIR}"
echo "  Output: ${OUTPUT_PEAK_FILE}"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================
if [[ ! -d "${SAMPLE_TAG_DIR}" ]]; then
    echo "ERROR: Sample tag directory not found: ${SAMPLE_TAG_DIR}"
    echo "Make sure 02_maketagdir_samples.sh completed successfully"
    exit 1
fi

if [[ ! -d "${CONTROL_TAG_DIR}" ]]; then
    echo "ERROR: Control tag directory not found: ${CONTROL_TAG_DIR}"
    echo "Make sure maketagdir_control.sh has been run (standalone, manual step)"
    exit 1
fi

# ============================================================================
# CALL PEAKS
# ============================================================================
echo "Running findPeaks with style=factor..."
echo "  Parameters: ${FINDPEAKS_PARAMS}"

findPeaks "${SAMPLE_TAG_DIR}" \
    -style factor \
    -i "${CONTROL_TAG_DIR}" \
    ${FINDPEAKS_PARAMS} \
    2>&1 | tee "${LOG_DIR}/findpeaks_${REGULATOR}_${REPLICATE}.log"

# ============================================================================
# VALIDATE OUTPUT
# ============================================================================
if [[ ! -f "${SAMPLE_TAG_DIR}/peaks.txt" ]]; then
    echo "ERROR: findPeaks did not produce output"
    exit 1
fi

# Copy peak file to output directory for organization
cp "${SAMPLE_TAG_DIR}/peaks.txt" "${OUTPUT_PEAK_FILE}"

# Extract summary statistics from peak file header
echo ""
echo "Peak calling summary for ${REGULATOR}_${REPLICATE}:"
grep "^# total peaks" "${OUTPUT_PEAK_FILE}" || true
grep "^# peak size" "${OUTPUT_PEAK_FILE}" || true
grep "^# fragment length" "${OUTPUT_PEAK_FILE}" || true
grep "^# Approximate IP efficiency" "${OUTPUT_PEAK_FILE}" || true
grep "^# number of putative peaks" "${OUTPUT_PEAK_FILE}" || true
grep "^# Putative peaks filtered by input" "${OUTPUT_PEAK_FILE}" || true
grep "^# Putative peaks filtered by local signal" "${OUTPUT_PEAK_FILE}" || true
grep "^# Fold over input required" "${OUTPUT_PEAK_FILE}" || true
grep "^# Fold over local region required" "${OUTPUT_PEAK_FILE}" || true

echo ""
echo "✓ Peak file ready: ${OUTPUT_PEAK_FILE}"
