#!/bin/bash
#SBATCH --job-name=chec_findpeaks
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=06:00:00
#SBATCH -o logs/findpeaks_%a.log
#SBATCH -e logs/findpeaks_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/homer_samtools:0e83b23821fcb7e6

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Usage: 03_findpeaks.sh <lookup_file> [bam_type] [--control-tag-dir=<path>] [--genome-size=<n>]
#   bam_type: "nuclear" (default) or "full" - informational only now (used in
#   the task summary echo below); no longer affects genome size, since that's
#   now passed explicitly via --genome-size (or omitted) rather than looked
#   up from a nuclear/full constant.
#   --control-tag-dir=<path>: path to the control tag directory built by
#   maketagdir_control.sh (default: results/tag_dirs/control_MNase, matching
#   that script's own default output location). Can appear anywhere in the
#   args, alongside the positional lookup_file/bam_type.
#   --genome-size=<n>: explicit genome size passed to findPeaks' -gsize flag
#   (the statistical background denominator for HOMER's Poisson model).
#   Optional - if omitted, -gsize is not passed to findPeaks at all, and
#   HOMER auto-estimates genome size from the tag directory itself instead
#   (this also matches the original Mahendrawada et al. findPeaks calls,
#   which never passed -gsize). Can appear anywhere in the args.
OUTPUT_DIR="results"
LOG_DIR="logs"
CONTROL_TAG_DIR="${OUTPUT_DIR}/tag_dirs/control_MNase"
GENOME_SIZE=""

POSITIONAL=()
for arg in "$@"; do
    case "${arg}" in
        --control-tag-dir=*)
            CONTROL_TAG_DIR="${arg#--control-tag-dir=}"
            ;;
        --genome-size=*)
            GENOME_SIZE="${arg#--genome-size=}"
            ;;
        *)
            POSITIONAL+=("${arg}")
            ;;
    esac
done

LOOKUP_FILE="${POSITIONAL[0]:?ERROR: lookup_file is required}"
BAM_TYPE="${POSITIONAL[1]:-nuclear}"

case "${BAM_TYPE}" in
    nuclear|full) ;;
    *)
        echo "ERROR: Invalid bam_type '${BAM_TYPE}' - must be 'nuclear' or 'full'"
        exit 1
        ;;
esac

# HOMER findPeaks parameters (from methods)
# -C 0: disable clonal filtering (appropriate for MNase-treated data)
# -L 6: 6-fold enrichment over local background (vs. default 4-fold)
# -F 10: 10-fold enrichment over control (vs. default 4-fold)
# -gsize: only included if --genome-size was provided; otherwise omitted
#         entirely, letting findPeaks auto-estimate genome size from the tag
#         directory instead.
if [[ -n "${GENOME_SIZE}" ]]; then
    FINDPEAKS_PARAMS="-o auto -C 0 -L 6 -F 10 -gsize ${GENOME_SIZE}"
else
    FINDPEAKS_PARAMS="-o auto -C 0 -L 6 -F 10"
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

SAMPLE_TAG_DIR="${OUTPUT_DIR}/tag_dirs/${REGULATOR}/${REPLICATE}"
OUTPUT_PEAK_FILE="${OUTPUT_DIR}/peaks/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_peaks.txt"

# Create sample-specific peak directory (must come after OUTPUT_PEAK_FILE is set)
PEAK_DIR="$(dirname "${OUTPUT_PEAK_FILE}")"
mkdir -p "${PEAK_DIR}"

echo "Task ${SLURM_ARRAY_TASK_ID}: Calling peaks for ${REGULATOR}_${REPLICATE}"
echo "  BAM type: ${BAM_TYPE}"
echo "  Genome size: $( [[ -n "${GENOME_SIZE}" ]] && echo "${GENOME_SIZE}" || echo "auto-estimated by HOMER (--genome-size not provided)" )"
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
