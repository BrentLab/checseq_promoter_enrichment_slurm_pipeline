#!/bin/bash
#SBATCH --job-name=chec_tagdir_samples
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH -o logs/tagdir_samples_%a.log
#SBATCH -e logs/tagdir_samples_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/homer_samtools:0e83b23821fcb7e6

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Usage: 02_maketagdir_samples.sh <lookup_file> [bam_type] [--authors_orig]
#   bam_type: "nuclear" (default) or "full"
#     nuclear -> {regulator}_{replicate}_nuclear.bam (chrM filtered out)
#     full    -> {regulator}_{replicate}.bam         (all chromosomes, incl. chrM)
#   --authors_orig: use -keepAll instead of -unique -mapq 10, matching the
#     original Mahendrawada et al. makeTagDirectory calls (which kept all
#     alignments, including multi-mappers/low-MAPQ reads, rather than HOMER's
#     default unique/MAPQ>=10 filtering). Can appear anywhere in the args.
#     ALSO skips passing -fragLength entirely (see below), matching the
#     original scripts, which never passed -fragLength and let HOMER's own
#     autocorrelation-based estimate run instead.
#   Fragment length: by default, derived per-sample at runtime from that
#   sample's own samtools_stats.txt "insert size average" (produced by
#   01_align.sh), rather than a single shared value -- each regulator's
#   ChEC-seq fragment size can genuinely differ from the free MNase control.
#   With --authors_orig, this derivation is skipped entirely and -fragLength
#   is not passed to makeTagDirectory at all, matching the original scripts'
#   behavior of relying on HOMER's own autocorrelation estimate instead.

# Pull --authors_orig out of the argument list wherever it appears, leaving
# the remaining positional args (lookup_file, bam_type) in order.
AUTHORS_ORIG=false
POSITIONAL=()
for arg in "$@"; do
    if [[ "${arg}" == "--authors_orig" ]]; then
        AUTHORS_ORIG=true
    else
        POSITIONAL+=("${arg}")
    fi
done

LOOKUP_FILE="${POSITIONAL[0]:?ERROR: lookup_file is required}"
BAM_TYPE="${POSITIONAL[1]:-nuclear}"
GENOME_FASTA="/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529_chr_normalized.fa"
OUTPUT_DIR="results"
LOG_DIR="logs"

if [[ "${AUTHORS_ORIG}" == "true" ]]; then
    READ_FILTER_FLAGS="-keepAll"
else
    READ_FILTER_FLAGS="-unique -mapq 10"
fi

case "${BAM_TYPE}" in
    nuclear) BAM_SUFFIX="_nuclear.bam" ;;
    full)    BAM_SUFFIX=".bam" ;;
    *)
        echo "ERROR: Invalid bam_type '${BAM_TYPE}' - must be 'nuclear' or 'full'"
        exit 1
        ;;
esac

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}/tag_dirs" "${LOG_DIR}" "${OUTPUT_DIR}/peaks"

# Get the task line from lookup file (skip header)
# Since SLURM_ARRAY_TASK_ID starts at 1, add 1 to skip the header line
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

if [[ -z "$LINE" ]]; then
    echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    exit 1
fi

# Parse TSV: regulator_symbol, replicate, fastq_1, fastq_2
read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

# Corresponding BAM from 01_align.sh (full or nuclear-only, per BAM_TYPE above)
INPUT_BAM="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}${BAM_SUFFIX}"
STATS_FILE="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_samtools_stats.txt"
TAG_DIR="${OUTPUT_DIR}/tag_dirs/${REGULATOR}/${REPLICATE}"

echo "Task ${SLURM_ARRAY_TASK_ID}: Creating tag directory for ${REGULATOR}_${REPLICATE}"
echo "  BAM type: ${BAM_TYPE}"
echo "  Read filtering: ${READ_FILTER_FLAGS}$( [[ "${AUTHORS_ORIG}" == "true" ]] && echo " (--authors_orig mode)" )"
echo "  Input BAM: ${INPUT_BAM}"
echo "  Output tag directory: ${TAG_DIR}"

# ============================================================================
# VALIDATE INPUT
# ============================================================================
if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM not found: ${INPUT_BAM}"
    echo "Make sure 01_align.sh completed successfully"
    exit 1
fi

if [[ ! -f "${INPUT_BAM}.bai" ]]; then
    echo "Indexing input BAM..."
    samtools index "${INPUT_BAM}"
fi

# ============================================================================
# DETERMINE FRAGMENT LENGTH FROM THIS SAMPLE'S OWN SAMTOOLS STATS
# (skipped entirely with --authors_orig - see header)
# ============================================================================
if [[ "${AUTHORS_ORIG}" == "true" ]]; then
    FRAGLENGTH_ARG=()
    echo "  Fragment length: letting HOMER auto-estimate (--authors_orig mode, skipping samtools stats)"
else
    if [[ ! -f "${STATS_FILE}" ]]; then
        echo "ERROR: Stats file not found: ${STATS_FILE}"
        echo "  This file is produced by 01_align.sh - make sure alignment completed first."
        exit 1
    fi

    FRAGMENT_LENGTH=$(grep "^SN.*insert size average" "${STATS_FILE}" | awk '{print int($NF + 0.5)}')

    if [[ -z "${FRAGMENT_LENGTH}" ]] || [[ "${FRAGMENT_LENGTH}" -lt 1 ]]; then
        echo "ERROR: Could not determine fragment length from ${STATS_FILE}"
        exit 1
    fi

    echo "  Fragment length (from this sample's insert size average): ${FRAGMENT_LENGTH} bp"
    FRAGLENGTH_ARG=(-fragLength "${FRAGMENT_LENGTH}")
fi

# ============================================================================
# MAKE TAG DIRECTORY
# ============================================================================
# Convert BAM to SAM (makeTagDirectory expects SAM format)
TEMP_SAM="${TAG_DIR}_temp.sam"
mkdir -p "$(dirname "${TAG_DIR}")"
samtools view -h "${INPUT_BAM}" > "${TEMP_SAM}"

echo "Running makeTagDirectory..."
# -format sam        : explicit format (avoids relying on autodetection)
# -unique -mapq 10   : HOMER's default filtering (or -keepAll if --authors_orig,
#                      matching the original paper's makeTagDirectory calls)
# -genome -checkGC   : run sequence bias (GC%) diagnostics for this sample
# -fragLength        : use this sample's own insert size average, rather than
#                      HOMER's autocorrelation estimate (omitted entirely
#                      with --authors_orig - see above)
makeTagDirectory "${TAG_DIR}" "${TEMP_SAM}" \
    -format sam ${READ_FILTER_FLAGS} -genome "${GENOME_FASTA}" -checkGC \
    "${FRAGLENGTH_ARG[@]}" \
    2>&1 | tee "${LOG_DIR}/tagdir_${REGULATOR}_${REPLICATE}_makeTagDirectory.log"

# Clean up temporary SAM
rm -f "${TEMP_SAM}"

# ============================================================================
# VALIDATE OUTPUT
# ============================================================================
if [[ ! -d "${TAG_DIR}" ]]; then
    echo "ERROR: Failed to create tag directory"
    exit 1
fi

echo ""
echo "Tag directory created successfully:"
ls -lh "${TAG_DIR}" | head -n 5

if [[ -f "${TAG_DIR}/tagInfo.txt" ]]; then
    echo ""
    echo "Tag directory statistics:"
    grep "^genome" "${TAG_DIR}/tagInfo.txt" || true
    grep "fragmentLengthEstimate=" "${TAG_DIR}/tagInfo.txt" || true
fi

echo ""
echo "✓ Tag directory ready: ${TAG_DIR}"
