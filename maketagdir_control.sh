#!/bin/bash
#SBATCH --job-name=chec_tagdir_control
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH -o logs/tagdir_control_%a.log
#SBATCH -e logs/tagdir_control_%a.log
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
# Two calling conventions:
#
# 1. DIRECT MODE (single job - the original, unchanged interface):
#      sbatch maketagdir_control.sh <control_bam> <control_stats_file> [--authors_orig]
#    Builds ONE combined control tag directory at results/tag_dirs/control_MNase
#    from an already-combined BAM (e.g. replicates merged via samtools merge -
#    see README: Building the control samples).
#
# 2. LOOKUP MODE (array job, per-replicate - NEW):
#      sbatch --array=1-N maketagdir_control.sh <lookup_file> [bam_type] [--authors_orig]
#    Auto-detected from --array being used (SLURM_ARRAY_TASK_ID is set) - no
#    extra flag needed, since direct mode is never run as an array job. An
#    explicit --lookup marker is also accepted (e.g. for local testing
#    outside SLURM) but isn't required: `--lookup <lookup_file> [bam_type]`.
#    Builds a SEPARATE tag directory per free-MNase replicate, one per lookup
#    row, instead of one combined directory - useful for inspecting/QC-ing
#    each replicate's own tag directory (tag counts, GC bias, autocorrelation)
#    before deciding how to combine them, or for any analysis that wants a
#    per-replicate control rather than the single merged one.
#    bam_type: "nuclear" (default) or "full", same convention as elsewhere.
#    Output: results/tag_dirs/{regulator}/{replicate}/ (nested, matching
#    02_maketagdir_samples.sh's convention for regular samples - this is
#    deliberately NOT named "control_MNase", so it never collides with the
#    single combined control tag directory from direct mode).
#
#   --authors_orig: use -keepAll instead of -unique -mapq 10, matching the
#     original Mahendrawada et al. makeTagDirectory calls (which kept all
#     alignments, including multi-mappers/low-MAPQ reads, for both the
#     control and sample tag directories). Can appear anywhere in the args,
#     in either mode.
#
# Run this after 01_align.sh has completed for the control sample(s), and
# before running submit_pipeline.sh (which expects results/tag_dirs/control_MNase
# to already exist -- see 03_findpeaks.sh, which errors clearly if it's missing).
GENOME_FASTA="/ref/mblab/data/KN99/KN99_genome_fungidb.fasta"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# MODE DETECTION AND ARGUMENT PARSING
# ============================================================================
# Lookup mode is always run as a SLURM array (--array=1-N); direct mode never
# is - so SLURM_ARRAY_TASK_ID being set is a reliable signal on its own,
# without needing an extra flag. --lookup is still accepted as an optional,
# explicit marker (e.g. for local testing without a real array job), but
# isn't required.
if [[ "${1:-}" == "--lookup" || -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    MODE="lookup"
else
    MODE="direct"
fi

# Pull --authors_orig out of the argument list wherever it appears, leaving
# the remaining args (mode-specific) in order.
AUTHORS_ORIG=false
POSITIONAL=()
for arg in "$@"; do
    if [[ "${arg}" == "--authors_orig" ]]; then
        AUTHORS_ORIG=true
    else
        POSITIONAL+=("${arg}")
    fi
done

if [[ "${AUTHORS_ORIG}" == "true" ]]; then
    READ_FILTER_FLAGS="-keepAll"
else
    READ_FILTER_FLAGS="-unique -mapq 10"
fi

mkdir -p "${OUTPUT_DIR}/tag_dirs" "${LOG_DIR}"

if [[ "${MODE}" == "lookup" ]]; then
    # Drop the "--lookup" marker if it was given explicitly - it's optional
    # now that SLURM_ARRAY_TASK_ID alone triggers this mode, but still
    # accepted for clarity/local-testing purposes.
    LOOKUP_ARGS=("${POSITIONAL[@]}")
    if [[ "${LOOKUP_ARGS[0]:-}" == "--lookup" ]]; then
        LOOKUP_ARGS=("${LOOKUP_ARGS[@]:1}")
    fi
    LOOKUP_FILE="${LOOKUP_ARGS[0]:?ERROR: must provide a lookup file}"
    BAM_TYPE="${LOOKUP_ARGS[1]:-nuclear}"

    case "${BAM_TYPE}" in
        nuclear) BAM_SUFFIX="_nuclear.bam" ;;
        full)    BAM_SUFFIX=".bam" ;;
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

    CONTROL_BAM="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}${BAM_SUFFIX}"
    CONTROL_STATS_FILE="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}/${REGULATOR}_${REPLICATE}_samtools_stats.txt"
    CONTROL_TAG_DIR="${OUTPUT_DIR}/tag_dirs/${REGULATOR}/${REPLICATE}"
    SAMPLE_LABEL="${REGULATOR}_${REPLICATE}"

else
    CONTROL_BAM="${POSITIONAL[0]:-/path/to/free_MNase_control.bam}"
    CONTROL_STATS_FILE="${POSITIONAL[1]:?ERROR: must provide path to control samtools_stats.txt as \$2}"
    CONTROL_TAG_DIR="${OUTPUT_DIR}/tag_dirs/control_MNase"
    SAMPLE_LABEL="control_MNase (combined)"
fi

echo "Creating HOMER tag directory for ${SAMPLE_LABEL}..."
echo "  Input BAM: ${CONTROL_BAM}"
echo "  Stats file: ${CONTROL_STATS_FILE}"
echo "  Read filtering: ${READ_FILTER_FLAGS}$( [[ "${AUTHORS_ORIG}" == "true" ]] && echo " (--authors_orig mode)" )"
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
# (skipped entirely with --authors_orig - matches the original scripts,
# which never passed -fragLength and let HOMER's own autocorrelation
# estimate run instead)
# ============================================================================
if [[ "${AUTHORS_ORIG}" == "true" ]]; then
    FRAGLENGTH_ARG=()
    echo "  Fragment length: letting HOMER auto-estimate (--authors_orig mode, skipping samtools stats)"
else
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
    FRAGLENGTH_ARG=(-fragLength "${FRAGMENT_LENGTH}")
fi

# ============================================================================
# MAKE TAG DIRECTORY
# ============================================================================
# Convert BAM to SAM (makeTagDirectory expects SAM format)
CONTROL_SAM="${CONTROL_TAG_DIR}_temp.sam"
mkdir -p "$(dirname "${CONTROL_SAM}")"
samtools view -h "${CONTROL_BAM}" > "${CONTROL_SAM}"

# -format sam        : explicit format (avoids relying on autodetection)
# -unique -mapq 10   : HOMER's default filtering (or -keepAll if --authors_orig,
#                      matching the original paper's makeTagDirectory calls)
# -fragLength        : override HOMER's autocorrelation estimate with the
#                      directly measured insert size average, since HOMER
#                      could not produce a reliable peak-size estimate on
#                      its own for this background control (omitted
#                      entirely with --authors_orig - see above)
# -genome -checkGC   : run sequence bias (GC%) diagnostics - important to
#                      check on a background control, since a GC-shifted
#                      library would look like fake "enrichment" at GC-rich
#                      regions (e.g. TSS/CpG-like regions) downstream
makeTagDirectory "${CONTROL_TAG_DIR}" "${CONTROL_SAM}" \
    -format sam ${READ_FILTER_FLAGS} \
    "${FRAGLENGTH_ARG[@]}" \
    -genome "${GENOME_FASTA}" -checkGC \
    2>&1 | tee "${LOG_DIR}/tagdir_$(basename "${CONTROL_TAG_DIR}")_makeTagDirectory.log"

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
echo "✓ Tag directory ready: ${CONTROL_TAG_DIR}"
