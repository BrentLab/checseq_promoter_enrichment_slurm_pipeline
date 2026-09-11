#!/bin/bash
#SBATCH --job-name=chec_promoter_scoring_lookup
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o logs/promoter_scoring_lookup_%a.log
#SBATCH -e logs/promoter_scoring_lookup_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/bioconductor-genomicranges_bioconductor-rtracklayer_r-optparse_r-tidyverse:8e8ae26758fa9e7f

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# STANDALONE script - NOT part of the automated submit_pipeline.sh chain.
# Runs promoter_scoring.R once per row of a lookup file, where EACH ROW
# defines its own tagged-sample/control pairing - unlike 07_promoter_scoring.sh,
# which uses one shared control for every regulator in a run. Built for
# designs where each regulator has its own matched background (e.g. the same
# regulator tagged in two genetic backgrounds), rather than one shared
# free-MNase-style control used across all regulators.
#
# Usage:
#   sbatch --array=1-N --export=ALL,PIPELINE_SCRIPT_DIR=/absolute/path/to/pipeline/dir \
#       promoter_scoring_standalone.sh <lookup_file> <promoter_bed> \
#       [--genomecov-dir=<path>] [--output-dir=<path>] [additional promoter_scoring.R flags...]
#
# <lookup_file>: TSV with header, one row per comparison, e.g.:
#   sample_name	tagged_regulator	control_bed
#   ACA1	ACA1∆msn2	ACA1
#   ADR1	ADR1∆msn2	ADR1
#   ...
#   - sample_name: label for this comparison; used to name the output subdir
#   - tagged_regulator: subdirectory name under --genomecov-dir to use as the
#     tagged/experiment sample (passed to promoter_scoring.R's --regulators)
#   - control_bed: subdirectory name under --genomecov-dir to use as the
#     control/background for THIS row. Default: the actual file path is
#     constructed as <genomecov-dir>/<control_bed>/<control_bed>_r1_5p.bed
#     (i.e. control_bed already has one pre-combined *_r1_5p.bed file, same
#     as a normal single-replicate regulator directory). If control_bed
#     instead has its OWN replicate subdirectories that haven't been
#     pre-combined (matching the same nested structure tagged_regulator
#     uses), pass --combine-control-replicates and this name is resolved
#     the same way tagged_regulator already is - replicates discovered and
#     summed automatically, rather than expecting one flat file.
#
# <promoter_bed>: path to the shared promoter regions BED (same for every
# row) - passed to promoter_scoring.R's --promoter-bed
#
# --genomecov-dir=<path>  (default: results/genomecov_5p)
# --combine-control-replicates  (default: off)
#     Off (default, unchanged): control_bed must already be one combined
#     *_r1_5p.bed file - existing behavior, nothing new required for the
#     normal single-shared-control case (e.g. free MNase, already merged).
#     On: control_bed is instead treated as having its own replicate
#     subdirectories to discover and combine, same as tagged_regulator
#     already does by default - this makes replicate-combining symmetric
#     between the tagged and control sides, which previously only happened
#     automatically on the tagged side.
# --output-dir=<path>     (default: <promoter_bed basename, extension
#                          stripped>/<sample_name> - e.g. for
#                          start_codon_500bp_upstream_promoters.bed and
#                          sample_name=ACA1:
#                          start_codon_500bp_upstream_promoters/ACA1
#                          promoter_scoring.R's own unchanged default
#                          behavior then nests a further
#                          {tagged_regulator}/ level inside that, same as it
#                          always does - see its header for why this isn't
#                          overridden here.)
#
# Any other flags (--pseudocount, --cores, etc.) are passed straight through
# to promoter_scoring.R.
#
# Since this always runs via `sbatch --array`, SLURM copies the submitted
# script into a per-job spool directory - BASH_SOURCE-based self-location
# won't find promoter_scoring.R there (see 07_promoter_scoring.sh's header
# for the full explanation). You MUST export PIPELINE_SCRIPT_DIR yourself
# when submitting (there's no submit_pipeline.sh wrapper doing it for you,
# since this script is standalone):
#   sbatch --array=1-N --export=ALL,PIPELINE_SCRIPT_DIR=/absolute/path/to/pipeline/dir ...

if [[ -n "${PIPELINE_SCRIPT_DIR:-}" ]]; then
    SCRIPT_DIR="${PIPELINE_SCRIPT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
R_SCRIPT="${SCRIPT_DIR}/promoter_scoring.R"

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    echo "  Expected promoter_scoring.R alongside this script"
    if [[ -n "${SLURM_JOB_ID:-}" && -z "${PIPELINE_SCRIPT_DIR:-}" ]]; then
        echo ""
        echo "  You're running under SLURM (SLURM_JOB_ID is set) without PIPELINE_SCRIPT_DIR"
        echo "  set. SLURM copies submitted scripts into a per-job spool directory, so this"
        echo "  script can't reliably find its own location to locate its sibling R script."
        echo "  Pass it explicitly:"
        echo "    sbatch --export=ALL,PIPELINE_SCRIPT_DIR=/absolute/path/to/pipeline/dir ..."
    fi
    exit 1
fi

# ============================================================================
# ARGUMENT PARSING
# ============================================================================
if [[ $# -lt 2 ]]; then
    echo "ERROR: Usage: promoter_scoring_lookup.sh <lookup_file> <promoter_bed> [--genomecov-dir=<path>] [--output-dir=<path>] [additional flags...]"
    exit 1
fi

LOOKUP_FILE="$1"
PROMOTER_BED="$2"
shift 2

GENOMECOV_DIR="results/genomecov_5p"
OUTPUT_DIR_OVERRIDE=""
COMBINE_CONTROL_REPLICATES=false
PASSTHROUGH_ARGS=()

for arg in "$@"; do
    case "${arg}" in
        --genomecov-dir=*)
            GENOMECOV_DIR="${arg#--genomecov-dir=}"
            ;;
        --output-dir=*)
            OUTPUT_DIR_OVERRIDE="${arg#--output-dir=}"
            ;;
        --combine-control-replicates)
            COMBINE_CONTROL_REPLICATES=true
            PASSTHROUGH_ARGS+=("${arg}")
            ;;
        *)
            PASSTHROUGH_ARGS+=("${arg}")
            ;;
    esac
done

if [[ ! -f "${LOOKUP_FILE}" ]]; then
    echo "ERROR: Lookup file not found: ${LOOKUP_FILE}"
    exit 1
fi

if [[ ! -f "${PROMOTER_BED}" ]]; then
    echo "ERROR: Promoter BED not found: ${PROMOTER_BED}"
    exit 1
fi

mkdir -p logs

# ============================================================================
# READ THIS TASK'S ROW (skip header, same convention as the rest of the pipeline)
# ============================================================================
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")
if [[ -z "${LINE}" ]]; then
    echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    exit 1
fi
read -r SAMPLE_NAME TAGGED_REGULATOR CONTROL_BED_NAME <<< "${LINE}"

# --combine-control-replicates changes what --control-bed means downstream:
# instead of a literal, already-combined file, promoter_scoring.R treats it
# as a regulator-style NAME and discovers+combines its own replicates the
# same way tagged regulators already are (see that script's header). Default
# (off) keeps the original literal-file behavior unchanged.
if [[ "${COMBINE_CONTROL_REPLICATES}" == "true" ]]; then
    CONTROL_BED_ARG="${CONTROL_BED_NAME}"
    if [[ ! -d "${GENOMECOV_DIR}/${CONTROL_BED_NAME}" ]]; then
        echo "ERROR: --combine-control-replicates is set, but no directory found for control_bed='${CONTROL_BED_NAME}' under --genomecov-dir=${GENOMECOV_DIR}"
        exit 1
    fi
else
    CONTROL_BED_ARG="${GENOMECOV_DIR}/${CONTROL_BED_NAME}/${CONTROL_BED_NAME}_r1_5p.bed"
    if [[ ! -f "${CONTROL_BED_ARG}" ]]; then
        echo "ERROR: Control bed file not found: ${CONTROL_BED_ARG}"
        echo "  (constructed from --genomecov-dir=${GENOMECOV_DIR} and control_bed='${CONTROL_BED_NAME}')"
        echo "  If '${CONTROL_BED_NAME}' actually has multiple replicate subdirectories rather"
        echo "  than one pre-combined file, pass --combine-control-replicates."
        exit 1
    fi
fi

if [[ ! -d "${GENOMECOV_DIR}/${TAGGED_REGULATOR}" ]]; then
    echo "ERROR: Tagged regulator directory not found: ${GENOMECOV_DIR}/${TAGGED_REGULATOR}"
    exit 1
fi

# ============================================================================
# OUTPUT DIRECTORY
# ============================================================================
if [[ -n "${OUTPUT_DIR_OVERRIDE}" ]]; then
    OUTPUT_DIR="${OUTPUT_DIR_OVERRIDE}"
else
    PROMOTER_BED_BASENAME="$(basename "${PROMOTER_BED}")"
    PROMOTER_BED_BASENAME="${PROMOTER_BED_BASENAME%.*}"
    OUTPUT_DIR="${PROMOTER_BED_BASENAME}/${SAMPLE_NAME}"
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_NAME}"
echo "  Tagged regulator: ${TAGGED_REGULATOR}"
echo "  Control bed: ${CONTROL_BED_ARG}"
echo "  Promoter bed: ${PROMOTER_BED}"
echo "  Output dir: ${OUTPUT_DIR}"

# ============================================================================
# RUN (with retry-with-backoff, same rationale as 06/07 - see their headers)
# ============================================================================
MAX_ATTEMPTS=3
ATTEMPT=1
until Rscript "${R_SCRIPT}" \
    --promoter-bed="${PROMOTER_BED}" \
    --control-bed="${CONTROL_BED_ARG}" \
    --genomecov-dir="${GENOMECOV_DIR}" \
    --regulators="${TAGGED_REGULATOR}" \
    --output-dir="${OUTPUT_DIR}" \
    "${PASSTHROUGH_ARGS[@]}"; do
    EXIT_CODE=$?
    if [[ "${ATTEMPT}" -ge "${MAX_ATTEMPTS}" ]]; then
        echo "ERROR: promoter_scoring.R failed after ${MAX_ATTEMPTS} attempts (exit code ${EXIT_CODE})"
        exit "${EXIT_CODE}"
    fi
    echo "WARNING: promoter_scoring.R failed (exit code ${EXIT_CODE}, attempt ${ATTEMPT}/${MAX_ATTEMPTS}) - retrying in 30s..."
    sleep 30
    ATTEMPT=$((ATTEMPT + 1))
done

echo ""
echo "✓ Promoter scoring complete for ${SAMPLE_NAME}"
