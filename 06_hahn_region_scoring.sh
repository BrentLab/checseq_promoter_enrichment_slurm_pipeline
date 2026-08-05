#!/bin/bash
#SBATCH --job-name=chec_hahn_region
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH -o logs/hahn_region_scoring_%a.log
#SBATCH -e logs/hahn_region_scoring_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/bioconductor-genomicranges_bioconductor-rtracklayer_r-optparse_r-tidyverse:8e8ae26758fa9e7f

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# ARRAY job - one task per UNIQUE regulator_symbol in <lookup_file>, not one
# task per lookup row/replicate (hahn_region_scoring.R already discovers and
# combines a regulator's own replicates internally - see that script's
# header). Depends on 01b_dmel_normalized_coverage.sh having completed for
# every sample, and on 04_pos2bed.sh's *_peaks_summits.txt files existing
# for peak assignment.
#
# Usage: sbatch --array=1-N_unique_regulators 06_hahn_region_scoring.sh <lookup_file> --tss-bed=<path> [additional hahn_region_scoring.R flags...]
#
# Array task N processes the Nth unique regulator_symbol found in
# <lookup_file> (alphabetically sorted, 1-indexed, so the index->regulator
# mapping is stable and deterministic across separate submissions/reruns of
# the same lookup file). --regulators is injected automatically from this -
# don't pass your own --regulators, it will be overridden.
#
# All other flags are passed straight through to the R script; run
# `Rscript hahn_region_scoring.R --help` (or see its own header
# comments) for the full parameter list: --peaks-dir, --coverage-dir,
# --output-dir, --promoter-upstream, --promoter-downstream,
# --peak-signal-window, --min-replicates-bound, --control-coverage,
# --pseudocount, --cores.
# NOTE: ${BASH_SOURCE[0]}-based self-location only works when this script
# runs directly via `bash`. Under `sbatch`, SLURM copies the submitted
# script into a per-job spool directory (/var/spool/slurmd/jobNNN/) and runs
# THAT COPY, so BASH_SOURCE[0] resolves to the spool path, not this script's
# real location - breaking any attempt to find sibling files (like the R
# script this wraps) relative to it. submit_pipeline.sh works around this by
# exporting PIPELINE_SCRIPT_DIR (via sbatch --export=) since IT can resolve
# its own location reliably (it's always run directly via bash, never
# sbatch). Prefer that if set; fall back to BASH_SOURCE for direct `bash`
# invocation (e.g. local testing without SLURM at all).
if [[ -n "${PIPELINE_SCRIPT_DIR:-}" ]]; then
    SCRIPT_DIR="${PIPELINE_SCRIPT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
R_SCRIPT="${SCRIPT_DIR}/hahn_region_scoring.R"

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    echo "  Expected hahn_region_scoring.R alongside this script"
    if [[ -n "${SLURM_JOB_ID:-}" && -z "${PIPELINE_SCRIPT_DIR:-}" ]]; then
        echo ""
        echo "  You're running under SLURM (SLURM_JOB_ID is set) without PIPELINE_SCRIPT_DIR"
        echo "  set. SLURM copies submitted scripts into a per-job spool directory, so this"
        echo "  script can't reliably find its own location to locate its sibling R script."
        echo "  If submitting directly with sbatch (not via submit_pipeline.sh), pass:"
        echo "    sbatch --export=ALL,PIPELINE_SCRIPT_DIR=/absolute/path/to/pipeline/dir ..."
    fi
    exit 1
fi

if [[ $# -lt 1 ]]; then
    echo "ERROR: Usage: 06_hahn_region_scoring.sh <lookup_file> --tss-bed=<path> [additional flags...]"
    exit 1
fi

LOOKUP_FILE="$1"
shift

if [[ ! -f "${LOOKUP_FILE}" ]]; then
    echo "ERROR: Lookup file not found: ${LOOKUP_FILE}"
    exit 1
fi

HAS_TSS_BED=false
for arg in "$@"; do
    [[ "${arg}" == --tss-bed=* ]] && HAS_TSS_BED=true
done

if [[ "${HAS_TSS_BED}" != "true" ]]; then
    echo "ERROR: --tss-bed=<path> is required"
    echo "Usage: 06_hahn_region_scoring.sh <lookup_file> --tss-bed=<path> [additional R script flags...]"
    exit 1
fi

# ============================================================================
# DETERMINE THIS TASK'S REGULATOR
# ============================================================================
# Unique regulator_symbol values (column 1), skipping the header, sorted
# alphabetically for a stable/deterministic index->regulator mapping.
REGULATOR=$(tail -n +2 "${LOOKUP_FILE}" | cut -f1 | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "${REGULATOR}" ]]; then
    echo "ERROR: Could not determine regulator for array task ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    echo "  (fewer unique regulators in the lookup file than the --array range submitted?)"
    exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: regulator = ${REGULATOR}"

mkdir -p logs

echo "Running hahn_region_scoring.R for ${REGULATOR} with args: $* --regulators=${REGULATOR}"

# Retry a couple of times with a short delay before giving up - guards
# against rare, transient failures observed when this step is chained
# immediately after a preceding array job via --dependency=afterok (e.g.
# momentary network filesystem inconsistency on a freshly-scheduled compute
# node), which didn't reproduce when the same command was run standalone.
MAX_ATTEMPTS=3
ATTEMPT=1
until Rscript "${R_SCRIPT}" "$@" --regulators="${REGULATOR}"; do
    EXIT_CODE=$?
    if [[ "${ATTEMPT}" -ge "${MAX_ATTEMPTS}" ]]; then
        echo "ERROR: hahn_region_scoring.R failed after ${MAX_ATTEMPTS} attempts (exit code ${EXIT_CODE})"
        exit "${EXIT_CODE}"
    fi
    echo "WARNING: hahn_region_scoring.R failed (exit code ${EXIT_CODE}, attempt ${ATTEMPT}/${MAX_ATTEMPTS}) - retrying in 30s..."
    sleep 30
    ATTEMPT=$((ATTEMPT + 1))
done

echo ""
echo "✓ Promoter scoring (dmel-normalized / Mahendrawada-Hahn method) complete for ${REGULATOR}"
