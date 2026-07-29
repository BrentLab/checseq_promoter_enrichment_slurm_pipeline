#!/bin/bash
# Convenience script to submit the automated ChEC-seq pipeline jobs with proper
# dependencies. Run after 00_prepare.sh and configuring parameterized paths in
# each script.
#
# PREREQUISITE: maketagdir_control.sh must already have been run manually to
# build results/tag_dirs/control_MNase before running this script. It is not
# part of the automated chain below, since the control tag directory is
# typically built once and reused across pipeline runs rather than
# regenerated every time.
#
# Usage: submit_pipeline.sh <lookup_file> [--bam-type=nuclear|full] [--start-at=STEP] [--authors-orig]
#
#   --bam-type=nuclear|full   (default: nuclear)
#       Which BAM 02_maketagdir_samples.sh builds tag directories from:
#         nuclear -> {regulator}_{replicate}_nuclear.bam (chrM filtered out)
#         full    -> {regulator}_{replicate}.bam         (all chromosomes)
#
#   --start-at=STEP           (default: 01_align)
#       Which step to begin submitting from - anything upstream of STEP is
#       assumed to have already completed successfully in a prior run, and
#       is not resubmitted. STEP can be any of the pipeline script names,
#       with or without the .sh extension, e.g. --start-at=03_findpeaks
#       Valid values: 01_align, 01a_map_to_dmel, 02_maketagdir_samples,
#                     03_findpeaks, 04_pos2bed, 05_annotatepeaks, 06_multiqc
#
#   --authors-orig            (default: off)
#       Passes --authors_orig to 02_maketagdir_samples.sh, which uses -keepAll
#       instead of -unique -mapq 10 when building sample tag directories -
#       matching the original Mahendrawada et al. makeTagDirectory calls.
#       NOTE: this only affects 02_maketagdir_samples.sh (the automated
#       chain). If you also want the control tag directory built the same
#       way, pass --authors_orig directly to maketagdir_control.sh yourself,
#       since that script is run manually and isn't part of this chain.
#
# Genome FASTA and GTF paths for the annotatePeaks step are hardcoded inside
# 05_annotatepeaks.sh itself, not passed as arguments here.

set -euo pipefail

# Directory this script itself lives in - lets you invoke submit_pipeline.sh
# from any working directory (e.g. `bash path/to/pipeline/submit_pipeline.sh
# test.lookup`) without first cd-ing into the pipeline directory. All the
# individual step scripts (01_align.sh, etc.) are expected to live alongside
# this one and are referenced via this path, NOT via the current directory.
#
# LOOKUP_FILE, by contrast, stays relative to wherever YOU run this from
# (your current directory) since that's where your sample sheet lives -
# results/ and logs/ are likewise created in your current directory when the
# submitted jobs actually run (SLURM jobs default to the submission dir).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# ARGUMENT PARSING
# ============================================================================
LOOKUP_FILE=""
BAM_TYPE="nuclear"
START_STEP="01_align"
AUTHORS_ORIG=false

for arg in "$@"; do
    case "${arg}" in
        --bam-type=*)
            BAM_TYPE="${arg#--bam-type=}"
            ;;
        --start-at=*)
            START_STEP="${arg#--start-at=}"
            START_STEP="${START_STEP%.sh}"   # tolerate a trailing .sh
            ;;
        --authors-orig)
            AUTHORS_ORIG=true
            ;;
        --*)
            echo "ERROR: Unrecognized option: ${arg}"
            exit 1
            ;;
        *)
            if [[ -z "${LOOKUP_FILE}" ]]; then
                LOOKUP_FILE="${arg}"
            else
                echo "ERROR: Unexpected extra argument: ${arg}"
                exit 1
            fi
            ;;
    esac
done

if [[ -z "${LOOKUP_FILE}" ]]; then
    echo "ERROR: lookup_file is required. Usage: submit_pipeline.sh <lookup_file> [--bam-type=nuclear|full] [--start-at=STEP] [--authors-orig]"
    exit 1
fi

if [[ "${BAM_TYPE}" != "nuclear" && "${BAM_TYPE}" != "full" ]]; then
    echo "ERROR: --bam-type must be 'nuclear' or 'full' (got '${BAM_TYPE}')"
    exit 1
fi

CONTROL_TAG_DIR="results/tag_dirs/control_MNase"

# ============================================================================
# PIPELINE DAG
# ============================================================================
# Ordered list of steps. Order matters for --start-at (anything with a lower
# index than START_STEP is skipped), but the actual submitted --dependency
# for each step comes from STEP_DEPS below, not from this ordering alone -
# so branches (01a/02 off 01_align, 04/05 off 03_findpeaks) are handled
# correctly regardless of where --start-at lands.
STEP_ORDER=(01_align 01a_map_to_dmel 02_maketagdir_samples 03_findpeaks 04_pos2bed 05_annotatepeaks 06_multiqc)

declare -A STEP_SCRIPT=(
    [01_align]="01_align.sh"
    [01a_map_to_dmel]="01a_map_to_dmel.sh"
    [02_maketagdir_samples]="02_maketagdir_samples.sh"
    [03_findpeaks]="03_findpeaks.sh"
    [04_pos2bed]="04_pos2bed.sh"
    [05_annotatepeaks]="05_annotatepeaks.sh"
    [06_multiqc]="06_multiqc.sh"
)

# Upstream dependencies for each step (space-separated step names, or empty)
declare -A STEP_DEPS=(
    [01_align]=""
    [01a_map_to_dmel]="01_align"
    [02_maketagdir_samples]="01_align"
    [03_findpeaks]="02_maketagdir_samples"
    [04_pos2bed]="03_findpeaks"
    [05_annotatepeaks]="03_findpeaks"
    [06_multiqc]="01a_map_to_dmel 04_pos2bed 05_annotatepeaks"
)

# 06_multiqc uses afterany (partial upstream failure still produces a report);
# every other step uses afterok.
declare -A STEP_DEP_TYPE=(
    [06_multiqc]="afterany"
)

# Steps that are single jobs, not SLURM arrays
declare -A STEP_IS_SINGLE=(
    [06_multiqc]=1
)

# ============================================================================
# VALIDATION
# ============================================================================
if [[ ! -f "${LOOKUP_FILE}" ]]; then
    echo "ERROR: ${LOOKUP_FILE} not found"
    exit 1
fi

if [[ ! -d "${CONTROL_TAG_DIR}" ]]; then
    echo "ERROR: Control tag directory not found: ${CONTROL_TAG_DIR}"
    echo "  Run maketagdir_control.sh first to build the control tag directory"
    echo "  (this is a manual, one-time step, not part of this automated pipeline)."
    exit 1
fi

if [[ ! -f "${CONTROL_TAG_DIR}/tagInfo.txt" ]]; then
    echo "ERROR: ${CONTROL_TAG_DIR}/tagInfo.txt not found - control tag directory looks incomplete"
    exit 1
fi

# Validate --start-at is a real step
START_INDEX=-1
for i in "${!STEP_ORDER[@]}"; do
    if [[ "${STEP_ORDER[$i]}" == "${START_STEP}" ]]; then
        START_INDEX=$i
        break
    fi
done
if [[ "${START_INDEX}" -eq -1 ]]; then
    echo "ERROR: Invalid --start-at value: '${START_STEP}'"
    echo "  Valid values: ${STEP_ORDER[*]}"
    exit 1
fi

# Get array size (total lines minus header)
ARRAY_SIZE=$(tail -n +2 "${LOOKUP_FILE}" | wc -l)

if [[ $ARRAY_SIZE -eq 0 ]]; then
    echo "ERROR: No samples found in ${LOOKUP_FILE}"
    exit 1
fi

echo "================================"
echo "ChEC-seq Pipeline Submission"
echo "================================"
echo "Lookup file: ${LOOKUP_FILE}"
echo "Format: Paired-end (regulator, replicate, fastq_1, fastq_2)"
echo "Array size (# of samples): ${ARRAY_SIZE}"
echo "Control tag directory: ${CONTROL_TAG_DIR}"
echo "BAM type for tag directories: ${BAM_TYPE}"
echo "Read filtering (02_maketagdir_samples): $( [[ "${AUTHORS_ORIG}" == "true" ]] && echo "-keepAll (--authors-orig)" || echo "-unique -mapq 10 (default)" )"
echo "Starting at step: ${START_STEP}"
if [[ "${START_INDEX}" -gt 0 ]]; then
    echo "  (steps before this are assumed already complete and will NOT be resubmitted)"
fi
echo ""

# ============================================================================
# VALIDATE SCRIPT FILES
# ============================================================================
for step in "${STEP_ORDER[@]}"; do
    script="${SCRIPT_DIR}/${STEP_SCRIPT[$step]}"
    if [[ ! -f "${script}" ]]; then
        echo "ERROR: Script not found: ${script}"
        exit 1
    fi
done

echo "✓ All scripts found in ${SCRIPT_DIR}"
echo ""

# ============================================================================
# HELPER: submit a job and return ONLY a validated numeric job ID
# ============================================================================
# Guards against the classic footgun where `sbatch ... | awk '{print $NF}'`
# silently returns an error message (not a job ID) when submission fails -
# which would then be spliced into a downstream --dependency and poison the
# chain. On any non-numeric result we abort immediately with the raw sbatch
# output.
submit_job() {
    local out jobid
    out=$(sbatch "$@" 2>&1)
    jobid=$(awk '{print $NF}' <<< "${out}")
    if ! [[ "${jobid}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: sbatch submission failed for: $*" >&2
        echo "  sbatch output: ${out}" >&2
        exit 1
    fi
    printf '%s' "${jobid}"
}

# ============================================================================
# SUBMIT JOBS
# ============================================================================
echo "Submitting pipeline jobs..."
echo ""

declare -A JOBIDS=()

for i in "${!STEP_ORDER[@]}"; do
    step="${STEP_ORDER[$i]}"

    if [[ "$i" -lt "${START_INDEX}" ]]; then
        echo "-- Skipping ${step} (before --start-at=${START_STEP})"
        continue
    fi

    script="${SCRIPT_DIR}/${STEP_SCRIPT[$step]}"

    # Build --dependency from whichever of this step's upstream deps were
    # ACTUALLY submitted in this run. Deps that were skipped (because they're
    # before START_STEP) are assumed already complete and simply omitted -
    # the step then submits with no dependency on that branch.
    dep_ids=()
    for d in ${STEP_DEPS[$step]}; do
        if [[ -n "${JOBIDS[$d]:-}" ]]; then
            dep_ids+=("${JOBIDS[$d]}")
        fi
    done

    dep_args=()
    if [[ "${#dep_ids[@]}" -gt 0 ]]; then
        dep_type="${STEP_DEP_TYPE[$step]:-afterok}"
        dep_str="${dep_type}"
        for jid in "${dep_ids[@]}"; do
            dep_str+=":${jid}"
        done
        dep_args=(--dependency="${dep_str}")
    fi

    # Extra positional args each script needs, beyond LOOKUP_FILE
    extra_args=("${LOOKUP_FILE}")
    if [[ "${step}" == "02_maketagdir_samples" || "${step}" == "03_findpeaks" ]]; then
        extra_args+=("${BAM_TYPE}")
    fi
    if [[ "${step}" == "02_maketagdir_samples" && "${AUTHORS_ORIG}" == "true" ]]; then
        extra_args+=("--authors_orig")
    fi

    if [[ -n "${STEP_IS_SINGLE[$step]:-}" ]]; then
        echo "Submitting ${step} (single job)..."
        echo "  Command: sbatch ${dep_args[*]:-} ${script}"
        JOBIDS[$step]=$(submit_job "${dep_args[@]}" "${script}")
    else
        echo "Submitting ${step} (array 1-${ARRAY_SIZE})..."
        echo "  Command: sbatch ${dep_args[*]:-} --array=1-${ARRAY_SIZE} ${script} ${extra_args[*]}"
        JOBIDS[$step]=$(submit_job "${dep_args[@]}" --array=1-${ARRAY_SIZE} "${script}" "${extra_args[@]}")
    fi
    echo "  Job ID: ${JOBIDS[$step]}"
    echo ""
done

# ============================================================================
# SUMMARY
# ============================================================================
echo "================================"
echo "Pipeline Submission Complete"
echo "================================"
echo ""
echo "Job dependencies:"
echo "  (control_MNase tag dir already built manually via maketagdir_control.sh)"
for step in "${STEP_ORDER[@]}"; do
    if [[ -n "${JOBIDS[$step]:-}" ]]; then
        echo "  ${step}: ${JOBIDS[$step]}"
    else
        echo "  ${step}: (skipped - assumed already complete)"
    fi
done
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  tail -f logs/align_1.log"
echo ""
if [[ -n "${JOBIDS[06_multiqc]:-}" ]]; then
    echo "Track final job:"
    echo "  squeue -j ${JOBIDS[06_multiqc]}"
    echo ""
    echo "When 06_multiqc finishes, the consolidated report is at:"
    echo "  results/multiqc/multiqc_report.html"
fi
echo ""
