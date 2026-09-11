#!/bin/bash
# Convenience script to submit the automated ChEC-seq pipeline jobs with proper
# dependencies. Run after 00_prepare.sh and configuring parameterized paths in
# each script.
#
# OPTIONAL PREREQUISITE: if you want HOMER peak calling (03_findpeaks.sh
# onward, plus 06_hahn_region_scoring.sh), maketagdir_control.sh must first
# be run manually to build a control tag directory, then passed via
# --control-tag-dir=<path>. It is not part of the automated chain below,
# since the control tag directory is typically built once and reused across
# pipeline runs rather than regenerated every time. Without
# --control-tag-dir, those steps are simply skipped - see its entry below.
#
# Usage: submit_pipeline.sh <lookup_file> \
#     --primary-bowtie-index=PATH --mito-chrom=NAME --primary-genome-fasta=PATH \
#     [--bam-type=nuclear|full] [--start-at=STEP] [--authors-orig] \
#     [--align_dmel --spikein-bowtie-index=PATH [--tss-bed=TSS.bed [--control-coverage=CONTROL.bedgraph]]] \
#     [--filter_genomecov --include-regions=REGIONS.bed [--promoter-bed=PROMOTERS.bed --control-bed=CONTROL.bed]] \
#     [--control-tag-dir=PATH --primary-gtf-file=PATH]
#
#   --primary-bowtie-index=<path>   (required - no default)
#       Bowtie2 index prefix for 01_align.sh, the main organism alignment.
#       01_align.sh always runs, so this is always required.
#
#   --mito-chrom=<name>   (required - no default)
#       Mitochondrial chromosome name, as it appears in the reference/BAM,
#       for 01_align.sh's nuclear/mito BAM split (e.g. "chrM" for many
#       S. cerevisiae assemblies - varies by organism/genome build).
#       01_align.sh always runs, so this is always required.
#
#   --primary-genome-fasta=<path>   (required - no default)
#       Reference genome FASTA. Used by 02_maketagdir_samples.sh (always
#       runs, so this is always required) and, if --control-tag-dir is also
#       set, by 05_annotatepeaks.sh too (same file, passed to both).
#
#   --spikein-bowtie-index=<path>   (required only if --align_dmel is set)
#       Bowtie2 index prefix for 01a_map_to_dmel.sh, the D. melanogaster
#       spike-in alignment. Only needed - and only required - when
#       --align_dmel is set, since 01a_map_to_dmel.sh isn't submitted at
#       all otherwise. Omitting both --align_dmel and this flag together is
#       fine and produces no error.
#
#   --primary-gtf-file=<path>   (required only if --control-tag-dir is set)
#       GTF for 05_annotatepeaks.sh's peak annotation. Only needed - and
#       only required - when --control-tag-dir is set, since
#       05_annotatepeaks.sh isn't submitted at all otherwise.
#
#   --bam-type=nuclear|full   (default: nuclear)
#       Which BAM 01b_dmel_normalized_coverage.sh, 01c_filter_bam.sh,
#       02_maketagdir_samples.sh, and 03_findpeaks.sh all use:
#         nuclear -> {regulator}_{replicate}_nuclear.bam (chrM filtered out)
#         full    -> {regulator}_{replicate}.bam         (all chromosomes)
#
#   --start-at=STEP           (default: 01_align)
#       Which step to begin submitting from - anything upstream of STEP is
#       assumed to have already completed successfully in a prior run, and
#       is not resubmitted. STEP can be any of the pipeline script names,
#       with or without the .sh extension, e.g. --start-at=03_findpeaks
#       Valid values: 01_align, 01a_map_to_dmel, 01b_dmel_coverage,
#                     01c_filter_bam, 01d_genomecov_5p, 02_maketagdir_samples,
#                     03_findpeaks, 04_pos2bed, 05_annotatepeaks,
#                     06_hahn_region_scoring, 07_promoter_scoring,
#                     08_multiqc
#       NOTE: --start-at=01a_map_to_dmel or 01b_dmel_coverage also requires
#       --align_dmel; --start-at=01c_filter_bam or 01d_genomecov_5p also
#       requires --filter_genomecov; --start-at=03_findpeaks, 04_pos2bed, or
#       05_annotatepeaks also requires --control-tag-dir;
#       --start-at=06_hahn_region_scoring requires both --tss-bed and
#       --control-tag-dir - these steps are otherwise disabled entirely.
#
#   --authors-orig            (default: off)
#       Passes --authors_orig to both 02_maketagdir_samples.sh and
#       03_findpeaks.sh, matching the original Mahendrawada et al. scripts
#       more closely in three ways at once:
#         - 02_maketagdir_samples.sh: -keepAll instead of -unique -mapq 10
#           when building sample tag directories
#         - 02_maketagdir_samples.sh: skips -fragLength entirely (no longer
#           derives it from this sample's samtools stats), letting HOMER's
#           own autocorrelation estimate run instead - the original scripts
#           never passed -fragLength
#         - 03_findpeaks.sh: skips -gsize entirely (no longer passes the
#           hardcoded nuclear/full genome size constant), letting findPeaks
#           auto-estimate genome size from the tag directory instead - the
#           original scripts never passed -gsize either
#       NOTE: this only affects the automated chain (02/03). If you also
#       want the control tag directory built the same way (-keepAll, no
#       -fragLength), pass --authors_orig directly to maketagdir_control.sh
#       yourself, since that script is run manually and isn't part of this
#       chain.
#
#   --align_dmel              (default: off)
#       Enables the D. melanogaster spike-in branch: 01a_map_to_dmel.sh and
#       01b_dmel_normalized_coverage.sh. Off by default since not all sample
#       sets have a dmel spike-in. When off, these two steps are not
#       submitted at all (not "assumed already done" - genuinely skipped),
#       and 08_multiqc's dependency on them is dropped automatically.
#
#   --tss-bed=TSS.bed
#       Enables 06_hahn_region_scoring.sh (the Mahendrawada 2025/Donczew &
#       Hahn 2020 promoter-scoring method), once --align_dmel has produced
#       coverage for every sample. ALSO requires --control-tag-dir=<path>,
#       since this step reads 04_pos2bed's peak-assignment output too, not
#       just dmel coverage. Optional even with --align_dmel set - if
#       omitted (or --control-tag-dir is omitted), 01a/01b still run but
#       scoring is skipped. Runs as an array
#       job, one task per unique regulator_symbol in <lookup_file> (not one
#       per lookup row/replicate). See hahn_region_scoring.R's own header
#       for additional tunable parameters (promoter window, signal window,
#       min replicates bound); pass extras straight through, e.g. by editing
#       06_hahn_region_scoring.sh or invoking it directly once its
#       dependencies exist.
#
#   --control-coverage=CONTROL.bedgraph
#       Optional, only used alongside --tss-bed. Path to the combined
#       control's dmel-normalized bedgraph (see README: Building the control
#       samples). If provided, 06_hahn_region_scoring.sh additionally
#       reports background_tag_count and enrichment (tagged/background
#       ratio) - see promoter_enrichment_coverage.R's header for why this is
#       an addition on top of, not a replacement for, the paper's own
#       tagged_tag_count ("signal per promoter") metric.
#
#   --filter_genomecov        (default: off, requires --include-regions=)
#       Enables an independent, dmel-free quantification branch:
#       01c_filter_bam.sh (region-restricted, MAPQ/pairing-filtered BAM) ->
#       01d_genomecov_5p.sh (per-base, per-strand 5' cut-site coverage).
#       Off by default.
#
#   --include-regions=REGIONS.bed
#       Required when --filter_genomecov is set. BED file of regions to
#       restrict 01c_filter_bam.sh to (e.g. promoters).
#
#   --promoter-bed=PROMOTERS.bed / --control-bed=CONTROL.bed
#       Both enable 07_promoter_scoring.sh (the calling-cards-
#       style promoter enrichment method), once --filter_genomecov has
#       produced 5' cut-site coverage for every sample. Both optional even
#       with --filter_genomecov set - if either is omitted, 01c/01d still
#       run but scoring is skipped. Runs as an array job, one task per
#       unique regulator_symbol in <lookup_file> (not one per lookup
#       row/replicate). --control-bed points at a combined control 5'
#       cut-site BED (see README: Building the control samples).
#
#   --control-tag-dir=<path>  (default: none - HOMER peak-calling steps skipped)
#       Path to the HOMER control tag directory built by maketagdir_control.sh
#       (manual, standalone step - see README: Building the control samples).
#       Optional - if omitted, 03_findpeaks.sh, 04_pos2bed.sh,
#       05_annotatepeaks.sh, and 06_hahn_region_scoring.sh (which reads
#       04_pos2bed's output too) are all skipped entirely, rather than
#       requiring a control tag directory to exist at some default location.
#       02_maketagdir_samples.sh still runs either way. If provided, it's
#       passed through to 03_findpeaks.sh, and must already exist (built via
#       maketagdir_control.sh) with a valid tagInfo.txt.
#
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
ALIGN_DMEL=false
FILTER_GENOMECOV=false
INCLUDE_REGIONS_BED=""
TSS_BED=""
CONTROL_COVERAGE=""
PROMOTER_BED=""
CONTROL_BED=""
CONTROL_TAG_DIR=""
PRIMARY_BOWTIE_INDEX=""
MITO_CHROM=""
SPIKEIN_BOWTIE_INDEX=""
PRIMARY_GENOME_FASTA=""
PRIMARY_GTF_FILE=""

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
        --align_dmel)
            ALIGN_DMEL=true
            ;;
        --filter_genomecov)
            FILTER_GENOMECOV=true
            ;;
        --include-regions=*)
            INCLUDE_REGIONS_BED="${arg#--include-regions=}"
            ;;
        --tss-bed=*)
            TSS_BED="${arg#--tss-bed=}"
            ;;
        --control-coverage=*)
            CONTROL_COVERAGE="${arg#--control-coverage=}"
            ;;
        --promoter-bed=*)
            PROMOTER_BED="${arg#--promoter-bed=}"
            ;;
        --control-bed=*)
            CONTROL_BED="${arg#--control-bed=}"
            ;;
        --control-tag-dir=*)
            CONTROL_TAG_DIR="${arg#--control-tag-dir=}"
            ;;
        --primary-bowtie-index=*)
            PRIMARY_BOWTIE_INDEX="${arg#--primary-bowtie-index=}"
            ;;
        --mito-chrom=*)
            MITO_CHROM="${arg#--mito-chrom=}"
            ;;
        --spikein-bowtie-index=*)
            SPIKEIN_BOWTIE_INDEX="${arg#--spikein-bowtie-index=}"
            ;;
        --primary-genome-fasta=*)
            PRIMARY_GENOME_FASTA="${arg#--primary-genome-fasta=}"
            ;;
        --primary-gtf-file=*)
            PRIMARY_GTF_FILE="${arg#--primary-gtf-file=}"
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
    echo "ERROR: lookup_file is required. Usage: submit_pipeline.sh <lookup_file> --primary-bowtie-index=PATH --mito-chrom=NAME --primary-genome-fasta=PATH [--bam-type=nuclear|full] [--start-at=STEP] [--authors-orig] [--align_dmel --spikein-bowtie-index=PATH [--tss-bed=TSS.bed [--control-coverage=CONTROL.bedgraph]]] [--filter_genomecov --include-regions=REGIONS.bed [--promoter-bed=PROMOTERS.bed --control-bed=CONTROL.bed]] [--control-tag-dir=PATH --primary-gtf-file=PATH]"
    exit 1
fi


if [[ "${BAM_TYPE}" != "nuclear" && "${BAM_TYPE}" != "full" ]]; then
    echo "ERROR: --bam-type must be 'nuclear' or 'full' (got '${BAM_TYPE}')"
    exit 1
fi

if [[ "${FILTER_GENOMECOV}" == "true" && -z "${INCLUDE_REGIONS_BED}" ]]; then
    echo "ERROR: --filter_genomecov requires --include-regions=<path/to/regions.bed>"
    exit 1
fi

if [[ -n "${INCLUDE_REGIONS_BED}" && ! -f "${INCLUDE_REGIONS_BED}" ]]; then
    echo "ERROR: --include-regions file not found: ${INCLUDE_REGIONS_BED}"
    exit 1
fi

# --tss-bed / --promoter-bed / --control-bed are optional - if omitted, the
# corresponding scoring step is simply skipped (01a/01b or 01c/01d still run
# on their own via --align_dmel / --filter_genomecov). If provided, though,
# they should actually exist.
if [[ -n "${TSS_BED}" && ! -f "${TSS_BED}" ]]; then
    echo "ERROR: --tss-bed file not found: ${TSS_BED}"
    exit 1
fi
if [[ -n "${CONTROL_COVERAGE}" && ! -f "${CONTROL_COVERAGE}" ]]; then
    echo "ERROR: --control-coverage file not found: ${CONTROL_COVERAGE}"
    exit 1
fi
if [[ -n "${PROMOTER_BED}" && ! -f "${PROMOTER_BED}" ]]; then
    echo "ERROR: --promoter-bed file not found: ${PROMOTER_BED}"
    exit 1
fi
if [[ -n "${CONTROL_BED}" && ! -f "${CONTROL_BED}" ]]; then
    echo "ERROR: --control-bed file not found: ${CONTROL_BED}"
    exit 1
fi

# ============================================================================
# PIPELINE DAG
# ============================================================================
# Ordered list of steps. Order matters for --start-at (anything with a lower
# index than START_STEP is skipped), but the actual submitted --dependency
# for each step comes from STEP_DEPS below, not from this ordering alone -
# so branches (01a/02 off 01_align, 04/05 off 03_findpeaks) are handled
# correctly regardless of where --start-at lands.
STEP_ORDER=(01_align 01a_map_to_dmel 01b_dmel_coverage 01c_filter_bam 01d_genomecov_5p 02_maketagdir_samples 03_findpeaks 04_pos2bed 05_annotatepeaks 06_hahn_region_scoring 07_promoter_scoring 08_multiqc)

declare -A STEP_SCRIPT=(
    [01_align]="01_align.sh"
    [01a_map_to_dmel]="01a_map_to_dmel.sh"
    [01b_dmel_coverage]="01b_dmel_normalized_coverage.sh"
    [01c_filter_bam]="01c_filter_bam.sh"
    [01d_genomecov_5p]="01d_genomecov_5p.sh"
    [02_maketagdir_samples]="02_maketagdir_samples.sh"
    [03_findpeaks]="03_findpeaks.sh"
    [04_pos2bed]="04_pos2bed.sh"
    [05_annotatepeaks]="05_annotatepeaks.sh"
    [06_hahn_region_scoring]="06_hahn_region_scoring.sh"
    [07_promoter_scoring]="07_promoter_scoring.sh"
    [08_multiqc]="08_multiqc.sh"
)

# Upstream dependencies for each step (space-separated step names, or empty)
declare -A STEP_DEPS=(
    [01_align]=""
    [01a_map_to_dmel]="01_align"
    [01b_dmel_coverage]="01_align 01a_map_to_dmel"
    [01c_filter_bam]="01_align"
    [01d_genomecov_5p]="01c_filter_bam"
    [02_maketagdir_samples]="01_align"
    [03_findpeaks]="02_maketagdir_samples"
    [04_pos2bed]="03_findpeaks"
    [05_annotatepeaks]="03_findpeaks"
    [06_hahn_region_scoring]="01b_dmel_coverage 04_pos2bed"
    [07_promoter_scoring]="01d_genomecov_5p"
    [08_multiqc]="01b_dmel_coverage 01d_genomecov_5p 06_hahn_region_scoring 07_promoter_scoring 04_pos2bed 05_annotatepeaks"
)

# 08_multiqc uses afterany (partial upstream failure still produces a report);
# every other step uses afterok.
declare -A STEP_DEP_TYPE=(
    [08_multiqc]="afterany"
)

# Steps that are single jobs, not SLURM arrays
declare -A STEP_IS_SINGLE=(
    [08_multiqc]=1
)

# ============================================================================
# VALIDATION
# ============================================================================
if [[ ! -f "${LOOKUP_FILE}" ]]; then
    echo "ERROR: ${LOOKUP_FILE} not found"
    exit 1
fi

# --control-tag-dir is optional - if omitted, 03_findpeaks/04_pos2bed/
# 05_annotatepeaks (and 06_hahn_region_scoring, which reads 04_pos2bed's
# output too) are simply skipped rather than requiring a control tag
# directory to exist at some default location. If provided, though, it
# should actually exist and look complete.
if [[ -n "${CONTROL_TAG_DIR}" ]]; then
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
fi

# --primary-bowtie-index, --mito-chrom, and --primary-genome-fasta are
# required unconditionally: 01_align.sh and 02_maketagdir_samples.sh both
# always run, and both need them (no default in either script - see their
# own headers).
if [[ -z "${PRIMARY_BOWTIE_INDEX}" ]]; then
    echo "ERROR: --primary-bowtie-index=<path> is required (01_align.sh always runs and needs a bowtie2 index)"
    exit 1
fi
if [[ -z "${MITO_CHROM}" ]]; then
    echo "ERROR: --mito-chrom=<name> is required (01_align.sh always runs and needs the mitochondrial chromosome name for its nuclear/mito BAM split - this varies by organism/genome build, e.g. 'chrM' for many S. cerevisiae assemblies)"
    exit 1
fi
if [[ -z "${PRIMARY_GENOME_FASTA}" ]]; then
    echo "ERROR: --primary-genome-fasta=<path> is required (02_maketagdir_samples.sh always runs and needs it)"
    exit 1
fi
if [[ ! -f "${PRIMARY_GENOME_FASTA}" ]]; then
    echo "ERROR: --primary-genome-fasta file not found: ${PRIMARY_GENOME_FASTA}"
    exit 1
fi

# --spikein-bowtie-index is only required if --align_dmel is set, since
# 01a_map_to_dmel.sh (which needs it) is only submitted in that case - see
# the user's explicit request: no error for a missing dmel index when the
# dmel branch (and everything downstream of it, incl. 06_hahn_region_scoring)
# isn't even running.
if [[ "${ALIGN_DMEL}" == "true" && -z "${SPIKEIN_BOWTIE_INDEX}" ]]; then
    echo "ERROR: --align_dmel requires --spikein-bowtie-index=<path> (01a_map_to_dmel.sh needs a D. melanogaster bowtie2 index)"
    exit 1
fi

# --primary-gtf-file is only required if --control-tag-dir is set, since
# 05_annotatepeaks.sh (which needs it) is only submitted in that case.
if [[ -n "${CONTROL_TAG_DIR}" && -z "${PRIMARY_GTF_FILE}" ]]; then
    echo "ERROR: --control-tag-dir requires --primary-gtf-file=<path> (05_annotatepeaks.sh needs a GTF for peak annotation)"
    exit 1
fi
if [[ -n "${PRIMARY_GTF_FILE}" && ! -f "${PRIMARY_GTF_FILE}" ]]; then
    echo "ERROR: --primary-gtf-file file not found: ${PRIMARY_GTF_FILE}"
    exit 1
fi

# Bowtie2 indexes are a filename PREFIX, not a single real file (the actual
# files are <prefix>.1.bt2, <prefix>.2.bt2, etc., or .bt2l for large
# genomes) - a strict -f check on the prefix itself would always fail, so
# this is a soft warning rather than a hard requirement.
if [[ -n "${PRIMARY_BOWTIE_INDEX}" && ! -f "${PRIMARY_BOWTIE_INDEX}.1.bt2" && ! -f "${PRIMARY_BOWTIE_INDEX}.1.bt2l" ]]; then
    echo "WARNING: --primary-bowtie-index=${PRIMARY_BOWTIE_INDEX} - no ${PRIMARY_BOWTIE_INDEX}.1.bt2(l) found; double check this is the correct bowtie2 index prefix"
fi
if [[ -n "${SPIKEIN_BOWTIE_INDEX}" && ! -f "${SPIKEIN_BOWTIE_INDEX}.1.bt2" && ! -f "${SPIKEIN_BOWTIE_INDEX}.1.bt2l" ]]; then
    echo "WARNING: --spikein-bowtie-index=${SPIKEIN_BOWTIE_INDEX} - no ${SPIKEIN_BOWTIE_INDEX}.1.bt2(l) found; double check this is the correct bowtie2 index prefix"
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

if [[ ( "${START_STEP}" == "01a_map_to_dmel" || "${START_STEP}" == "01b_dmel_coverage" ) && "${ALIGN_DMEL}" != "true" ]]; then
    echo "ERROR: --start-at=${START_STEP} requires --align_dmel (dmel steps are disabled by default)"
    exit 1
fi

if [[ ( "${START_STEP}" == "01c_filter_bam" || "${START_STEP}" == "01d_genomecov_5p" ) && "${FILTER_GENOMECOV}" != "true" ]]; then
    echo "ERROR: --start-at=${START_STEP} requires --filter_genomecov (disabled by default)"
    exit 1
fi

if [[ ( "${START_STEP}" == "03_findpeaks" || "${START_STEP}" == "04_pos2bed" || "${START_STEP}" == "05_annotatepeaks" ) && -z "${CONTROL_TAG_DIR}" ]]; then
    echo "ERROR: --start-at=${START_STEP} requires --control-tag-dir=<path> (these steps are skipped without it)"
    exit 1
fi

if [[ "${START_STEP}" == "06_hahn_region_scoring" && ( -z "${TSS_BED}" || -z "${CONTROL_TAG_DIR}" ) ]]; then
    echo "ERROR: --start-at=06_hahn_region_scoring requires both --tss-bed=<path> and --control-tag-dir=<path> (step is skipped without either - it reads 04_pos2bed's output, which needs a control tag directory to exist)"
    exit 1
fi

if [[ "${START_STEP}" == "07_promoter_scoring" && ( -z "${PROMOTER_BED}" || -z "${CONTROL_BED}" ) ]]; then
    echo "ERROR: --start-at=07_promoter_scoring requires --promoter-bed=<path> and --control-bed=<path> (step is skipped without both)"
    exit 1
fi

# Get array size (total lines minus header)
ARRAY_SIZE=$(tail -n +2 "${LOOKUP_FILE}" | wc -l)

if [[ $ARRAY_SIZE -eq 0 ]]; then
    echo "ERROR: No samples found in ${LOOKUP_FILE}"
    exit 1
fi

# 06_hahn_region_scoring and 07_promoter_scoring are array jobs indexed by
# UNIQUE regulator (one task per regulator, not one per lookup row/replicate
# - each task's own R script invocation already discovers and combines that
# regulator's replicates internally). This is a different array size than
# ARRAY_SIZE above (which is per-row, used by 01-05).
UNIQUE_REGULATOR_COUNT=$(tail -n +2 "${LOOKUP_FILE}" | cut -f1 | sort -u | wc -l)

echo "================================"
echo "ChEC-seq Pipeline Submission"
echo "================================"
echo "Lookup file: ${LOOKUP_FILE}"
echo "Format: Paired-end (regulator, replicate, fastq_1, fastq_2)"
echo "Array size (# of samples): ${ARRAY_SIZE}"
echo "Unique regulators (for 06/07 array size): ${UNIQUE_REGULATOR_COUNT}"
echo "Control tag directory: $( [[ -n "${CONTROL_TAG_DIR}" ]] && echo "${CONTROL_TAG_DIR}" || echo "not provided - 03_findpeaks/04_pos2bed/05_annotatepeaks/06_hahn_region_scoring will be skipped (pass --control-tag-dir=<path> to enable)" )"
echo "Primary bowtie2 index: ${PRIMARY_BOWTIE_INDEX}"
echo "Mitochondrial chromosome name: ${MITO_CHROM}"
echo "Spike-in (dmel) bowtie2 index: $( [[ -n "${SPIKEIN_BOWTIE_INDEX}" ]] && echo "${SPIKEIN_BOWTIE_INDEX}" || echo "not provided (only needed with --align_dmel)" )"
echo "Genome FASTA: ${PRIMARY_GENOME_FASTA}"
echo "GTF file: $( [[ -n "${PRIMARY_GTF_FILE}" ]] && echo "${PRIMARY_GTF_FILE}" || echo "not provided (only needed with --control-tag-dir)" )"
echo "BAM type for tag directories: ${BAM_TYPE}"
echo "D. melanogaster spike-in steps: $( [[ "${ALIGN_DMEL}" == "true" ]] && echo "enabled (--align_dmel)" || echo "disabled (default - pass --align_dmel to enable)" )"
echo "Filter+genomecov 5' cut-site steps: $( [[ "${FILTER_GENOMECOV}" == "true" ]] && echo "enabled (--filter_genomecov, regions=${INCLUDE_REGIONS_BED})" || echo "disabled (default - pass --filter_genomecov --include-regions=... to enable)" )"
echo "Promoter scoring (dmel/Mahendrawada-Hahn): $( [[ -n "${TSS_BED}" ]] && echo "enabled (--tss-bed=${TSS_BED})" || echo "skipped (no --tss-bed provided)" )"
echo "Promoter scoring (calling-cards style): $( [[ -n "${PROMOTER_BED}" && -n "${CONTROL_BED}" ]] && echo "enabled (--promoter-bed=${PROMOTER_BED}, --control-bed=${CONTROL_BED})" || echo "skipped (--promoter-bed/--control-bed not both provided)" )"
echo "Read filtering (02_maketagdir_samples): $( [[ "${AUTHORS_ORIG}" == "true" ]] && echo "-keepAll (--authors-orig)" || echo "-unique -mapq 10 (default)" )"
echo "Fragment length (02_maketagdir_samples): $( [[ "${AUTHORS_ORIG}" == "true" ]] && echo "HOMER auto-estimate (--authors-orig)" || echo "from samtools stats (default)" )"
echo "Genome size (03_findpeaks): $( [[ "${AUTHORS_ORIG}" == "true" ]] && echo "HOMER auto-estimate (--authors-orig)" || echo "${BAM_TYPE}-based constant (default)" )"
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

    if [[ ( "${step}" == "01a_map_to_dmel" || "${step}" == "01b_dmel_coverage" ) && "${ALIGN_DMEL}" != "true" ]]; then
        echo "-- Skipping ${step} (dmel alignment disabled by default; pass --align_dmel to enable)"
        continue
    fi

    if [[ ( "${step}" == "01c_filter_bam" || "${step}" == "01d_genomecov_5p" ) && "${FILTER_GENOMECOV}" != "true" ]]; then
        echo "-- Skipping ${step} (filter+genomecov branch disabled by default; pass --filter_genomecov --include-regions=... to enable)"
        continue
    fi

    if [[ ( "${step}" == "03_findpeaks" || "${step}" == "04_pos2bed" || "${step}" == "05_annotatepeaks" ) && -z "${CONTROL_TAG_DIR}" ]]; then
        echo "-- Skipping ${step} (no --control-tag-dir provided; 02_maketagdir_samples still ran)"
        continue
    fi

    if [[ "${step}" == "06_hahn_region_scoring" && ( -z "${TSS_BED}" || -z "${CONTROL_TAG_DIR}" ) ]]; then
        echo "-- Skipping ${step} (requires both --tss-bed and --control-tag-dir; 01a/01b still ran if --align_dmel was set)"
        continue
    fi

    if [[ "${step}" == "07_promoter_scoring" && ( -z "${PROMOTER_BED}" || -z "${CONTROL_BED}" ) ]]; then
        echo "-- Skipping ${step} (--promoter-bed and --control-bed both required; 01c/01d still ran if --filter_genomecov was set)"
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

    # Extra positional/flag args each script needs. 01c_filter_bam has a
    # different positional order (lookup, include_regions, bam_type) than
    # everything else (lookup, bam_type). 06/07 take LOOKUP_FILE as their
    # first positional (to derive the array task's regulator from), then
    # flag-style args pointing at their required BED files.
    if [[ "${step}" == "01c_filter_bam" ]]; then
        extra_args=("${LOOKUP_FILE}" "${INCLUDE_REGIONS_BED}" "${BAM_TYPE}")
    elif [[ "${step}" == "06_hahn_region_scoring" ]]; then
        extra_args=("${LOOKUP_FILE}" "--tss-bed=${TSS_BED}")
        if [[ -n "${CONTROL_COVERAGE}" ]]; then
            extra_args+=("--control-coverage=${CONTROL_COVERAGE}")
        fi
    elif [[ "${step}" == "07_promoter_scoring" ]]; then
        extra_args=("${LOOKUP_FILE}" "--promoter-bed=${PROMOTER_BED}" "--control-bed=${CONTROL_BED}")
    elif [[ "${step}" == "08_multiqc" ]]; then
        extra_args=()
    else
        extra_args=("${LOOKUP_FILE}")
        if [[ "${step}" == "01b_dmel_coverage" || "${step}" == "01d_genomecov_5p" || "${step}" == "02_maketagdir_samples" || "${step}" == "03_findpeaks" ]]; then
            extra_args+=("${BAM_TYPE}")
        fi
        if [[ ( "${step}" == "02_maketagdir_samples" || "${step}" == "03_findpeaks" ) && "${AUTHORS_ORIG}" == "true" ]]; then
            extra_args+=("--authors_orig")
        fi
        if [[ "${step}" == "03_findpeaks" ]]; then
            extra_args+=("--control-tag-dir=${CONTROL_TAG_DIR}")
        fi
        if [[ "${step}" == "01_align" ]]; then
            extra_args+=("--bowtie-index=${PRIMARY_BOWTIE_INDEX}" "--mito-chrom=${MITO_CHROM}")
        fi
        if [[ "${step}" == "01a_map_to_dmel" ]]; then
            extra_args+=("--bowtie-index=${SPIKEIN_BOWTIE_INDEX}")
        fi
        if [[ "${step}" == "02_maketagdir_samples" ]]; then
            extra_args+=("--genome-fasta=${PRIMARY_GENOME_FASTA}")
        fi
        if [[ "${step}" == "05_annotatepeaks" ]]; then
            extra_args+=("--genome-fasta=${PRIMARY_GENOME_FASTA}" "--gtf-file=${PRIMARY_GTF_FILE}")
        fi
    fi

    # 06/07 wrap an R script sitting alongside them and need to locate it at
    # runtime. Under sbatch, SLURM copies the submitted script into a
    # per-job spool directory and runs THAT COPY, so the script's own
    # BASH_SOURCE-based self-location resolves to the spool path, not here -
    # export the (correctly-resolved, since we're running via bash, not
    # sbatch) SCRIPT_DIR so the wrapper can find its sibling R script.
    if [[ "${step}" == "06_hahn_region_scoring" || "${step}" == "07_promoter_scoring" ]]; then
        dep_args+=("--export=ALL,PIPELINE_SCRIPT_DIR=${SCRIPT_DIR}")
    fi

    if [[ -n "${STEP_IS_SINGLE[$step]:-}" ]]; then
        echo "Submitting ${step} (single job)..."
        echo "  Command: sbatch ${dep_args[*]:-} ${script} ${extra_args[*]:-}"
        JOBIDS[$step]=$(submit_job "${dep_args[@]}" "${script}" "${extra_args[@]}")
    else
        if [[ "${step}" == "06_hahn_region_scoring" || "${step}" == "07_promoter_scoring" ]]; then
            step_array_size="${UNIQUE_REGULATOR_COUNT}"
        else
            step_array_size="${ARRAY_SIZE}"
        fi
        echo "Submitting ${step} (array 1-${step_array_size})..."
        echo "  Command: sbatch ${dep_args[*]:-} --array=1-${step_array_size} ${script} ${extra_args[*]}"
        JOBIDS[$step]=$(submit_job "${dep_args[@]}" --array=1-${step_array_size} "${script}" "${extra_args[@]}")
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
echo "  ($( [[ -n "${CONTROL_TAG_DIR}" ]] && echo "control tag dir: ${CONTROL_TAG_DIR} (already built manually via maketagdir_control.sh)" || echo "no --control-tag-dir provided - peak-calling steps were skipped" ))"
for step in "${STEP_ORDER[@]}"; do
    if [[ -n "${JOBIDS[$step]:-}" ]]; then
        echo "  ${step}: ${JOBIDS[$step]}"
    elif [[ ( "${step}" == "01a_map_to_dmel" || "${step}" == "01b_dmel_coverage" ) && "${ALIGN_DMEL}" != "true" ]]; then
        echo "  ${step}: (disabled - pass --align_dmel to enable)"
    elif [[ ( "${step}" == "01c_filter_bam" || "${step}" == "01d_genomecov_5p" ) && "${FILTER_GENOMECOV}" != "true" ]]; then
        echo "  ${step}: (disabled - pass --filter_genomecov --include-regions=... to enable)"
    elif [[ ( "${step}" == "03_findpeaks" || "${step}" == "04_pos2bed" || "${step}" == "05_annotatepeaks" ) && -z "${CONTROL_TAG_DIR}" ]]; then
        echo "  ${step}: (skipped - no --control-tag-dir provided)"
    elif [[ "${step}" == "06_hahn_region_scoring" && ( -z "${TSS_BED}" || -z "${CONTROL_TAG_DIR}" ) ]]; then
        echo "  ${step}: (skipped - requires both --tss-bed and --control-tag-dir)"
    elif [[ "${step}" == "07_promoter_scoring" && ( -z "${PROMOTER_BED}" || -z "${CONTROL_BED}" ) ]]; then
        echo "  ${step}: (skipped - --promoter-bed/--control-bed not both provided)"
    else
        echo "  ${step}: (skipped - assumed already complete)"
    fi
done
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  tail -f logs/align_1.log"
echo ""
if [[ -n "${JOBIDS[08_multiqc]:-}" ]]; then
    echo "Track final job:"
    echo "  squeue -j ${JOBIDS[08_multiqc]}"
    echo ""
    echo "When 08_multiqc finishes, the consolidated report is at:"
    echo "  results/multiqc/multiqc_report.html"
fi
echo ""
