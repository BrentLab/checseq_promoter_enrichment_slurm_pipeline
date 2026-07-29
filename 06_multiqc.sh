#!/bin/bash
#SBATCH --job-name=chec_multiqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH -o logs/multiqc.log
#SBATCH -e logs/multiqc.log
#SBATCH --container=docker://ghcr.io/multiqc/multiqc:pdf-dev

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================
# Scans the results/ and logs/ directories for QC-relevant outputs produced by
# the pipeline (samtools stats/flagstats/idxstats, bowtie2 logs, HOMER tag
# directory QC) and builds a single consolidated MultiQC report.
OUTPUT_DIR="results"
LOG_DIR="logs"
MULTIQC_OUT="${OUTPUT_DIR}/multiqc"

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${MULTIQC_OUT}"

echo "================================"
echo "MultiQC Report Generation"
echo "================================"
echo "  Scanning: ${OUTPUT_DIR}/ and ${LOG_DIR}/"
echo "  Output:   ${MULTIQC_OUT}/"
echo ""

# ============================================================================
# RUN MULTIQC
# ============================================================================
# MultiQC auto-detects and parses (among others):
#   - samtools stats / flagstats / idxstats  (*.txt in results/bams/**)
#   - bowtie2 alignment logs                 (logs/*_bowtie.log)
#   - HOMER tag directory QC                 (results/tag_dirs/**/tagInfo.txt etc.)
# We point it at both results/ and logs/ so it picks up everything in one pass.
multiqc \
    "${OUTPUT_DIR}" \
    "${LOG_DIR}" \
    --outdir "${MULTIQC_OUT}" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc_run.log"

# ============================================================================
# VALIDATE OUTPUT
# ============================================================================
REPORT="${MULTIQC_OUT}/multiqc_report.html"
if [[ ! -f "${REPORT}" ]]; then
    echo "ERROR: MultiQC did not produce a report at ${REPORT}"
    exit 1
fi

echo ""
echo "✓ MultiQC report complete"
echo "  Report: ${REPORT}"
echo "  Data:   ${MULTIQC_OUT}/multiqc_data/"
