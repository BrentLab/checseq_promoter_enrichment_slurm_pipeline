#!/bin/bash
# Helper script to prepare the pipeline and compute SLURM array sizes
# Input: TSV file with columns: regulator_symbol, replicate, fastq_1, fastq_2

# Don't use set -e so we can see all output even if something fails
set -u

# ============================================================================
# CONFIGURATION
# ============================================================================
LOOKUP_FILE="${1:-samples.tsv}"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}/bams" "${OUTPUT_DIR}/tag_dirs" "${OUTPUT_DIR}/peaks" "${LOG_DIR}" || true

echo "================================"
echo "ChEC-seq Pipeline Preparation"
echo "================================"
echo ""

# Validate lookup file exists
if [[ ! -f "$LOOKUP_FILE" ]]; then
    echo "ERROR: Lookup file '$LOOKUP_FILE' not found"
    exit 1
fi

echo "Lookup file: $LOOKUP_FILE"
echo ""

# ============================================================================
# VALIDATE PAIRED-END FORMAT
# ============================================================================
echo "Validating lookup file format..."
echo "Expected columns: regulator_symbol, replicate, fastq_1, fastq_2"
echo ""

# Count lines (skip header if present)
TOTAL_SAMPLES=$(tail -n +2 "$LOOKUP_FILE" 2>/dev/null | wc -l)
echo "Total data rows found: $TOTAL_SAMPLES"
echo ""

if [[ $TOTAL_SAMPLES -eq 0 ]]; then
    echo "ERROR: No data rows found in lookup file"
    echo "File contents:"
    head -n 5 "$LOOKUP_FILE"
    exit 1
fi

# Validate each row
VALID_SAMPLES=0
INVALID_SAMPLES=0

echo "Validating each sample..."
while IFS=$'\t' read -r regulator_symbol replicate fastq_1 fastq_2; do
    # Skip header row
    if [[ "$regulator_symbol" == "regulator_symbol" ]]; then
        continue
    fi

    # Skip empty lines and comments
    if [[ -z "$regulator_symbol" ]]; then
        continue
    fi

    if [[ "$regulator_symbol" =~ ^# ]]; then
        continue
    fi

    # Check if both FASTQ files exist (handle gzipped files)
    if [[ -f "$fastq_1" ]] && [[ -f "$fastq_2" ]]; then
        ((VALID_SAMPLES++))
        echo "  ✓ ${regulator_symbol}_${replicate}"
    else
        echo "  ✗ ${regulator_symbol}_${replicate}"
        if [[ ! -f "$fastq_1" ]]; then
            echo "      Missing: $fastq_1"
        fi
        if [[ ! -f "$fastq_2" ]]; then
            echo "      Missing: $fastq_2"
        fi
        ((INVALID_SAMPLES++))
    fi
done < "$LOOKUP_FILE"

echo ""

# ============================================================================
# REPORT SUMMARY
# ============================================================================
echo "Pipeline Setup Summary:"
echo "  Lookup file: $LOOKUP_FILE"
echo "  Format: Paired-end (regulator, replicate, fastq_1, fastq_2)"
echo "  Output directory: $OUTPUT_DIR"
echo "  Log directory: $LOG_DIR"
echo "  Valid samples: $VALID_SAMPLES"
echo "  Invalid samples: $INVALID_SAMPLES"
echo ""

if [[ $VALID_SAMPLES -eq 0 ]]; then
    echo "ERROR: No valid sample pairs found in lookup"
    exit 1
fi

echo "✓ Validation complete!"
echo ""
echo "Array size for SLURM jobs: $VALID_SAMPLES"
echo ""
echo "Next steps:"
echo "  1. Verify SLURM configuration in each script"
echo "  2. Fill in parameterized paths in:"
echo "     - 01_align.sh: SCER_BOWTIE_INDEX"
echo "     - maketagdir_control.sh: CONTROL_BAM (manual prereq)"
echo "     - 01a_map_to_dmel.sh: DMEL_BOWTIE_INDEX"
echo ""
echo "Then run:"
echo "  bash submit_pipeline.sh"
echo ""
echo "✓ Pipeline ready for submission"
