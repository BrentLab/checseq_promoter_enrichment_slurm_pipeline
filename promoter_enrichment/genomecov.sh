#!/bin/bash
#SBATCH --time=0-00:30:00 #days-hh:mm:ss
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/cov.out
#SBATCH --error=logs/cov.err

# Enable safe script execution
set -euo pipefail

# Usage:
# sbatch --array=1-12 genomecov.sh lookup.txt output_dirname

# Load necessary modules
eval $(spack load --sh bedtools2)

# Read BAM file from lookup.txt based on the array task ID
if ! read -r bam_file < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$1"); then
    echo "Error: Unable to read BAM file from lookup.txt at line ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

if [[ -z "${bam_file// }" ]]; then
  echo "Error: BAM file path is empty or whitespace for task ID $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

if [[ ! -f "$bam_file" ]]; then
  echo "Error: BAM file does not exist at: $bam_file" >&2
  exit 1
fi

OUTPUT_DIR=$2

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define output filename (final sorted BED file)
BASENAME="$(basename "$bam_file" .bam)_r1_5p"
OUTPUT_BED="$OUTPUT_DIR/${BASENAME}.bed"

# Create temporary directory
TMP_DIR=$(mktemp -d -t genomecov_XXXXXX)

# Ensure cleanup on exit
trap "rm -rf $TMP_DIR" EXIT

echo "Processing: $bam_file"
echo "Temporary directory: $TMP_DIR"

# Generate coverage for plus strand
echo "Generating plus strand coverage..."
bedtools genomecov -ibam "$bam_file" -5 -dz -strand + > "$TMP_DIR/plus.cov"

# Generate coverage for minus strand
echo "Generating minus strand coverage..."
bedtools genomecov -ibam "$bam_file" -5 -dz -strand - > "$TMP_DIR/minus.cov"

# Convert plus strand to BED format
echo "Converting plus strand to BED..."
awk 'OFS="\t" {print $1, $2, $2+1, ".", $3, "+"}' "$TMP_DIR/plus.cov" > "$TMP_DIR/plus.bed"

# Convert minus strand to BED format
echo "Converting minus strand to BED..."
awk 'OFS="\t" {print $1, $2, $2+1, ".", $3, "-"}' "$TMP_DIR/minus.cov" > "$TMP_DIR/minus.bed"

# Combine both strands (unsorted)
echo "Combining strands..."
cat "$TMP_DIR/plus.bed" "$TMP_DIR/minus.bed" > "$TMP_DIR/combined_unsorted.bed"

# Sort the combined BED file
echo "Sorting combined BED file..."
bedtools sort -i "$TMP_DIR/combined_unsorted.bed" > "$OUTPUT_BED"

# Report results
PLUS_COUNT=$(wc -l < "$TMP_DIR/plus.bed")
MINUS_COUNT=$(wc -l < "$TMP_DIR/minus.bed")
TOTAL_COUNT=$(wc -l < "$OUTPUT_BED")

echo "Complete!"
echo "Plus strand positions: $PLUS_COUNT"
echo "Minus strand positions: $MINUS_COUNT"
echo "Total positions: $TOTAL_COUNT"
echo "Output file: $OUTPUT_BED"

