#!/bin/bash
#SBATCH --job-name=chec_map_dmel
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o logs/map_dmel_%a.log
#SBATCH -e logs/map_dmel_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/bowtie2_samtools_gzip:0c7fd2d5085a1394

set -uo pipefail

# ============================================================================
# CONFIGURATION - Fill in these paths before running
# ============================================================================
LOOKUP_FILE="${1:-samples.tsv}"
DMEL_BOWTIE_INDEX="/ref/mblab/data/dmelanogaster/bowtie2_index/dmel-all-chromosome-r6.65"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# BOWTIE PARAMETERS (same as S. cerevisiae)
# ============================================================================
BOWTIE_PARAMS="-I 10 -X 700 --local --very-sensitive-local --no-unal --no-mixed --no-discordant"

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}/bams" "${LOG_DIR}"

# Get the task line from lookup file (skip header)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

if [[ -z "$LINE" ]]; then
    echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    exit 1
fi

# Parse TSV: regulator_symbol, replicate, fastq_1, fastq_2
read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

# All dmel outputs live alongside the other BAMs in bams/{reg}/{rep}/
SAMPLE_DIR="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}"
UNMAPPED_R1="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_unmapped_R1.fastq.gz"
UNMAPPED_R2="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_unmapped_R2.fastq.gz"
DMEL_BAM="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel.bam"
DMEL_COUNTS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel_counts.txt"

echo "Task ${SLURM_ARRAY_TASK_ID}: Mapping unmapped reads to D. melanogaster for ${REGULATOR}_${REPLICATE}"
echo "  Unmapped R1: ${UNMAPPED_R1}"
echo "  Unmapped R2: ${UNMAPPED_R2}"

# ============================================================================
# VALIDATE INPUT
# ============================================================================
mkdir -p "${SAMPLE_DIR}"

if [[ ! -f "${UNMAPPED_R1}" ]] || [[ ! -f "${UNMAPPED_R2}" ]]; then
    echo "WARNING: Unmapped FASTQ files not found"
    [[ ! -f "${UNMAPPED_R1}" ]] && echo "  Missing: ${UNMAPPED_R1}"
    [[ ! -f "${UNMAPPED_R2}" ]] && echo "  Missing: ${UNMAPPED_R2}"
    echo "  Skipping D. melanogaster mapping for this sample"
    echo "0" > "${DMEL_COUNTS_FILE}"
    exit 0
fi

# Count unmapped reads
UNMAPPED_COUNT=$(gunzip -c "${UNMAPPED_R1}" 2>/dev/null | grep -c "^@" || echo "0")

if [[ "$UNMAPPED_COUNT" -eq 0 ]]; then
    echo "No unmapped reads found, skipping D. melanogaster alignment"
    echo "0" > "${DMEL_COUNTS_FILE}"
    exit 0
fi

echo "  Unmapped read pairs: ${UNMAPPED_COUNT}"

# ============================================================================
# ALIGN TO D. MELANOGASTER
# ============================================================================
echo "Aligning unmapped reads to D. melanogaster..."
bowtie2 \
    -p 8 \
    -q \
    -x "${DMEL_BOWTIE_INDEX}" \
    ${BOWTIE_PARAMS} \
    -1 <(gunzip -c "${UNMAPPED_R1}") \
    -2 <(gunzip -c "${UNMAPPED_R2}") \
    2> "${LOG_DIR}/${REGULATOR}_${REPLICATE}_dmel_bowtie.log" | \
    samtools view -b -h -F 4 - | \
    samtools sort -@ 4 -o "${DMEL_BAM}" -

samtools index "${DMEL_BAM}"

# ============================================================================
# GENERATE STATISTICS
# ============================================================================
DMEL_STATS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel_samtools_stats.txt"
DMEL_FLAGSTATS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel_samtools_flagstats.txt"
DMEL_COVERAGE_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel_coverage.txt"
DMEL_IDXSTATS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_dmel_idxstats.txt"

echo "Generating D. melanogaster alignment statistics..."
samtools stats "${DMEL_BAM}" > "${DMEL_STATS_FILE}"
samtools flagstats "${DMEL_BAM}" > "${DMEL_FLAGSTATS_FILE}"
samtools coverage "${DMEL_BAM}" > "${DMEL_COVERAGE_FILE}"
samtools idxstats "${DMEL_BAM}" > "${DMEL_IDXSTATS_FILE}"

# ============================================================================
# COUNT D. MELANOGASTER READS
# ============================================================================
DMEL_MAPPED_COUNT=$(samtools view -c "${DMEL_BAM}")

echo ""
echo "D. melanogaster mapping results for ${REGULATOR}_${REPLICATE}:"
echo "  Unmapped S. cerevisiae reads: ${UNMAPPED_COUNT}"
echo "  Mapped to D. melanogaster: ${DMEL_MAPPED_COUNT}"
echo ""

# Write count to file for use in downstream spike-in normalization
echo -e "${REGULATOR}\t${REPLICATE}\t${DMEL_MAPPED_COUNT}" > "${DMEL_COUNTS_FILE}"

echo "✓ D. melanogaster mapping complete"
echo "  BAM: ${DMEL_BAM}"
echo "  Stats: ${DMEL_STATS_FILE}"
echo "  Flagstats: ${DMEL_FLAGSTATS_FILE}"
echo "  Coverage: ${DMEL_COVERAGE_FILE}"
echo "  Idxstats: ${DMEL_IDXSTATS_FILE}"
echo "  Counts file: ${DMEL_COUNTS_FILE}"
