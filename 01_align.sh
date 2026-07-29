#!/bin/bash
#SBATCH --job-name=chec_align
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o logs/align_%a.log
#SBATCH -e logs/align_%a.log
#SBATCH --container=oras://community.wave.seqera.io/library/bowtie2_samtools_gzip:0c7fd2d5085a1394

set -uo pipefail

# ============================================================================
# CONFIGURATION - Fill in these paths before running
# ============================================================================
LOOKUP_FILE="$1"
SCER_BOWTIE_INDEX="/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/bowtie2_index/S288C_reference_sequence_R64-5-1_20240529_chr_normalized"
OUTPUT_DIR="results"
LOG_DIR="logs"

# ============================================================================
# BOWTIE PARAMETERS (from methods) - paired-end
# ============================================================================
BOWTIE_PARAMS="-I 10 -X 700 --local --very-sensitive-local --no-unal --no-mixed --no-discordant"

# ============================================================================
# SETUP
# ============================================================================
mkdir -p "${OUTPUT_DIR}/bams" "${LOG_DIR}"

# Create temporary directory for intermediate files
TMPDIR="${TMPDIR:-/tmp}"
WORK_TMP=$(mktemp -d "${TMPDIR}/chec_align_${SLURM_ARRAY_TASK_ID}.XXXXXX")
trap "rm -rf ${WORK_TMP}" EXIT

echo "Using temporary directory: ${WORK_TMP}"

# Get the task line from lookup file (skip header)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${LOOKUP_FILE}")

if [[ -z "$LINE" ]]; then
    echo "ERROR: Could not read line ${SLURM_ARRAY_TASK_ID} from ${LOOKUP_FILE}"
    exit 1
fi

# Parse TSV columns: regulator_symbol, replicate, fastq_1, fastq_2
read -r REGULATOR REPLICATE FASTQ_R1 FASTQ_R2 <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: Aligning ${REGULATOR}_${REPLICATE}"
echo "  R1: ${FASTQ_R1}"
echo "  R2: ${FASTQ_R2}"

# ============================================================================
# VALIDATE INPUT FILES
# ============================================================================
if [[ ! -f "$FASTQ_R1" ]]; then
    echo "ERROR: FASTQ R1 not found: $FASTQ_R1"
    exit 1
fi

if [[ ! -f "$FASTQ_R2" ]]; then
    echo "ERROR: FASTQ R2 not found: $FASTQ_R2"
    exit 1
fi

# ============================================================================
# ESTIMATE READ LENGTH (from R1)
# ============================================================================
echo "Estimating read length from R1..."

if [[ "$FASTQ_R1" == *.gz ]]; then
    echo "  (gzipped file, using gunzip)"
    READ_LENGTH=$(gunzip -c "$FASTQ_R1" 2>/dev/null | head -n 1000 | \
        awk 'NR % 4 == 2 { sum += length($0); count++ } END { if (count > 0) print int(sum / count) }')
else
    echo "  (uncompressed file)"
    READ_LENGTH=$(head -n 1000 "$FASTQ_R1" | \
        awk 'NR % 4 == 2 { sum += length($0); count++ } END { if (count > 0) print int(sum / count) }')
fi

if [[ -z "$READ_LENGTH" ]] || [[ "$READ_LENGTH" -lt 1 ]]; then
    echo "ERROR: Could not determine read length from FASTQ"
    exit 1
fi

echo "Detected read length: ${READ_LENGTH} bp"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY STRUCTURE
# ============================================================================
SAMPLE_DIR="${OUTPUT_DIR}/bams/${REGULATOR}/${REPLICATE}"
mkdir -p "${SAMPLE_DIR}"

echo "Output directory: ${SAMPLE_DIR}"

# ============================================================================
# ALIGNMENT WITH BOWTIE2 (PAIRED-END)
# ============================================================================
TEMP_SAM="${WORK_TMP}/${REGULATOR}_${REPLICATE}_all.sam"
FULL_BAM="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}.bam"
NUCLEAR_BAM="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_nuclear.bam"
MITO_BAM="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_mito.bam"
STATS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_samtools_stats.txt"
FLAGSTATS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_samtools_flagstats.txt"

# bowtie2's --un-conc-gz writes pairs that fail to align concordantly directly
# to gzipped FASTQ. This is necessary because --no-unal/--no-mixed/--no-discordant
# mean the SAM output itself never contains unmapped/discordant reads to extract
# after the fact.
UNMAPPED_PREFIX="${WORK_TMP}/${REGULATOR}_${REPLICATE}_unmapped.fastq.gz"
UNMAPPED_R1="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_unmapped_R1.fastq.gz"
UNMAPPED_R2="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_unmapped_R2.fastq.gz"

echo "Running bowtie2 alignment (paired-end)..."
echo "  Index: ${SCER_BOWTIE_INDEX}"
echo "  Params: ${BOWTIE_PARAMS}"

bowtie2 \
    -p 8 \
    -q \
    -x "${SCER_BOWTIE_INDEX}" \
    ${BOWTIE_PARAMS} \
    --un-conc-gz "${UNMAPPED_PREFIX}" \
    -1 "$FASTQ_R1" \
    -2 "$FASTQ_R2" \
    2> "${LOG_DIR}/${REGULATOR}_${REPLICATE}_bowtie.log" > "${TEMP_SAM}"

BOWTIE_EXIT=$?
if [[ $BOWTIE_EXIT -ne 0 ]]; then
    echo "ERROR: Bowtie2 failed with exit code $BOWTIE_EXIT"
    echo "Bowtie2 stderr:"
    cat "${LOG_DIR}/${REGULATOR}_${REPLICATE}_bowtie.log"
    exit 1
fi

echo "  Bowtie2 completed successfully"

echo "DEBUG: files in WORK_TMP matching unmapped pattern:"
ls -la "${WORK_TMP}"/*unmapped* 2>/dev/null || echo "  (none found)"

# bowtie2 inserts .1/.2 before the LAST extension in the path, not after it.
# So a prefix of "foo.fastq.gz" produces "foo.fastq.1.gz" and "foo.fastq.2.gz"
# (not "foo.fastq.gz.1.gz"). Use globbing to find them reliably rather than
# assuming the exact suffix.
UNMAPPED_PREFIX_BASE="${UNMAPPED_PREFIX%.gz}"

FOUND_R1=$(ls "${UNMAPPED_PREFIX_BASE}".1.gz 2>/dev/null | head -n 1)
FOUND_R2=$(ls "${UNMAPPED_PREFIX_BASE}".2.gz 2>/dev/null | head -n 1)

if [[ -n "${FOUND_R1}" && -f "${FOUND_R1}" ]]; then
    mv "${FOUND_R1}" "${UNMAPPED_R1}"
    R1_COUNT=$(( $(gunzip -c "${UNMAPPED_R1}" | wc -l) / 4 ))
    echo "  Unmapped R1: ${UNMAPPED_R1} (${R1_COUNT} reads)"
else
    echo "  WARNING: No unmapped R1 file found (looked for ${UNMAPPED_PREFIX_BASE}.1.gz)"
fi

if [[ -n "${FOUND_R2}" && -f "${FOUND_R2}" ]]; then
    mv "${FOUND_R2}" "${UNMAPPED_R2}"
    R2_COUNT=$(( $(gunzip -c "${UNMAPPED_R2}" | wc -l) / 4 ))
    echo "  Unmapped R2: ${UNMAPPED_R2} (${R2_COUNT} reads)"
else
    echo "  WARNING: No unmapped R2 file found (looked for ${UNMAPPED_PREFIX_BASE}.2.gz)"
fi

# ============================================================================
# PROCESS ALIGNMENT RESULTS
# ============================================================================
echo "Converting SAM to sorted BAM (all mapped reads, including chrM)..."
samtools view -b -h -F 4 "${TEMP_SAM}" | \
    samtools sort -@ 4 -o "${FULL_BAM}" -

echo "Indexing full BAM..."
samtools index "${FULL_BAM}"

# ============================================================================
# GENERATE STATISTICS ON FULL BAM (for QC - includes chrM)
# ============================================================================
echo "Generating samtools statistics (full BAM, includes chrM)..."
samtools stats "${FULL_BAM}" > "${STATS_FILE}"
samtools flagstats "${FULL_BAM}" > "${FLAGSTATS_FILE}"

COVERAGE_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_coverage.txt"
samtools coverage "${FULL_BAM}" > "${COVERAGE_FILE}"

IDXSTATS_FILE="${SAMPLE_DIR}/${REGULATOR}_${REPLICATE}_idxstats.txt"
samtools idxstats "${FULL_BAM}" > "${IDXSTATS_FILE}"

# ============================================================================
# SPLIT INTO NUCLEAR AND MITOCHONDRIAL BAMs
# ============================================================================
echo "Splitting into nuclear and mitochondrial BAMs..."

# Mitochondrial reads only (uses index, since BAM is indexed and coordinate-sorted)
samtools view -b -h "${FULL_BAM}" chrM > "${MITO_BAM}"
samtools index "${MITO_BAM}"

# Nuclear reads only (everything except chrM)
samtools view -h "${FULL_BAM}" | \
    awk '$0 ~ /^@/ || $3 != "chrM" {print}' | \
    samtools view -b -h - > "${NUCLEAR_BAM}"
samtools index "${NUCLEAR_BAM}"

MITO_COUNT=$(samtools view -c "${MITO_BAM}")
NUCLEAR_COUNT=$(samtools view -c "${NUCLEAR_BAM}")
echo "  Mitochondrial reads: ${MITO_COUNT} -> ${MITO_BAM}"
echo "  Nuclear reads: ${NUCLEAR_COUNT} -> ${NUCLEAR_BAM}"

# ============================================================================
# GATHER MAPPING STATS FROM SAMTOOLS (FULL BAM, INCLUDES chrM)
# ============================================================================
echo ""
echo "Mapping stats for ${REGULATOR}_${REPLICATE} (full BAM, includes chrM):"
grep "^SN" "${STATS_FILE}" | cut -f 2- | head -n 5

echo ""
echo "Flagstats:"
head -n 5 "${FLAGSTATS_FILE}"

TOTAL_READS=$(grep "^SN.*raw total sequences" "${STATS_FILE}" | awk -F'\t' '{print $3}')
TOTAL_READS="${TOTAL_READS:-0}"

# Guard against a non-numeric result (defensive, in case the stats format
# ever changes again) rather than letting it hit an unbound-variable/arithmetic
# error further down.
if ! [[ "${TOTAL_READS}" =~ ^[0-9]+$ ]]; then
    echo "WARNING: Could not parse a numeric total read count from ${STATS_FILE} (got: '${TOTAL_READS}'); defaulting to 0"
    TOTAL_READS=0
fi

MAPPED_READS=$(samtools view -c -F 4 "${FULL_BAM}")

if [[ "$TOTAL_READS" -gt 0 ]]; then
    MITO_PCT=$(awk "BEGIN {printf \"%.2f\", ($MITO_COUNT / $TOTAL_READS) * 100}")
else
    MITO_PCT="0.00"
fi

echo ""
echo "Summary:"
echo "  Total mapped reads: ${TOTAL_READS}"
echo "  Nuclear reads: ${NUCLEAR_COUNT}"
echo "  Mitochondrial reads: ${MITO_COUNT} (${MITO_PCT}%)"

echo ""
echo "✓ Alignment complete"
echo "  Full BAM (all chromosomes): ${FULL_BAM}"
echo "  Nuclear BAM: ${NUCLEAR_BAM}"
echo "  Mitochondrial BAM: ${MITO_BAM}"
echo "  Unmapped R1: ${UNMAPPED_R1}"
echo "  Unmapped R2: ${UNMAPPED_R2}"
echo "  Stats (full BAM): ${STATS_FILE}"
echo "  Flagstats (full BAM): ${FLAGSTATS_FILE}"
echo "  Coverage (full BAM): ${COVERAGE_FILE}"
echo "  Idxstats (full BAM): ${IDXSTATS_FILE}"
