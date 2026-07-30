# 5' Cut-Site Quantification (Region-Filtered)

An alternative, downstream-of-alignment quantification method: instead of
HOMER peak calling (tag pileups across a fragment-length window), this
computes single-nucleotide-resolution **5' cut-site density**, per strand,
restricted to a set of regions of interest. This is the classic
cut-site/footprint style of quantification used for ChEC-seq/MNase-type
data, where the informative signal is *where the cut happened* (the read's
5' end), not the shape of the full fragment pileup.

Two scripts, run in sequence:

| Script | Purpose |
|---|---|
| `filter_bam_improved.sh` | Filters a BAM to properly-paired, high-quality, region-of-interest reads |
| `genomecov.sh` | Computes per-base, per-strand 5' end coverage from a filtered BAM, output as BED |

## Requirements

- SLURM with `spack` module support
- `samtools` (loaded via `spack load samtools` inside `filter_bam_improved.sh`)
- `bedtools2` (loaded via `spack load bedtools2` inside `genomecov.sh`)

## Step 1: `filter_bam_improved.sh`

Filters an aligned BAM down to the reads you actually want to quantify:
properly-paired, mapping quality ≥10, overlapping a supplied BED of regions
of interest, with unmapped/secondary/supplementary alignments excluded.

```
sbatch --array=1-N filter_bam_improved.sh <lookup.txt> <output_dir> <include_regions.bed> [--single_end] [--keep_both_reads]
```

**Arguments:**
- `lookup.txt` — plain text, **no header**, one BAM file path per line.
  `SLURM_ARRAY_TASK_ID` picks the line (task 1 → line 1, etc.), so
  `--array=1-N` must match the number of lines in the file.
- `output_dir` — where filtered BAMs are written (created if it doesn't exist)
- `include_regions.bed` — a BED file of regions to restrict to (e.g.
  promoters, a specific locus set); reads outside these regions are dropped

**Optional flags** (mutually relevant only to paired-end data):
- *(default, no flag)* — paired-end, keeps **R1 only**, requires the pair to
  be properly paired. This is the typical mode for cut-site analysis, since
  R1's 5' end is the position of interest and R2 would just double-count it.
- `--keep_both_reads` — paired-end, keeps both R1 and R2 (still requires
  properly paired)
- `--single_end` — single-end mode: basic quality/mapping filtering only, no
  pairing requirement

**What gets filtered, concretely (via `samtools view -F/-f/-q/-L`):**

| Mode | Excludes (`-F`) | Requires (`-f`) | MAPQ | Regions |
|---|---|---|---|---|
| Single-end | unmapped, secondary, supplementary | — | ≥10 | `-L include_regions.bed` |
| Paired, R1 only (default) | unmapped, mate unmapped, secondary, supplementary | properly paired + read1 | ≥10 | `-L include_regions.bed` |
| Paired, both reads (`--keep_both_reads`) | unmapped, mate unmapped, secondary, supplementary | properly paired | ≥10 | `-L include_regions.bed` |

**Output:** `<output_dir>/<input_basename>.filtered.bam` (sorted + indexed),
plus a summary of total vs. filtered read counts and % retained printed to
the SLURM log.

## Step 2: `genomecov.sh`

Takes filtered BAMs (typically the output of step 1) and computes per-base
coverage of read 5' ends, separately for each strand, as a BED file.

```
sbatch --array=1-N genomecov.sh <lookup.txt> <output_dirname>
```

**Arguments:**
- `lookup.txt` — plain text, **no header**, one BAM path per line (point
  this at the `*.filtered.bam` files from step 1)
- `output_dirname` — where output BED files are written (created if it
  doesn't exist)

**What it does:**
1. Runs `bedtools genomecov -ibam <bam> -5 -dz -strand +` and `-strand -`
   separately, to get per-base counts of read 5' ends on each strand
2. Converts each strand's output to BED6 format: `chrom  pos  pos+1  .  count  strand`
3. Concatenates both strands and sorts the combined result

**Output:** `<output_dirname>/<bam_basename>_r1_5p.bed` — a sorted BED file
where each line is a single-nucleotide position with its 5'-end read count
and strand. The SLURM log reports plus/minus/total position counts.

## Typical workflow

```bash
# 1. Filter aligned BAMs to your regions of interest (R1 only, paired-end default)
sbatch --array=1-12 filter_bam_improved.sh bam_lookup.txt filtered_bams/ regions_of_interest.bed

# 2. Point a new lookup at the filtered BAMs
ls filtered_bams/*.filtered.bam > filtered_bam_lookup.txt

# 3. Compute per-base, per-strand 5' cut-site coverage
sbatch --array=1-12 genomecov.sh filtered_bam_lookup.txt cutsite_coverage/
```

The resulting `_r1_5p.bed` files give you exact cut-site positions and
counts per strand within your regions of interest — useful for footprinting,
comparing cut-site density across conditions/replicates at specific loci, or
as an input to other position-resolution downstream analyses, as a
complement to (not a replacement for) the HOMER-based peak calling used
elsewhere in the main pipeline.
