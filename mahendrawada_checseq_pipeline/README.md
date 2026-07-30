# ChEC-seq Analysis Pipeline

SLURM-based pipeline for processing paired-end ChEC-seq data in *S. cerevisiae*
from FASTQ to annotated peaks, following the methods from Mahendrawada et al.
2025, with D. melanogaster spike-in mapping and a consolidated MultiQC report.

**NOTE**: The  option `--authors_orig` sets the same settings used in Mahendrawada
et al 2025 according to both the Methods section, and scripts shared upon request.

## Pipeline steps

| Script | Purpose |
|---|---|
| `00_prepare.sh` | Validates the sample lookup file and FASTQ paths |
| `maketagdir_control.sh` | **Manual, one-time step.** Builds the HOMER tag directory for the free MNase control |
| `01_align.sh` | Aligns paired-end reads to *S. cerevisiae* (bowtie2); splits output into full/nuclear/mitochondrial BAMs; writes unmapped reads to FASTQ |
| `01a_map_to_dmel.sh` | Aligns unmapped reads from `01_align.sh` to *D. melanogaster* for spike-in normalization |
| `02_maketagdir_samples.sh` | Builds a HOMER tag directory per sample replicate |
| `03_findpeaks.sh` | Calls peaks per replicate with HOMER `findPeaks`, using the control tag directory as background |
| `04_pos2bed.sh` | Converts HOMER peak files to BED format and extracts peak summits |
| `05_annotatepeaks.sh` | Annotates peaks with HOMER `annotatePeaks.pl` against a custom GTF |
| `06_multiqc.sh` | Scans `results/` and `logs/` and builds a single consolidated MultiQC report |
| `submit_pipeline.sh` | Submits all of the above (except `maketagdir_control.sh`) with the correct SLURM dependencies |

## Pipeline dependency graph

```
maketagdir_control.sh   (manual, run once, before submit_pipeline.sh)

01_align
  ├─→ 01a_map_to_dmel
  └─→ 02_maketagdir_samples
       └─→ 03_findpeaks        (also reads the control tag directory)
            ├─→ 04_pos2bed
            └─→ 05_annotatepeaks

06_multiqc   (runs after 01a, 04, and 05 all finish; scans the full results/ tree)
```

## Prerequisites

- SLURM job scheduler with container support (`#SBATCH --container=...`)
- bowtie2, samtools, HOMER, MultiQC (each pulled via container image per script)
- Reference genomes indexed for bowtie2:
  - *S. cerevisiae* (sacCer3 / R64-5-1)
  - *D. melanogaster* (dm6, release 6.65)
- A custom GTF for the sacCer3 assembly (used by `05_annotatepeaks.sh`)

## Input format

A tab-separated lookup file with a header row and one line per sample replicate:

```
regulator_symbol	replicate	fastq_1	fastq_2
AFT2	rep1	/path/to/AFT2_rep1_R1.fastq.gz	/path/to/AFT2_rep1_R2.fastq.gz
AFT2	rep2	/path/to/AFT2_rep2_R1.fastq.gz	/path/to/AFT2_rep2_R2.fastq.gz
MED1	rep1	/path/to/MED1_rep1_R1.fastq.gz	/path/to/MED1_rep1_R2.fastq.gz
```

- Column 1: regulator/transcription factor symbol
- Column 2: replicate label (e.g. `rep1`, `A`, `B`, `C`)
- Column 3/4: paths to R1/R2 FASTQ (gzipped or uncompressed)
- The free MNase control is also listed as a row in this file, aligned by
  `01_align.sh` like any other sample; its resulting BAM and stats file are
  what you point `maketagdir_control.sh` at.

## Hardcoded reference paths

These paths are set directly inside the scripts (not passed as arguments).
Update them in-place if your reference locations change:

| Script | Variable | Current value |
|---|---|---|
| `01_align.sh` | `SCER_BOWTIE_INDEX` | `/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/bowtie2_index/S288C_reference_sequence_R64-5-1_20240529_chr_normalized` |
| `01a_map_to_dmel.sh` | `DMEL_BOWTIE_INDEX` | `/ref/mblab/data/dmelanogaster/bowtie2_index/dmel-all-chromosome-r6.65` |
| `03_findpeaks.sh` | `GENOME_SIZE` | `12071326` (nuclear) or `12157105` (full) - picked automatically from `bam_type`, see below |
| `maketagdir_control.sh` / `02_maketagdir_samples.sh` | `GENOME_FASTA` | `/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529_chr_normalized.fa` |
| `05_annotatepeaks.sh` | `GENOME_FASTA` / `GTF_FILE` | sacCer3 FASTA (above) / `/ref/mblab/data/yeast_data/reprocess_mahendrawada/mahendrawada_slurm_pipeline/sacCer3.ensGene.gtf` |

## Options

`submit_pipeline.sh` accepts three optional flags, in addition to the
required `<lookup_file>`:

```bash
bash submit_pipeline.sh <lookup_file> [--bam-type=nuclear|full] [--start-at=STEP] [--authors-orig]
```

**`--bam-type=nuclear|full`** (default: `nuclear`)
Which BAM `02_maketagdir_samples.sh` and `03_findpeaks.sh` use:
- `nuclear` — `{regulator}_{replicate}_nuclear.bam` (chrM filtered out), genome size `12071326`
- `full` — `{regulator}_{replicate}.bam` (all chromosomes, incl. chrM), genome size `12157105`

Both scripts must agree on this, which is why the flag threads through to
both automatically rather than being set independently.

**`--start-at=STEP`** (default: `01_align`)
Resume the pipeline partway through instead of resubmitting everything.
Anything upstream of `STEP` is assumed to have already completed
successfully; it is not resubmitted, and the step you start at is submitted
with no dependency. Accepts any step name, with or without `.sh`:
`01_align`, `01a_map_to_dmel`, `02_maketagdir_samples`, `03_findpeaks`,
`04_pos2bed`, `05_annotatepeaks`, `06_multiqc`.

**`--authors-orig`** (default: off)
Passes `--authors_orig` to `02_maketagdir_samples.sh`, which uses `-keepAll`
instead of `-unique -mapq 10` when building sample tag directories —
matching the original Mahendrawada et al. `makeTagDirectory` calls (which
kept all alignments, including multi-mappers/low-MAPQ reads). This only
affects the automated `02_maketagdir_samples.sh` step; if you also want the
manually-run `maketagdir_control.sh` built the same way, pass
`--authors_orig` to it directly.

## Output layout

```
results/
├── bams/{regulator}/{replicate}/
│   ├── *.bam, *_nuclear.bam, *_mito.bam        (+ .bai indexes)
│   ├── *_unmapped_R1.fastq.gz, *_unmapped_R2.fastq.gz
│   ├── *_dmel.bam                              (D. melanogaster spike-in)
│   ├── *_dmel_counts.txt
│   └── *_samtools_stats.txt, *_flagstats.txt, *_idxstats.txt, *_coverage.txt
│       (generated for the S. cerevisiae BAM and, separately, the dmel BAM)
├── tag_dirs/
│   ├── control_MNase/               (built by maketagdir_control.sh)
│   └── {regulator}_{replicate}/     (flat naming - NOT nested - so every
│                                     sample gets a unique name in MultiQC)
├── peaks/{regulator}/{replicate}/
│   ├── *_peaks.txt              (HOMER native format)
│   ├── *_peaks.bed, *_peaks_summits.txt
│   └── *_annotated_peaks.txt, *_annotatePeaks.err
└── multiqc/
    ├── multiqc_report.html
    └── multiqc_data/
```

## Building the control samples

This is a one-time, manual process per `bam_type` (`nuclear` and/or `full`),
done before `submit_pipeline.sh` is ever run. Since the free MNase control is
sequenced as two replicates, they need to be merged into a single combined
BAM before building the control tag directory.

**1. Align both control replicates**, same as any other sample:
```bash
sbatch --array=1-2 01_align.sh freemnase_lookup.txt
```
This produces, per replicate, the usual full/nuclear/mito BAM trio under
`results/bams/free_mnase/{rep1,rep2}/`.

**2. Merge the two replicates into a combined BAM** - once for whichever
`bam_type`(s) you plan to use downstream (nuclear-only, full, or both):
```bash
# Nuclear-only combined control
samtools merge combined_freemnase.bam \
    results/bams/free_mnase/rep1/free_mnase_rep1_nuclear.bam \
    results/bams/free_mnase/rep2/free_mnase_rep2_nuclear.bam
samtools index combined_freemnase.bam
samtools stats combined_freemnase.bam > combined_freemnase.stats

# Full (chrM-included) combined control
samtools merge combined_freemnase_full.bam \
    results/bams/free_mnase/rep1/free_mnase_rep1.bam \
    results/bams/free_mnase/rep2/free_mnase_rep2.bam
samtools index combined_freemnase_full.bam
samtools stats combined_freemnase_full.bam > combined_freemnase_full.stats
```

**3. Build the tag directory** with `maketagdir_control.sh`, pointing it at
whichever combined BAM/stats pair matches the `bam_type` you're targeting:
```bash
sbatch maketagdir_control.sh combined_freemnase.bam combined_freemnase.stats
```
This writes to `results/tag_dirs/control_MNase` (a path hardcoded relative
to wherever you run it from - not to the archive locations below).

**4. Archive each `bam_type`'s BAM and tag directory separately**, since
you'll want to switch between them later without rebuilding. A layout like:
```
control_data/
├── nuclear/
│   ├── bams/
│   │   ├── combined_freemnase.bam(.bai)
│   │   └── combined_freemnase.stats
│   └── tag_dir/
│       └── control_MNase/       (copied from results/tag_dirs/control_MNase
│                                  after step 3, for the nuclear build)
└── full/
    ├── bams/
    │   ├── combined_freemnase_full.bam(.bai)
    │   └── combined_freemnase_full.stats
    └── tag_dir/
        └── control_MNase/       (same, for the full build)
```

**5. Before running `submit_pipeline.sh` with a given `--bam-type`**, stage
the matching archived tag directory into place, since that's the path
`03_findpeaks.sh` actually reads from:
```bash
# example: about to run with --bam-type=nuclear (the default)
rm -rf results/tag_dirs/control_MNase
cp -r control_data/nuclear/tag_dir/control_MNase results/tag_dirs/control_MNase
```
Swap in `control_data/full/tag_dir/control_MNase` instead if running with
`--bam-type=full`. Mismatching the staged control against the `--bam-type`
you're running samples with will silently produce inconsistent results
(different genome size, different chrM composition in the background model)
rather than an obvious error - so this is worth double-checking each time
you switch.

## Running it

```bash
# One-time per bam_type: build + archive the control tag directory
# (see "Building the control samples" above)

# Before each run: stage the matching archived control tag directory
cp -r control_data/nuclear/tag_dir/control_MNase results/tag_dirs/control_MNase

# Validate your lookup file
bash 00_prepare.sh samples.tsv

# Submit the rest of the pipeline (can be run from any directory)
bash submit_pipeline.sh samples.tsv
```

Monitor progress:
```bash
squeue -u $USER
tail -f logs/align_1.log
```

When it's done:
```
results/multiqc/multiqc_report.html   # start here
results/peaks/{regulator}/{replicate}/*_annotated_peaks.txt
```

## Troubleshooting

### A step is stuck in `PD` with reason `DependencyNeverSatisfied`

```bash
squeue -u $USER
#   44106823_[1-534]   general chec_fin   chasem PD  0:00  1 (DependencyNeverSatisfied)
```

This means an array job it depends on had at least one failed task. SLURM's
`afterok` dependency on an array requires **every single task** in that
array to exit 0 - one failure permanently blocks everything downstream,
even if 99% of the array succeeded.

**1. Find out what it's actually waiting on:**
```bash
scontrol show job 44106823 | grep -i depend
#   Dependency=afterok:44106822_*(failed)
```
The `(failed)` tag confirms it, and gives you the upstream job ID.

**2. Find which specific task(s) failed:**
```bash
sacct -j 44106822 --format=JobID,JobName,State,ExitCode | grep -v COMPLETED
```
This shows the exact array index(es) and their state (`FAILED`, `TIMEOUT`,
`CANCELLED`, etc.) - note this state, since it tells you what kind of fix is
needed (see below).

**3. Check the log for each failed task to find the root cause:**
```bash
cat logs/tagdir_samples_57.log   # substitute the actual step + task number
```

**4. Address the root cause, then either:**
- **Fix and re-run** just the failed array indices, if the problem was
  something environmental (e.g. `TIMEOUT` - bump `--time`; a transient
  cluster issue; a missing dependency) rather than a real data problem:
  ```bash
  sbatch --time=01:00:00 --array=57,181,203 02_maketagdir_samples.sh samples.tsv
  ```
- **OR, if the problem is with the input data itself** (e.g. a genuinely bad
  or mispaired FASTQ pair that can't be fixed by re-running), make a new
  lookup file that omits those specific rows, and proceed with the rest of
  the pipeline on everything else. Revisit the excluded samples separately
  once the underlying data issue is resolved.

**Important:** re-running the same array indices under a *new* job ID does
**not** retroactively fix the original job. `44106822` stays permanently
`(failed)` no matter how many times you successfully rerun those indices
elsewhere - and anything still depending on `44106822` specifically will
never unstick.

**5. Cancel the stuck downstream jobs** - they can never resolve themselves:
```bash
scancel 44106823 44106825 44106826 44106827
```

**6. Resubmit from wherever you actually need to resume**, using
`--start-at` so you don't redo everything that already succeeded:
```bash
bash submit_pipeline.sh samples.tsv --start-at=03_findpeaks
```
This builds a fresh dependency chain rooted in new job IDs, with no
reference to the dead job from step 4.

Remember to pass along any other flags you were using (`--bam-type`,
`--authors-orig`) so the resumed steps stay consistent with the ones that
already ran.
