# ChEC-seq Analysis Pipeline

SLURM-based pipeline for processing paired-end ChEC-seq data in *S. cerevisiae*
from FASTQ to annotated peaks, following the methods from Mahendrawada et al.
2025, with D. melanogaster spike-in mapping and a consolidated MultiQC report.

## Pipeline steps

| Script | Purpose |
|---|---|
| `00_prepare.sh` | Validates the sample lookup file and FASTQ paths |
| `maketagdir_control.sh` | **Manual, one-time step.** Builds the HOMER tag directory for the free MNase control |
| `01_align.sh` | Aligns paired-end reads to *S. cerevisiae* (bowtie2); splits output into full/nuclear/mitochondrial BAMs; writes unmapped reads to FASTQ |
| `01a_map_to_dmel.sh` | Aligns unmapped reads from `01_align.sh` to *D. melanogaster* for spike-in normalization |
| `01b_dmel_normalized_coverage.sh` | Computes genome-wide coverage from the *S. cerevisiae* alignment, normalized to the dmel spike-in read count (per Mahendrawada et al. methods) |
| `01c_filter_bam.sh` | Filters the *S. cerevisiae* alignment to properly-paired, MAPQ≥10 reads overlapping a supplied regions-of-interest BED |
| `01d_genomecov_5p.sh` | Computes per-base, per-strand 5' read-end ("cut site") coverage from `01c_filter_bam.sh`'s output |
| `02_maketagdir_samples.sh` | Builds a HOMER tag directory per sample replicate |
| `03_findpeaks.sh` | Calls peaks per replicate with HOMER `findPeaks`, using the control tag directory as background |
| `04_pos2bed.sh` | Converts HOMER peak files to BED format and extracts peak summits |
| `05_annotatepeaks.sh` | Annotates peaks with HOMER `annotatePeaks.pl` against a custom GTF |
| `06_hahn_region_scoring.sh` | Runs `hahn_region_scoring.R` (Mahendrawada 2025/Donczew & Hahn 2020 method) - array job, one task per unique regulator |
| `07_promoter_scoring.sh` | Runs `promoter_scoring.R` (calling-cards-style method) - array job, one task per unique regulator |
| `08_multiqc.sh` | Scans `results/` and `logs/` and builds a single consolidated MultiQC report |
| `submit_pipeline.sh` | Submits all of the above (except `maketagdir_control.sh`) with the correct SLURM dependencies |

## Pipeline dependency graph

```
maketagdir_control.sh   (manual, run once, before submit_pipeline.sh)

01_align
  ├─→ 01a_map_to_dmel              (opt-in: only submitted with --align_dmel)
  │    └─→ 01b_dmel_coverage       (also depends on 01_align directly)
  │         └─→ 06_hahn_region_scoring     (opt-in: also requires --tss-bed=;
  │              ↑                          ALSO depends on 04_pos2bed below,
  │              │                          since it reads peak-assignment
  │              │                          output too - see next branch)
  ├─→ 01c_filter_bam                (opt-in: only submitted with --filter_genomecov)
  │    └─→ 01d_genomecov_5p
  │         └─→ 07_promoter_scoring  (opt-in: also requires --promoter-bed= and --control-bed=)
  └─→ 02_maketagdir_samples
       └─→ 03_findpeaks        (also reads the control tag directory)
            ├─→ 04_pos2bed ───────────────→ 06_hahn_region_scoring (see above)
            └─→ 05_annotatepeaks

08_multiqc   (runs after 04, 05, and whichever of 01b/06/01d/07 actually ran; scans the full results/ tree)
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
| `05_annotatepeaks.sh` | `GENOME_FASTA` / `GTF_FILE` | sacCer3 FASTA (above) / `sacCer3.ensGene.gtf` |

## Options

`submit_pipeline.sh` accepts several optional flags, in addition to the
required `<lookup_file>`:

```bash
bash submit_pipeline.sh <lookup_file> [--bam-type=nuclear|full] [--start-at=STEP] [--authors-orig] \
    [--align_dmel [--tss-bed=TSS.bed]] \
    [--filter_genomecov --include-regions=REGIONS.bed [--promoter-bed=PROMOTERS.bed --control-bed=CONTROL.bed]]
```

**`--bam-type=nuclear|full`** (default: `nuclear`)
Which BAM `01b_dmel_normalized_coverage.sh`, `01c_filter_bam.sh`,
`02_maketagdir_samples.sh`, and `03_findpeaks.sh` all use:
- `nuclear` — `{regulator}_{replicate}_nuclear.bam` (chrM filtered out), genome size `12071326`
- `full` — `{regulator}_{replicate}.bam` (all chromosomes, incl. chrM), genome size `12157105`

All four scripts must agree on this, which is why the flag threads through
to all of them automatically rather than being set independently.

**`--start-at=STEP`** (default: `01_align`)
Resume the pipeline partway through instead of resubmitting everything.
Anything upstream of `STEP` is assumed to have already completed
successfully; it is not resubmitted, and the step you start at is submitted
with no dependency. Accepts any step name, with or without `.sh`:
`01_align`, `01a_map_to_dmel`, `01b_dmel_coverage`, `01c_filter_bam`,
`01d_genomecov_5p`, `02_maketagdir_samples`, `03_findpeaks`, `04_pos2bed`,
`05_annotatepeaks`, `06_hahn_region_scoring`,
`07_promoter_scoring`, `08_multiqc`.
Starting at `01a_map_to_dmel` or `01b_dmel_coverage` also requires
`--align_dmel`; starting at `01c_filter_bam` or `01d_genomecov_5p` also
requires `--filter_genomecov`; starting at `06_hahn_region_scoring` also
requires `--tss-bed=`; starting at `07_promoter_scoring` also
requires both `--promoter-bed=` and `--control-bed=` - these steps are
otherwise disabled/skipped entirely rather than assumed-already-complete.

**`--authors-orig`** (default: off)
Passes `--authors_orig` to both `02_maketagdir_samples.sh` and
`03_findpeaks.sh`, matching the original Mahendrawada et al. scripts more
closely in three ways at once:
- `02_maketagdir_samples.sh`: `-keepAll` instead of `-unique -mapq 10` when
  building sample tag directories (kept all alignments, including
  multi-mappers/low-MAPQ reads)
- `02_maketagdir_samples.sh`: skips `-fragLength` entirely — no longer
  derives it from this sample's own `samtools stats`, letting HOMER's own
  autocorrelation estimate run instead (the original scripts never passed
  `-fragLength`)
- `03_findpeaks.sh`: skips `-gsize` entirely — no longer passes the
  hardcoded nuclear/full genome size constant, letting `findPeaks`
  auto-estimate genome size from the tag directory instead (the original
  scripts never passed `-gsize` either)

This only affects the automated `02`/`03` steps; if you also want the
manually-run `maketagdir_control.sh` built the same way (`-keepAll`, no
`-fragLength`), pass `--authors_orig` to it directly.

**`--align_dmel`** (default: off)
Enables the D. melanogaster spike-in branch: `01a_map_to_dmel.sh` and
`01b_dmel_normalized_coverage.sh`. **Off by default**, since not every
sample set has a dmel spike-in. When off, these two steps aren't submitted
at all — not treated as already-complete, genuinely skipped — and
`08_multiqc`'s dependency on them is dropped automatically so it doesn't
wait on jobs that were never submitted.

**`--tss-bed=TSS.bed`**
Enables `06_hahn_region_scoring.sh` (the Mahendrawada 2025/Donczew & Hahn
2020 promoter-scoring method) once `--align_dmel` has produced coverage for
every sample. **Optional even with `--align_dmel` set** — if omitted, `01a`/
`01b` still run but scoring is skipped. Runs as an array job, one task per
unique regulator in the lookup file (not one per row/replicate) — each
task's own R invocation discovers and combines that regulator's own
replicates internally. Additional tunable parameters (promoter window,
signal window, min replicates bound) are set inside
`06_hahn_region_scoring.sh`/`hahn_region_scoring.R` — see that
script's own header for the full list; run it directly (once its
dependencies exist) to override them without going through
`submit_pipeline.sh`.

**`--filter_genomecov`** (default: off, requires `--include-regions=`)
Enables an independent, dmel-free quantification branch: `01c_filter_bam.sh`
(region-restricted, MAPQ≥10, properly-paired BAM) → `01d_genomecov_5p.sh`
(per-base, per-strand 5' cut-site coverage). This only needs the
*S. cerevisiae* alignment — it works whether or not `--align_dmel` is used,
and is meant for sample sets that don't have a dmel spike-in at all.
Starting at `01c_filter_bam` or `01d_genomecov_5p` also requires this flag.

**`--include-regions=REGIONS.bed`**
Required when `--filter_genomecov` is set. BED file of regions
`01c_filter_bam.sh` restricts reads to (e.g. promoters).

**`--promoter-bed=PROMOTERS.bed` / `--control-bed=CONTROL.bed`**
Both enable `07_promoter_scoring.sh` (the calling-cards-style
promoter enrichment method) once `--filter_genomecov` has produced 5'
cut-site coverage for every sample. **Both optional even with
`--filter_genomecov` set** — if either is omitted, `01c`/`01d` still run but
scoring is skipped. Runs as an array job, one task per unique regulator in
the lookup file (not one per row/replicate) — each task's own R invocation
discovers and combines that regulator's own replicates internally.
`--control-bed` points at a combined control 5' cut-site
BED (see "Building the control samples" below). Additional tunable
parameters (pseudocount) are set inside
`07_promoter_scoring.sh`/`promoter_scoring.R`.

**`--control-tag-dir=<path>`** (default: `results/tag_dirs/control_MNase`)
Path to the HOMER control tag directory `03_findpeaks.sh` uses as its `-i`
background — the one built manually by `maketagdir_control.sh` (see
"Building the control samples" below). Change this if you'd rather point
directly at an archived control tag directory than stage/copy it into the
default `results/tag_dirs/control_MNase` location before each run.

## Two independent promoter-quantification pathways

This pipeline produces two different coverage tracks, each scored by its
own automated step (`06`/`07`) once its required inputs are provided:

| | dmel-normalized (`--align_dmel`) | filter+genomecov (`--filter_genomecov`) |
|---|---|---|
| Coverage script | `01b_dmel_normalized_coverage.sh` | `01c_filter_bam.sh` → `01d_genomecov_5p.sh` |
| Coverage type | Whole-fragment depth (`-pc`), dmel-normalized | Per-strand 5' cut-site depth, unnormalized |
| Requires dmel spike-in | Yes | No |
| Scoring step | `06_hahn_region_scoring.sh` (needs `--tss-bed=`) | `07_promoter_scoring.sh` (needs `--promoter-bed=` and `--control-bed=`) |
| Underlying R script | `hahn_region_scoring.R` | `promoter_scoring.R` |
| Method | Mahendrawada 2025 / Donczew & Hahn 2020: per-replicate peak-summit-anchored signal, bound-in-≥N-replicates filter, averaged across replicates | Calling-cards-style: `sum_overlap_scores` per promoter region, replicates combined via `Reduce`, Poisson/hypergeometric enrichment vs. a control coverage track |
| Background | Control genome-wide coverage (implicitly, via the paper's TSS/peak-summit window) | An explicit control 5' cut-site BED, built the same way as the control tag directory (see below) |
| Output | `results/hahn_region_scoring/{regulator}/{regulator}_summary.tsv` + `_replicate.tsv` | `results/promoter_scoring/{regulator}/{regulator}_combined.tsv` + `_replicate.tsv` |

Both R scripts are also runnable directly (not just via `submit_pipeline.sh`)
once their coverage inputs exist, with the rest of their parameters
(promoter/signal window sizes, min replicates bound, pseudocount, core
count) exposed as additional CLI flags — see each script's own header
comments and `option_list` for the full set.

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
│   └── {regulator}/{replicate}/     (nested, matching every other output
│                                     subdirectory - see 08_multiqc.sh's
│                                     --dirs --dirs-depth 2 for how sample
│                                     names stay unique in MultiQC despite
│                                     every regulator sharing leaf names
│                                     like "A"/"B"/"C")
├── peaks/{regulator}/{replicate}/
│   ├── *_peaks.txt              (HOMER native format)
│   ├── *_peaks.bed, *_peaks_summits.txt
│   └── *_annotated_peaks.txt, *_annotatePeaks.err
├── coverage/{regulator}/{replicate}/
│   └── *_dmel_norm.bedgraph     (genome-wide coverage / dmel_reads * 10000)
├── filtered_bams/{regulator}/{replicate}/
│   └── *_filtered.bam           (+.bai; region-restricted, MAPQ>=10, R1-only;
│                                  filename derived from the input BAM's own
│                                  basename, e.g. {reg}_{rep}_nuclear_filtered.bam
│                                  vs {reg}_{rep}_filtered.bam, per --bam-type)
├── genomecov_5p/{regulator}/{replicate}/
│   └── *_r1_5p.bed              (per-base, per-strand 5' cut-site coverage)
├── hahn_region_scoring/{regulator}/
│   ├── {regulator}_summary.tsv    (per-promoter)
│   └── {regulator}_replicate.tsv  (per-promoter, per-replicate)
├── promoter_scoring/{regulator}/
│   ├── {regulator}_combined.tsv   (per-promoter)
│   └── {regulator}_replicate.tsv  (per-promoter, per-replicate)
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

**Inspecting replicates individually before merging.** If you want to check
each free-MNase replicate's own tag directory (tag counts, GC bias,
autocorrelation) before deciding how to combine them, `maketagdir_control.sh`
also supports a lookup-driven array mode that builds one tag directory per
replicate instead of a single combined one:
```bash
sbatch --array=1-2 maketagdir_control.sh --lookup freemnase_lookup.txt nuclear
```
This writes to `results/tag_dirs/{regulator}/{replicate}/` (e.g.
`results/tag_dirs/free_mnase/A`) - nested, matching
`02_maketagdir_samples.sh`'s convention, deliberately distinct from the
singular `control_MNase` directory step 3 produces, so the two never
collide. This is independent of, and doesn't replace, steps 1-3 above -
it's for inspection/QC, not for use as `03_findpeaks.sh`'s actual control.

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

**5. Before running `submit_pipeline.sh` with a given `--bam-type`**, either
stage the matching archived tag directory into the default location, or
point `--control-tag-dir=` directly at it instead - your choice:
```bash
# Option A: stage into the default location (example: --bam-type=nuclear)
rm -rf results/tag_dirs/control_MNase
cp -r control_data/nuclear/tag_dir/control_MNase results/tag_dirs/control_MNase

# Option B: skip staging, point directly at the archived copy instead
bash submit_pipeline.sh samples.tsv --control-tag-dir=control_data/nuclear/tag_dir/control_MNase ...
```
Swap in `control_data/full/tag_dir/control_MNase` instead if running with
`--bam-type=full`. Mismatching the control against the `--bam-type`
you're running samples with will silently produce inconsistent results
(different genome size, different chrM composition in the background model)
rather than an obvious error - so this is worth double-checking each time
you switch, regardless of which option you use.

**If you're also using `--filter_genomecov`**, build the control's 5'
cut-site file the same way, using `01c_filter_bam.sh`/`01d_genomecov_5p.sh`
in **direct mode** on the same combined control BAM from step 2. Note the
output filename is derived from the *input* BAM's own basename (not
`--output-name`, which only controls the output *directory*) - so
`combined_freemnase.bam` produces `combined_freemnase_filtered.bam`:
```bash
sbatch 01c_filter_bam.sh --bam control_data/nuclear/bams/combined_freemnase.bam \
    --include-regions promoters.bed --output-name control_combined
sbatch 01d_genomecov_5p.sh \
    --bam results/filtered_bams/control_combined/combined_freemnase_filtered.bam \
    --output-name control_combined
```
This produces `results/genomecov_5p/control_combined/control_combined_r1_5p.bed`,
which you then pass as `--control-bed=` (to `submit_pipeline.sh`, or directly
to `07_promoter_scoring.sh`/`promoter_scoring.R`).

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

# With both scoring pathways enabled:
bash submit_pipeline.sh samples.tsv \
    --align_dmel --tss-bed=tss.bed \
    --filter_genomecov --include-regions=promoters.bed \
    --promoter-bed=promoters.bed --control-bed=results/genomecov_5p/control_combined/control_combined_r1_5p.bed
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
results/hahn_region_scoring/{regulator}/{regulator}_summary.tsv (if --tss-bed was used)
results/promoter_scoring/{regulator}/{regulator}_combined.tsv   (if --promoter-bed/--control-bed were used)
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
