# Mahendrawada ChEC-seq Reprocessing

Reprocessing pipeline for the ChEC-seq dataset from Mahendrawada et al.
2024, plus a complementary region-based cut-site quantification method.

## Contents

**[`mahendrawada_checseq_pipeline/`](mahendrawada_checseq_pipeline/README.md)**
— The main SLURM pipeline: paired-end alignment, D. melanogaster spike-in
mapping, HOMER tag directories and peak calling, and a consolidated MultiQC
report. Takes raw FASTQs through to annotated peaks.

**[`promoter_enrichment/`](promoter_enrichment/README.md)**
— A downstream, region-restricted alternative to peak calling: filters BAMs
to a set of regions of interest (e.g. promoters) and computes single-base,
per-strand 5' cut-site density as BED files. Complements the peak-calling
pipeline with nucleotide-resolution footprinting at specific loci.

**`sacCer3.ensGene.gtf`** — Shared sacCer3 GTF annotation, used by
`05_annotatepeaks.sh` in the main pipeline.

See each subdirectory's own README for setup, usage, and output details.
