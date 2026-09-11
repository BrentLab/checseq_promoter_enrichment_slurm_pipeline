bash checseq_promoter_enrichment_slurm_pipeline/submit_pipeline.sh \
    $1 \
    --primary-bowtie-index=/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/bowtie2_index/S288C_reference_sequence_R64-5-1_20240529_chr_normalized \
    --primary-genome-fasta=/ref/mblab/data/S288C_R64/S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529_chr_normalized.fa \
    --spikein-bowtie-index=/ref/mblab/data/dmelanogaster/bowtie2_index/dmel-all-chromosome-r6.65 \
    --primary-gtf-file=/ref/mblab/data/yeast_data/reprocess_mahendrawada/mahendrawada_slurm_pipeline/sacCer3.ensGene.gtf \
    --control-tag-dir=checseq_promoter_enrichment_slurm_pipeline/control_data/from_correspondance/tag_dir \
    --authors-orig \
    --align_dmel \
    --bam-type=full \
    --tss-bed=checseq_promoter_enrichment_slurm_pipeline/ypd_tss_sgd-5-1_verified_orf.bed \
    --control-coverage=checseq_promoter_enrichment_slurm_pipeline/control_data/from_correspondance/from_correspondance_A_dmel_norm.bedgraph \
    --filter_genomecov \
    --include-regions=checseq_promoter_enrichment_slurm_pipeline/nuclear_chroms.bed \
    --promoter-bed=checseq_promoter_enrichment_slurm_pipeline/start_codon_500bp_upstream_promoters.bed \
    --control-bed=checseq_promoter_enrichment_slurm_pipeline/control_data/from_correspondance/from_correspondance_A_r1_5p.bed \
    --start-at=02_maketagdir_samples
