bash checseq_promoter_enrichment_slurm_pipeline/submit_pipeline.sh \
    $1 \
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
