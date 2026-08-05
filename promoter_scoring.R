#!/usr/bin/env Rscript
# ==============================================================================
# Promoter enrichment from 5' cut-site coverage
# ==============================================================================
#
# For each promoter region, sums the read depth (not just presence/absence -
# a position covered by 3 reads contributes 3, not 1) of 5' cut-site reads
# that overlap it, from the pipeline's 5' cut-site coverage output
# (01c_filter_bam.sh -> 01d_genomecov_5p.sh). Replicates are combined by
# summing their per-promoter counts together. Enrichment/Poisson/
# hypergeometric statistics are then computed against a supplied control/
# background track. This is a DIFFERENT, independent quantification pathway
# from hahn_region_scoring.R (which mimics Mahendrawada 2025/Donczew & Hahn
# 2020's peak-summit + dmel-normalized-coverage method). This one:
#   - does not require the dmel spike-in branch (01a/01b) at all - it only
#     needs the S. cerevisiae alignment, filtered to regions of interest
#     and reduced to 5' cut-site positions
#   - uses a supplied control/background coverage track (analogous to how
#     the control tag directory is passed to findPeaks) instead of a dmel
#     spike-in count, to build a Poisson/hypergeometric background model
#
# CONTROL INPUT, two modes (--combine-control-replicates, default OFF):
#   Off (default, unchanged from before): --control-bed is a literal,
#     already-combined *_r1_5p.bed file, read directly as-is.
#   On: --control-bed is instead a regulator-style NAME (not a file path) -
#     its own replicate subdirectories under --genomecov-dir are discovered
#     and summed together, the exact same way tagged regulators already are
#     (this mode literally reuses combine_replicates() for the control too).
#     Useful when the control itself has multiple replicates that haven't
#     been pre-combined into a single file (e.g. promoter_scoring_lookup.sh's
#     paired tagged/control-per-row design).
#
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(optparse)
})

# ==============================================================================
# CLI ARGUMENTS
# ==============================================================================

option_list <- list(
    make_option("--promoter-bed", type = "character", default = NULL,
                help = "Path to promoter regions BED (chr, start, end, locus_tag, score, strand) [required]"),
    make_option("--control-bed", type = "character", default = NULL,
                help = "Path to the combined control *_r1_5p.bed (see README: Building the control samples) [required]. If --combine-control-replicates is set, this is instead a regulator-style NAME to resolve under --genomecov-dir, not a literal file path."),
    make_option("--combine-control-replicates", action = "store_true", default = FALSE,
                help = "Treat --control-bed as a regulator-style name with its own replicate subdirectories under --genomecov-dir, combining them the same way tagged regulators are already combined - rather than a literal, already-combined file. Default: off (unchanged, existing behavior: --control-bed is a literal file)."),
    make_option("--genomecov-dir", type = "character", default = "results/genomecov_5p",
                help = "Root of the pipeline's results/genomecov_5p/{regulator}/{replicate}/ tree [default: %default]"),
    make_option("--output-dir", type = "character", default = "results/promoter_scoring",
                help = "Output directory for the result TSVs [default: %default]"),
    make_option("--pseudocount", type = "double", default = 0.1,
                help = "Pseudocount for enrichment/Poisson/hypergeometric calculations [default: %default]"),
    make_option("--require-strand-match", action = "store_true", default = FALSE,
                help = "Only count a read toward a promoter if they're on the same strand. Default: off - these libraries are unstranded, so a read overlapping a promoter counts regardless of strand (overrides GenomicRanges::findOverlaps' own strand-specific default)."),
    make_option("--regulators", type = "character", default = NULL,
                help = "Comma-separated regulator list [default: auto-discover from --genomecov-dir]"),
    make_option("--cores", type = "integer", default = NA,
                help = "Parallel cores [default: SLURM_CPUS_PER_TASK env var, or 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt[["promoter-bed"]])) {
    stop("--promoter-bed is required")
}
if (is.null(opt[["control-bed"]])) {
    stop("--control-bed is required")
}
if (!file.exists(opt[["promoter-bed"]])) {
    stop(sprintf("--promoter-bed file not found: %s", opt[["promoter-bed"]]))
}
COMBINE_CONTROL_REPLICATES <- opt[["combine-control-replicates"]]

if (!COMBINE_CONTROL_REPLICATES && !file.exists(opt[["control-bed"]])) {
    stop(sprintf("--control-bed file not found: %s", opt[["control-bed"]]))
}
if (COMBINE_CONTROL_REPLICATES && !dir.exists(file.path(opt[["genomecov-dir"]], opt[["control-bed"]]))) {
    stop(sprintf("--combine-control-replicates is set, but no directory found for control-bed name '%s' under --genomecov-dir (%s)",
                 opt[["control-bed"]], opt[["genomecov-dir"]]))
}
if (!dir.exists(opt[["genomecov-dir"]])) {
    stop(sprintf("--genomecov-dir directory not found: %s", opt[["genomecov-dir"]]))
}
if (opt[["pseudocount"]] <= 0) {
    stop(sprintf("--pseudocount must be > 0 (got %s)", opt[["pseudocount"]]))
}
if (!is.na(opt[["cores"]]) && opt[["cores"]] < 1) {
    stop(sprintf("--cores must be >= 1 (got %s)", opt[["cores"]]))
}

PROMOTER_BED_PATH <- opt[["promoter-bed"]]
CONTROL_BED_PATH <- opt[["control-bed"]]
GENOMECOV_DIR <- opt[["genomecov-dir"]]
OUTPUT_DIR <- opt[["output-dir"]]
PSEUDOCOUNT <- opt[["pseudocount"]]
IGNORE_STRAND <- !opt[["require-strand-match"]]
N_CORES <- if (!is.na(opt[["cores"]])) opt[["cores"]] else {
    env_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    if (nzchar(env_cores)) as.integer(env_cores) else 4L
}


# ==============================================================================
# HELPERS
# ==============================================================================

#' Convert BED format data frame to GRanges
#'
#' Handles coordinate system conversion from 0-indexed half-open BED format
#' to 1-indexed closed GenomicRanges format
#'
#' @param bed_df Data frame with chr, start, end columns in BED format (0-indexed, half-open)
#' @param zero_indexed Logical, whether input is 0-indexed (default: TRUE)
#' @return GRanges object
bed_to_granges <- function(bed_df, zero_indexed = TRUE) {

    if (!all(c("chr", "start", "end") %in% names(bed_df))) {
        stop("bed_df must have columns: chr, start, end")
    }

    if (zero_indexed) {
        gr_start <- bed_df$start + 1
        gr_end <- bed_df$end
    } else {
        gr_start <- bed_df$start
        gr_end <- bed_df$end
    }

    strand_vals <- if ("strand" %in% names(bed_df)) bed_df$strand else "*"

    gr <- GenomicRanges::GRanges(
        seqnames = bed_df$chr,
        ranges = IRanges::IRanges(start = gr_start, end = gr_end),
        strand = strand_vals
    )

    extra_cols <- setdiff(names(bed_df), c("chr", "start", "end", "strand"))
    if (length(extra_cols) > 0) {
        GenomicRanges::mcols(gr) <- bed_df[, extra_cols, drop = FALSE]
    }

    gr
}

#' Read a 01d_genomecov_5p.sh output BED (chr, start, end, name, score, strand)
#'
#' @param path Path to the *_r1_5p.bed file
#' @return GRanges with mcols$score = read count at that position
read_genomecov_5p_bed <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("genomecov 5' BED not found: %s", path))
    }

    tryCatch({
        tbl <- read_tsv(
            path,
            col_names = c("chr", "start", "end", "name", "score", "strand"),
            col_types = cols(chr = col_character(), start = col_integer(),
                             end = col_integer(), name = col_character(),
                             score = col_double(), strand = col_character())
        )
        bed_to_granges(select(tbl, chr, start, end, score, strand))
    }, error = function(e) {
        stop(sprintf("Failed to read/parse genomecov 5' BED at %s: %s", path, conditionMessage(e)))
    })
}

#' Sum scores of overlapping insertions per region (adapted from reference
#' script - added ignore_strand, defaulting to TRUE)
#'
#' GenomicRanges::findOverlaps() is strand-specific by default
#' (ignore.strand=FALSE) - since both the promoter regions (real gene
#' strand, from the promoter BED's strand column) and the 5' cut-site reads
#' (real read strand, from the +/- split in 01d_genomecov_5p.sh) carry
#' actual +/- values rather than "*", a read overlapping a promoter on the
#' OPPOSITE strand would silently not be counted unless this is explicitly
#' overridden. These libraries are unstranded, so a read over a promoter
#' should count regardless of which strand it's on - hence ignore_strand
#' defaults to TRUE here, deliberately overriding findOverlaps' own default.
#'
#' @param insertions_gr GRanges object with insertions containing a 'score' metadata column
#' @param regions_gr GRanges object with regions
#' @param ignore_strand If TRUE (default), count overlaps regardless of
#'   strand - correct for unstranded libraries. Set FALSE to require the
#'   read and promoter to be on the same strand (findOverlaps' own default).
#' @return Numeric vector of summed scores per region
sum_overlap_scores <- function(insertions_gr, regions_gr, ignore_strand = TRUE) {
    overlaps <- GenomicRanges::findOverlaps(regions_gr, insertions_gr, ignore.strand = ignore_strand)

    if (length(overlaps) == 0) {
        return(rep(0, length(regions_gr)))
    }

    scores <- GenomicRanges::mcols(insertions_gr)$score[S4Vectors::subjectHits(overlaps)]
    summed_scores <- tapply(scores, S4Vectors::queryHits(overlaps), sum)

    result <- rep(0, length(regions_gr))
    result[as.integer(names(summed_scores))] <- summed_scores
    result
}

#' Discover replicate names for a regulator from the genomecov_5p directory
#'
#' @param regulator Regulator/TF symbol
#' @param genomecov_dir Root genomecov_5p directory
#' @return Character vector of replicate names (subdirectory names)
discover_replicates <- function(regulator, genomecov_dir) {
    reg_dir <- file.path(genomecov_dir, regulator)
    if (!dir.exists(reg_dir)) {
        stop(sprintf("No genomecov_5p directory found for regulator '%s': %s", regulator, reg_dir))
    }
    list.dirs(reg_dir, full.names = FALSE, recursive = FALSE)
}

# ==============================================================================
# STATISTICS FUNCTIONS - generic count-based enrichment/significance tests
# ==============================================================================

#' Calculate enrichment (experiment vs. background ratio)
#'
#' @param total_background_counts Total number of counts in background (scalar or vector)
#' @param total_experiment_counts Total number of counts in experiment (scalar or vector)
#' @param background_counts Number of counts in background per region (vector)
#' @param experiment_counts Number of counts in experiment per region (vector)
#' @param pseudocount Pseudocount to avoid division by zero (default: 0.1)
#' @return Enrichment values
calculate_enrichment <- function(total_background_counts,
                                 total_experiment_counts,
                                 background_counts,
                                 experiment_counts,
                                 pseudocount = 0.1) {

    if (!all(is.numeric(c(total_background_counts, total_experiment_counts,
                          background_counts, experiment_counts)))) {
        stop("All inputs must be numeric")
    }

    n_regions <- length(background_counts)

    if (length(experiment_counts) != n_regions) {
        stop("background_counts and experiment_counts must be the same length")
    }

    if (length(total_background_counts) == 1) {
        total_background_counts <- rep(total_background_counts, n_regions)
    }
    if (length(total_experiment_counts) == 1) {
        total_experiment_counts <- rep(total_experiment_counts, n_regions)
    }

    if (length(total_background_counts) != n_regions ||
        length(total_experiment_counts) != n_regions) {
        stop("All input vectors must be the same length or scalars")
    }

    numerator <- experiment_counts / total_experiment_counts
    denominator <- (background_counts + pseudocount) / total_background_counts
    enrichment <- numerator / denominator

    if (any(enrichment < 0, na.rm = TRUE)) {
        stop("Enrichment values must be non-negative")
    }
    if (any(is.na(enrichment))) {
        stop("Enrichment values must not be NA")
    }
    if (any(is.infinite(enrichment))) {
        stop("Enrichment values must not be infinite")
    }

    enrichment
}

#' Calculate Poisson p-values
#'
#' @param total_background_counts Total number of counts in background (scalar or vector)
#' @param total_experiment_counts Total number of counts in experiment (scalar or vector)
#' @param background_counts Number of counts in background per region (vector)
#' @param experiment_counts Number of counts in experiment per region (vector)
#' @param pseudocount Pseudocount for lambda calculation (default: 0.1)
#' @param ... additional arguments to `ppois`. lower tail is set to FALSE already
#' @return Poisson p-values
calculate_poisson_pval <- function(total_background_counts,
                                   total_experiment_counts,
                                   background_counts,
                                   experiment_counts,
                                   pseudocount = 0.1,
                                   ...) {

    n_regions <- length(background_counts)

    if (length(total_background_counts) == 1) {
        total_background_counts <- rep(total_background_counts, n_regions)
    }
    if (length(total_experiment_counts) == 1) {
        total_experiment_counts <- rep(total_experiment_counts, n_regions)
    }

    hop_ratio <- total_experiment_counts / total_background_counts
    mu <- (background_counts + pseudocount) * hop_ratio
    x <- experiment_counts

    ppois(x - 1, lambda = mu, lower.tail = FALSE, ...)
}

#' Calculate hypergeometric p-values
#'
#' @param total_background_counts Total number of counts in background (scalar or vector)
#' @param total_experiment_counts Total number of counts in experiment (scalar or vector)
#' @param background_counts Number of counts in background per region (vector)
#' @param experiment_counts Number of counts in experiment per region (vector)
#' @param ... additional arguments to phyper. lower tail is set to false already
#' @return Hypergeometric p-values
calculate_hypergeom_pval <- function(total_background_counts,
                                     total_experiment_counts,
                                     background_counts,
                                     experiment_counts,
                                     ...) {

    n_regions <- length(background_counts)

    if (length(total_background_counts) == 1) {
        total_background_counts <- rep(total_background_counts, n_regions)
    }
    if (length(total_experiment_counts) == 1) {
        total_experiment_counts <- rep(total_experiment_counts, n_regions)
    }

    M <- total_background_counts + total_experiment_counts
    n <- total_experiment_counts
    N <- background_counts + experiment_counts
    x <- experiment_counts - 1

    valid <- (M >= 1) & (N >= 1)
    pval <- rep(1, length(M))

    if (any(valid)) {
        pval[valid] <- phyper(x[valid], n[valid], M[valid] - n[valid], N[valid],
                              lower.tail = FALSE, ...)
    }

    pval
}

# ==============================================================================
# COMBINE REPLICATES
# ==============================================================================

#' Combine a regulator's replicates: per-replicate region counts, plus a
#' summed ("combined") version.
#'
#' @param regulator Regulator/TF symbol
#' @param regions_gr Promoter regions GRanges
#' @param genomecov_dir Root genomecov_5p directory
#' @return list(library_totals, replicates (named list of region-count
#'   vectors), combined (summed region-count vector))
combine_replicates <- function(regulator, regions_gr, genomecov_dir = GENOMECOV_DIR,
                                ignore_strand = IGNORE_STRAND) {

    message(sprintf("Working on regulator: %s", regulator))

    replicates <- discover_replicates(regulator, genomecov_dir)

    insertions_by_rep <- purrr::map(replicates, function(rep_name) {
        path <- file.path(genomecov_dir, regulator, rep_name,
                           sprintf("%s_%s_r1_5p.bed", regulator, rep_name))
        read_genomecov_5p_bed(path)
    })
    names(insertions_by_rep) <- replicates

    # Library total = total 5' cut-site reads (sum of the score column - each
    # read contributes exactly 1 to the depth at its own 5' position, so this
    # is exactly the total read count, not just "distinct positions covered")
    library_totals <- tibble(
        replicate = replicates,
        n = purrr::map_dbl(insertions_by_rep, ~sum(.x$score))
    )

    replicate_region_counts <- purrr::map(insertions_by_rep, ~sum_overlap_scores(.x, regions_gr, ignore_strand = ignore_strand))
    names(replicate_region_counts) <- replicates

    list(
        library_totals = library_totals,
        replicates = replicate_region_counts,
        combined = Reduce(`+`, replicate_region_counts)
    )
}

#' Load the control/background 5' cut-site coverage and compute per-region counts
#'
#' Two modes, selected by combine_replicates_flag:
#'   FALSE (default, unchanged from before): control_bed_path is a literal,
#'     already-combined *_r1_5p.bed file - read directly, no replicate
#'     discovery/combining performed here.
#'   TRUE: control_bed_path is instead a regulator-style NAME (not a file
#'     path) - its replicates are discovered under genomecov_dir the exact
#'     same way any tagged regulator's are, and summed via the same
#'     Reduce(`+`, ...) logic already used for tagged regulators (this
#'     literally reuses combine_replicates() - the control is treated as
#'     just another regulator to combine, once this mode is on).
#'
#' @param regions_gr Promoter regions GRanges
#' @param control_bed_path Literal file path (default mode) or regulator-
#'   style name (combine_replicates_flag = TRUE mode)
#' @param combine_replicates_flag If TRUE, control_bed_path is a name to
#'   combine replicates for, not a literal file path
#' @param genomecov_dir Root genomecov_5p directory (only used when
#'   combine_replicates_flag = TRUE)
#' @return list(total, counts) - total background reads, and per-region counts
load_control <- function(regions_gr, control_bed_path = CONTROL_BED_PATH,
                          combine_replicates_flag = COMBINE_CONTROL_REPLICATES,
                          genomecov_dir = GENOMECOV_DIR,
                          ignore_strand = IGNORE_STRAND) {
    if (combine_replicates_flag) {
        combined <- combine_replicates(control_bed_path, regions_gr, genomecov_dir, ignore_strand = ignore_strand)
        list(
            total = sum(combined$library_totals$n),
            counts = combined$combined
        )
    } else {
        control_gr <- read_genomecov_5p_bed(control_bed_path)
        list(
            total = sum(control_gr$score),
            counts = sum_overlap_scores(control_gr, regions_gr, ignore_strand = ignore_strand)
        )
    }
}

# ==============================================================================
# MAIN DRIVER: enrichment analysis for one regulator
# ==============================================================================

#' Compute promoter enrichment statistics for one regulator: per-replicate
#' AND combined (Reduce-summed) versions.
#'
#' @param regulator Regulator/TF symbol
#' @param regions_gr Promoter regions GRanges (shared across all regulators)
#' @param control_counts Per-region background counts (from load_control())
#' @param control_total Total background reads (from load_control())
#' @param genomecov_dir Root genomecov_5p directory
#' @param pseudocount Pseudocount for enrichment/Poisson calculations
#' @return list(replicates = named list of GRanges with stats, combined = GRanges with stats)
enrichment_analysis <- function(regulator,
                                  regions_gr,
                                  control_counts,
                                  control_total,
                                  genomecov_dir = GENOMECOV_DIR,
                                  pseudocount = PSEUDOCOUNT) {

    counts_regulator <- combine_replicates(regulator, regions_gr, genomecov_dir)

    replicate_quants <- purrr::map(names(counts_regulator$replicates), function(rep_name) {
        message(sprintf("  Working on replicate: %s", rep_name))

        af <- regions_gr
        experiment_counts <- counts_regulator$replicates[[rep_name]]
        total_experiment_counts <- counts_regulator$library_totals %>%
            filter(replicate == rep_name) %>%
            pull(n)

        # Raw counts the enrichment/statistics below are computed from -
        # included directly so the numbers can be recomputed/audited without
        # re-reading the source bed files. tagged_total_tags and
        # background_total_tags are scalars, recycled across every promoter
        # row for this replicate (background_total_tags is additionally the
        # same across every regulator/replicate, since there's only one
        # control).
        GenomicRanges::mcols(af)$tagged_tag_count <- experiment_counts
        GenomicRanges::mcols(af)$background_tag_count <- control_counts
        GenomicRanges::mcols(af)$tagged_total_tags <- total_experiment_counts
        GenomicRanges::mcols(af)$background_total_tags <- control_total

        GenomicRanges::mcols(af)$enrichment <- calculate_enrichment(
            control_total, total_experiment_counts, control_counts, experiment_counts, pseudocount
        )
        GenomicRanges::mcols(af)$poisson_pval <- calculate_poisson_pval(
            control_total, total_experiment_counts, control_counts, experiment_counts, pseudocount
        )
        GenomicRanges::mcols(af)$log_poisson_pval <- calculate_poisson_pval(
            control_total, total_experiment_counts, control_counts, experiment_counts, pseudocount, log.p = TRUE
        )
        GenomicRanges::mcols(af)$hypergeometric_pval <- calculate_hypergeom_pval(
            control_total, total_experiment_counts, control_counts, experiment_counts
        )
        GenomicRanges::mcols(af)$log_hypergeometric_pval <- calculate_hypergeom_pval(
            control_total, total_experiment_counts, control_counts, experiment_counts, log.p = TRUE
        )
        GenomicRanges::mcols(af)$poisson_qval <- p.adjust(GenomicRanges::mcols(af)$poisson_pval, method = "fdr")
        GenomicRanges::mcols(af)$hypergeometric_qval <- p.adjust(GenomicRanges::mcols(af)$hypergeometric_pval, method = "fdr")

        af
    })
    names(replicate_quants) <- names(counts_regulator$replicates)

    message(sprintf("Working on the combined for regulator %s", regulator))

    combined_gr <- regions_gr
    combined_experiment_counts <- counts_regulator$combined
    combined_total_experiment_counts <- sum(counts_regulator$library_totals$n)

    GenomicRanges::mcols(combined_gr)$tagged_tag_count <- combined_experiment_counts
    GenomicRanges::mcols(combined_gr)$background_tag_count <- control_counts
    GenomicRanges::mcols(combined_gr)$tagged_total_tags <- combined_total_experiment_counts
    GenomicRanges::mcols(combined_gr)$background_total_tags <- control_total

    GenomicRanges::mcols(combined_gr)$enrichment <- calculate_enrichment(
        control_total, combined_total_experiment_counts, control_counts, combined_experiment_counts, pseudocount
    )
    GenomicRanges::mcols(combined_gr)$poisson_pval <- calculate_poisson_pval(
        control_total, combined_total_experiment_counts, control_counts, combined_experiment_counts, pseudocount
    )
    GenomicRanges::mcols(combined_gr)$log_poisson_pval <- calculate_poisson_pval(
        control_total, combined_total_experiment_counts, control_counts, combined_experiment_counts, pseudocount, log.p = TRUE
    )
    GenomicRanges::mcols(combined_gr)$hypergeometric_pval <- calculate_hypergeom_pval(
        control_total, combined_total_experiment_counts, control_counts, combined_experiment_counts
    )
    GenomicRanges::mcols(combined_gr)$log_hypergeometric_pval <- calculate_hypergeom_pval(
        control_total, combined_total_experiment_counts, control_counts, combined_experiment_counts, log.p = TRUE
    )
    GenomicRanges::mcols(combined_gr)$poisson_qval <- p.adjust(GenomicRanges::mcols(combined_gr)$poisson_pval, method = "fdr")
    GenomicRanges::mcols(combined_gr)$hypergeometric_qval <- p.adjust(GenomicRanges::mcols(combined_gr)$hypergeometric_pval, method = "fdr")

    list(replicates = replicate_quants, combined = combined_gr)
}

# ==============================================================================
# MAIN DRIVER
# ==============================================================================

message(sprintf("Loading promoter regions from %s...", PROMOTER_BED_PATH))
promoter_tbl <- tryCatch({
    read_tsv(
        PROMOTER_BED_PATH,
        col_names = c("chr", "start", "end", "locus_tag", "score", "strand"),
        col_types = cols(chr = col_character(), start = col_integer(),
                         end = col_integer(), locus_tag = col_character(),
                         score = col_double(), strand = col_character())
    )
}, error = function(e) {
    stop(sprintf("Failed to read/parse --promoter-bed at %s: %s", PROMOTER_BED_PATH, conditionMessage(e)))
})
regions_gr <- bed_to_granges(select(promoter_tbl, chr, start, end, locus_tag, strand))

message(sprintf("Loading control coverage from %s...", CONTROL_BED_PATH))
control <- tryCatch({
    load_control(regions_gr, CONTROL_BED_PATH)
}, error = function(e) {
    stop(sprintf("Failed to load --control-bed at %s: %s", CONTROL_BED_PATH, conditionMessage(e)))
})
message(sprintf("Control total reads: %s", format(control$total, big.mark = ",")))

regulator_list <- if (!is.null(opt$regulators)) {
    strsplit(opt$regulators, ",")[[1]]
} else {
    reg_dirs <- list.dirs(GENOMECOV_DIR, full.names = FALSE, recursive = FALSE)
    # Exclude the control's own directory from being auto-discovered as a
    # regulator. In combine-replicates mode CONTROL_BED_PATH is already a
    # bare name (not a path) - exclude it directly. In the default (literal
    # file) mode, derive the name from the file's own parent directory
    # (e.g. a "control_combined" direct-mode output living under genomecov_dir).
    control_dir_to_exclude <- if (COMBINE_CONTROL_REPLICATES) {
        CONTROL_BED_PATH
    } else {
        basename(dirname(CONTROL_BED_PATH))
    }
    setdiff(reg_dirs, control_dir_to_exclude)
}

message(sprintf("Scoring %d regulator(s) using %d core(s)...", length(regulator_list), N_CORES))

all_results <- parallel::mclapply(
    regulator_list,
    enrichment_analysis,
    regions_gr = regions_gr,
    control_counts = control$counts,
    control_total = control$total,
    genomecov_dir = GENOMECOV_DIR,
    pseudocount = PSEUDOCOUNT,
    mc.cores = N_CORES
)
names(all_results) <- regulator_list

errored <- purrr::keep(all_results, ~inherits(.x, "try-error"))
if (length(errored) > 0) {
    for (reg in names(errored)) {
        message(sprintf("ERROR processing regulator '%s': %s", reg, as.character(errored[[reg]])))
    }
    all_results <- all_results[!names(all_results) %in% names(errored)]
}

if (length(all_results) == 0) {
    stop("No regulators were successfully processed - see errors above.")
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Write per-regulator, matching the {output-dir}/{regulator}/{regulator}_...
# convention used everywhere else in the pipeline (results/bams/{regulator}/,
# results/peaks/{regulator}/, etc.), rather than one flat file combining
# every regulator together.
for (regulator in names(all_results)) {
    reg_dir <- file.path(OUTPUT_DIR, regulator)
    dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)

    res <- all_results[[regulator]]

    per_replicate_tbl <- purrr::imap_dfr(res$replicates, function(gr, rep_name) {
        as_tibble(gr) %>% mutate(regulator = regulator, replicate = rep_name)
    }) %>% relocate(regulator, replicate)

    combined_tbl <- as_tibble(res$combined) %>%
        mutate(regulator = regulator) %>%
        relocate(regulator)

    summary_path <- file.path(reg_dir, sprintf("%s_combined.tsv", regulator))
    per_replicate_path <- file.path(reg_dir, sprintf("%s_replicate.tsv", regulator))

    readr::write_tsv(combined_tbl, summary_path)
    readr::write_tsv(per_replicate_tbl, per_replicate_path)

    message(sprintf("Wrote %s: combined (%d rows), per-replicate (%d rows)",
                     regulator, nrow(combined_tbl), nrow(per_replicate_tbl)))
}

message("Done.")
