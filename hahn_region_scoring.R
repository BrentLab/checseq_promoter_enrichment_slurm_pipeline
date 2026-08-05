#!/usr/bin/env Rscript
# ==============================================================================
# Promoter enrichment (signal per promoter) from dmel-normalized coverage
# ==============================================================================
#
# Reimplements the "signal per promoter" quantification described in
# Mahendrawada et al. 2025 / Donczew & Hahn 2020 (elifesciences.org/articles/50109),
# following a LITERAL, per-replicate reading of both methods texts:
#
#   1. PER REPLICATE, independently: a promoter is "bound" if a called peak's
#      summit falls within a window around its TSS (paper: -400/+200bp). If
#      more than one peak's summit falls in that window, the one closest to
#      the TSS is used. This is done separately for each replicate - no
#      pooling of peaks across replicates at this stage.
#
#   2. A promoter is included in the final bound-promoter list if it was
#      independently bound (step 1) in at least N replicates (Mahendrawada
#      2025: 2 of 3; Donczew 2020: 3 of 4).
#
#   3. For promoters on the final bound list, signal is computed PER
#      REPLICATE, using that replicate's own normalized coverage track:
#        - if this replicate had an assigned peak (step 1), its summit is
#          used as the anchor
#        - if this replicate specifically lacked a peak for this promoter
#          (even though the promoter cleared the group threshold via other
#          replicates), the position of the strongest coverage signal within
#          the promoter window in THIS replicate's own track is used as a
#          substitute summit ("we used the position of the strongest signal
#          around the TSS as the peak summit" - Donczew 2020)
#        - signal = sum of THIS replicate's own normalized coverage in a
#          fixed window (paper: 300bp) around that summit
#
#   4. The promoter's overall signal is the arithmetic mean of the per-
#      replicate signal values from step 3.
#
# WHY THIS READING (not the combine-first version from an earlier draft of
# this script): the missing-peak fallback rule in step 3 only makes sense if
# signal is computed independently per replicate. If replicates were pooled/
# combined before computing signal, a single replicate lacking its own peak
# wouldn't need a special substitute at all - the pooled data from the other
# replicates would already cover for it. The existence of a per-replicate
# fallback is the clearest evidence the paper's actual pipeline kept
# replicates separate all the way through step 3, only averaging at the very
# end (step 4).
#
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(rtracklayer)
    library(optparse)
})

# ==============================================================================
# CLI ARGUMENTS
# ==============================================================================

option_list <- list(
    make_option("--tss-bed", type = "character", default = NULL,
                help = "Path to TSS BED file (chr, start, end, locus_tag, score, strand) [required]"),
    make_option("--peaks-dir", type = "character", default = "results/peaks",
                help = "Root of the pipeline's results/peaks/{regulator}/{replicate}/ tree [default: %default]"),
    make_option("--coverage-dir", type = "character", default = "results/coverage",
                help = "Root of the pipeline's results/coverage/{regulator}/{replicate}/ tree [default: %default]"),
    make_option("--output-dir", type = "character", default = "results/hahn_region_scoring",
                help = "Output directory for the two result TSVs [default: %default]"),
    make_option("--promoter-upstream", type = "integer", default = 400,
                help = "bp upstream of TSS for the promoter window [default: %default]"),
    make_option("--promoter-downstream", type = "integer", default = 200,
                help = "bp downstream of TSS for the promoter window [default: %default]"),
    make_option("--peak-signal-window", type = "integer", default = 300,
                help = "Window width around the peak summit for signal sum [default: %default]"),
    make_option("--min-replicates-bound", type = "integer", default = 2,
                help = "Minimum replicates independently bound to call a promoter 'bound' [default: %default]"),
    make_option("--control-coverage", type = "character", default = NULL,
                help = "Path to the combined control's dmel-normalized bedgraph (e.g. results/coverage/combined_freemnase/combined_freemnase_dmel_norm.bedgraph). Optional - if provided, adds background_tag_count and enrichment (tagged/background ratio) columns; if omitted, only the raw tagged_tag_count (paper's literal 'signal') is reported."),
    make_option("--pseudocount", type = "double", default = 0.1,
                help = "Pseudocount added to background_tag_count when computing enrichment, to avoid division by zero [default: %default]"),
    make_option("--regulators", type = "character", default = NULL,
                help = "Comma-separated regulator list [default: auto-discover from --peaks-dir]"),
    make_option("--cores", type = "integer", default = NA,
                help = "Parallel cores [default: SLURM_CPUS_PER_TASK env var, or 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt[["tss-bed"]])) {
    stop("--tss-bed is required")
}

TSS_BED_PATH <- opt[["tss-bed"]]
PEAKS_DIR <- opt[["peaks-dir"]]
COVERAGE_DIR <- opt[["coverage-dir"]]
OUTPUT_DIR <- opt[["output-dir"]]
PROMOTER_UPSTREAM <- opt[["promoter-upstream"]]
PROMOTER_DOWNSTREAM <- opt[["promoter-downstream"]]
PEAK_SIGNAL_WINDOW <- opt[["peak-signal-window"]]
MIN_REPLICATES_BOUND <- opt[["min-replicates-bound"]]
CONTROL_COVERAGE_PATH <- opt[["control-coverage"]]
PSEUDOCOUNT <- opt[["pseudocount"]]
N_CORES <- if (!is.na(opt[["cores"]])) opt[["cores"]] else {
    env_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    if (nzchar(env_cores)) as.integer(env_cores) else 4L
}

if (!is.null(CONTROL_COVERAGE_PATH) && !file.exists(CONTROL_COVERAGE_PATH)) {
    stop(sprintf("--control-coverage file not found: %s", CONTROL_COVERAGE_PATH))
}

# NOTE on the missing-peak fallback window: neither methods text specifies
# whether "the position of the strongest signal around the TSS" means
# searching within the same promoter window used for peak assignment, or
# some other/smaller window specifically centered on the TSS. This script
# assumes the SAME promoter window (TSS -upstream/+downstream) for both,
# since that's the only window either text actually defines.

# ==============================================================================
# HELPERS
# ==============================================================================

#' Read the TSS BED file and build strand-aware promoter GRanges
#'
#' @param tss_bed_path Path to TSS BED file (chr, start, end, locus_tag, score, strand)
#' @param upstream bp upstream of TSS (relative to gene strand)
#' @param downstream bp downstream of TSS (relative to gene strand)
#' @return GRanges of promoter windows, one per TSS row, with mcols$locus_tag
build_promoter_granges <- function(tss_bed_path, upstream, downstream) {

    tss_tbl <- read_tsv(
        tss_bed_path,
        col_names = c("chr", "start", "end", "locus_tag", "score", "strand"),
        col_types = cols(chr = col_character(), start = col_integer(),
                         end = col_integer(), locus_tag = col_character(),
                         score = col_double(), strand = col_character())
    )

    # TSS BED rows are single-bp, 0-indexed half-open -> 1-indexed TSS position
    tss_pos <- tss_tbl$start + 1L

    # Strand-aware promoter window: "upstream"/"downstream" flip direction
    # depending on gene strand, since a BED file's start/end are always
    # genomic-coordinate order regardless of strand.
    is_plus <- tss_tbl$strand == "+"

    promoter_start <- ifelse(is_plus, tss_pos - upstream, tss_pos - downstream)
    promoter_end   <- ifelse(is_plus, tss_pos + downstream, tss_pos + upstream)

    GRanges(
        seqnames = tss_tbl$chr,
        ranges = IRanges(start = pmax(promoter_start, 1L), end = promoter_end),
        strand = tss_tbl$strand,
        locus_tag = tss_tbl$locus_tag,
        tss_pos = tss_pos
    )
}

#' Discover replicate names for a regulator from the peaks directory
#'
#' @param regulator Regulator/TF symbol
#' @param peaks_dir Root peaks directory
#' @return Character vector of replicate names (subdirectory names)
discover_replicates <- function(regulator, peaks_dir) {
    reg_dir <- file.path(peaks_dir, regulator)
    if (!dir.exists(reg_dir)) {
        stop(sprintf("No peaks directory found for regulator '%s': %s", regulator, reg_dir))
    }
    list.dirs(reg_dir, full.names = FALSE, recursive = FALSE)
}

#' Read a peaks_summits.txt file (04_pos2bed.sh output) as GRanges
#'
#' @param regulator Regulator/TF symbol
#' @param replicate Replicate name
#' @param peaks_dir Root peaks directory
#' @return GRanges of peak summits (1bp each), or an empty GRanges if the
#'   file is missing/empty (treated as "no peaks called for this replicate")
read_summits <- function(regulator, replicate, peaks_dir) {
    path <- file.path(peaks_dir, regulator, replicate,
                       sprintf("%s_%s_peaks_summits.txt", regulator, replicate))

    if (!file.exists(path)) {
        # A missing summits file is ambiguous on its own: it could mean this
        # replicate genuinely had zero peaks pass HOMER's thresholds (a real,
        # valid outcome), OR it could mean 03_findpeaks.sh/04_pos2bed.sh
        # hasn't actually finished for this replicate yet - e.g. a transient
        # network filesystem propagation delay right after the upstream
        # array job completes (the same class of issue the retry-with-
        # backoff in 06_hahn_region_scoring.sh handles for hard failures,
        # but which a silent "treat missing as zero peaks" here would mask
        # entirely, since the script would still exit 0 with silently
        # incomplete data).
        #
        # 03_findpeaks.sh's raw *_peaks.txt is always written once that step
        # genuinely completes, even when zero peaks pass thresholds (it's
        # the raw HOMER output with header/metadata, never literally absent
        # for a completed run) - use it as a sentinel to tell the two cases
        # apart: if it's ALSO missing, this replicate hasn't actually been
        # processed yet, and silently treating it as "zero peaks" would be
        # wrong - error instead of masking that.
        raw_peaks_path <- file.path(peaks_dir, regulator, replicate,
                                     sprintf("%s_%s_peaks.txt", regulator, replicate))
        if (!file.exists(raw_peaks_path)) {
            stop(sprintf(
                "Neither %s nor %s exist for %s_%s. This looks like this replicate hasn't actually finished processing yet (or a filesystem propagation delay) rather than a genuine zero-peaks result - re-run once 03_findpeaks.sh/04_pos2bed.sh have fully completed and their output is visible.",
                path, raw_peaks_path, regulator, replicate))
        }

        message(sprintf("  No summits (0 peaks passed thresholds) for %s_%s", regulator, replicate))
        return(GRanges())
    }

    if (file.size(path) == 0) {
        message(sprintf("  Empty summits file (0 peaks passed thresholds) for %s_%s", regulator, replicate))
        return(GRanges())
    }

    summits_tbl <- read_tsv(
        path,
        col_names = c("chr", "start", "end", "peak_id", "score", "strand"),
        col_types = cols(chr = col_character(), start = col_integer(),
                         end = col_integer(), peak_id = col_character(),
                         score = col_double(), strand = col_character())
    )

    if (nrow(summits_tbl) == 0) return(GRanges())

    GRanges(
        seqnames = summits_tbl$chr,
        ranges = IRanges(start = summits_tbl$start + 1L, end = summits_tbl$end),
        strand = "*",   # summit itself isn't meaningfully stranded
        peak_id = summits_tbl$peak_id
    )
}

#' Read a dmel-normalized bedgraph and return genome-wide coverage as an RleList
#'
#' @param regulator Regulator/TF symbol
#' @param replicate Replicate name
#' @param coverage_dir Root coverage directory
#' @return RleList of per-base normalized coverage (weighted by bedgraph score)
read_normalized_coverage <- function(regulator, replicate, coverage_dir) {
    path <- file.path(coverage_dir, regulator, replicate,
                       sprintf("%s_%s_dmel_norm.bedgraph", regulator, replicate))

    if (!file.exists(path)) {
        stop(sprintf("Coverage file not found: %s", path))
    }

    gr <- rtracklayer::import(path, format = "bedGraph")
    GenomicRanges::coverage(gr, weight = gr$score)
}

#' Read a dmel-normalized bedgraph from an explicit path (used for the
#' control coverage track, which isn't keyed by regulator/replicate)
#'
#' @param path Path to a *_dmel_norm.bedgraph file
#' @return RleList of per-base normalized coverage (weighted by bedgraph score)
read_normalized_coverage_from_path <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("Coverage file not found: %s", path))
    }

    gr <- rtracklayer::import(path, format = "bedGraph")
    GenomicRanges::coverage(gr, weight = gr$score)
}

#' Assign ONE replicate's own peak summits to promoters (no pooling)
#'
#' For each promoter, finds this replicate's summits within the promoter
#' window and picks the one closest to the TSS if there's more than one.
#'
#' @param promoters_gr Promoter GRanges (from build_promoter_granges)
#' @param summits_gr This replicate's summit GRanges (from read_summits)
#' @return Tibble: locus_tag, tss_pos, peak_summit_pos (NA if none),
#'   peak_assigned (logical)
assign_replicate_peaks_to_promoters <- function(promoters_gr, summits_gr) {

    n_promoters <- length(promoters_gr)
    peak_summit_pos <- rep(NA_integer_, n_promoters)

    if (length(summits_gr) > 0) {
        hits <- GenomicRanges::findOverlaps(promoters_gr, summits_gr)

        if (length(hits) > 0) {
            for (i in unique(S4Vectors::queryHits(hits))) {
                summit_idx <- S4Vectors::subjectHits(hits)[S4Vectors::queryHits(hits) == i]
                summit_positions <- GenomicRanges::start(summits_gr)[summit_idx]

                # Closest-to-TSS tie-break (paper's rule for multiple peaks
                # assigned to the same promoter)
                tss_pos_i <- promoters_gr$tss_pos[i]
                closest_idx <- which.min(abs(summit_positions - tss_pos_i))
                peak_summit_pos[i] <- summit_positions[closest_idx]
            }
        }
    }

    tibble(
        locus_tag = promoters_gr$locus_tag,
        tss_pos = promoters_gr$tss_pos,
        peak_summit_pos = peak_summit_pos,
        peak_assigned = !is.na(peak_summit_pos)
    )
}

#' Find the position of maximum coverage within each promoter window, using
#' one replicate's own coverage track. Used as the fallback substitute
#' summit for promoters where this replicate lacks its own peak.
#'
#' @param promoters_gr Promoter GRanges
#' @param coverage_rle This replicate's coverage RleList (from read_normalized_coverage)
#' @return Integer vector of argmax positions (1-indexed genomic coordinate),
#'   one per promoter (NA if the chromosome isn't present in the coverage track)
find_strongest_signal_position <- function(promoters_gr, coverage_rle) {

    chroms <- as.character(GenomicRanges::seqnames(promoters_gr))
    starts <- IRanges::start(promoters_gr)
    ends   <- IRanges::end(promoters_gr)

    purrr::pmap_int(list(chroms, starts, ends), function(chr, s, e) {
        if (!(chr %in% names(coverage_rle))) return(NA_integer_)

        chr_len <- length(coverage_rle[[chr]])
        win_start <- max(s, 1L)
        win_end <- min(e, chr_len)
        if (win_start > win_end) return(NA_integer_)

        window_vals <- as.numeric(coverage_rle[[chr]][win_start:win_end])
        win_start + which.max(window_vals) - 1L
    })
}

#' Sum one replicate's own normalized coverage in a window around a summit
#'
#' @param chroms Character vector of chromosomes, one per promoter
#' @param summit_pos Integer vector of summit positions, one per promoter (NA allowed)
#' @param coverage_rle This replicate's coverage RleList
#' @param window Total window width around the summit (paper default: 300)
#' @return Numeric vector of signal values, one per promoter (NA where summit_pos is NA)
compute_replicate_signal <- function(chroms, summit_pos, coverage_rle, window) {

    half_window_down <- window %/% 2
    half_window_up   <- window - half_window_down

    purrr::map2_dbl(chroms, summit_pos, function(chr, pos) {
        if (is.na(pos) || !(chr %in% names(coverage_rle))) return(NA_real_)

        chr_len <- length(coverage_rle[[chr]])
        win_start <- max(pos - half_window_down, 1L)
        win_end   <- min(pos + half_window_up - 1L, chr_len)
        if (win_start > win_end) return(NA_real_)

        sum(as.numeric(coverage_rle[[chr]][win_start:win_end]))
    })
}

# ==============================================================================
# MAIN DRIVER: per-regulator enrichment analysis
# ==============================================================================

#' Compute promoter binding/signal for one regulator, following the paper's
#' literal per-replicate-then-average approach (see header note).
#'
#' @param regulator Regulator/TF symbol
#' @param promoters_gr Promoter GRanges (shared across all regulators)
#' @param peaks_dir Root peaks directory
#' @param coverage_dir Root coverage directory
#' @param peak_signal_window Window width around summit for signal sum
#' @param min_replicates_bound Minimum replicates independently bound to call "bound"
#' @param control_coverage_rle Optional control coverage RleList (from
#'   read_normalized_coverage_from_path()). If supplied, adds
#'   background_tag_count (control's coverage summed at the SAME window as
#'   tagged_tag_count) and enrichment (their ratio) to the output.
#' @param pseudocount Added to background_tag_count when computing enrichment
#' @return list(per_replicate = tibble of per-replicate values (long format),
#'   summary = tibble of one row per promoter: n_replicates_bound, bound,
#'   tagged_tag_count (mean across replicates), and, if control_coverage_rle
#'   was supplied, background_tag_count and enrichment (also averaged))
promoter_enrichment_for_regulator <- function(regulator,
                                                promoters_gr,
                                                peaks_dir = PEAKS_DIR,
                                                coverage_dir = COVERAGE_DIR,
                                                peak_signal_window = PEAK_SIGNAL_WINDOW,
                                                min_replicates_bound = MIN_REPLICATES_BOUND,
                                                control_coverage_rle = NULL,
                                                pseudocount = PSEUDOCOUNT) {

    message(sprintf("Working on regulator: %s", regulator))

    replicates <- discover_replicates(regulator, peaks_dir)
    message(sprintf("  Found %d replicate(s): %s", length(replicates), paste(replicates, collapse = ", ")))

    chroms <- as.character(GenomicRanges::seqnames(promoters_gr))

    # ---- Step 1: per-replicate, independent peak assignment ----
    peak_assignment_by_rep <- purrr::map(
        replicates, ~assign_replicate_peaks_to_promoters(promoters_gr, read_summits(regulator, .x, peaks_dir))
    )
    names(peak_assignment_by_rep) <- replicates

    # ---- Step 2: bound/unbound (group-level threshold on independent peak calls) ----
    n_replicates_bound <- Reduce(`+`, purrr::map(peak_assignment_by_rep, "peak_assigned"))
    bound <- n_replicates_bound >= min_replicates_bound

    # ---- Step 3: per-replicate signal, only for bound promoters ----
    message("  Computing per-replicate signal (with missing-peak fallback)...")
    coverage_by_rep <- purrr::map(replicates, ~read_normalized_coverage(regulator, .x, coverage_dir))
    names(coverage_by_rep) <- replicates

    per_replicate_tbl <- purrr::map_dfr(replicates, function(rep_name) {
        assignment <- peak_assignment_by_rep[[rep_name]]
        cov_rle <- coverage_by_rep[[rep_name]]

        # Final summit for this replicate: its own peak summit if assigned;
        # otherwise the strongest-signal fallback from its own coverage
        # track - but ONLY for promoters that cleared the group bound
        # threshold (no point computing a fallback for promoters that will
        # be excluded anyway).
        final_summit_pos <- assignment$peak_summit_pos
        needs_fallback <- bound & !assignment$peak_assigned

        if (any(needs_fallback)) {
            fallback_pos <- find_strongest_signal_position(
                promoters_gr[needs_fallback], cov_rle
            )
            final_summit_pos[needs_fallback] <- fallback_pos
        }

        # Signal only computed for bound promoters
        final_summit_pos[!bound] <- NA_integer_

        # tagged_tag_count: this replicate's own dmel-normalized coverage,
        # summed in the window around its final summit - the paper's
        # literal "signal per promoter" quantity, on its own.
        tagged_tag_count <- compute_replicate_signal(chroms, final_summit_pos, cov_rle, peak_signal_window)

        result <- tibble(
            replicate = rep_name,
            locus_tag = promoters_gr$locus_tag,
            tss_pos = promoters_gr$tss_pos,
            peak_assigned = assignment$peak_assigned,
            used_fallback_summit = needs_fallback,
            summit_pos = final_summit_pos,
            tagged_tag_count = tagged_tag_count
        )

        # background_tag_count / enrichment are only added if a control
        # coverage track was supplied - the paper's own "signal" metric
        # doesn't require a control at all (see script header), so this is
        # an addition on top of, not a replacement for, tagged_tag_count.
        # Computed with the SAME window-sum function, at the SAME summit
        # position, just reading from the control's coverage instead - so
        # tagged and background are directly comparable, same window.
        if (!is.null(control_coverage_rle)) {
            background_tag_count <- compute_replicate_signal(
                chroms, final_summit_pos, control_coverage_rle, peak_signal_window
            )
            result$background_tag_count <- background_tag_count
            result$enrichment <- tagged_tag_count / (background_tag_count + pseudocount)
        }

        result
    })

    # ---- Step 4: overall values per promoter = mean across replicates ----
    avg_cols <- intersect(c("tagged_tag_count", "background_tag_count", "enrichment"), names(per_replicate_tbl))

    summary_tbl <- per_replicate_tbl %>%
        group_by(locus_tag, tss_pos) %>%
        summarise(across(all_of(avg_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
        mutate(across(all_of(avg_cols), ~ifelse(is.nan(.x), NA_real_, .x)))

    summary_tbl <- tibble(
        locus_tag = promoters_gr$locus_tag,
        tss_pos = promoters_gr$tss_pos,
        n_replicates_bound = n_replicates_bound,
        bound = bound
    ) %>%
        left_join(summary_tbl, by = c("locus_tag", "tss_pos")) %>%
        mutate(regulator = regulator) %>%
        relocate(regulator, locus_tag)

    per_replicate_tbl <- per_replicate_tbl %>%
        mutate(regulator = regulator) %>%
        relocate(regulator, replicate, locus_tag)

    list(per_replicate = per_replicate_tbl, summary = summary_tbl)
}

# ==============================================================================
# MAIN DRIVER
# =============================================================================

message(sprintf("Building promoter windows from %s (upstream=%d, downstream=%d)...",
                 TSS_BED_PATH, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM))
promoters_gr <- build_promoter_granges(TSS_BED_PATH, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM)

regulator_list <- if (!is.null(opt$regulators)) {
    strsplit(opt$regulators, ",")[[1]]
} else {
    list.dirs(PEAKS_DIR, full.names = FALSE, recursive = FALSE)
}

# Load the control coverage ONCE (not per-regulator) if provided, since it's
# a single shared background track regardless of which regulator is being
# scored - same rationale as the calling-cards script loading its control
# once up front.
control_coverage_rle <- NULL
if (!is.null(CONTROL_COVERAGE_PATH)) {
    message(sprintf("Loading control coverage from %s...", CONTROL_COVERAGE_PATH))
    control_coverage_rle <- tryCatch({
        read_normalized_coverage_from_path(CONTROL_COVERAGE_PATH)
    }, error = function(e) {
        stop(sprintf("Failed to load --control-coverage at %s: %s", CONTROL_COVERAGE_PATH, conditionMessage(e)))
    })
} else {
    message("No --control-coverage provided - only tagged_tag_count (paper's literal 'signal') will be reported, no background/enrichment.")
}

message(sprintf("Scoring %d regulator(s) using %d core(s)...", length(regulator_list), N_CORES))

all_results <- parallel::mclapply(
    regulator_list,
    promoter_enrichment_for_regulator,
    promoters_gr = promoters_gr,
    peaks_dir = PEAKS_DIR,
    coverage_dir = COVERAGE_DIR,
    peak_signal_window = PEAK_SIGNAL_WINDOW,
    min_replicates_bound = MIN_REPLICATES_BOUND,
    control_coverage_rle = control_coverage_rle,
    pseudocount = PSEUDOCOUNT,
    mc.cores = N_CORES
)
names(all_results) <- regulator_list

# Surface any per-regulator errors clearly rather than silently dropping them
# (mclapply returns a try-error object for a failed element instead of stopping)
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

    summary_tbl <- all_results[[regulator]]$summary
    per_replicate_tbl <- all_results[[regulator]]$per_replicate

    summary_path <- file.path(reg_dir, sprintf("%s_summary.tsv", regulator))
    per_replicate_path <- file.path(reg_dir, sprintf("%s_replicate.tsv", regulator))

    readr::write_tsv(summary_tbl, summary_path)
    readr::write_tsv(per_replicate_tbl, per_replicate_path)

    message(sprintf("Wrote %s: summary (%d rows), per-replicate (%d rows)",
                     regulator, nrow(summary_tbl), nrow(per_replicate_tbl)))
}

message("Done.")
