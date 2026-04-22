#!/usr/bin/env Rscript
# inst/scripts/convert_sample.R
#
# CLI script for converting a single BAM/CRAM into a binned BED.gz or NPZ file.
#
# Supports RWisecondorX, upstream WisecondorX, and NIPTeR modes.
#
# Usage:
#   Rscript convert_sample.R --mode rwisecondorx --bam sample.bam --out sample.bed.gz
#   Rscript convert_sample.R --mode rwisecondorx --bam sample.bam --out sample.npz --npz
#   Rscript convert_sample.R --mode wisecondorx --bam sample.bam --out sample.npz
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz --gc-table hg38_gc.tsv.bgz
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz --fasta hg38.fa
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz \
#     --compute-seqff --compute-y-unique --sample-qc-out sample.sample_qc.tsv
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz \
#     --seqff-tsv cohort.seqff.tsv --y-unique-tsv cohort.y_unique.tsv --y-regions-file grch37_Y_chrom_blacklist.bed

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with: install.packages('optparse')",
       call. = FALSE)
}
if (!requireNamespace("this.path", quietly = TRUE)) {
  stop("this.path is required. Install it with: install.packages('this.path')",
       call. = FALSE)
}

.script_path <- normalizePath(this.path::this.path(), winslash = "/", mustWork = TRUE)
.script_dir <- dirname(.script_path)

library(optparse)

.compact_metadata <- function(x) {
  keep <- !vapply(x, is.null, logical(1L))
  x[keep]
}

option_list <- list(
  make_option("--mode", type = "character", default = "rwisecondorx",
              help = "Conversion mode: 'rwisecondorx', 'wisecondorx', or 'nipter' [default: %default]"),

  # --- Input / output ---
  make_option("--bam", type = "character", default = NULL,
              help = "Input BAM or CRAM file [required]"),
  make_option("--out", type = "character", default = NULL,
              help = "Output path (.bed.gz or .npz) [required]"),

  # --- Output format ---
  make_option("--npz", action = "store_true", default = FALSE,
              help = paste0("Write WisecondorX-compatible NPZ instead of BED.gz. ",
                            "Only valid in rwisecondorx mode [default: %default]")),

  # --- Binning parameters ---
  make_option("--binsize", type = "integer", default = NULL,
              help = paste0("Bin size in bp. Defaults: 100000 (wisecondorx), ",
                            "50000 (nipter) [default: mode-dependent]")),
  make_option("--mapq", type = "integer", default = NULL,
              help = "Minimum MAPQ [default: 1 wisecondorx, 0 nipter]"),
  make_option("--rmdup", type = "character", default = NULL,
              help = paste0("Dedup mode: 'streaming', 'flag', 'none'. ",
                            "Defaults: 'streaming' (wisecondorx), 'none' (nipter)")),
  make_option("--exclude-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask to exclude (samtools -F) [default: %default]"),
  make_option("--require-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask to require (samtools -f) [default: %default]"),
  make_option("--reference", type = "character", default = NULL,
              help = "FASTA reference for CRAM decoding"),

  # --- NIPTeR-specific ---
  make_option("--separate-strands", action = "store_true", default = FALSE,
              help = "NIPTeR mode: produce SeparatedStrands 9-column BED [default: %default]"),
  make_option("--gc-table", type = "character", default = NULL,
              help = paste0("Precomputed GC table (TSV.bgz from nipter_gc_precompute). ",
                            "NIPTeR mode: GC-correct and embed corrected counts in BED")),
  make_option("--fasta", type = "character", default = NULL,
              help = paste0("FASTA reference for on-the-fly GC correction (alternative to --gc-table). ",
                            "NIPTeR mode only. Slower than --gc-table for multiple samples")),
  make_option("--nipter-gc-include-sex", action = "store", type = "logical", default = TRUE,
              help = paste0("NIPTeR mode: also GC-correct X/Y bins when writing corrected BED columns. ",
                            "Pass --nipter-gc-include-sex=FALSE to disable [default: %default]")),
  make_option("--nipter-gc-loess-span", type = "double", default = 0.75,
              help = "NIPTeR mode: LOESS span used for GC correction [default: %default]"),
  make_option("--gc-curve-plot-out", type = "character", default = NULL,
              help = paste0("Optional PNG path for a per-sample NIPTeR GC-curve plot. ",
                            "Only used in nipter mode when GC correction is applied.")),
  make_option("--corrected-bin-ratio-genome-plot-out", type = "character", default = NULL,
              help = paste0(
                "Optional PNG path for a genome-order plot of post-GC corrected bin ratios ",
                "excluding out-of-range bins. Only used in nipter mode when GC correction is applied."
              )),
  make_option("--qc-plot-theme", type = "character", default = "minimal",
              help = "QC plot theme: minimal, light, bw, classic [default: %default]"),
  make_option("--qc-plot-base-size", type = "integer", default = 11L,
              help = "QC plot base font size [default: %default]"),
  # --- Optional per-sample FF / Y-unique metrics ---
  make_option("--compute-seqff", action = "store_true", default = FALSE,
              help = paste0("Compute SeqFF FF metrics twice (nonfiltered + filtered) for this sample. ",
                            "Can be written to --sample-qc-out and/or BED metadata [default: %default]")),
  make_option("--compute-y-unique", action = "store_true", default = FALSE,
              help = paste0("Compute Y-unique ratios twice (nonfiltered + filtered) for this sample. ",
                            "Can be written to --sample-qc-out and/or BED metadata [default: %default]")),
  make_option("--seqff-tsv", type = "character", default = NULL,
              help = paste0("Optional TSV/CSV with per-sample SeqFF metrics. ",
                            "Must contain the explicit *_pre/*_post columns written by RWisecondorX.")),
  make_option("--y-unique-tsv", type = "character", default = NULL,
              help = paste0("Optional TSV/CSV with per-sample Y-unique metrics. ",
                            "Must contain the explicit *_pre/*_post columns written by RWisecondorX.")),
  make_option("--sample-qc-out", type = "character", default = NULL,
              help = "Optional one-row TSV path to write the canonical sample QC row for this conversion."),
  make_option("--y-regions-file", type = "character", default = NULL,
              help = paste0("Optional Y-unique regions file used for --compute-y-unique and metadata provenance. ",
                            "Accepts headered TSV or headerless BED-like files.")),
  make_option("--metrics-pre-mapq", type = "integer", default = 0L,
              help = "MAPQ used for nonfiltered SeqFF/Y-unique metrics [default: %default]"),
  make_option("--metrics-pre-require-flags", type = "integer", default = 0L,
              help = "Required flags used for nonfiltered SeqFF/Y-unique metrics [default: %default]"),
  make_option("--metrics-pre-exclude-flags", type = "integer", default = 0L,
              help = "Excluded flags used for nonfiltered SeqFF/Y-unique metrics [default: %default]"),
  make_option("--metrics-post-mapq", type = "integer", default = NULL,
              help = paste0(
                "MAPQ used for filtered SeqFF/Y-unique metrics. ",
                "Defaults to 40 in nipter mode and to the main conversion --mapq otherwise."
              )),
  make_option("--metrics-post-require-flags", type = "integer", default = NULL,
              help = paste0(
                "Required flags used for filtered SeqFF/Y-unique metrics. ",
                "Defaults to 0 in nipter mode and to the main conversion --require-flags otherwise."
              )),
  make_option("--metrics-post-exclude-flags", type = "integer", default = NULL,
              help = paste0(
                "Excluded flags used for filtered SeqFF/Y-unique metrics. ",
                "Defaults to 1024 in nipter mode and to the main conversion --exclude-flags otherwise."
              ))
)

parser <- OptionParser(
  usage = "usage: %prog [options]",
  description = paste(
    "Convert a single BAM/CRAM file to a binned BED.gz or NPZ file.",
    "",
    "Supports RWisecondorX, upstream WisecondorX, and NIPTeR conversion modes.",
    "",
    "Examples:",
    "  # Native RWisecondorX BED output",
    "  %prog --mode rwisecondorx --bam sample.bam --out sample.bed.gz",
    "",
    "  # Native RWisecondorX NPZ output",
    "  %prog --mode rwisecondorx --bam sample.bam --out sample.npz --npz",
    "",
    "  # Upstream WisecondorX NPZ output via condathis wrapper",
    "  %prog --mode wisecondorx --bam sample.bam --out sample.npz",
    "",
    "  # NIPTeR BED output with GC correction (precomputed table)",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz --gc-table hg38_gc.tsv.bgz",
    "",
    "  # NIPTeR BED output with GC correction (on-the-fly from FASTA)",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz --fasta hg38.fa",
    "",
    "  # NIPTeR BED output with GC-corrected sex chromosomes too",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz --gc-table hg38_gc.tsv.bgz --nipter-gc-include-sex",
    "",
    "  # NIPTeR with pre-filtering (MAPQ 40, exclude duplicates)",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz --mapq 40 --exclude-flags 1024",
    "",
    "  # NIPTeR BED with FF/Y-unique computed on the fly and written both to header metadata and a sidecar TSV",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz \\",
    "    --compute-seqff --compute-y-unique --sample-qc-out sample.sample_qc.tsv",
    "",
    "  # NIPTeR BED with FF/Y-unique loaded from cohort summary tables and explicit Y-region provenance",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz \\",
    "    --seqff-tsv cohort.seqff.tsv --y-unique-tsv cohort.y_unique.tsv \\",
    "    --y-regions-file grch37_Y_chrom_blacklist.bed",
    sep = "\n"
  ),
  option_list = option_list
)

opts <- parse_args(parser)

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

mode <- tolower(opts$mode)
if (!mode %in% c("rwisecondorx", "wisecondorx", "nipter")) {
  stop("--mode must be 'rwisecondorx', 'wisecondorx', or 'nipter', got: ", opts$mode, call. = FALSE)
}

if (is.null(opts$bam)) {
  stop("--bam is required (path to input BAM/CRAM file)", call. = FALSE)
}

if (!file.exists(opts$bam)) {
  stop("BAM file not found: ", opts$bam, call. = FALSE)
}

if (is.null(opts$out)) {
  stop("--out is required (output path for BED.gz or NPZ file)", call. = FALSE)
}

if (opts$npz && mode != "rwisecondorx") {
  stop("--npz is only valid in rwisecondorx mode", call. = FALSE)
}

if (!is.null(opts$`gc-table`) && mode != "nipter") {
  stop("--gc-table is only valid in nipter mode", call. = FALSE)
}

if (!is.null(opts$fasta) && mode != "nipter") {
  stop("--fasta GC correction is only valid in nipter mode", call. = FALSE)
}

if (is.na(opts$`nipter-gc-include-sex`)) {
  stop("--nipter-gc-include-sex must be TRUE or FALSE.", call. = FALSE)
}

if (isTRUE(opts$`nipter-gc-include-sex`) && mode != "nipter") {
  stop("--nipter-gc-include-sex is only valid in nipter mode", call. = FALSE)
}
if (!is.null(opts$`gc-curve-plot-out`) && mode != "nipter") {
  stop("--gc-curve-plot-out is only valid in nipter mode", call. = FALSE)
}
if (!is.null(opts$`corrected-bin-ratio-genome-plot-out`) && mode != "nipter") {
  stop("--corrected-bin-ratio-genome-plot-out is only valid in nipter mode", call. = FALSE)
}

if (!is.null(opts$`gc-table`) && !is.null(opts$fasta)) {
  stop("--gc-table and --fasta are mutually exclusive; provide one or the other", call. = FALSE)
}
if (!is.numeric(opts$`nipter-gc-loess-span`) ||
    length(opts$`nipter-gc-loess-span`) != 1L ||
    !is.finite(opts$`nipter-gc-loess-span`) ||
    opts$`nipter-gc-loess-span` <= 0 ||
    opts$`nipter-gc-loess-span` > 1) {
  stop("--nipter-gc-loess-span must be a finite number in (0, 1].", call. = FALSE)
}

if (opts$`separate-strands` && mode != "nipter") {
  stop("--separate-strands is only valid in nipter mode", call. = FALSE)
}

if (mode == "wisecondorx" && !grepl("\\.npz$", opts$out, ignore.case = TRUE)) {
  stop("wisecondorx mode writes upstream NPZ output; --out must end in .npz", call. = FALSE)
}

if (mode == "wisecondorx" && opts$npz) {
  stop("wisecondorx mode already uses the upstream NPZ wrapper; do not pass --npz", call. = FALSE)
}

for (path_opt in c("seqff-tsv", "y-unique-tsv", "y-regions-file")) {
  path_val <- opts[[path_opt]]
  if (!is.null(path_val) && !file.exists(path_val)) {
    stop("--", path_opt, " does not exist: ", path_val, call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Set mode-dependent defaults
# ---------------------------------------------------------------------------

if (mode == "rwisecondorx") {
  binsize <- if (is.null(opts$binsize)) 100000L else opts$binsize
  mapq    <- if (is.null(opts$mapq)) 1L else opts$mapq
  rmdup   <- if (is.null(opts$rmdup)) "streaming" else opts$rmdup
} else if (mode == "wisecondorx") {
  binsize <- if (is.null(opts$binsize)) 100000L else opts$binsize
  mapq    <- if (is.null(opts$mapq)) 1L else opts$mapq
  rmdup   <- if (is.null(opts$rmdup)) "streaming" else opts$rmdup
} else {
  binsize <- if (is.null(opts$binsize)) 50000L else opts$binsize
  mapq    <- if (is.null(opts$mapq)) 0L else opts$mapq
  rmdup   <- if (is.null(opts$rmdup)) "none" else opts$rmdup
}

library(RWisecondorX)

.match_metric_row <- getFromNamespace(".match_metric_row", "RWisecondorX")
.seqff_metric_record <- getFromNamespace(".seqff_metric_record", "RWisecondorX")
.y_unique_metric_record <- getFromNamespace(".y_unique_metric_record", "RWisecondorX")
.sample_qc_row <- getFromNamespace(".sample_qc_row", "RWisecondorX")
.sample_qc_row_with_processing_status <- getFromNamespace(".sample_qc_row_with_processing_status", "RWisecondorX")
.sample_qc_row_with_nipter_status <- getFromNamespace(".sample_qc_row_with_nipter_status", "RWisecondorX")
.sample_qc_tabix_metadata <- getFromNamespace(".sample_qc_tabix_metadata", "RWisecondorX")
.merge_tabix_metadata <- getFromNamespace(".merge_tabix_metadata", "RWisecondorX")
.samtools_stats_metric_record <- getFromNamespace(".samtools_stats_metric_record", "RWisecondorX")
.read_counts_metric_record <- getFromNamespace(".read_counts_metric_record", "RWisecondorX")
.bam_bin_count_rows <- getFromNamespace(".bam_bin_count_rows", "RWisecondorX")
.bam_chr_lengths <- getFromNamespace(".bam_chr_lengths", "RWisecondorX")
.nipter_rows_to_sample <- getFromNamespace(".nipter_rows_to_sample", "RWisecondorX")
.bam_alignment_count_metadata <- getFromNamespace(".bam_alignment_count_metadata", "RWisecondorX")
.bam_read_counts_metadata <- getFromNamespace(".bam_read_counts_metadata", "RWisecondorX")
.nipter_preprocess_qc <- getFromNamespace(".nipter_preprocess_qc", "RWisecondorX")
.write_nipter_gc_curve_plot <- getFromNamespace(".write_nipter_gc_curve_plot", "RWisecondorX")
.write_nipter_corrected_bin_ratio_genome_plot <- getFromNamespace(".write_nipter_corrected_bin_ratio_genome_plot", "RWisecondorX")

opts$`qc-plot-theme` <- match.arg(
  tolower(opts$`qc-plot-theme`),
  c("minimal", "light", "bw", "classic")
)

.resolve_y_regions_file <- function(path = NULL) {
  if (!is.null(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  normalizePath(system.file(
    "extdata",
    "grch37_Y_chrom_blacklist.bed",
    package = "RWisecondorX",
    mustWork = TRUE
  ), winslash = "/", mustWork = TRUE)
}

.compute_read_counts_metrics <- function(bam,
                                         binsize,
                                         mapq,
                                         require_flags,
                                         exclude_flags,
                                         rmdup,
                                         reference) {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  Rduckhts::rduckhts_load(con)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  read_count_rows <- .bam_bin_count_rows(
    con = con,
    bam = bam,
    binsize = binsize,
    mapq = mapq,
    require_flags = as.integer(require_flags),
    exclude_flags = as.integer(exclude_flags),
    rmdup = rmdup,
    reference = reference,
    stats = "gc,mq",
    include_unmapped = TRUE
  )

  .read_counts_metric_record(
    record = .merge_tabix_metadata(
      .bam_alignment_count_metadata(con, bam),
      .bam_read_counts_metadata(read_count_rows)
    ),
    source = "samtools_idxstats+bam_bin_counts(gc,mq)"
  )
}

.compute_bam_samtools_stats_metrics <- function(bam,
                                                reference,
                                                filters_pre,
                                                filters_post) {
  .samtools_stats_metric_record(
    pre = bam_samtools_stats_summary(
      path = bam,
      reference = reference,
      min_mapq = as.integer(filters_pre$mapq),
      require_flags = as.integer(filters_pre$require_flags),
      exclude_flags = as.integer(filters_pre$exclude_flags)
    ),
    post = bam_samtools_stats_summary(
      path = bam,
      reference = reference,
      min_mapq = as.integer(filters_post$mapq),
      require_flags = as.integer(filters_post$require_flags),
      exclude_flags = as.integer(filters_post$exclude_flags),
      report_filtered_stream_bookkeeping = TRUE
    ),
    source = "native_htslib_samtools_stats"
  )
}

# ---------------------------------------------------------------------------
# Convert
# ---------------------------------------------------------------------------

bam <- opts$bam
out <- opts$out
sample_name <- sub("\\.cram$|\\.bam$", "", basename(bam), ignore.case = TRUE)

gc_fasta <- if (mode == "nipter") opts$fasta else NULL
metrics_reference <- if (is.null(opts$reference)) gc_fasta else opts$reference
y_regions_file <- if (!is.null(opts$`y-unique-tsv`) || isTRUE(opts$`compute-y-unique`)) {
  .resolve_y_regions_file(opts$`y-regions-file`)
} else {
  NULL
}
metrics_post_defaults <- if (identical(mode, "nipter")) {
  list(mapq = 40L, require_flags = 0L, exclude_flags = 1024L)
} else {
  list(
    mapq = as.integer(mapq),
    require_flags = as.integer(opts$`require-flags`),
    exclude_flags = as.integer(opts$`exclude-flags`)
  )
}
metrics_pre_filters <- list(
  mapq = as.integer(opts$`metrics-pre-mapq`),
  require_flags = as.integer(opts$`metrics-pre-require-flags`),
  exclude_flags = as.integer(opts$`metrics-pre-exclude-flags`)
)
metrics_post_filters <- list(
  mapq = as.integer(if (is.null(opts$`metrics-post-mapq`)) metrics_post_defaults$mapq else opts$`metrics-post-mapq`),
  require_flags = as.integer(if (is.null(opts$`metrics-post-require-flags`)) metrics_post_defaults$require_flags else opts$`metrics-post-require-flags`),
  exclude_flags = as.integer(if (is.null(opts$`metrics-post-exclude-flags`)) metrics_post_defaults$exclude_flags else opts$`metrics-post-exclude-flags`)
)

needs_sample_qc <- !is.null(opts$`sample-qc-out`) ||
  !is.null(opts$`seqff-tsv`) || !is.null(opts$`y-unique-tsv`) ||
  isTRUE(opts$`compute-seqff`) || isTRUE(opts$`compute-y-unique`)

qc_state <- new.env(parent = emptyenv())
qc_state$seqff_metrics <- NULL
qc_state$y_unique_metrics <- NULL
qc_state$nipter_qc_metrics <- NULL
qc_state$read_counts_metrics <- NULL
qc_state$bam_samtools_stats_metrics <- NULL
qc_state$gc_correction <- NULL

conversion_status <- new.env(parent = emptyenv())
conversion_status$ok <- FALSE
conversion_status$stage <- "initialization"
conversion_status$step <- "startup"
conversion_status$message <- "Conversion did not complete."

.write_sample_qc_on_exit <- function() {
  if (is.null(opts$`sample-qc-out`)) {
    return(invisible(NULL))
  }

  row <- .sample_qc_row(
    sample_name = sample_name,
    bam = bam,
    seqff = qc_state$seqff_metrics,
    y_unique = qc_state$y_unique_metrics,
    read_counts = qc_state$read_counts_metrics,
    bam_stats = qc_state$bam_samtools_stats_metrics,
    nipter_qc = qc_state$nipter_qc_metrics,
    gc_correction = qc_state$gc_correction,
    filters_pre = metrics_pre_filters,
    filters_post = metrics_post_filters
  )
  row <- .sample_qc_row_with_processing_status(
    row,
    status = if (isTRUE(conversion_status$ok)) "ok" else "failed",
    failure_stage = if (isTRUE(conversion_status$ok)) NULL else conversion_status$stage,
    failure_step = if (isTRUE(conversion_status$ok)) NULL else conversion_status$step,
    failure_message = if (isTRUE(conversion_status$ok)) NULL else conversion_status$message
  )
  row <- .sample_qc_row_with_nipter_status(
    row,
    status = if (identical(mode, "nipter")) {
      if (isTRUE(conversion_status$ok)) "ok" else "failed"
    } else {
      "not_run"
    },
    error = if (identical(mode, "nipter") && !isTRUE(conversion_status$ok)) {
      conversion_status$message
    } else {
      NULL
    }
  )

  utils::write.table(
    row,
    file = opts$`sample-qc-out`,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  invisible(NULL)
}

main <- function() {
if (needs_sample_qc) {
  conversion_status$stage <- "samtools_stats"
  conversion_status$step <- "summary"
  qc_state$bam_samtools_stats_metrics <- .compute_bam_samtools_stats_metrics(
    bam = bam,
    reference = metrics_reference,
    filters_pre = metrics_pre_filters,
    filters_post = metrics_post_filters
  )
  if (!is.null(opts$`seqff-tsv`)) {
    conversion_status$stage <- "seqff"
    conversion_status$step <- "load_from_file"
    qc_state$seqff_metrics <- .seqff_metric_record(
      record = .match_metric_row(opts$`seqff-tsv`, sample_name),
      source = paste0("file:", basename(opts$`seqff-tsv`))
    )
  } else if (isTRUE(opts$`compute-seqff`)) {
    conversion_status$stage <- "seqff"
    conversion_status$step <- "compute"
    qc_state$seqff_metrics <- .seqff_metric_record(
      pre = seqff_predict(
        input = bam,
        input_type = "bam",
        mapq = metrics_pre_filters$mapq,
        require_flags = metrics_pre_filters$require_flags,
        exclude_flags = metrics_pre_filters$exclude_flags,
        reference = metrics_reference
      ),
      post = seqff_predict(
        input = bam,
        input_type = "bam",
        mapq = metrics_post_filters$mapq,
        require_flags = metrics_post_filters$require_flags,
        exclude_flags = metrics_post_filters$exclude_flags,
        reference = metrics_reference
      ),
      source = "computed"
    )
  }

  if (!is.null(opts$`y-unique-tsv`)) {
    conversion_status$stage <- "y_unique"
    conversion_status$step <- "load_from_file"
    qc_state$y_unique_metrics <- .y_unique_metric_record(
      record = .match_metric_row(opts$`y-unique-tsv`, sample_name),
      regions_file = y_regions_file,
      source = paste0("file:", basename(opts$`y-unique-tsv`))
    )
  } else if (isTRUE(opts$`compute-y-unique`)) {
    conversion_status$stage <- "y_unique"
    conversion_status$step <- "compute"
    qc_state$y_unique_metrics <- .y_unique_metric_record(
      pre = nipter_y_unique_ratio(
        bam = bam,
        mapq = metrics_pre_filters$mapq,
        require_flags = metrics_pre_filters$require_flags,
        exclude_flags = metrics_pre_filters$exclude_flags,
        regions_file = y_regions_file,
        reference = metrics_reference
      ),
      post = nipter_y_unique_ratio(
        bam = bam,
        mapq = metrics_post_filters$mapq,
        require_flags = metrics_post_filters$require_flags,
        exclude_flags = metrics_post_filters$exclude_flags,
        regions_file = y_regions_file,
        reference = metrics_reference
      ),
      regions_file = y_regions_file,
      source = "computed"
    )
  }

  if (mode != "nipter") {
    conversion_status$stage <- "read_counts"
    conversion_status$step <- "summarize"
    qc_state$read_counts_metrics <- .compute_read_counts_metrics(
      bam = bam,
      binsize = binsize,
      mapq = mapq,
      require_flags = opts$`require-flags`,
      exclude_flags = opts$`exclude-flags`,
      rmdup = rmdup,
      reference = opts$reference
    )
  }
}

if (mode == "rwisecondorx") {
  conversion_status$stage <- if (isTRUE(opts$npz)) "rwisecondorx_npz" else "rwisecondorx_bed"
  conversion_status$step <- "metadata"
  bed_metadata <- if (!isTRUE(opts$npz)) {
    .compact_metadata(list(
      format = "rwisecondorx_bed",
      schema = "count_v1",
      binsize = as.integer(binsize),
      mapq = as.integer(mapq),
      require_flags = as.integer(opts$`require-flags`),
      exclude_flags = as.integer(opts$`exclude-flags`),
      rmdup = rmdup,
      reference = if (!is.null(opts$reference)) normalizePath(opts$reference, winslash = "/", mustWork = TRUE) else NULL
    ))
  } else NULL
  if (!is.null(bed_metadata) && needs_sample_qc) {
    bed_metadata <- .merge_tabix_metadata(
      bed_metadata,
      .sample_qc_tabix_metadata(
        seqff = qc_state$seqff_metrics,
        y_unique = qc_state$y_unique_metrics,
        bam_stats = qc_state$bam_samtools_stats_metrics,
        filters_pre = metrics_pre_filters,
        filters_post = metrics_post_filters
      )
    )
  }

  if (opts$npz) {
    conversion_status$step <- "convert"
    cat(sprintf("Converting %s → %s (native rwisecondorx NPZ, binsize=%d, rmdup=%s)\n",
                basename(bam), basename(out), binsize, rmdup))
    bam_convert_npz(
      bam       = bam,
      npz       = out,
      binsize   = binsize,
      mapq      = mapq,
      require_flags = opts$`require-flags`,
      exclude_flags = opts$`exclude-flags`,
      rmdup     = rmdup,
      reference = opts$reference
    )
  } else {
    conversion_status$step <- "convert"
    cat(sprintf("Converting %s → %s (native rwisecondorx BED, binsize=%d, mapq=%d, rmdup=%s)\n",
                basename(bam), basename(out), binsize, mapq, rmdup))
    bam_convert_bed(
      bam           = bam,
      bed           = out,
      binsize       = binsize,
      mapq          = mapq,
      require_flags = opts$`require-flags`,
      exclude_flags = opts$`exclude-flags`,
      rmdup         = rmdup,
      reference     = opts$reference,
      metadata      = bed_metadata
    )
    if (needs_sample_qc) {
      conversion_status$step <- "read_tabix_metadata"
      qc_state$read_counts_metrics <- .read_counts_metric_record(
        record = as.list(tabix_metadata(out)),
        source = "samtools_idxstats+bam_bin_counts(gc,mq)"
      )
    }
  }
} else if (mode == "wisecondorx") {
  conversion_status$stage <- "wisecondorx_npz"
  conversion_status$step <- "convert"
  cat(sprintf("Converting %s → %s (upstream wisecondorx NPZ wrapper, binsize=%d)\n",
              basename(bam), basename(out), binsize))
  wisecondorx_convert(
    bam = bam,
    npz = out,
    reference = opts$reference,
    binsize = binsize,
    normdup = identical(rmdup, "none"),
    extra_args = character(0)
  )
} else {
  # NIPTeR mode
  conversion_status$stage <- "nipter_bed"
  conversion_status$step <- "prepare_gc_correction"
  gc_table <- opts$`gc-table`
  gc_fasta <- opts$fasta
  qc_state$gc_correction <- list(
    applied = !is.null(gc_table) || !is.null(gc_fasta),
    method = if (!is.null(gc_table) || !is.null(gc_fasta)) "loess" else "none",
    loess_span = if (!is.null(gc_table) || !is.null(gc_fasta)) as.numeric(opts$`nipter-gc-loess-span`) else NA_real_,
    include_sex = if (!is.null(gc_table) || !is.null(gc_fasta)) isTRUE(opts$`nipter-gc-include-sex`) else FALSE,
    binsize = as.integer(binsize),
    table_bgz = if (!is.null(gc_table)) normalizePath(gc_table, winslash = "/", mustWork = TRUE) else NA_character_,
    fasta = if (!is.null(gc_fasta)) normalizePath(gc_fasta, winslash = "/", mustWork = TRUE) else NA_character_
  )
  bed_metadata <- .compact_metadata(list(
      format = "nipter_bed",
      schema = if (isTRUE(opts$`separate-strands`)) "separated_v1" else "combined_v1",
      binsize = as.integer(binsize),
      mapq = as.integer(mapq),
      require_flags = as.integer(opts$`require-flags`),
      exclude_flags = as.integer(opts$`exclude-flags`),
      rmdup = rmdup,
      corrected_columns = !is.null(gc_table) || !is.null(gc_fasta),
      reference = if (!is.null(opts$reference)) normalizePath(opts$reference, winslash = "/", mustWork = TRUE) else NULL
    ))

  cat(sprintf("Converting %s → %s (NIPTeR BED%s, binsize=%d, mapq=%d, rmdup=%s)\n",
              basename(bam), basename(out),
              if (opts$`separate-strands`) " 9-col" else " 5-col",
              binsize, mapq, rmdup))

  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  Rduckhts::rduckhts_load(con)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  conversion_status$step <- "read_counts"
  read_count_rows <- .bam_bin_count_rows(
    con = con,
    bam = bam,
    binsize = binsize,
    mapq = mapq,
    require_flags = as.integer(opts$`require-flags`),
    exclude_flags = as.integer(opts$`exclude-flags`),
    rmdup = rmdup,
    reference = opts$reference,
    stats = "gc,mq",
    include_unmapped = TRUE
  )
  qc_state$read_counts_metrics <- .read_counts_metric_record(
    record = .merge_tabix_metadata(
      .bam_alignment_count_metadata(con, bam),
      .bam_read_counts_metadata(read_count_rows)
    ),
    source = "samtools_idxstats+bam_bin_counts(gc,mq)"
  )
  chr_lengths <- .bam_chr_lengths(con, bam)
  raw_sample <- .nipter_rows_to_sample(
    rows = read_count_rows,
    chr_lengths = chr_lengths,
    binsize = binsize,
    name = sub("\\.cram$|\\.bam$", "", basename(bam), ignore.case = TRUE),
    separate_strands = isTRUE(opts$`separate-strands`)
  )
  corrected_sample <- NULL
  if (!is.null(gc_table) || !is.null(gc_fasta)) {
    conversion_status$step <- "gc_correction"
    cat("Applying GC correction ...\n")
    corrected_sample <- if (!is.null(gc_table)) {
      nipter_gc_correct(
        raw_sample,
        gc_table = gc_table,
        span = as.numeric(opts$`nipter-gc-loess-span`),
        include_sex = isTRUE(opts$`nipter-gc-include-sex`),
        con = con
      )
    } else {
      nipter_gc_correct(
        raw_sample,
        fasta = gc_fasta,
        span = as.numeric(opts$`nipter-gc-loess-span`),
        binsize = binsize,
        include_sex = isTRUE(opts$`nipter-gc-include-sex`),
        con = con
      )
    }
  }

  if (!is.null(corrected_sample)) {
    conversion_status$step <- "qc_metrics"
    nipter_qc <- .nipter_preprocess_qc(
      sample = raw_sample,
      corrected = corrected_sample,
      gc_table = if (!is.null(gc_table)) gc_table else NULL,
      fasta = if (is.null(gc_table)) gc_fasta else NULL,
      include_sex = isTRUE(opts$`nipter-gc-include-sex`),
      binsize = binsize,
      con = con
    )
    if (!is.null(opts$`gc-curve-plot-out`)) {
      conversion_status$step <- "qc_gc_curve_plot"
      .write_nipter_gc_curve_plot(
        curve_data = nipter_qc$curve_data,
        sample_name = sample_name,
        out_png = opts$`gc-curve-plot-out`,
        theme = opts$`qc-plot-theme`,
        base_size = as.integer(opts$`qc-plot-base-size`),
        loess_span = as.numeric(opts$`nipter-gc-loess-span`)
      )
      nipter_qc$metrics$gc_curve_plot <- opts$`gc-curve-plot-out`
    }
    if (!is.null(opts$`corrected-bin-ratio-genome-plot-out`)) {
      conversion_status$step <- "qc_corrected_bin_ratio_plot"
      .write_nipter_corrected_bin_ratio_genome_plot(
        curve_data = nipter_qc$curve_data,
        sample_name = sample_name,
        out_png = opts$`corrected-bin-ratio-genome-plot-out`,
        theme = opts$`qc-plot-theme`,
        base_size = as.integer(opts$`qc-plot-base-size`)
      )
      nipter_qc$metrics$corrected_bin_ratio_genome_plot <- opts$`corrected-bin-ratio-genome-plot-out`
    }
    qc_state$nipter_qc_metrics <- nipter_qc$metrics
  }

  if (!is.null(bed_metadata)) {
    bed_metadata <- .merge_tabix_metadata(
      bed_metadata,
      .sample_qc_tabix_metadata(
        seqff = qc_state$seqff_metrics,
        y_unique = qc_state$y_unique_metrics,
        read_counts = qc_state$read_counts_metrics,
        bam_stats = qc_state$bam_samtools_stats_metrics,
        nipter_qc = qc_state$nipter_qc_metrics,
        gc_correction = qc_state$gc_correction,
        filters_pre = metrics_pre_filters,
        filters_post = metrics_post_filters
      )
    )
  }

  conversion_status$step <- "write_bed"
  nipter_sample_to_bed(
    sample = raw_sample,
    bed = out,
    binsize = binsize,
    corrected = corrected_sample,
    con = con,
    metadata = bed_metadata
  )
}

conversion_status$ok <- TRUE
conversion_status$stage <- NA_character_
conversion_status$step <- NA_character_
conversion_status$message <- NA_character_

.write_sample_qc_on_exit()

cat("Done.\n")
  invisible(NULL)
}

tryCatch(
  main(),
  error = function(e) {
    conversion_status$message <- conditionMessage(e)
    if (!is.null(opts$`sample-qc-out`)) {
      try(.write_sample_qc_on_exit(), silent = TRUE)
    }
    stop(e)
  }
)
