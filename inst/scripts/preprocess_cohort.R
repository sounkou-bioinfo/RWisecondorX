#!/usr/bin/env Rscript
#
# preprocess_cohort.R
#
# Manifest-driven preprocessing for real BAM/CRAM cohorts. This script:
#   1. stages a BAM manifest under --out-root/manifests
#   2. precomputes a NIPTeR GC table once
#   3. optionally computes SeqFF fetal-fraction estimates from BAMs
#   4. computes Y-unique sex-predictor ratios for every BAM
#   5. converts every BAM to RWisecondorX BED.gz
#   6. optionally writes upstream WisecondorX NPZ files via the Python CLI
#   7. converts every BAM to NIPTeR BED.gz with separated strands and GC-corrected columns
#   8. writes native mosdepth-compatible coverage outputs for pre-filter and
#      post-filter NIPT coverage using an explicit coverage bin size
#
# It intentionally does not build references or score samples.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with install.packages('optparse').",
       call. = FALSE)
}
library(optparse)

.read_file_list <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) >= 1L && identical(lines[[1L]], "bam")) {
    warning(
      "Manifest ", path, " starts with a 'bam' header line; stripping it.",
      call. = FALSE
    )
    lines <- lines[-1L]
  }
  missing <- lines[!file.exists(lines)]
  if (length(missing) > 0L) {
    stop("Manifest contains missing BAM/CRAM paths:\n  ",
         paste(head(missing, 10L), collapse = "\n  "),
         if (length(missing) > 10L) {
           paste0("\n  ... and ", length(missing) - 10L, " more")
         },
         call. = FALSE)
  }
  lines
}

.write_manifest <- function(paths, out_path) {
  utils::write.table(
    data.frame(bam = paths),
    file = out_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

.compact_metadata <- function(x) {
  keep <- !vapply(x, is.null, logical(1L))
  x[keep]
}

.ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("Failed to create directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.normalize_legacy_wisecondorx_npz <- function(path) {
  legacy_path <- paste0(path, ".npz")
  if (!file.exists(path) && file.exists(legacy_path)) {
    ok <- file.rename(legacy_path, path)
    if (!isTRUE(ok) || !file.exists(path)) {
      stop(
        "Failed to rename legacy WisecondorX NPZ output:\n  ",
        legacy_path, "\n  -> ", path,
        call. = FALSE
      )
    }
  }
  invisible(path)
}

.run_one <- function(expr, label) {
  cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), label))
  force(expr)
  invisible(NULL)
}

.log_sample <- function(i, n, label, stem) {
  cat(sprintf("[%s] [%d/%d] %s %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              i, n, label, stem))
}

.run_stage_samples <- function(bams, jobs, worker) {
  jobs <- as.integer(jobs)
  stopifnot(length(jobs) == 1L, jobs >= 1L)
  if (length(bams) == 0L) {
    return(list())
  }
  res <- if (.Platform$OS.type == "unix" && jobs > 1L) {
    parallel::mclapply(
      X = seq_along(bams),
      FUN = function(i) worker(i, bams[[i]]),
      mc.cores = jobs,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_along(bams), function(i) worker(i, bams[[i]]))
  }

  bad <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(bad)) {
    msgs <- vapply(res[bad], as.character, character(1))
    stop(
      "Parallel stage failed:\n",
      paste(sprintf("[[%d]] %s", which(bad), msgs), collapse = "\n"),
      call. = FALSE
    )
  }

  res
}

.with_duckhts_con <- function(fun) {
  stopifnot(is.function(fun))
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  Rduckhts::rduckhts_load(con)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fun(con)
}

.per_worker_threads <- function(total_threads, jobs) {
  total_threads <- as.integer(total_threads)
  jobs <- as.integer(jobs)
  max(1L, total_threads %/% max(1L, jobs))
}

.mosdepth_thread_args <- function(total_threads, jobs) {
  worker_threads <- .per_worker_threads(total_threads, jobs)
  if (worker_threads <= 1L) {
    return(list(threads = 1L, processing_threads = 1L))
  }
  list(
    threads = worker_threads - 1L,
    processing_threads = 1L
  )
}

.SAM_FLAG_UNMAPPED <- 0x4L
.SAM_FLAG_SECONDARY <- 0x100L
.SAM_FLAG_QCFAIL <- 0x200L
.SAM_FLAG_DUPLICATE <- 0x400L
.MOSDEPTH_EXCLUDE_FLAGS_BASE <- bitwOr(
  .SAM_FLAG_UNMAPPED,
  bitwOr(
    .SAM_FLAG_SECONDARY,
    bitwOr(.SAM_FLAG_QCFAIL, .SAM_FLAG_DUPLICATE)
  )
)

option_list <- list(
  make_option("--bam-list", type = "character", default = NULL,
              help = "Text file with one BAM/CRAM path per line [required]"),
  make_option("--out-root", type = "character", default = NULL,
              help = "Working root for manifests/BEDs/NPZs/logs [required]"),
  make_option("--fasta", type = "character", default = NULL,
              help = "Indexed reference FASTA for GC precompute [required]"),
  make_option("--threads", type = "integer", default = 20L,
              help = "Reserved thread budget for downstream work [default: %default]"),
  make_option("--jobs", type = "integer", default = 4L,
              help = "Number of samples to process in parallel per stage [default: %default]"),
  make_option("--wcx-binsize", type = "integer", default = 100000L,
              help = "WisecondorX convert binsize [default: %default]"),
  make_option("--wcx-mapq", type = "integer", default = 1L,
              help = "WisecondorX BED MAPQ filter [default: %default]"),
  make_option("--wcx-require-flags", type = "integer", default = 0L,
              help = "WisecondorX BED required SAM flags [default: %default]"),
  make_option("--wcx-exclude-flags", type = "integer", default = 0L,
              help = "WisecondorX BED excluded SAM flags [default: %default]"),
  make_option("--wcx-rmdup", type = "character", default = "streaming",
              help = "WisecondorX BED duplicate mode: streaming, flag, none [default: %default]"),
  make_option("--nipter-binsize", type = "integer", default = 50000L,
              help = "NIPTeR binsize [default: %default]"),
  make_option("--nipter-mapq", type = "integer", default = 40L,
              help = "NIPTeR MAPQ filter [default: %default]"),
  make_option("--nipter-exclude-flags", type = "integer", default = 1024L,
              help = "NIPTeR exclude flags [default: %default]"),
  make_option("--nipter-gc-include-sex", action = "store", type = "logical", default = TRUE,
              help = "Also GC-correct NIPTeR X/Y bins before BED export. Pass --nipter-gc-include-sex=FALSE to disable [default: %default]"),
  make_option("--nipter-gc-loess-span", type = "double", default = 0.75,
              help = "LOESS span used for NIPTeR GC correction [default: %default]"),
  make_option("--nipter-gc-curves", action = "store", type = "logical", default = TRUE,
              help = "Write per-sample NIPTeR GC-curve PNGs. Pass --nipter-gc-curves=FALSE to disable [default: %default]"),
  make_option("--nipter-corrected-bin-ratio-genome-plots", action = "store", type = "logical", default = TRUE,
              help = paste0(
                "Write per-sample genome-order plots of post-GC corrected bin ratios ",
                "excluding out-of-range bins. Pass --nipter-corrected-bin-ratio-genome-plots=FALSE to disable [default: %default]"
              )),
  make_option("--nipter-separate-strands", action = "store", type = "logical", default = TRUE,
              help = "Write NIPTeR 9-column separated-strand BEDs. Pass --nipter-separate-strands=FALSE to disable [default: %default]"),
  make_option("--qc-plot-theme", type = "character", default = "minimal",
              help = "QC plot theme: minimal, light, bw, classic [default: %default]"),
  make_option("--qc-plot-base-size", type = "integer", default = 11L,
              help = "QC plot base font size [default: %default]"),
  make_option("--seqff", action = "store", type = "logical", default = TRUE,
              help = "Compute SeqFF fetal-fraction estimates from BAMs. Pass --seqff=FALSE to disable [default: %default]"),
  make_option("--seqff-tsv", type = "character", default = NULL,
              help = paste0("Optional TSV/CSV with per-sample SeqFF metrics. ",
                            "Must contain the explicit *_pre/*_post columns written by RWisecondorX.")),
  make_option("--y-unique-tsv", type = "character", default = NULL,
              help = paste0("Optional TSV/CSV with per-sample Y-unique metrics. ",
                            "Must contain the explicit *_pre/*_post columns written by RWisecondorX.")),
  make_option("--y-regions-file", type = "character", default = NULL,
              help = paste0("Optional Y-unique regions file used for on-the-fly Y-unique computation ",
                            "and written into BED metadata. Accepts headered TSV or headerless BED-like files.")),
  make_option("--metrics-pre-mapq", type = "integer", default = 0L,
              help = "MAPQ used for nonfiltered SeqFF/Y-unique metrics [default: %default]"),
  make_option("--metrics-pre-require-flags", type = "integer", default = 0L,
              help = "Required flags used for nonfiltered SeqFF/Y-unique metrics [default: %default]"),
  make_option("--metrics-pre-exclude-flags", type = "integer", default = 0L,
              help = "Excluded flags used for nonfiltered SeqFF/Y-unique metrics [default: %default]"),
  make_option("--metrics-post-mapq", type = "integer", default = NULL,
              help = "MAPQ used for filtered SeqFF/Y-unique metrics. Defaults to --nipter-mapq."),
  make_option("--metrics-post-require-flags", type = "integer", default = 0L,
              help = "Required flags used for filtered SeqFF/Y-unique metrics [default: %default]"),
  make_option("--metrics-post-exclude-flags", type = "integer", default = NULL,
              help = "Excluded flags used for filtered SeqFF/Y-unique metrics. Defaults to --nipter-exclude-flags."),
  make_option("--wcx-write-npz", action = "store", type = "logical", default = FALSE,
              help = "Also write upstream WisecondorX NPZ files via the Python CLI. Pass --wcx-write-npz=TRUE to enable [default: %default]"),
  make_option("--coverage-binsize", type = "integer", default = 50000L,
              help = "Bin size used for pre/post coverage summaries [default: %default]"),
  make_option("--mosdepth-exclude-flags-base", type = "integer", default = .MOSDEPTH_EXCLUDE_FLAGS_BASE,
              help = paste0(
                "SAM flags excluded in mosdepth coverage outputs. Default excludes unmapped (0x4), ",
                "secondary (0x100), QC-fail (0x200), and duplicate (0x400) reads [default: %default]"
              )),
  make_option("--mosdepth-filtered-mapq", type = "integer", default = NULL,
              help = "MAPQ used for the filtered mosdepth output. Defaults to --nipter-mapq."),
  make_option("--mosdepth-filtered-exclude-flags", type = "integer", default = NULL,
              help = "Excluded SAM flags used for the filtered mosdepth output. Defaults to --mosdepth-exclude-flags-base."),
  make_option("--overwrite", action = "store", type = "logical", default = FALSE,
              help = "Overwrite existing outputs. Pass --overwrite=TRUE to enable [default: %default]")
)

parser <- OptionParser(
  usage = "%prog --bam-list cohort.txt --out-root /mnt/data/niptseq/testdata --fasta ref.fa [options]",
  option_list = option_list
)
opts <- parse_args(parser)

if (is.null(opts$`bam-list`) || !file.exists(opts$`bam-list`)) {
  stop("--bam-list is required and must exist.", call. = FALSE)
}
if (is.null(opts$`out-root`) || !nzchar(opts$`out-root`)) {
  stop("--out-root is required.", call. = FALSE)
}
if (is.null(opts$fasta) || !file.exists(opts$fasta)) {
  stop("--fasta is required and must exist.", call. = FALSE)
}
if (is.na(opts$jobs) || as.integer(opts$jobs) < 1L) {
  stop("--jobs must be >= 1.", call. = FALSE)
}
if (is.na(opts$seqff)) {
  stop("--seqff must be TRUE or FALSE.", call. = FALSE)
}
if (is.na(opts$`nipter-gc-include-sex`)) {
  stop("--nipter-gc-include-sex must be TRUE or FALSE.", call. = FALSE)
}
if (is.na(opts$`nipter-gc-curves`)) {
  stop("--nipter-gc-curves must be TRUE or FALSE.", call. = FALSE)
}
if (is.na(opts$`nipter-corrected-bin-ratio-genome-plots`)) {
  stop("--nipter-corrected-bin-ratio-genome-plots must be TRUE or FALSE.", call. = FALSE)
}
if (is.na(opts$`nipter-separate-strands`)) {
  stop("--nipter-separate-strands must be TRUE or FALSE.", call. = FALSE)
}
if (is.na(opts$`wcx-write-npz`)) {
  stop("--wcx-write-npz must be TRUE or FALSE.", call. = FALSE)
}
if (is.na(opts$overwrite)) {
  stop("--overwrite must be TRUE or FALSE.", call. = FALSE)
}
if (!is.numeric(opts$`nipter-gc-loess-span`) ||
    length(opts$`nipter-gc-loess-span`) != 1L ||
    !is.finite(opts$`nipter-gc-loess-span`) ||
    opts$`nipter-gc-loess-span` <= 0 ||
    opts$`nipter-gc-loess-span` > 1) {
  stop("--nipter-gc-loess-span must be a finite number in (0, 1].", call. = FALSE)
}
if (is.na(opts$`coverage-binsize`) || as.integer(opts$`coverage-binsize`) < 1L) {
  stop("--coverage-binsize must be >= 1.", call. = FALSE)
}
opts$`qc-plot-theme` <- match.arg(
  tolower(opts$`qc-plot-theme`),
  c("minimal", "light", "bw", "classic")
)
opts$`wcx-rmdup` <- match.arg(
  tolower(opts$`wcx-rmdup`),
  c("streaming", "flag", "none")
)
for (path_opt in c("seqff-tsv", "y-unique-tsv", "y-regions-file")) {
  path_val <- opts[[path_opt]]
  if (!is.null(path_val) && !file.exists(path_val)) {
    stop("--", path_opt, " does not exist: ", path_val, call. = FALSE)
  }
}

bams <- .read_file_list(opts$`bam-list`)

out_root <- .ensure_dir(opts$`out-root`)
dirs <- list(
  manifests = .ensure_dir(file.path(out_root, "manifests")),
  gc = .ensure_dir(file.path(out_root, "gc")),
  seqff = .ensure_dir(file.path(out_root, "seqff")),
  y_unique = .ensure_dir(file.path(out_root, "y_unique")),
  sample_qc = .ensure_dir(file.path(out_root, "sample_qc")),
  rwcx_beds = .ensure_dir(file.path(out_root, "rwcx_beds")),
  wisecondorx_npz = .ensure_dir(file.path(out_root, "wisecondorx_npz")),
  nipter_beds = .ensure_dir(file.path(out_root, "nipter_beds")),
  nipter_gc_curves = .ensure_dir(file.path(out_root, "nipter_gc_curves")),
  nipter_gc_curve_data = .ensure_dir(file.path(out_root, "nipter_gc_curve_data")),
  nipter_corrected_bin_ratio_genome_plots = .ensure_dir(file.path(out_root, "nipter_corrected_bin_ratio_genome_plots")),
  coverage_pre = .ensure_dir(file.path(out_root, "coverage_pre")),
  coverage_post = .ensure_dir(file.path(out_root, "coverage_post")),
  logs = .ensure_dir(file.path(out_root, "logs"))
)

staged_manifest <- file.path(dirs$manifests, basename(opts$`bam-list`))
if (!file.exists(staged_manifest) || isTRUE(opts$overwrite)) {
  .write_manifest(bams, staged_manifest)
}

gc_table <- file.path(
  dirs$gc,
  sprintf("gc_%s.tsv.bgz", as.integer(opts$`nipter-binsize`))
)

library(RWisecondorX)

.match_metric_row <- getFromNamespace(".match_metric_row", "RWisecondorX")
.seqff_metric_record <- getFromNamespace(".seqff_metric_record", "RWisecondorX")
.y_unique_metric_record <- getFromNamespace(".y_unique_metric_record", "RWisecondorX")
.sample_qc_row <- getFromNamespace(".sample_qc_row", "RWisecondorX")
.sample_qc_row_with_processing_status <- getFromNamespace(".sample_qc_row_with_processing_status", "RWisecondorX")
.sample_qc_row_with_nipter_status <- getFromNamespace(".sample_qc_row_with_nipter_status", "RWisecondorX")
.sample_qc_tabix_metadata <- getFromNamespace(".sample_qc_tabix_metadata", "RWisecondorX")
.samtools_stats_metric_record <- getFromNamespace(".samtools_stats_metric_record", "RWisecondorX")
.seqff_summary_schema <- getFromNamespace(".seqff_summary_schema", "RWisecondorX")
.y_unique_summary_schema <- getFromNamespace(".y_unique_summary_schema", "RWisecondorX")
.sample_qc_schema <- getFromNamespace(".sample_qc_schema", "RWisecondorX")
.read_counts_metric_record <- getFromNamespace(".read_counts_metric_record", "RWisecondorX")
.coverage_metric_record <- getFromNamespace(".coverage_metric_record", "RWisecondorX")
.nipter_preprocess_qc_metric_record <- getFromNamespace(".nipter_preprocess_qc_metric_record", "RWisecondorX")
.gc_correction_record <- getFromNamespace(".gc_correction_record", "RWisecondorX")
.metric_filter_record <- getFromNamespace(".metric_filter_record", "RWisecondorX")
.bam_bin_count_rows <- getFromNamespace(".bam_bin_count_rows", "RWisecondorX")
.bam_chr_lengths <- getFromNamespace(".bam_chr_lengths", "RWisecondorX")
.nipter_rows_to_sample <- getFromNamespace(".nipter_rows_to_sample", "RWisecondorX")
.merge_tabix_metadata <- getFromNamespace(".merge_tabix_metadata", "RWisecondorX")
.bam_alignment_count_metadata <- getFromNamespace(".bam_alignment_count_metadata", "RWisecondorX")
.bam_read_counts_metadata <- getFromNamespace(".bam_read_counts_metadata", "RWisecondorX")
.nipter_preprocess_qc <- getFromNamespace(".nipter_preprocess_qc", "RWisecondorX")
.write_nipter_gc_curve_plot <- getFromNamespace(".write_nipter_gc_curve_plot", "RWisecondorX")
.write_nipter_corrected_bin_ratio_genome_plot <- getFromNamespace(".write_nipter_corrected_bin_ratio_genome_plot", "RWisecondorX")
.write_nipter_gc_curve_data_bgz <- getFromNamespace(".write_nipter_gc_curve_data_bgz", "RWisecondorX")
.tabix_metadata_frame <- getFromNamespace(".tabix_metadata_frame", "RWisecondorX")
.ensure_wisecondorx_env <- getFromNamespace(".ensure_wisecondorx_env", "RWisecondorX")

.nipter_gc_correction_record <- function() {
  .gc_correction_record(
    applied = TRUE,
    method = "loess",
    loess_span = as.numeric(opts$`nipter-gc-loess-span`),
    include_sex = isTRUE(opts$`nipter-gc-include-sex`),
    binsize = as.integer(opts$`nipter-binsize`),
    table_bgz = normalizePath(gc_table, winslash = "/", mustWork = TRUE),
    fasta = normalizePath(opts$fasta, winslash = "/", mustWork = TRUE)
  )
}

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

y_regions_file <- .resolve_y_regions_file(opts$`y-regions-file`)

.write_bed_header_qc_summary <- function(bed_dir, out_tsv) {
  bed_files <- sort(Sys.glob(file.path(bed_dir, "*.bed.gz")))
  header_df <- if (length(bed_files)) {
    .with_duckhts_con(function(con) {
      .tabix_metadata_frame(bed_files, con = con)
    })
  } else {
    data.frame(
      sample_name = character(0),
      artifact_path = character(0),
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(
    header_df,
    file = out_tsv,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  invisible(out_tsv)
}

metrics_pre_filters <- list(
  mapq = as.integer(opts$`metrics-pre-mapq`),
  require_flags = as.integer(opts$`metrics-pre-require-flags`),
  exclude_flags = as.integer(opts$`metrics-pre-exclude-flags`)
)
metrics_post_filters <- list(
  mapq = as.integer(if (is.null(opts$`metrics-post-mapq`)) opts$`nipter-mapq` else opts$`metrics-post-mapq`),
  require_flags = as.integer(opts$`metrics-post-require-flags`),
  exclude_flags = as.integer(if (is.null(opts$`metrics-post-exclude-flags`)) opts$`nipter-exclude-flags` else opts$`metrics-post-exclude-flags`)
)
coverage_binsize <- as.integer(opts$`coverage-binsize`)
coverage_pre_filters <- list(
  mapq = 0L,
  exclude_flags = as.integer(opts$`mosdepth-exclude-flags-base`)
)
coverage_post_filters <- list(
  mapq = as.integer(if (is.null(opts$`mosdepth-filtered-mapq`)) {
    opts$`nipter-mapq`
  } else {
    opts$`mosdepth-filtered-mapq`
  }),
  exclude_flags = as.integer(if (is.null(opts$`mosdepth-filtered-exclude-flags`)) {
    opts$`mosdepth-exclude-flags-base`
  } else {
    opts$`mosdepth-filtered-exclude-flags`
  })
)

.seqff_summary_cols <- .seqff_summary_schema()
.y_unique_summary_cols <- .y_unique_summary_schema()
.sample_qc_summary_cols <- .sample_qc_schema()

.sample_stem <- function(path) {
  sub("\\.(bam|cram)$", "", basename(path), ignore.case = TRUE)
}

.stage_success <- function(sample_name, ...) {
  c(list(ok = TRUE, sample_name = sample_name), list(...))
}

.stage_failure <- function(sample_name, stage, step, message) {
  list(
    ok = FALSE,
    sample_name = sample_name,
    failure_stage = stage,
    failure_step = step,
    failure_message = message
  )
}

sample_names <- vapply(bams, .sample_stem, character(1L))
names(bams) <- sample_names

seqff_metrics_by_sample <- setNames(vector("list", length(sample_names)), sample_names)
y_unique_metrics_by_sample <- setNames(vector("list", length(sample_names)), sample_names)
read_counts_by_sample <- setNames(vector("list", length(sample_names)), sample_names)
bam_stats_by_sample <- setNames(vector("list", length(sample_names)), sample_names)
coverage_metrics_by_sample <- setNames(vector("list", length(sample_names)), sample_names)
sample_failures_by_sample <- setNames(vector("list", length(sample_names)), sample_names)
final_sample_qc_rows <- setNames(vector("list", length(sample_names)), sample_names)

.record_sample_failure <- function(item, i = NULL, n = NULL) {
  stopifnot(is.list(item), !is.null(item$sample_name), !is.null(item$failure_stage))
  stem <- as.character(item$sample_name)[[1L]]
  if (is.null(sample_failures_by_sample[[stem]])) {
    sample_failures_by_sample[[stem]] <<- list(
      stage = as.character(item$failure_stage)[[1L]],
      step = if (is.null(item$failure_step)) NA_character_ else as.character(item$failure_step)[[1L]],
      message = as.character(item$failure_message)[[1L]]
    )
  }
  if (!is.null(i) && !is.null(n)) {
    cat(sprintf(
      "[%s] [%d/%d] %s failed %s: %s\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      as.integer(i), as.integer(n),
      as.character(item$failure_stage)[[1L]],
      stem,
      as.character(item$failure_message)[[1L]]
    ))
  }
  invisible(NULL)
}

.active_bams <- function() {
  bams[vapply(names(bams), function(stem) is.null(sample_failures_by_sample[[stem]]), logical(1L))]
}

.build_sample_qc_row <- function(stem, bam) {
  row <- .sample_qc_row(
    sample_name = stem,
    bam = bam,
    seqff = seqff_metrics_by_sample[[stem]],
    y_unique = y_unique_metrics_by_sample[[stem]],
    read_counts = read_counts_by_sample[[stem]],
    bam_stats = bam_stats_by_sample[[stem]],
    coverage = coverage_metrics_by_sample[[stem]],
    filters_pre = metrics_pre_filters,
    filters_post = metrics_post_filters
  )

  failure <- sample_failures_by_sample[[stem]]
  if (is.null(failure)) {
    row <- .sample_qc_row_with_processing_status(row, status = "ok")
    row <- .sample_qc_row_with_nipter_status(row, status = "not_run")
  } else {
    row <- .sample_qc_row_with_processing_status(
      row,
      status = "failed",
      failure_stage = failure$stage,
      failure_step = failure$step,
      failure_message = failure$message
    )
    row <- .sample_qc_row_with_nipter_status(
      row,
      status = if (identical(failure$stage, "nipter_bed")) "failed" else "not_run",
      error = if (identical(failure$stage, "nipter_bed")) failure$message else NULL
    )
  }

  row[, .sample_qc_summary_cols, drop = FALSE]
}

.run_one({
  nipter_gc_precompute(
    fasta = opts$fasta,
    out = gc_table,
    binsize = as.integer(opts$`nipter-binsize`)
  )
}, sprintf("Precompute GC table → %s", gc_table))

if (isTRUE(opts$seqff) || !is.null(opts$`seqff-tsv`)) {
  .run_one({
    stage_bams <- .active_bams()
    if (!length(stage_bams)) {
      return(invisible(NULL))
    }
    seqff_rows <- .run_stage_samples(stage_bams, opts$jobs, function(i, bam) {
      stem <- .sample_stem(bam)
      tryCatch({
        out_tsv <- file.path(dirs$seqff, paste0(stem, ".seqff.tsv"))
        if (file.exists(out_tsv) && !isTRUE(opts$overwrite)) {
          row <- utils::read.delim(
            out_tsv,
            sep = "\t",
            header = TRUE,
            check.names = FALSE
          )
          metrics <- .seqff_metric_record(
            record = row,
            source = if (!is.null(opts$`seqff-tsv`)) {
              paste0("file:", basename(opts$`seqff-tsv`))
            } else {
              "computed"
            }
          )
          return(.stage_success(stem, row = row, metrics = metrics))
        }
        .log_sample(i, length(stage_bams), "SeqFF", stem)
        seqff_metrics <- if (!is.null(opts$`seqff-tsv`)) {
          .seqff_metric_record(
            record = .match_metric_row(opts$`seqff-tsv`, stem),
            source = paste0("file:", basename(opts$`seqff-tsv`))
          )
        } else {
          .seqff_metric_record(
            pre = seqff_predict(
              input = bam,
              input_type = "bam",
              mapq = metrics_pre_filters$mapq,
              require_flags = metrics_pre_filters$require_flags,
              exclude_flags = metrics_pre_filters$exclude_flags,
              reference = opts$fasta
            ),
            post = seqff_predict(
              input = bam,
              input_type = "bam",
              mapq = metrics_post_filters$mapq,
              require_flags = metrics_post_filters$require_flags,
              exclude_flags = metrics_post_filters$exclude_flags,
              reference = opts$fasta
            ),
            source = "computed"
          )
        }
        row <- .sample_qc_row(
          sample_name = stem,
          bam = bam,
          seqff = seqff_metrics,
          filters_pre = metrics_pre_filters,
          filters_post = metrics_post_filters
        )[, .seqff_summary_cols, drop = FALSE]
        utils::write.table(
          row,
          file = out_tsv,
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE,
          quote = FALSE
        )
        .stage_success(stem, row = row, metrics = seqff_metrics)
      }, error = function(e) {
        .stage_failure(stem, "seqff", "compute_or_load", conditionMessage(e))
      })
    })

    seqff_summary_rows <- list()
    for (item in seqff_rows) {
      if (isTRUE(item$ok)) {
        seqff_metrics_by_sample[[item$sample_name]] <- item$metrics
        seqff_summary_rows[[length(seqff_summary_rows) + 1L]] <- item$row
      } else {
        .record_sample_failure(item)
      }
    }
    seqff_summary <- if (length(seqff_summary_rows)) {
      do.call(rbind, seqff_summary_rows)
    } else {
      data.frame(matrix(nrow = 0L, ncol = length(.seqff_summary_cols)))
    }
    names(seqff_summary) <- .seqff_summary_cols
    utils::write.table(
      seqff_summary,
      file = file.path(dirs$seqff, "seqff_summary.tsv"),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }, sprintf("Compute SeqFF fetal-fraction estimates for %d BAMs", length(.active_bams())))
}

 .run_one({
   stage_bams <- .active_bams()
   if (!length(stage_bams)) {
     return(invisible(NULL))
   }
   y_unique_rows <- .run_stage_samples(stage_bams, opts$jobs, function(i, bam) {
     stem <- .sample_stem(bam)
     tryCatch({
       out_tsv <- file.path(dirs$y_unique, paste0(stem, ".y_unique.tsv"))
       if (file.exists(out_tsv) && !isTRUE(opts$overwrite)) {
         row <- utils::read.delim(
           out_tsv,
           sep = "\t",
           header = TRUE,
           check.names = FALSE
         )
         metrics <- .y_unique_metric_record(
           record = row,
           regions_file = y_regions_file,
           source = if (!is.null(opts$`y-unique-tsv`)) {
             paste0("file:", basename(opts$`y-unique-tsv`))
           } else {
             "computed"
           }
         )
         return(.stage_success(stem, row = row, metrics = metrics))
       }
       .log_sample(i, length(stage_bams), "Y-unique ratio", stem)
       y_unique_metrics <- if (!is.null(opts$`y-unique-tsv`)) {
         .y_unique_metric_record(
           record = .match_metric_row(opts$`y-unique-tsv`, stem),
           regions_file = y_regions_file,
           source = paste0("file:", basename(opts$`y-unique-tsv`))
         )
       } else {
         .y_unique_metric_record(
           pre = nipter_y_unique_ratio(
             bam = bam,
             mapq = metrics_pre_filters$mapq,
             require_flags = metrics_pre_filters$require_flags,
             exclude_flags = metrics_pre_filters$exclude_flags,
             regions_file = y_regions_file,
             reference = opts$fasta
           ),
           post = nipter_y_unique_ratio(
             bam = bam,
             mapq = metrics_post_filters$mapq,
             require_flags = metrics_post_filters$require_flags,
             exclude_flags = metrics_post_filters$exclude_flags,
             regions_file = y_regions_file,
             reference = opts$fasta
           ),
           regions_file = y_regions_file,
           source = "computed"
         )
       }
       row <- .sample_qc_row(
         sample_name = stem,
         bam = bam,
         y_unique = y_unique_metrics,
         filters_pre = metrics_pre_filters,
         filters_post = metrics_post_filters
       )[, .y_unique_summary_cols, drop = FALSE]
       utils::write.table(
         row,
         file = out_tsv,
         sep = "\t",
         row.names = FALSE,
         col.names = TRUE,
         quote = FALSE
       )
       .stage_success(stem, row = row, metrics = y_unique_metrics)
     }, error = function(e) {
       .stage_failure(stem, "y_unique", "compute_or_load", conditionMessage(e))
     })
   })
   y_unique_summary_rows <- list()
   for (item in y_unique_rows) {
     if (isTRUE(item$ok)) {
       y_unique_metrics_by_sample[[item$sample_name]] <- item$metrics
       y_unique_summary_rows[[length(y_unique_summary_rows) + 1L]] <- item$row
     } else {
       .record_sample_failure(item)
     }
   }
   y_unique_summary <- if (length(y_unique_summary_rows)) {
     do.call(rbind, y_unique_summary_rows)
   } else {
     data.frame(matrix(nrow = 0L, ncol = length(.y_unique_summary_cols)))
   }
   names(y_unique_summary) <- .y_unique_summary_cols
   utils::write.table(
     y_unique_summary,
     file = file.path(dirs$y_unique, "y_unique_summary.tsv"),
     sep = "\t",
     row.names = FALSE,
     col.names = TRUE,
     quote = FALSE
   )
 }, sprintf("Compute Y-unique sex-predictor ratios for %d BAMs", length(.active_bams())))

.run_one({
  stage_bams <- .active_bams()
  if (!length(stage_bams)) {
    return(invisible(NULL))
  }
  mosdepth_threads <- .mosdepth_thread_args(opts$threads, opts$jobs)
  coverage_rows <- .run_stage_samples(stage_bams, opts$jobs, function(i, bam) {
    stem <- .sample_stem(bam)
    tryCatch({
      pre_prefix <- file.path(dirs$coverage_pre, stem)
      pre_summary <- paste0(pre_prefix, ".mosdepth.summary.txt")
      if (!file.exists(pre_summary) || isTRUE(opts$overwrite)) {
        .log_sample(i, length(stage_bams), "Coverage pre", stem)
        .with_duckhts_con(function(con) {
          Rduckhts::rduckhts_mosdepth(
            con = con,
            prefix = pre_prefix,
            path = bam,
            by = as.character(coverage_binsize),
            fasta = opts$fasta,
            no_per_base = TRUE,
            threads = mosdepth_threads$threads,
            processing_threads = mosdepth_threads$processing_threads,
            flag = coverage_pre_filters$exclude_flags,
            include_flag = 0L,
            fast_mode = TRUE,
            mapq = coverage_pre_filters$mapq,
            overwrite = isTRUE(opts$overwrite)
          )
        })
      }

      post_prefix <- file.path(dirs$coverage_post, stem)
      post_summary <- paste0(post_prefix, ".mosdepth.summary.txt")
      if (!file.exists(post_summary) || isTRUE(opts$overwrite)) {
        .log_sample(i, length(stage_bams), "Coverage post", stem)
        .with_duckhts_con(function(con) {
          Rduckhts::rduckhts_mosdepth(
            con = con,
            prefix = post_prefix,
            path = bam,
            by = as.character(coverage_binsize),
            fasta = opts$fasta,
            no_per_base = TRUE,
            threads = mosdepth_threads$threads,
            processing_threads = mosdepth_threads$processing_threads,
            flag = coverage_post_filters$exclude_flags,
            include_flag = 0L,
            fast_mode = TRUE,
            mapq = coverage_post_filters$mapq,
            overwrite = isTRUE(opts$overwrite)
          )
        })
      }

      .stage_success(
        stem,
        coverage_metrics = .coverage_metric_record(
          pre_summary = pre_summary,
          post_summary = post_summary,
          binsize = coverage_binsize,
          pre_mapq = coverage_pre_filters$mapq,
          pre_exclude_flags = coverage_pre_filters$exclude_flags,
          post_mapq = coverage_post_filters$mapq,
          post_exclude_flags = coverage_post_filters$exclude_flags,
          source = "native_mosdepth_dense"
        )
      )
    }, error = function(e) {
      .stage_failure(stem, "coverage", "pre_or_post_mosdepth", conditionMessage(e))
    })
  })
  for (item in coverage_rows) {
    if (isTRUE(item$ok)) {
      coverage_metrics_by_sample[[item$sample_name]] <- item$coverage_metrics
    } else {
      .record_sample_failure(item)
    }
  }
}, sprintf("Write pre/post coverage outputs for %d BAMs", length(.active_bams())))

.run_one({
  stage_bams <- .active_bams()
  if (!length(stage_bams)) {
    return(invisible(NULL))
  }
  bam_stats_rows <- .run_stage_samples(stage_bams, opts$jobs, function(i, bam) {
    stem <- .sample_stem(bam)
    tryCatch({
      .log_sample(i, length(stage_bams), "samtools-stats summary", stem)
      .stage_success(
        stem,
        bam_stats = .samtools_stats_metric_record(
          pre = bam_samtools_stats_summary(
            path = bam,
            reference = opts$fasta,
            min_mapq = metrics_pre_filters$mapq,
            require_flags = metrics_pre_filters$require_flags,
            exclude_flags = metrics_pre_filters$exclude_flags
          ),
          post = bam_samtools_stats_summary(
            path = bam,
            reference = opts$fasta,
            min_mapq = metrics_post_filters$mapq,
            require_flags = metrics_post_filters$require_flags,
            exclude_flags = metrics_post_filters$exclude_flags,
            report_filtered_stream_bookkeeping = TRUE
          ),
          source = "native_htslib_samtools_stats"
        )
      )
    }, error = function(e) {
      .stage_failure(stem, "samtools_stats", "summary", conditionMessage(e))
    })
  })
  for (item in bam_stats_rows) {
    if (isTRUE(item$ok)) {
      bam_stats_by_sample[[item$sample_name]] <- item$bam_stats
    } else {
      .record_sample_failure(item)
    }
  }
}, sprintf("Compute pre/post samtools-stats summaries for %d BAMs", length(.active_bams())))

.run_one({
  stage_bams <- .active_bams()
  if (!length(stage_bams)) {
    return(invisible(NULL))
  }
  rwcx_rows <- .run_stage_samples(stage_bams, opts$jobs, function(i, bam) {
    stem <- .sample_stem(bam)
    out_bed <- file.path(dirs$rwcx_beds, paste0(stem, ".bed.gz"))
    tryCatch({
      if (file.exists(out_bed) && !isTRUE(opts$overwrite)) {
        return(.stage_success(stem))
      }
      .log_sample(i, length(stage_bams), "RWisecondorX BED", stem)
      sample_qc_md <- .sample_qc_tabix_metadata(
        seqff = seqff_metrics_by_sample[[stem]],
        y_unique = y_unique_metrics_by_sample[[stem]],
        read_counts = read_counts_by_sample[[stem]],
        bam_stats = bam_stats_by_sample[[stem]],
        coverage = coverage_metrics_by_sample[[stem]],
        filters_pre = metrics_pre_filters,
        filters_post = metrics_post_filters
      )
      bam_convert_bed(
        bam = bam,
        bed = out_bed,
        binsize = as.integer(opts$`wcx-binsize`),
        mapq = as.integer(opts$`wcx-mapq`),
        require_flags = as.integer(opts$`wcx-require-flags`),
        exclude_flags = as.integer(opts$`wcx-exclude-flags`),
        rmdup = opts$`wcx-rmdup`,
        reference = opts$fasta,
        metadata = .merge_tabix_metadata(
          .compact_metadata(list(
            format = "rwisecondorx_bed",
            schema = "count_v1",
            binsize = as.integer(opts$`wcx-binsize`),
            mapq = as.integer(opts$`wcx-mapq`),
            require_flags = as.integer(opts$`wcx-require-flags`),
            exclude_flags = as.integer(opts$`wcx-exclude-flags`),
            rmdup = opts$`wcx-rmdup`,
            reference = normalizePath(opts$fasta, winslash = "/", mustWork = TRUE)
          )),
          sample_qc_md
        )
      )
      .stage_success(stem)
    }, error = function(e) {
      partial_paths <- c(out_bed, paste0(out_bed, ".tbi"))
      partial_paths <- partial_paths[file.exists(partial_paths)]
      if (length(partial_paths)) {
        unlink(partial_paths, force = TRUE)
      }
      .stage_failure(stem, "rwcx_bed", "convert", conditionMessage(e))
    })
  })
  for (item in rwcx_rows) {
    if (!isTRUE(item$ok)) {
      .record_sample_failure(item)
    }
  }
}, sprintf("Convert %d BAMs to RWisecondorX BED.gz", length(.active_bams())))

.run_one({
  .write_bed_header_qc_summary(
    bed_dir = dirs$rwcx_beds,
    out_tsv = file.path(dirs$sample_qc, "rwcx_bed_header_qc.tsv")
  )
}, sprintf("Summarize RWisecondorX BED headers for %d samples", length(bams)))

if (isTRUE(opts$`wcx-write-npz`)) {
  .run_one({
    stage_bams <- .active_bams()
    if (!length(stage_bams)) {
      return(invisible(NULL))
    }
    .ensure_wisecondorx_env("wisecondorx")
    npz_rows <- .run_stage_samples(stage_bams, opts$jobs, function(i, bam) {
      stem <- .sample_stem(bam)
      out_npz <- file.path(dirs$wisecondorx_npz, paste0(stem, ".npz"))
      tryCatch({
        .normalize_legacy_wisecondorx_npz(out_npz)
        if (file.exists(out_npz) && !isTRUE(opts$overwrite)) {
          return(.stage_success(stem))
        }
        .log_sample(i, length(stage_bams), "WisecondorX CLI NPZ", stem)
        wisecondorx_convert(
          bam = bam,
          npz = out_npz,
          binsize = as.integer(opts$`wcx-binsize`),
          reference = opts$fasta,
          normdup = FALSE
        )
        .stage_success(stem)
      }, error = function(e) {
        partial_paths <- c(out_npz, paste0(out_npz, ".npz"))
        partial_paths <- partial_paths[file.exists(partial_paths)]
        if (length(partial_paths)) {
          unlink(partial_paths, force = TRUE)
        }
        .stage_failure(stem, "wisecondorx_npz", "convert", conditionMessage(e))
      })
    })
    for (item in npz_rows) {
      if (!isTRUE(item$ok)) {
        .record_sample_failure(item)
      }
    }
  }, sprintf("Convert %d BAMs to upstream WisecondorX NPZ via Python CLI", length(.active_bams())))
}

.run_one({
  stage_bams <- .active_bams()
  sample_qc_rows <- .run_stage_samples(stage_bams, opts$jobs, function(i, bam) {
    stem <- .sample_stem(bam)
    out_bed <- file.path(dirs$nipter_beds, paste0(stem, ".bed.gz"))
    out_qc <- file.path(dirs$sample_qc, paste0(stem, ".sample_qc.tsv"))
    if (file.exists(out_bed) && file.exists(out_qc) && !isTRUE(opts$overwrite)) {
      cached_row <- utils::read.delim(
        out_qc,
        sep = "\t",
        header = TRUE,
        check.names = FALSE
      )
      read_counts_metrics <- .read_counts_metric_record(
        record = cached_row,
        source = "samtools_idxstats+bam_bin_counts(gc,mq)"
      )
      coverage_metrics <- coverage_metrics_by_sample[[stem]]
      if (is.null(coverage_metrics) || !length(coverage_metrics)) {
        coverage_metrics <- .coverage_metric_record(record = cached_row)
      }
      row <- tryCatch(
        .sample_qc_row_with_nipter_status(
          cached_row,
          written = file.exists(out_bed),
          status = "ok"
        ),
        error = function(e) NULL
      )
      if (is.null(row) || !all(.sample_qc_summary_cols %in% names(row))) {
        row <- .sample_qc_row(
          sample_name = stem,
          bam = bam,
          seqff = seqff_metrics_by_sample[[stem]],
          y_unique = y_unique_metrics_by_sample[[stem]],
          read_counts = read_counts_metrics,
          bam_stats = bam_stats_by_sample[[stem]],
          coverage = coverage_metrics,
          nipter_qc = .nipter_preprocess_qc_metric_record(record = cached_row),
          gc_correction = .gc_correction_record(record = cached_row),
          filters_pre = metrics_pre_filters,
          filters_post = metrics_post_filters
        )
        row <- .sample_qc_row_with_processing_status(
          row,
          status = "ok"
        )
        row <- .sample_qc_row_with_nipter_status(
          row,
          written = file.exists(out_bed),
          status = "ok"
        )
        utils::write.table(
          row[, .sample_qc_summary_cols, drop = FALSE],
          file = out_qc,
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE,
          quote = FALSE
        )
      }
      row <- row[, .sample_qc_summary_cols, drop = FALSE]
      return(list(row = row, read_counts_metrics = read_counts_metrics))
    }
    .log_sample(i, length(stage_bams), "NIPTeR BED", stem)
    nipter_failure_step <- "read_counts"
    result <- tryCatch(
      .with_duckhts_con(function(con) {
        read_count_rows <- .bam_bin_count_rows(
          con = con,
          bam = bam,
          binsize = as.integer(opts$`nipter-binsize`),
          mapq = as.integer(opts$`nipter-mapq`),
          require_flags = 0L,
          exclude_flags = as.integer(opts$`nipter-exclude-flags`),
          rmdup = "none",
          reference = opts$fasta,
          stats = "gc,mq",
          include_unmapped = TRUE
        )
        read_counts_metrics <- .read_counts_metric_record(
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
          binsize = as.integer(opts$`nipter-binsize`),
          name = stem,
          separate_strands = isTRUE(opts$`nipter-separate-strands`)
        )
        nipter_failure_step <<- "gc_correction"
        corrected_sample <- nipter_gc_correct(
          raw_sample,
          gc_table = gc_table,
          span = as.numeric(opts$`nipter-gc-loess-span`),
          include_sex = isTRUE(opts$`nipter-gc-include-sex`),
          con = con
        )
        nipter_failure_step <<- "qc_metrics"
        nipter_qc <- .nipter_preprocess_qc(
          sample = raw_sample,
          corrected = corrected_sample,
          gc_table = gc_table,
          include_sex = isTRUE(opts$`nipter-gc-include-sex`),
          binsize = as.integer(opts$`nipter-binsize`),
          con = con
        )
        curve_bgz <- file.path(dirs$nipter_gc_curve_data, paste0(stem, ".gc_curve_data.tsv.bgz"))
        .write_nipter_gc_curve_data_bgz(
          curve_data = nipter_qc$curve_data,
          out_bgz = curve_bgz,
          con = con,
          metadata = .merge_tabix_metadata(
            list(
              artifact_type = "nipter_gc_curve_data",
              sample_name = stem
            ),
            .sample_qc_tabix_metadata(
              seqff = seqff_metrics_by_sample[[stem]],
              y_unique = y_unique_metrics_by_sample[[stem]],
              read_counts = read_counts_metrics,
              bam_stats = bam_stats_by_sample[[stem]],
              coverage = coverage_metrics_by_sample[[stem]],
              nipter_qc = nipter_qc$metrics,
              gc_correction = .nipter_gc_correction_record(),
              filters_pre = metrics_pre_filters,
              filters_post = metrics_post_filters
            )
          )
        )
        nipter_qc$metrics$gc_curve_data_bgz <- curve_bgz
        if (isTRUE(opts$`nipter-gc-curves`)) {
          nipter_failure_step <<- "qc_gc_curve_plot"
          curve_png <- file.path(dirs$nipter_gc_curves, paste0(stem, ".gc_curve.png"))
          .write_nipter_gc_curve_plot(
            curve_data = nipter_qc$curve_data,
            sample_name = stem,
            out_png = curve_png,
            theme = opts$`qc-plot-theme`,
            base_size = as.integer(opts$`qc-plot-base-size`),
            loess_span = as.numeric(opts$`nipter-gc-loess-span`)
          )
          nipter_qc$metrics$gc_curve_plot <- curve_png
        }
        if (isTRUE(opts$`nipter-corrected-bin-ratio-genome-plots`)) {
          nipter_failure_step <<- "qc_corrected_bin_ratio_plot"
          genome_plot_png <- file.path(
            dirs$nipter_corrected_bin_ratio_genome_plots,
            paste0(stem, ".corrected_bin_ratio_genome.png")
          )
          .write_nipter_corrected_bin_ratio_genome_plot(
            curve_data = nipter_qc$curve_data,
            sample_name = stem,
            out_png = genome_plot_png,
            theme = opts$`qc-plot-theme`,
            base_size = as.integer(opts$`qc-plot-base-size`)
          )
          nipter_qc$metrics$corrected_bin_ratio_genome_plot <- genome_plot_png
        }
        nipter_failure_step <<- "write_bed"
        nipter_sample_to_bed(
          sample = raw_sample,
          corrected = corrected_sample,
          bed = out_bed,
          binsize = as.integer(opts$`nipter-binsize`),
          con = con,
          metadata = .merge_tabix_metadata(
            .compact_metadata(list(
              format = "nipter_bed",
              schema = if (isTRUE(opts$`nipter-separate-strands`)) "separated_v1" else "combined_v1",
              binsize = as.integer(opts$`nipter-binsize`),
              mapq = as.integer(opts$`nipter-mapq`),
              require_flags = 0L,
              exclude_flags = as.integer(opts$`nipter-exclude-flags`),
              rmdup = "none",
              corrected_columns = TRUE,
              reference = normalizePath(opts$fasta, winslash = "/", mustWork = TRUE)
            )),
            .merge_tabix_metadata(
              .sample_qc_tabix_metadata(
                seqff = seqff_metrics_by_sample[[stem]],
                y_unique = y_unique_metrics_by_sample[[stem]],
                read_counts = read_counts_metrics,
                bam_stats = bam_stats_by_sample[[stem]],
                coverage = coverage_metrics_by_sample[[stem]],
                nipter_qc = nipter_qc$metrics,
                gc_correction = .nipter_gc_correction_record(),
                filters_pre = metrics_pre_filters,
                filters_post = metrics_post_filters
              )
            )
          )
        )
        list(
          ok = TRUE,
          read_counts_metrics = read_counts_metrics,
          nipter_qc = nipter_qc$metrics,
          failure_step = NULL,
          failure_message = NULL
        )
      }),
      error = function(e) {
        list(
          ok = FALSE,
          read_counts_metrics = NULL,
          nipter_qc = NULL,
          failure_step = nipter_failure_step,
          failure_message = conditionMessage(e)
        )
      }
    )
    if (!isTRUE(result$ok)) {
      partial_paths <- c(out_bed, paste0(out_bed, ".tbi"))
      partial_paths <- partial_paths[file.exists(partial_paths)]
      if (length(partial_paths)) {
        unlink(partial_paths, force = TRUE)
      }
    }
    read_counts_metrics <- result$read_counts_metrics
    failure_message <- result$failure_message
    row <- .sample_qc_row(
      sample_name = stem,
      bam = bam,
      seqff = seqff_metrics_by_sample[[stem]],
      y_unique = y_unique_metrics_by_sample[[stem]],
      read_counts = read_counts_metrics,
      bam_stats = bam_stats_by_sample[[stem]],
      coverage = coverage_metrics_by_sample[[stem]],
      nipter_qc = result$nipter_qc,
      gc_correction = .nipter_gc_correction_record(),
      filters_pre = metrics_pre_filters,
      filters_post = metrics_post_filters
    )
    row <- .sample_qc_row_with_processing_status(
      row,
      status = if (isTRUE(result$ok)) "ok" else "failed",
      failure_stage = if (isTRUE(result$ok)) NULL else "nipter_bed",
      failure_step = if (isTRUE(result$ok)) NULL else result$failure_step,
      failure_message = if (isTRUE(result$ok)) NULL else failure_message
    )
    row <- .sample_qc_row_with_nipter_status(
      row,
      written = isTRUE(result$ok) && file.exists(out_bed),
      error = if (isTRUE(result$ok)) NULL else failure_message,
      status = if (isTRUE(result$ok)) "ok" else "failed"
    )[, .sample_qc_summary_cols, drop = FALSE]
    utils::write.table(
      row,
      file = out_qc,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
    if (!isTRUE(result$ok)) {
      cat(sprintf(
        "[%s] [%d/%d] NIPTeR BED failed %s: %s\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        i, length(bams), stem, failure_message
      ))
    }
    list(row = row, read_counts_metrics = read_counts_metrics)
  })
  for (item in sample_qc_rows) {
    row <- item$row[, .sample_qc_summary_cols, drop = FALSE]
    stem <- as.character(row$sample_name[[1L]])
    read_counts_by_sample[[stem]] <- item$read_counts_metrics
    final_sample_qc_rows[[stem]] <- row
    if (identical(row$sample_processing_status[[1L]], "failed")) {
      sample_failures_by_sample[[stem]] <- list(
        stage = as.character(row$sample_failure_stage[[1L]]),
        step = as.character(row$sample_failure_step[[1L]]),
        message = as.character(row$sample_failure_message[[1L]])
      )
    }
  }
  for (stem in sample_names) {
    if (is.null(final_sample_qc_rows[[stem]])) {
      final_sample_qc_rows[[stem]] <- .build_sample_qc_row(stem, bams[[stem]])
    }
    utils::write.table(
      final_sample_qc_rows[[stem]][, .sample_qc_summary_cols, drop = FALSE],
      file = file.path(dirs$sample_qc, paste0(stem, ".sample_qc.tsv")),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }
  sample_qc_summary <- do.call(rbind, final_sample_qc_rows[sample_names])
  sample_qc_summary <- sample_qc_summary[, .sample_qc_summary_cols, drop = FALSE]
  utils::write.table(
    sample_qc_summary,
    file = file.path(dirs$sample_qc, "sample_qc.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  sample_failures <- sample_qc_summary[
    sample_qc_summary$sample_processing_status %in% "failed",
    c(
      "sample_name",
      "bam",
      "sample_processing_status",
      "sample_failure_stage",
      "sample_failure_step",
      "sample_failure_message",
      "nipter_bed_status",
      "nipter_bed_written",
      "nipter_bed_error"
    ),
    drop = FALSE
  ]
  utils::write.table(
    sample_failures,
    file = file.path(dirs$sample_qc, "sample_failures.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

}, sprintf("Convert %d BAMs to NIPTeR BED.gz", length(.active_bams())))

.run_one({
  .write_bed_header_qc_summary(
    bed_dir = dirs$nipter_beds,
    out_tsv = file.path(dirs$sample_qc, "nipter_bed_header_qc.tsv")
  )
}, sprintf("Summarize NIPTeR BED headers for %d samples", length(bams)))

cat("\nPreprocessing layout ready under:\n")
cat("  out_root:    ", out_root, "\n", sep = "")
cat("  manifest:    ", staged_manifest, "\n", sep = "")
cat("  jobs:        ", as.integer(opts$jobs), "\n", sep = "")
cat("  gc_table:    ", gc_table, "\n", sep = "")
cat("  y_regions:   ", y_regions_file, "\n", sep = "")
cat("  seqff:       ", dirs$seqff, "\n", sep = "")
cat("  y_unique:    ", dirs$y_unique, "\n", sep = "")
cat("  sample_qc:   ", file.path(dirs$sample_qc, "sample_qc.tsv"), "\n", sep = "")
cat("  sample_failures: ", file.path(dirs$sample_qc, "sample_failures.tsv"), "\n", sep = "")
cat("  rwcx_header_qc: ", file.path(dirs$sample_qc, "rwcx_bed_header_qc.tsv"), "\n", sep = "")
cat("  nipter_header_qc: ", file.path(dirs$sample_qc, "nipter_bed_header_qc.tsv"), "\n", sep = "")
cat("  rwcx_beds:   ", dirs$rwcx_beds, "\n", sep = "")
cat("  wisecondorx_npz: ", dirs$wisecondorx_npz, "\n", sep = "")
cat("  nipter_beds: ", dirs$nipter_beds, "\n", sep = "")
cat("  nipter_gc_curves: ", dirs$nipter_gc_curves, "\n", sep = "")
cat("  nipter_gc_curve_data: ", dirs$nipter_gc_curve_data, "\n", sep = "")
cat("  nipter_corrected_bin_ratio_genome_plots: ", dirs$nipter_corrected_bin_ratio_genome_plots, "\n", sep = "")
cat("  coverage_pre: ", dirs$coverage_pre, "\n", sep = "")
cat("  coverage_post: ", dirs$coverage_post, "\n", sep = "")
cat("  coverage_binsize: ", coverage_binsize, "\n", sep = "")
