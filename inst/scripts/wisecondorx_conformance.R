#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(RWisecondorX)
})

option_list <- list(
  make_option("--rwisecondorx-ref", type = "character",
              help = "Path to native RWisecondorX reference RDS."),
  make_option("--wisecondorx-ref", type = "character",
              help = "Path to upstream WisecondorX reference NPZ."),
  make_option("--case-manifest", type = "character", default = NULL,
              help = "TSV with columns sample, bed, npz."),
  make_option("--bed-dir", type = "character", default = NULL,
              help = "Directory containing native rwisecondorx BED.gz test cases."),
  make_option("--npz-dir", type = "character", default = NULL,
              help = "Directory containing upstream wisecondorx NPZ test cases."),
  make_option("--bed-pattern", type = "character", default = "\\.bed\\.gz$",
              help = "Regex pattern used with --bed-dir [default %default]."),
  make_option("--npz-pattern", type = "character", default = "\\.npz$",
              help = "Regex pattern used with --npz-dir [default %default]."),
  make_option("--out-dir", type = "character", default = "wisecondorx_conformance",
              help = "Output directory [default %default]."),
  make_option("--sample-binsize", type = "integer", default = 100000L,
              help = "Sample binsize used to build the prepared test cases [default %default]."),
  make_option("--cpus", type = "integer", default = 4L,
              help = "Thread count for native predict [default %default]."),
  make_option("--minrefbins", type = "integer", default = 150L,
              help = "Predict minrefbins [default %default]."),
  make_option("--maskrepeats", type = "integer", default = 5L,
              help = "Predict maskrepeats [default %default]."),
  make_option("--zscore", type = "double", default = 5,
              help = "Predict z-score threshold [default %default]."),
  make_option("--alpha", type = "double", default = 1e-4,
              help = "Predict CBS alpha [default %default]."),
  make_option("--seed", type = "integer", default = 42L,
              help = "Predict seed [default %default]."),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Overwrite per-case output directories.")
)

parser <- OptionParser(
  usage = "%prog --rwisecondorx-ref ref.rds --wisecondorx-ref ref.npz --case-manifest cases.tsv [options]",
  option_list = option_list
)
opt <- parse_args(parser)

.fail <- function(...) stop(sprintf(...), call. = FALSE)

.validate_inputs <- function(opt) {
  if (is.null(opt$rwisecondorx_ref) || !nzchar(opt$rwisecondorx_ref)) {
    .fail("--rwisecondorx-ref is required.")
  }
  if (is.null(opt$wisecondorx_ref) || !nzchar(opt$wisecondorx_ref)) {
    .fail("--wisecondorx-ref is required.")
  }
  if (!file.exists(opt$rwisecondorx_ref)) {
    .fail("Native reference not found: %s", opt$rwisecondorx_ref)
  }
  if (!file.exists(opt$wisecondorx_ref)) {
    .fail("Upstream reference not found: %s", opt$wisecondorx_ref)
  }
  has_manifest <- !is.null(opt$case_manifest)
  has_dirs <- !is.null(opt$bed_dir) || !is.null(opt$npz_dir)
  if (has_manifest == has_dirs) {
    .fail("Provide either --case-manifest or the pair --bed-dir/--npz-dir.")
  }
  if (has_dirs && (is.null(opt$bed_dir) || is.null(opt$npz_dir))) {
    .fail("--bed-dir and --npz-dir must be provided together.")
  }
}

.read_case_manifest <- function(path) {
  if (!file.exists(path)) {
    .fail("Case manifest not found: %s", path)
  }
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  needed <- c("sample", "bed", "npz")
  if (!all(needed %in% names(df))) {
    .fail("Case manifest must contain columns: %s", paste(needed, collapse = ", "))
  }
  if (nrow(df) == 0L) {
    .fail("Case manifest is empty: %s", path)
  }
  df <- df[, needed]
  df$sample <- trimws(df$sample)
  df$bed <- trimws(df$bed)
  df$npz <- trimws(df$npz)
  bad <- !nzchar(df$sample) | !nzchar(df$bed) | !nzchar(df$npz)
  if (any(bad)) {
    .fail("Case manifest contains blank sample/bed/npz fields.")
  }
  df
}

.resolve_cases <- function(opt) {
  cases <- if (!is.null(opt$case_manifest)) {
    .read_case_manifest(opt$case_manifest)
  } else {
    if (!dir.exists(opt$bed_dir)) {
      .fail("BED directory not found: %s", opt$bed_dir)
    }
    if (!dir.exists(opt$npz_dir)) {
      .fail("NPZ directory not found: %s", opt$npz_dir)
    }
    beds <- list.files(opt$bed_dir, pattern = opt$bed_pattern, full.names = TRUE)
    npzs <- list.files(opt$npz_dir, pattern = opt$npz_pattern, full.names = TRUE)
    if (length(beds) == 0L) {
      .fail("No BED cases found in %s", opt$bed_dir)
    }
    if (length(npzs) == 0L) {
      .fail("No NPZ cases found in %s", opt$npz_dir)
    }
    bed_stems <- sub("\\.bed\\.gz$", "", basename(beds), ignore.case = TRUE)
    npz_stems <- sub("\\.npz$", "", basename(npzs), ignore.case = TRUE)
    shared <- intersect(bed_stems, npz_stems)
    if (length(shared) == 0L) {
      .fail("No shared sample stems between %s and %s", opt$bed_dir, opt$npz_dir)
    }
    data.frame(
      sample = shared,
      bed = beds[match(shared, bed_stems)],
      npz = npzs[match(shared, npz_stems)],
      stringsAsFactors = FALSE
    )
  }
  cases$bed <- normalizePath(cases$bed, winslash = "/", mustWork = FALSE)
  cases$npz <- normalizePath(cases$npz, winslash = "/", mustWork = FALSE)
  missing <- unique(c(cases$bed[!file.exists(cases$bed)], cases$npz[!file.exists(cases$npz)]))
  if (length(missing) > 0L) {
    .fail("Missing case files:\n  %s", paste(missing, collapse = "\n  "))
  }
  cases
}

.load_bins_bed <- function(path) {
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (!all(c("chr", "start", "end", "id", "ratio", "zscore") %in% names(df))) {
    .fail("Unexpected bins BED columns in %s", path)
  }
  df
}

.compare_bins <- function(native_bed, upstream_bed, tol = 1e-6) {
  native <- .load_bins_bed(native_bed)
  upstream <- .load_bins_bed(upstream_bed)

  shared_cols <- c("chr", "start", "end", "id")
  same_shape <- nrow(native) == nrow(upstream)
  same_keys <- same_shape && identical(native[shared_cols], upstream[shared_cols])

  if (!same_keys) {
    return(list(
      ok = FALSE,
      native_rows = nrow(native),
      upstream_rows = nrow(upstream),
      key_mismatches = NA_integer_,
      ratio_mismatches = NA_integer_,
      zscore_mismatches = NA_integer_,
      max_abs_ratio_diff = NA_real_,
      max_abs_zscore_diff = NA_real_
    ))
  }

  ratio_diff <- abs(native$ratio - upstream$ratio)
  zscore_diff <- abs(native$zscore - upstream$zscore)
  ratio_bad <- !(is.na(ratio_diff) | ratio_diff <= tol)
  zscore_bad <- !(is.na(zscore_diff) | zscore_diff <= tol)

  list(
    ok = !any(ratio_bad) && !any(zscore_bad),
    native_rows = nrow(native),
    upstream_rows = nrow(upstream),
    key_mismatches = 0L,
    ratio_mismatches = sum(ratio_bad),
    zscore_mismatches = sum(zscore_bad),
    max_abs_ratio_diff = if (length(ratio_diff)) max(ratio_diff, na.rm = TRUE) else 0,
    max_abs_zscore_diff = if (length(zscore_diff)) max(zscore_diff, na.rm = TRUE) else 0
  )
}

.validate_inputs(opt)
cases <- .resolve_cases(opt)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

native_ref <- readRDS(opt$rwisecondorx_ref)
native_ref_binsize <- as.integer(native_ref$binsize)
if (is.na(native_ref_binsize) || native_ref_binsize < 1L) {
  .fail("Could not determine binsize from native reference: %s", opt$rwisecondorx_ref)
}

summary_rows <- vector("list", nrow(cases))

for (i in seq_len(nrow(cases))) {
  case <- cases[i, , drop = FALSE]
  sample_name <- case$sample[[1L]]
  bed <- case$bed[[1L]]
  npz <- case$npz[[1L]]
  case_dir <- file.path(opt$out_dir, sample_name)

  if (dir.exists(case_dir) && !isTRUE(opt$overwrite)) {
    .fail("Case output exists, use --overwrite to replace: %s", case_dir)
  }
  unlink(case_dir, recursive = TRUE, force = TRUE)
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("[%d/%d] %s", i, length(cases), sample_name))

  native_prefix <- file.path(case_dir, "native")
  upstream_prefix <- file.path(case_dir, "upstream")
  sample <- bed_to_sample(bed, binsize = opt$sample_binsize)

  pred_native <- rwisecondorx_predict(
    sample = sample,
    reference = native_ref,
    sample_binsize = opt$sample_binsize,
    outprefix = native_prefix,
    minrefbins = opt$minrefbins,
    maskrepeats = opt$maskrepeats,
    zscore = opt$zscore,
    alpha = opt$alpha,
    seed = opt$seed,
    cpus = opt$cpus
  )

  wisecondorx_predict(
    npz = npz,
    ref = opt$wisecondorx_ref,
    output_prefix = upstream_prefix,
    minrefbins = opt$minrefbins,
    maskrepeats = opt$maskrepeats,
    zscore = opt$zscore,
    alpha = opt$alpha,
    bed = TRUE,
    plot = FALSE,
    seed = opt$seed
  )

  native_bins <- paste0(native_prefix, "_bins.bed")
  upstream_bins <- paste0(upstream_prefix, "_bins.bed")
  if (!file.exists(native_bins)) {
    .fail("Native bins BED not produced: %s", native_bins)
  }
  if (!file.exists(upstream_bins)) {
    .fail("Upstream bins BED not produced: %s", upstream_bins)
  }

  cmp <- .compare_bins(native_bins, upstream_bins)
  summary_rows[[i]] <- data.frame(
    sample = sample_name,
    bed = bed,
    npz = npz,
    native_bins = native_bins,
    upstream_bins = upstream_bins,
    native_ref_binsize = native_ref_binsize,
    sample_binsize = opt$sample_binsize,
    ok = cmp$ok,
    native_rows = cmp$native_rows,
    upstream_rows = cmp$upstream_rows,
    key_mismatches = cmp$key_mismatches,
    ratio_mismatches = cmp$ratio_mismatches,
    zscore_mismatches = cmp$zscore_mismatches,
    max_abs_ratio_diff = cmp$max_abs_ratio_diff,
    max_abs_zscore_diff = cmp$max_abs_zscore_diff,
    native_aberrations = nrow(pred_native$aberrations),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(opt$out_dir, "wisecondorx_conformance_summary.tsv")
utils::write.table(
  summary_df,
  file = summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("Wrote summary: ", summary_path)
message("Cases: ", nrow(summary_df))
message("Exact bin matches: ", sum(summary_df$ok), "/", nrow(summary_df))

if (!all(summary_df$ok)) {
  quit(status = 1L)
}
