#!/usr/bin/env Rscript
#
# predict_cohort.R
#
# Score a preprocessed cohort against RWisecondorX and/or NIPTeR
# references. This script expects BED-based cohort artifacts that were already
# produced by preprocess_cohort.R plus reference objects from build_reference.R.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with install.packages('optparse').",
       call. = FALSE)
}
library(optparse)

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
}

.read_file_list <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) >= 1L && identical(lines[[1L]], "bam")) {
    warning(
      "List ", path, " starts with a 'bam' header line; stripping it.",
      call. = FALSE
    )
    lines <- lines[-1L]
  }
  lines
}

.sample_stem <- function(path) {
  sub("\\.(bam|cram)$", "", basename(path), ignore.case = TRUE)
}

.ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.match_bed_paths <- function(sample_names, paths, label) {
  stems <- vapply(paths, .sample_stem, character(1L))
  dup <- duplicated(stems)
  if (any(dup)) {
    stop(
      label, " contains duplicate sample stems: ",
      paste(unique(stems[dup]), collapse = ", "),
      call. = FALSE
    )
  }
  idx <- match(sample_names, stems)
  missing <- sample_names[is.na(idx)]
  if (length(missing)) {
    stop(
      label, " is missing samples: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  normalizePath(paths[idx], winslash = "/", mustWork = FALSE)
}

.write_tsv <- function(df, path) {
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

.run_stage_samples <- function(items, jobs, worker) {
  jobs <- as.integer(jobs)
  stopifnot(length(jobs) == 1L, jobs >= 1L)
  if (!length(items)) {
    return(list())
  }
  if (.Platform$OS.type == "unix" && jobs > 1L) {
    parallel::mclapply(
      X = seq_along(items),
      FUN = function(i) worker(i, items[[i]]),
      mc.cores = jobs,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_along(items), function(i) worker(i, items[[i]]))
  }
}

.per_worker_threads <- function(total_threads, jobs) {
  total_threads <- as.integer(total_threads)
  jobs <- as.integer(jobs)
  max(1L, total_threads %/% max(1L, jobs))
}

.log_sample <- function(i, n, label, stem) {
  cat(sprintf("[%s] [%d/%d] %s %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              i, n, label, stem))
}

.match_sample_qc_row <- function(sample_qc, sample_name) {
  hit <- which(as.character(sample_qc$sample_name) == sample_name)
  if (!length(hit)) {
    return(NULL)
  }
  sample_qc[hit[[1L]], , drop = FALSE]
}

.finite_numeric_predictors <- function(row) {
  if (is.null(row) || !nrow(row)) {
    return(NULL)
  }
  vals <- lapply(row[1L, , drop = FALSE], function(x) suppressWarnings(as.numeric(x[[1L]])))
  keep <- vapply(vals, function(x) length(x) == 1L && is.finite(x), logical(1L))
  if (!any(keep)) {
    return(NULL)
  }
  as.list(stats::setNames(unlist(vals[keep], use.names = FALSE), names(vals)[keep]))
}

.qc_scalar <- function(row, col) {
  if (is.null(row) || !nrow(row) || !(col %in% names(row))) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(row[[col]][[1L]]))
}

.qc_logical <- function(row, col) {
  if (is.null(row) || !nrow(row) || !(col %in% names(row))) {
    return(NA)
  }
  as.logical(row[[col]][[1L]])
}

.qc_text <- function(row, col) {
  if (is.null(row) || !nrow(row) || !(col %in% names(row))) {
    return(NA_character_)
  }
  val <- as.character(row[[col]][[1L]])
  if (is.na(val) || !nzchar(val) || identical(val, "NA")) NA_character_ else val
}

.rwx_stat_scalar <- function(prediction, chrom, field = "zscore") {
  stats_df <- prediction$statistics
  hit <- which(as.character(stats_df$chr) == as.character(chrom))
  if (!length(hit) || !(field %in% names(stats_df))) {
    return(NA_real_)
  }
  as.numeric(stats_df[[field]][hit[[1L]]])
}

.rwx_aberration_label <- function(prediction) {
  ab <- prediction$aberrations
  if (is.null(ab) || !nrow(ab)) {
    return(NA_character_)
  }
  paste(
    sprintf("%s:%s", as.character(ab$chr), as.character(ab$type)),
    collapse = ","
  )
}

.normalize_sex_label <- function(x) {
  x <- trimws(as.character(x))
  if (!nzchar(x)) {
    return(NA_character_)
  }
  if (x %in% c("M", "male", "Male")) {
    return("male")
  }
  if (x %in% c("F", "female", "Female")) {
    return("female")
  }
  x
}

.read_wisecondorx_statistics <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) {
    stop("Empty WisecondorX statistics file: ", path, call. = FALSE)
  }
  footer_idx <- which(!grepl("\t", lines, fixed = TRUE))
  table_end <- if (length(footer_idx)) footer_idx[[1L]] - 1L else length(lines)
  if (table_end < 1L) {
    stop("Could not locate the chromosome statistics table in: ", path, call. = FALSE)
  }
  stats_df <- utils::read.delim(
    text = paste(lines[seq_len(table_end)], collapse = "\n"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  footer_lines <- if (table_end < length(lines)) lines[(table_end + 1L):length(lines)] else character()
  list(
    statistics = stats_df,
    predicted_gender = .normalize_sex_label(sub("^.*: ", "", grep("^Gender based on", footer_lines, value = TRUE)[1] %||% "")),
    n_reads = suppressWarnings(as.numeric(sub("^.*: ", "", grep("^Number of reads:", footer_lines, value = TRUE)[1] %||% NA_character_))),
    ratio_sd = suppressWarnings(as.numeric(sub("^.*: ", "", grep("^Standard deviation of the ratios per chromosome:", footer_lines, value = TRUE)[1] %||% NA_character_))),
    median_segment_variance = suppressWarnings(as.numeric(sub("^.*: ", "", grep("^Median segment variance per bin", footer_lines, value = TRUE)[1] %||% NA_character_))),
    cpa_score = suppressWarnings(as.numeric(sub("^.*: ", "", grep("^Copy number profile abnormality \\(CPA\\) score", footer_lines, value = TRUE)[1] %||% NA_character_)))
  )
}

.read_wisecondorx_aberrations <- function(path) {
  if (!file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      ratio = numeric(),
      zscore = numeric(),
      type = character(),
      stringsAsFactors = FALSE
    ))
  }
  utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.wcx_upstream_stat_scalar <- function(stats_df, chrom, field = "zscore") {
  hit <- which(as.character(stats_df$chr) == as.character(chrom))
  if (!length(hit) || !(field %in% names(stats_df))) {
    return(NA_real_)
  }
  as.numeric(stats_df[[field]][hit[[1L]]])
}

.predict_wisecondorx_one <- function(sample_name,
                                     npz_path,
                                     reference,
                                     out_dir,
                                     minrefbins,
                                     maskrepeats,
                                     zscore,
                                     alpha,
                                     beta,
                                     plot,
                                     seed,
                                     overwrite) {
  sample_dir <- .ensure_dir(file.path(out_dir, sample_name))
  outprefix <- file.path(sample_dir, sample_name)
  stat_path <- paste0(outprefix, "_statistics.txt")
  aberr_path <- paste0(outprefix, "_aberrations.bed")

  if (isTRUE(overwrite) || !file.exists(stat_path)) {
    wisecondorx_predict(
      npz = npz_path,
      ref = reference,
      output_prefix = outprefix,
      minrefbins = as.integer(minrefbins),
      maskrepeats = as.integer(maskrepeats),
      zscore = as.numeric(zscore),
      alpha = as.numeric(alpha),
      beta = if (is.null(beta) || !is.finite(beta)) NULL else as.numeric(beta),
      bed = TRUE,
      plot = isTRUE(plot),
      seed = if (is.null(seed)) NULL else as.integer(seed)
    )
  }

  stat_obj <- .read_wisecondorx_statistics(stat_path)
  aberr <- .read_wisecondorx_aberrations(aberr_path)

  data.frame(
    sample_name = sample_name,
    status = "ok",
    error = NA_character_,
    predicted_gender = stat_obj$predicted_gender,
    n_reads = stat_obj$n_reads,
    ratio_sd = stat_obj$ratio_sd,
    median_segment_variance = stat_obj$median_segment_variance,
    cpa_score = stat_obj$cpa_score,
    n_aberrations = nrow(aberr),
    aberrations = if (nrow(aberr)) paste(sprintf("%s:%s", as.character(aberr$chr), as.character(aberr$type)), collapse = ",") else NA_character_,
    z_13 = .wcx_upstream_stat_scalar(stat_obj$statistics, 13, "zscore"),
    z_18 = .wcx_upstream_stat_scalar(stat_obj$statistics, 18, "zscore"),
    z_21 = .wcx_upstream_stat_scalar(stat_obj$statistics, 21, "zscore"),
    z_X = .wcx_upstream_stat_scalar(stat_obj$statistics, "X", "zscore"),
    out_dir = sample_dir,
    stringsAsFactors = FALSE
  )
}

.suffix_nonkey_cols <- function(df, suffix, key = "sample_name") {
  nms <- names(df)
  idx <- nms != key
  nms[idx] <- paste0(nms[idx], "_", suffix)
  names(df) <- nms
  df
}

.sex_summary_row <- function(sample_name, sex_score, sample_qc_row = NULL) {
  out <- data.frame(
    sample_name = sample_name,
    predicted_sex = sex_score$predicted_sex,
    y_unique_ratio_post = .qc_scalar(sample_qc_row, "y_unique_ratio_post"),
    gc_read_perc_post = .qc_scalar(sample_qc_row, "gc_read_perc_post"),
    stringsAsFactors = FALSE
  )
  for (nm in names(sex_score$z_scores)) {
    out[[nm]] <- unname(sex_score$z_scores[[nm]])
  }
  for (nm in names(sex_score$cv)) {
    out[[nm]] <- unname(sex_score$cv[[nm]])
  }
  for (nm in names(sex_score$reference_sizes)) {
    out[[paste0("reference_size_", nm)]] <- unname(sex_score$reference_sizes[[nm]])
  }
  out
}

.autosomal_prediction_tables <- function(sample,
                                         control_group,
                                         focus_chromosomes,
                                         regression_n_models,
                                         regression_n_predictors,
                                         regression_train_fraction,
                                         regression_exclude_chromosomes,
                                         regression_seed) {
  summary_rows <- vector("list", length(focus_chromosomes))
  regression_rows <- list()

  for (i in seq_along(focus_chromosomes)) {
    chrom <- as.integer(focus_chromosomes[[i]])
    z_obj <- nipter_z_score(sample, control_group, chromo_focus = chrom)
    ncv_obj <- nipter_ncv_score(sample, control_group, chromo_focus = chrom)
    reg_obj <- nipter_regression(
      sample = sample,
      control_group = control_group,
      chromo_focus = chrom,
      n_models = as.integer(regression_n_models),
      n_predictors = as.integer(regression_n_predictors),
      train_fraction = regression_train_fraction,
      exclude_chromosomes = as.integer(regression_exclude_chromosomes),
      force_practical_cv = TRUE,
      seed = if (is.null(regression_seed)) NULL else as.integer(regression_seed)
    )

    reg_scores <- vapply(reg_obj$models, function(m) m$z_score, numeric(1L))
    summary_rows[[i]] <- data.frame(
      chromosome = chrom,
      z_score = unname(z_obj$sample_z_score),
      z_cv = 100 * unname(z_obj$control_statistics[["sd"]]) /
        unname(z_obj$control_statistics[["mean"]]),
      z_shapiro_p = unname(z_obj$control_statistics[["shapiro_p_value"]]),
      ncv_score = unname(ncv_obj$sample_score),
      ncv_cv = unname(ncv_obj$best_cv) * 100,
      ncv_shapiro_p = unname(ncv_obj$control_statistics[["shapiro_p_value"]]),
      ncv_denominators = paste(ncv_obj$denominators, collapse = ","),
      rbz_mean = mean(reg_scores),
      rbz_median = stats::median(reg_scores),
      rbz_train_fraction_requested = reg_obj$train_fraction_requested %||% regression_train_fraction,
      rbz_train_fraction_effective = reg_obj$train_fraction_effective %||% NA_real_,
      rbz_split_mode = reg_obj$split_mode %||% NA_character_,
      stringsAsFactors = FALSE
    )

    regression_rows[[length(regression_rows) + 1L]] <- do.call(
      rbind,
      lapply(seq_along(reg_obj$models), function(j) {
        mdl <- reg_obj$models[[j]]
        data.frame(
          chromosome = chrom,
          model_id = paste0("predictor_set_", j),
          regression_z = mdl$z_score,
          regression_cv = 100 * mdl$cv,
          regression_cv_type = mdl$cv_type,
          regression_shapiro_p = mdl$shapiro_p_value,
          predictors = paste(mdl$predictors, collapse = ","),
          n_reference_samples = mdl$n_reference_samples %||% NA_integer_,
          n_training_samples = mdl$n_training_samples %||% NA_integer_,
          n_stat_samples = mdl$n_stat_samples %||% NA_integer_,
          train_fraction_requested = mdl$train_fraction_requested %||% regression_train_fraction,
          train_fraction_effective = mdl$train_fraction_effective %||% NA_real_,
          split_mode = mdl$split_mode %||% NA_character_,
          split_seed = mdl$split_seed %||% NA_integer_,
          stringsAsFactors = FALSE
        )
      })
    )
  }

  list(
    summary = do.call(rbind, summary_rows),
    regression = do.call(rbind, regression_rows)
  )
}

.predict_rwisecondorx_one <- function(sample_name,
                                      bed_path,
                                      reference,
                                      out_dir,
                                      sample_binsize,
                                      minrefbins,
                                      maskrepeats,
                                      zscore,
                                      alpha,
                                      beta,
                                      optimal_cutoff_sd_multiplier,
                                      within_sample_mask_iterations,
                                      within_sample_mask_quantile,
                                      cbs_split_min_gap_bp,
                                      segment_zscore_cap,
                                      cpus,
                                      overwrite,
                                      seed,
                                      sample_qc_row = NULL) {
  sample_dir <- .ensure_dir(file.path(out_dir, sample_name))
  outprefix <- file.path(sample_dir, sample_name)
  pred_rds <- file.path(sample_dir, "prediction.rds")

  if (!isTRUE(overwrite) && file.exists(pred_rds)) {
    prediction <- readRDS(pred_rds)
  } else {
    sample <- bed_to_sample(bed_path)
    prediction <- rwisecondorx_predict(
      sample = sample,
      reference = reference,
      sample_binsize = as.integer(sample_binsize),
      outprefix = outprefix,
      minrefbins = as.integer(minrefbins),
      maskrepeats = as.integer(maskrepeats),
      zscore = as.numeric(zscore),
      alpha = as.numeric(alpha),
      beta = if (is.finite(beta)) as.numeric(beta) else NULL,
      optimal_cutoff_sd_multiplier = as.numeric(optimal_cutoff_sd_multiplier),
      within_sample_mask_iterations = as.integer(within_sample_mask_iterations),
      within_sample_mask_quantile = as.numeric(within_sample_mask_quantile),
      cbs_split_min_gap_bp = as.integer(cbs_split_min_gap_bp),
      segment_zscore_cap = as.numeric(segment_zscore_cap),
      seed = if (is.null(seed)) NULL else as.integer(seed),
      cpus = as.integer(cpus)
    )
    saveRDS(prediction, pred_rds)
  }

  data.frame(
    sample_name = sample_name,
    status = "ok",
    error = NA_character_,
    predicted_gender = prediction$gender,
    n_reads = prediction$n_reads,
    n_aberrations = nrow(prediction$aberrations),
    aberrations = .rwx_aberration_label(prediction),
    z_13 = .rwx_stat_scalar(prediction, 13, "zscore"),
    z_18 = .rwx_stat_scalar(prediction, 18, "zscore"),
    z_21 = .rwx_stat_scalar(prediction, 21, "zscore"),
    z_X = .rwx_stat_scalar(prediction, "X", "zscore"),
    read_counts_input_total = .qc_scalar(sample_qc_row, "read_counts_input_total"),
    read_counts_mapped_total = .qc_scalar(sample_qc_row, "read_counts_mapped_total"),
    read_counts_binned_post_sum = .qc_scalar(sample_qc_row, "read_counts_binned_post_sum"),
    fetal_fraction_post = .qc_scalar(sample_qc_row, "fetal_fraction_post"),
    y_unique_ratio_post = .qc_scalar(sample_qc_row, "y_unique_ratio_post"),
    gc_curve_has_loess_support = .qc_logical(sample_qc_row, "gc_curve_has_loess_support"),
    nipter_bed_status = .qc_text(sample_qc_row, "nipter_bed_status"),
    out_dir = sample_dir,
    stringsAsFactors = FALSE
  )
}

.predict_nipter_one <- function(sample_name,
                                bed_path,
                                reference,
                                sample_qc_row,
                                out_dir,
                                autosomal_source,
                                sex_counts,
                                chi_cutoff,
                                focus_chromosomes,
                                regression_n_models,
                                regression_n_predictors,
                                regression_train_fraction,
                                regression_exclude_chromosomes,
                                regression_seed,
                                overwrite) {
  sample_dir <- .ensure_dir(file.path(out_dir, sample_name))
  out_rds <- file.path(sample_dir, "prediction.rds")

  y_unique_ratio <- .qc_scalar(sample_qc_row, "y_unique_ratio_post")
  if ("y_unique" %in% names(reference$sex_models) && !is.finite(y_unique_ratio)) {
    stop(
      "Reference includes a y_unique sex model but the sample lacks y_unique_ratio_post in sample_qc.tsv.",
      call. = FALSE
    )
  }

  if (!isTRUE(overwrite) && file.exists(out_rds)) {
    payload <- readRDS(out_rds)
    sex_summary <- payload$sex_summary
    gauno_summary <- payload$gaunosome_summary
    autosomal_summary <- payload$autosomal_summary
    regression_summary <- payload$autosomal_regression
  } else {
    sample_raw <- bed_to_nipter_sample(
      bed_path,
      autosomal_source = autosomal_source,
      sex_source = sex_counts
    )
    chi <- nipter_chi_correct(
      sample = sample_raw,
      control_group = reference$control_group,
      chi_cutoff = as.numeric(chi_cutoff)
    )
    sample_chi <- chi$sample
    control_chi <- chi$control_group
    sample_predictors <- .finite_numeric_predictors(sample_qc_row)

    sex_score <- nipter_sex_score(
      sample = sample_chi,
      reference = reference,
      y_unique_ratio = if (is.finite(y_unique_ratio)) y_unique_ratio else NULL
    )
    gauno <- nipter_gaunosome_score(
      sample = sample_chi,
      reference = reference,
      y_unique_ratio = if (is.finite(y_unique_ratio)) y_unique_ratio else NULL,
      sample_predictors = sample_predictors
    )
    autosomal <- .autosomal_prediction_tables(
      sample = sample_chi,
      control_group = control_chi,
      focus_chromosomes = focus_chromosomes,
      regression_n_models = regression_n_models,
      regression_n_predictors = regression_n_predictors,
      regression_train_fraction = regression_train_fraction,
      regression_exclude_chromosomes = regression_exclude_chromosomes,
      regression_seed = regression_seed
    )

    sex_summary <- .sex_summary_row(sample_name, sex_score, sample_qc_row = sample_qc_row)
    gauno_summary <- data.frame(sample_name = sample_name, gauno$summary, stringsAsFactors = FALSE)
    autosomal_summary <- autosomal$summary
    regression_summary <- autosomal$regression

    .write_tsv(sex_summary, file.path(sample_dir, "sex_summary.tsv"))
    .write_tsv(gauno_summary, file.path(sample_dir, "gaunosome_summary.tsv"))
    .write_tsv(autosomal_summary, file.path(sample_dir, "autosomal_summary.tsv"))
    .write_tsv(regression_summary, file.path(sample_dir, "autosomal_regression.tsv"))

    payload <- list(
      sample_name = sample_name,
      sex_summary = sex_summary,
      gaunosome_summary = gauno_summary,
      autosomal_summary = autosomal_summary,
      autosomal_regression = regression_summary
    )
    saveRDS(payload, out_rds)
  }

  chr21 <- autosomal_summary[autosomal_summary$chromosome == 21L, , drop = FALSE]
  chrX <- gauno_summary[gauno_summary$chromosome == "X", , drop = FALSE]
  chrY <- gauno_summary[gauno_summary$chromosome == "Y", , drop = FALSE]
  data.frame(
    sample_name = sample_name,
    status = "ok",
    error = NA_character_,
    predicted_sex = sex_summary$predicted_sex[[1L]],
    read_counts_input_total = .qc_scalar(sample_qc_row, "read_counts_input_total"),
    read_counts_mapped_total = .qc_scalar(sample_qc_row, "read_counts_mapped_total"),
    read_counts_binned_post_sum = .qc_scalar(sample_qc_row, "read_counts_binned_post_sum"),
    fetal_fraction_post = .qc_scalar(sample_qc_row, "fetal_fraction_post"),
    y_unique_ratio_post = .qc_scalar(sample_qc_row, "y_unique_ratio_post"),
    gc_curve_has_valid_bins = .qc_logical(sample_qc_row, "gc_curve_has_valid_bins"),
    gc_curve_has_loess_support = .qc_logical(sample_qc_row, "gc_curve_has_loess_support"),
    nipter_bed_status = .qc_text(sample_qc_row, "nipter_bed_status"),
    z_13 = autosomal_summary$z_score[autosomal_summary$chromosome == 13L][[1L]],
    z_18 = autosomal_summary$z_score[autosomal_summary$chromosome == 18L][[1L]],
    z_21 = if (nrow(chr21)) chr21$z_score[[1L]] else NA_real_,
    ncv_21 = if (nrow(chr21)) chr21$ncv_score[[1L]] else NA_real_,
    rbz_21 = if (nrow(chr21)) chr21$rbz_median[[1L]] else NA_real_,
    z_X = if (nrow(chrX)) chrX$z_score[[1L]] else NA_real_,
    z_Y = if (nrow(chrY)) chrY$z_score[[1L]] else NA_real_,
    ncv_X = if (nrow(chrX)) chrX$ncv_score_selected[[1L]] else NA_real_,
    ncv_Y = if (nrow(chrY)) chrY$ncv_score_selected[[1L]] else NA_real_,
    reg_X = if (nrow(chrX)) chrX$regression_score_selected[[1L]] else NA_real_,
    reg_Y = if (nrow(chrY)) chrY$regression_score_selected[[1L]] else NA_real_,
    out_dir = sample_dir,
    stringsAsFactors = FALSE
  )
}

option_list <- list(
  make_option("--bam-list", type = "character", default = NULL,
              help = "Text file with one BAM/CRAM path per line [required]"),
  make_option("--out-root", type = "character", default = NULL,
              help = "Prediction output root [required]"),
  make_option("--rwcx-bed-dir", type = "character", default = NULL,
              help = "Directory of RWisecondorX BED.gz files. Expected naming contract: <sample_name>.bed.gz unless --rwcx-bed-list is supplied."),
  make_option("--rwcx-bed-list", type = "character", default = NULL,
              help = "Explicit file list with one RWisecondorX BED.gz path per line. Overrides the implicit --rwcx-bed-dir naming contract."),
  make_option("--rwcx-reference", type = "character", default = NULL,
              help = "Native RWisecondorX reference RDS [required]"),
  make_option("--rwcx-binsize", type = "integer", default = 100000L,
              help = "Input RWisecondorX BED binsize [default: %default]"),
  make_option("--rwcx-minrefbins", type = "integer", default = 150L,
              help = "Native RWisecondorX minimum reference bins per target [default: %default]"),
  make_option("--rwcx-maskrepeats", type = "integer", default = 5L,
              help = "Native RWisecondorX iterative distance cutoff passes [default: %default]"),
  make_option("--rwcx-zscore", type = "double", default = 5,
              help = "Native RWisecondorX aberration z-score cutoff [default: %default]"),
  make_option("--rwcx-alpha", type = "double", default = 1e-4,
              help = "Native RWisecondorX CBS alpha [default: %default]"),
  make_option("--rwcx-beta", type = "double", default = NA_real_,
              help = "Native RWisecondorX beta-mode ratio cutoff; leave NA to disable [default: %default]"),
  make_option("--rwcx-optimal-cutoff-sd-multiplier", type = "double", default = 3,
              help = "Native RWisecondorX multiplier applied to the population SD when deriving the optimal distance cutoff [default: %default]"),
  make_option("--rwcx-within-sample-mask-iterations", type = "integer", default = 3L,
              help = "Native RWisecondorX within-sample masking passes [default: %default]"),
  make_option("--rwcx-within-sample-mask-quantile", type = "double", default = 0.99,
              help = "Native RWisecondorX z-quantile used to mask aberrant bins between passes [default: %default]"),
  make_option("--rwcx-cbs-split-min-gap-bp", type = "integer", default = 2000000L,
              help = "Native RWisecondorX minimum NA gap size in bp used to split CBS segments [default: %default]"),
  make_option("--rwcx-segment-zscore-cap", type = "double", default = 1000,
              help = "Native RWisecondorX absolute cap applied to segment z-scores [default: %default]"),
  make_option("--nipter-bed-dir", type = "character", default = NULL,
              help = "Directory of NIPTeR BED.gz files. Expected naming contract: <sample_name>.bed.gz unless --nipter-bed-list is supplied."),
  make_option("--nipter-bed-list", type = "character", default = NULL,
              help = "Explicit file list with one NIPTeR BED.gz path per line. Overrides the implicit --nipter-bed-dir naming contract."),
  make_option("--nipter-reference", type = "character", default = NULL,
              help = "NIPTeR reference-model RDS [required]"),
  make_option("--sample-qc-tsv", type = "character", default = NULL,
              help = "Sample QC TSV from preprocess_cohort.R [required]"),
  make_option("--threads", type = "integer", default = 20L,
              help = "Reserved thread budget for prediction [default: %default]"),
  make_option("--jobs", type = "integer", default = 4L,
              help = "Number of samples to score in parallel per stage [default: %default]"),
  make_option("--nipter-autosomal-source", type = "character", default = "auto",
              help = "BED import source for NIPTeR autosomes: auto, raw, corrected [default: %default]"),
  make_option("--nipter-sex-counts", type = "character", default = "match",
              help = "BED import source for NIPTeR sex counts: match, raw, corrected [default: %default]"),
  make_option("--nipter-chi-cutoff", type = "double", default = 3.5,
              help = "Sample-level chi correction cutoff for NIPTeR scoring [default: %default]"),
  make_option("--autosomal-focus", type = "character", default = "13,18,21",
              help = "Comma-separated autosomal chromosomes to summarise [default: %default]"),
  make_option("--regression-n-models", type = "integer", default = 4L,
              help = "NIPTeR autosomal RBZ models per chromosome [default: %default]"),
  make_option("--regression-n-predictors", type = "integer", default = 4L,
              help = "NIPTeR autosomal RBZ max predictors per model [default: %default]"),
  make_option("--regression-train-fraction", type = "double", default = 1,
              help = "Fraction of matched controls used to fit autosomal RBZ models; use 1 to fit and score on all controls [default: %default]"),
  make_option("--regression-exclude-chromosomes", type = "character", default = "13,18,21",
              help = "Comma-separated chromosomes excluded from autosomal RBZ predictors [default: %default]"),
  make_option("--seed", type = "integer", default = 1995L,
              help = "Random seed for reproducible prediction-side splits [default: %default]"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Overwrite existing prediction outputs [default: %default]")
)

parser <- OptionParser(
  usage = paste(
    "%prog --bam-list cohort.txt --out-root predict/ \\",
    "  --rwcx-bed-dir rwcx_beds --rwcx-reference rwisecondorx_ref.rds \\",
    "  --nipter-bed-dir nipter_beds --nipter-reference nipter_reference_model.rds \\",
    "  --sample-qc-tsv sample_qc/sample_qc.tsv",
    sep = "\n"
  ),
  option_list = option_list
)
opts <- parse_args(parser)

required_paths <- c("bam-list", "rwcx-reference", "nipter-reference", "sample-qc-tsv")
for (nm in required_paths) {
  path <- opts[[nm]]
  if (is.null(path) || !file.exists(path)) {
    stop("--", nm, " is required and must exist.", call. = FALSE)
  }
}
if (is.null(opts$`out-root`) || !nzchar(opts$`out-root`)) {
  stop("--out-root is required.", call. = FALSE)
}
if (!xor(is.null(opts$`rwcx-bed-dir`), is.null(opts$`rwcx-bed-list`))) {
  stop("Provide exactly one of --rwcx-bed-dir or --rwcx-bed-list.", call. = FALSE)
}
if (!xor(is.null(opts$`nipter-bed-dir`), is.null(opts$`nipter-bed-list`))) {
  stop("Provide exactly one of --nipter-bed-dir or --nipter-bed-list.", call. = FALSE)
}
for (nm in c("rwcx-bed-list", "nipter-bed-list")) {
  path <- opts[[nm]]
  if (!is.null(path) && !file.exists(path)) {
    stop("--", nm, " must exist: ", path, call. = FALSE)
  }
}
for (nm in c("rwcx-bed-dir", "nipter-bed-dir")) {
  path <- opts[[nm]]
  if (!is.null(path) && !dir.exists(path)) {
    stop("--", nm, " must exist: ", path, call. = FALSE)
  }
}

autosomal_source <- match.arg(tolower(opts$`nipter-autosomal-source`), c("auto", "raw", "corrected"))
sex_counts <- match.arg(tolower(opts$`nipter-sex-counts`), c("match", "raw", "corrected"))
focus_chromosomes <- as.integer(strsplit(opts$`autosomal-focus`, ",", fixed = TRUE)[[1L]])
focus_chromosomes <- unique(focus_chromosomes[is.finite(focus_chromosomes)])
regression_exclude <- as.integer(strsplit(opts$`regression-exclude-chromosomes`, ",", fixed = TRUE)[[1L]])
regression_exclude <- unique(regression_exclude[is.finite(regression_exclude)])
if (!length(focus_chromosomes)) {
  stop("--autosomal-focus resolved to an empty set.", call. = FALSE)
}

library(RWisecondorX)

bam_paths <- .read_file_list(opts$`bam-list`)
sample_names <- vapply(bam_paths, .sample_stem, character(1L))
sample_qc <- utils::read.delim(
  opts$`sample-qc-tsv`,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!("sample_name" %in% names(sample_qc))) {
  stop("--sample-qc-tsv must contain a sample_name column.", call. = FALSE)
}

rwx_ref <- readRDS(opts$`rwcx-reference`)
nipter_ref <- readRDS(opts$`nipter-reference`)

out_root <- .ensure_dir(opts$`out-root`)
dirs <- list(
  rwisecondorx = .ensure_dir(file.path(out_root, "rwisecondorx")),
  nipter = .ensure_dir(file.path(out_root, "nipter"))
)

worker_threads <- .per_worker_threads(opts$threads, opts$jobs)
rwcx_beds <- if (!is.null(opts$`rwcx-bed-list`)) {
  .match_bed_paths(sample_names, .read_file_list(opts$`rwcx-bed-list`), "--rwcx-bed-list")
} else {
  normalizePath(
    file.path(opts$`rwcx-bed-dir`, paste0(sample_names, ".bed.gz")),
    winslash = "/",
    mustWork = FALSE
  )
}
nipter_beds <- if (!is.null(opts$`nipter-bed-list`)) {
  .match_bed_paths(sample_names, .read_file_list(opts$`nipter-bed-list`), "--nipter-bed-list")
} else {
  normalizePath(
    file.path(opts$`nipter-bed-dir`, paste0(sample_names, ".bed.gz")),
    winslash = "/",
    mustWork = FALSE
  )
}
sample_manifest <- data.frame(
  sample_name = sample_names,
  bam = bam_paths,
  rwcx_bed = rwcx_beds,
  nipter_bed = nipter_beds,
  stringsAsFactors = FALSE
)

cat(sprintf("Predicting %d samples with %d parallel job(s) and %d worker thread(s).\n",
            nrow(sample_manifest), as.integer(opts$jobs), worker_threads))

rwx_rows <- .run_stage_samples(
  split(sample_manifest, seq_len(nrow(sample_manifest))),
  opts$jobs,
  function(i, item) {
    row <- item[[1L]]
    sample_name <- row$sample_name[[1L]]
    bed_path <- row$rwcx_bed[[1L]]
    sample_qc_row <- .match_sample_qc_row(sample_qc, sample_name)
    .log_sample(i, nrow(sample_manifest), "RWisecondorX predict", sample_name)
    if (!file.exists(bed_path)) {
      return(data.frame(
        sample_name = sample_name,
        status = "missing_input",
        error = paste("Missing RWisecondorX BED:", bed_path),
        predicted_gender = NA_character_,
        n_reads = NA_real_,
        n_aberrations = NA_real_,
        aberrations = NA_character_,
        z_13 = NA_real_,
        z_18 = NA_real_,
        z_21 = NA_real_,
        z_X = NA_real_,
        read_counts_input_total = .qc_scalar(sample_qc_row, "read_counts_input_total"),
        read_counts_mapped_total = .qc_scalar(sample_qc_row, "read_counts_mapped_total"),
        read_counts_binned_post_sum = .qc_scalar(sample_qc_row, "read_counts_binned_post_sum"),
        fetal_fraction_post = .qc_scalar(sample_qc_row, "fetal_fraction_post"),
        y_unique_ratio_post = .qc_scalar(sample_qc_row, "y_unique_ratio_post"),
        gc_curve_has_loess_support = .qc_logical(sample_qc_row, "gc_curve_has_loess_support"),
        nipter_bed_status = .qc_text(sample_qc_row, "nipter_bed_status"),
        out_dir = file.path(dirs$rwisecondorx, sample_name),
        stringsAsFactors = FALSE
      ))
    }
    tryCatch(
      .predict_rwisecondorx_one(
        sample_name = sample_name,
        bed_path = bed_path,
        reference = rwx_ref,
        out_dir = dirs$rwisecondorx,
        sample_binsize = as.integer(opts$`rwcx-binsize`),
        minrefbins = as.integer(opts$`rwcx-minrefbins`),
        maskrepeats = as.integer(opts$`rwcx-maskrepeats`),
        zscore = as.numeric(opts$`rwcx-zscore`),
        alpha = as.numeric(opts$`rwcx-alpha`),
        beta = as.numeric(opts$`rwcx-beta`),
        optimal_cutoff_sd_multiplier = as.numeric(opts$`rwcx-optimal-cutoff-sd-multiplier`),
        within_sample_mask_iterations = as.integer(opts$`rwcx-within-sample-mask-iterations`),
        within_sample_mask_quantile = as.numeric(opts$`rwcx-within-sample-mask-quantile`),
        cbs_split_min_gap_bp = as.integer(opts$`rwcx-cbs-split-min-gap-bp`),
        segment_zscore_cap = as.numeric(opts$`rwcx-segment-zscore-cap`),
        cpus = worker_threads,
        overwrite = isTRUE(opts$overwrite),
        seed = opts$seed,
        sample_qc_row = sample_qc_row
      ),
      error = function(e) {
        data.frame(
          sample_name = sample_name,
          status = "failed",
          error = conditionMessage(e),
          predicted_gender = NA_character_,
          n_reads = NA_real_,
          n_aberrations = NA_real_,
          aberrations = NA_character_,
          z_13 = NA_real_,
          z_18 = NA_real_,
          z_21 = NA_real_,
          z_X = NA_real_,
          read_counts_input_total = .qc_scalar(sample_qc_row, "read_counts_input_total"),
          read_counts_mapped_total = .qc_scalar(sample_qc_row, "read_counts_mapped_total"),
          read_counts_binned_post_sum = .qc_scalar(sample_qc_row, "read_counts_binned_post_sum"),
          fetal_fraction_post = .qc_scalar(sample_qc_row, "fetal_fraction_post"),
          y_unique_ratio_post = .qc_scalar(sample_qc_row, "y_unique_ratio_post"),
          gc_curve_has_loess_support = .qc_logical(sample_qc_row, "gc_curve_has_loess_support"),
          nipter_bed_status = .qc_text(sample_qc_row, "nipter_bed_status"),
          out_dir = file.path(dirs$rwisecondorx, sample_name),
          stringsAsFactors = FALSE
        )
      }
    )
  }
)
rwx_summary <- do.call(rbind, rwx_rows)
.write_tsv(rwx_summary, file.path(out_root, "rwisecondorx_summary.tsv"))

nipter_rows <- .run_stage_samples(
  split(sample_manifest, seq_len(nrow(sample_manifest))),
  opts$jobs,
  function(i, item) {
    row <- item[[1L]]
    sample_name <- row$sample_name[[1L]]
    bed_path <- row$nipter_bed[[1L]]
    sample_qc_row <- .match_sample_qc_row(sample_qc, sample_name)
    .log_sample(i, nrow(sample_manifest), "NIPTeR predict", sample_name)
    if (!file.exists(bed_path)) {
      return(data.frame(
        sample_name = sample_name,
        status = "missing_input",
        error = paste("Missing NIPTeR BED:", bed_path),
        predicted_sex = NA_character_,
        read_counts_input_total = .qc_scalar(sample_qc_row, "read_counts_input_total"),
        read_counts_mapped_total = .qc_scalar(sample_qc_row, "read_counts_mapped_total"),
        read_counts_binned_post_sum = .qc_scalar(sample_qc_row, "read_counts_binned_post_sum"),
        fetal_fraction_post = .qc_scalar(sample_qc_row, "fetal_fraction_post"),
        y_unique_ratio_post = .qc_scalar(sample_qc_row, "y_unique_ratio_post"),
        gc_curve_has_valid_bins = .qc_logical(sample_qc_row, "gc_curve_has_valid_bins"),
        gc_curve_has_loess_support = .qc_logical(sample_qc_row, "gc_curve_has_loess_support"),
        nipter_bed_status = .qc_text(sample_qc_row, "nipter_bed_status"),
        z_13 = NA_real_,
        z_18 = NA_real_,
        z_21 = NA_real_,
        ncv_21 = NA_real_,
        rbz_21 = NA_real_,
        z_X = NA_real_,
        z_Y = NA_real_,
        ncv_X = NA_real_,
        ncv_Y = NA_real_,
        reg_X = NA_real_,
        reg_Y = NA_real_,
        out_dir = file.path(dirs$nipter, sample_name),
        stringsAsFactors = FALSE
      ))
    }
    tryCatch(
      .predict_nipter_one(
        sample_name = sample_name,
        bed_path = bed_path,
        reference = nipter_ref,
        sample_qc_row = sample_qc_row,
        out_dir = dirs$nipter,
        autosomal_source = autosomal_source,
        sex_counts = sex_counts,
        chi_cutoff = opts$`nipter-chi-cutoff`,
        focus_chromosomes = focus_chromosomes,
        regression_n_models = opts$`regression-n-models`,
        regression_n_predictors = opts$`regression-n-predictors`,
        regression_train_fraction = opts$`regression-train-fraction`,
        regression_exclude_chromosomes = regression_exclude,
        regression_seed = opts$seed,
        overwrite = isTRUE(opts$overwrite)
      ),
      error = function(e) {
        data.frame(
          sample_name = sample_name,
          status = "failed",
          error = conditionMessage(e),
          predicted_sex = NA_character_,
          read_counts_input_total = .qc_scalar(sample_qc_row, "read_counts_input_total"),
          read_counts_mapped_total = .qc_scalar(sample_qc_row, "read_counts_mapped_total"),
          read_counts_binned_post_sum = .qc_scalar(sample_qc_row, "read_counts_binned_post_sum"),
          fetal_fraction_post = .qc_scalar(sample_qc_row, "fetal_fraction_post"),
          y_unique_ratio_post = .qc_scalar(sample_qc_row, "y_unique_ratio_post"),
          gc_curve_has_valid_bins = .qc_logical(sample_qc_row, "gc_curve_has_valid_bins"),
          gc_curve_has_loess_support = .qc_logical(sample_qc_row, "gc_curve_has_loess_support"),
          nipter_bed_status = .qc_text(sample_qc_row, "nipter_bed_status"),
          z_13 = NA_real_,
          z_18 = NA_real_,
          z_21 = NA_real_,
          ncv_21 = NA_real_,
          rbz_21 = NA_real_,
          z_X = NA_real_,
          z_Y = NA_real_,
          ncv_X = NA_real_,
          ncv_Y = NA_real_,
          reg_X = NA_real_,
          reg_Y = NA_real_,
          out_dir = file.path(dirs$nipter, sample_name),
          stringsAsFactors = FALSE
        )
      }
    )
  }
)
nipter_summary <- do.call(rbind, nipter_rows)
.write_tsv(nipter_summary, file.path(out_root, "nipter_summary.tsv"))

cohort_summary <- Reduce(
  function(x, y) merge(x, y, by = "sample_name", all = TRUE, suffixes = c("_rwisecondorx", "_nipter")),
  list(
    sample_manifest[, c("sample_name", "bam"), drop = FALSE],
    rwx_summary,
    nipter_summary
  )
)
.write_tsv(cohort_summary, file.path(out_root, "cohort_summary.tsv"))

cat(sprintf("RWisecondorX summary: %s\n", file.path(out_root, "rwisecondorx_summary.tsv")))
cat(sprintf("NIPTeR summary: %s\n", file.path(out_root, "nipter_summary.tsv")))
cat(sprintf("Combined cohort summary: %s\n", file.path(out_root, "cohort_summary.tsv")))
