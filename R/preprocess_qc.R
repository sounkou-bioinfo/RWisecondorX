.cv_percent <- function(x) {
  vals <- suppressWarnings(as.numeric(x))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2L) {
    return(NA_real_)
  }
  mu <- mean(vals)
  if (!is.finite(mu) || abs(mu) <= .Machine$double.eps) {
    return(NA_real_)
  }
  100 * stats::sd(vals) / mu
}

.expand_gc_vector <- function(gc_table, chrom_names, n_bins) {
  stopifnot(
    is.list(gc_table),
    is.character(chrom_names),
    length(n_bins) == 1L,
    is.numeric(n_bins),
    n_bins >= 1L
  )
  unlist(lapply(chrom_names, function(chr) {
    gc <- gc_table[[chr]]
    if (length(gc) < n_bins) {
      gc <- c(gc, rep(NA_real_, n_bins - length(gc)))
    } else if (length(gc) > n_bins) {
      gc <- gc[seq_len(n_bins)]
    }
    gc
  }))
}

.nipter_preprocess_qc <- function(sample,
                                  corrected = NULL,
                                  gc_table = NULL,
                                  fasta = NULL,
                                  include_sex = FALSE,
                                  binsize = 50000L,
                                  con = NULL) {
  stopifnot(!is.null(gc_table) || !is.null(fasta))
  stopifnot(is.logical(include_sex), length(include_sex) == 1L)

  gc_tbl <- .resolve_gc_table(gc_table, fasta, binsize, con)
  raw_auto <- autosomal_matrix(sample)
  corrected_auto <- if (is.null(corrected)) raw_auto else autosomal_matrix(corrected)
  n_bins <- ncol(raw_auto)

  chromosomes <- rep(as.character(seq_len(22L)), each = n_bins)
  bin_index <- rep(seq_len(n_bins), times = 22L)
  bins_per_chromosome <- rep(
    vapply(as.character(seq_len(22L)), function(chr) length(gc_tbl[[chr]]), integer(1L)),
    each = n_bins
  )
  gc_auto <- .expand_gc_vector(gc_tbl, chrom_names = as.character(seq_len(22L)), n_bins = n_bins)
  raw_flat <- as.numeric(t(raw_auto))
  corrected_flat <- as.numeric(t(corrected_auto))
  in_chromosome_range <- bin_index <= bins_per_chromosome
  gc_is_finite <- is.finite(gc_auto)
  gc_is_positive <- gc_is_finite & gc_auto > 0
  invalid_raw_nonfinite <- !is.finite(raw_flat)
  invalid_raw_nonpositive <- is.finite(raw_flat) & raw_flat <= 0
  corrected_count_is_finite <- is.finite(corrected_flat)
  valid_for_gc_correction_fit <- in_chromosome_range & gc_is_positive & !invalid_raw_nonfinite & !invalid_raw_nonpositive
  valid_for_post_gc_summary <- valid_for_gc_correction_fit & corrected_count_is_finite
  valid_for_corrected_ratio_genome_plot <- in_chromosome_range & corrected_count_is_finite

  fit_exclusion_reason <- rep("valid", length(valid_for_gc_correction_fit))
  fit_exclusion_reason[!in_chromosome_range] <- "out_of_chromosome_range"
  fit_exclusion_reason[
    in_chromosome_range &
      !gc_is_finite &
      fit_exclusion_reason == "valid"
  ] <- "gc_missing_or_nonfinite"
  fit_exclusion_reason[
    in_chromosome_range &
      gc_is_finite &
      !gc_is_positive &
      fit_exclusion_reason == "valid"
  ] <- "gc_nonpositive"
  fit_exclusion_reason[
    in_chromosome_range &
      gc_is_positive &
      invalid_raw_nonfinite &
      fit_exclusion_reason == "valid"
  ] <- "raw_count_nonfinite"
  fit_exclusion_reason[
    in_chromosome_range &
      gc_is_positive &
      !invalid_raw_nonfinite &
      invalid_raw_nonpositive &
      fit_exclusion_reason == "valid"
  ] <- "raw_count_nonpositive"

  valid_n <- sum(valid_for_gc_correction_fit)
  total_n <- length(valid_for_gc_correction_fit)
  invalid <- !valid_for_gc_correction_fit

  out_metrics <- list(
    gc_loess_valid_bins = valid_n,
    gc_loess_total_bins = total_n,
    gc_loess_invalid_bins = sum(invalid),
    gc_loess_valid_bin_fraction_pct = if (total_n > 0L) 100 * valid_n / total_n else NA_real_,
    gc_loess_invalid_out_of_chromosome_range = sum(!in_chromosome_range),
    gc_loess_invalid_gc_missing_or_nonfinite = sum(invalid & in_chromosome_range & !gc_is_finite),
    gc_loess_invalid_gc_nonpositive = sum(invalid & in_chromosome_range & gc_is_finite & !gc_is_positive),
    gc_loess_invalid_raw_nonfinite = sum(invalid & invalid_raw_nonfinite),
    gc_loess_invalid_raw_nonpositive = sum(invalid & invalid_raw_nonpositive),
    gc_loess_invalid_corrected_nonfinite = sum(valid_for_gc_correction_fit & !corrected_count_is_finite),
    gc_curve_has_valid_bins = isTRUE(valid_n > 0L),
    gc_curve_has_loess_support = isTRUE(valid_n >= 10L && length(unique(gc_auto[valid_for_gc_correction_fit])) >= 5L),
    nipter_autosomal_bin_cv_pre_gc_correction = NA_real_,
    nipter_autosomal_bin_cv_post_gc_correction = NA_real_,
    nipter_gc_correction_bin_scale_mean = NA_real_,
    nipter_gc_correction_bin_scale_sd = NA_real_,
    nipter_gc_correction_bin_scale_cv = NA_real_,
    nipter_gc_correlation_pre_gc_correction = NA_real_,
    nipter_gc_correlation_post_gc_correction = NA_real_
  )

  curve_df <- data.frame(
    chrom = chromosomes,
    start = as.integer((bin_index - 1L) * binsize),
    end = as.integer(bin_index * binsize),
    bin = as.integer(bin_index),
    gc = as.numeric(gc_auto),
    raw_count = as.numeric(raw_flat),
    corrected_count = as.numeric(corrected_flat),
    in_chromosome_range = as.logical(in_chromosome_range),
    gc_is_finite = as.logical(gc_is_finite),
    gc_is_positive = as.logical(gc_is_positive),
    raw_count_is_finite = as.logical(!invalid_raw_nonfinite),
    raw_count_is_positive = as.logical(!invalid_raw_nonpositive),
    corrected_count_is_finite = as.logical(corrected_count_is_finite),
    valid_for_gc_correction_fit = as.logical(valid_for_gc_correction_fit),
    valid_for_corrected_ratio_genome_plot = as.logical(valid_for_corrected_ratio_genome_plot),
    fit_exclusion_reason = as.character(fit_exclusion_reason),
    raw_ratio_to_fit_eligible_median = NA_real_,
    corrected_ratio_to_fit_eligible_median = NA_real_,
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_sex)) {
    raw_sex <- sex_matrix(sample)
    corrected_sex <- if (is.null(corrected)) raw_sex else sex_matrix(corrected)
    sex_n_bins <- ncol(raw_sex)
    sex_chroms <- c("X", "Y")
    sex_bins_per_chromosome <- rep(
      vapply(sex_chroms, function(chr) length(gc_tbl[[chr]]), integer(1L)),
      each = sex_n_bins
    )
    sex_bin_index <- rep(seq_len(sex_n_bins), times = length(sex_chroms))
    sex_gc <- .expand_gc_vector(gc_tbl, chrom_names = sex_chroms, n_bins = sex_n_bins)
    raw_sex_flat <- as.numeric(t(raw_sex))
    corrected_sex_flat <- as.numeric(t(corrected_sex))
    sex_in_chromosome_range <- sex_bin_index <= sex_bins_per_chromosome
    sex_gc_is_finite <- is.finite(sex_gc)
    sex_gc_is_positive <- sex_gc_is_finite & sex_gc > 0
    sex_invalid_raw_nonfinite <- !is.finite(raw_sex_flat)
    sex_invalid_raw_nonpositive <- is.finite(raw_sex_flat) & raw_sex_flat <= 0
    sex_corrected_count_is_finite <- is.finite(corrected_sex_flat)
    sex_valid_for_corrected_ratio_genome_plot <- sex_in_chromosome_range & sex_corrected_count_is_finite
    sex_fit_exclusion_reason <- rep("sex_chromosome_not_fit_eligible", length(sex_bin_index))
    sex_fit_exclusion_reason[!sex_in_chromosome_range] <- "out_of_chromosome_range"

    sex_curve_df <- data.frame(
      chrom = rep(sex_chroms, each = sex_n_bins),
      start = as.integer((sex_bin_index - 1L) * binsize),
      end = as.integer(sex_bin_index * binsize),
      bin = as.integer(sex_bin_index),
      gc = as.numeric(sex_gc),
      raw_count = as.numeric(raw_sex_flat),
      corrected_count = as.numeric(corrected_sex_flat),
      in_chromosome_range = as.logical(sex_in_chromosome_range),
      gc_is_finite = as.logical(sex_gc_is_finite),
      gc_is_positive = as.logical(sex_gc_is_positive),
      raw_count_is_finite = as.logical(!sex_invalid_raw_nonfinite),
      raw_count_is_positive = as.logical(!sex_invalid_raw_nonpositive),
      corrected_count_is_finite = as.logical(sex_corrected_count_is_finite),
      valid_for_gc_correction_fit = FALSE,
      valid_for_corrected_ratio_genome_plot = as.logical(sex_valid_for_corrected_ratio_genome_plot),
      fit_exclusion_reason = as.character(sex_fit_exclusion_reason),
      raw_ratio_to_fit_eligible_median = NA_real_,
      corrected_ratio_to_fit_eligible_median = NA_real_,
      stringsAsFactors = FALSE
    )
    curve_df <- rbind(curve_df, sex_curve_df)
  }

  if (any(valid_for_gc_correction_fit)) {
    raw_med <- stats::median(raw_flat[valid_for_gc_correction_fit])
    corrected_med <- stats::median(corrected_flat[valid_for_post_gc_summary])
    auto_fit_idx <- curve_df$valid_for_gc_correction_fit
    auto_plot_idx <- curve_df$chrom %in% as.character(seq_len(22L)) &
      curve_df$valid_for_corrected_ratio_genome_plot
    sex_plot_idx <- curve_df$chrom %in% c("X", "Y") &
      curve_df$valid_for_corrected_ratio_genome_plot

    if (is.finite(raw_med) && raw_med > 0) {
      curve_df$raw_ratio_to_fit_eligible_median[auto_fit_idx] <-
        curve_df$raw_count[auto_fit_idx] / raw_med
    }
    if (is.finite(corrected_med) && corrected_med > 0) {
      curve_df$corrected_ratio_to_fit_eligible_median[auto_plot_idx] <-
        curve_df$corrected_count[auto_plot_idx] / corrected_med
      if (isTRUE(include_sex)) {
        curve_df$corrected_ratio_to_fit_eligible_median[sex_plot_idx] <-
          curve_df$corrected_count[sex_plot_idx] / corrected_med
      }
    }

    out_metrics$nipter_autosomal_bin_cv_pre_gc_correction <- .cv_percent(raw_flat[valid_for_gc_correction_fit])
    out_metrics$nipter_autosomal_bin_cv_post_gc_correction <- .cv_percent(corrected_flat[valid_for_post_gc_summary])

    corrected_ratio_valid <- curve_df$corrected_ratio_to_fit_eligible_median[auto_fit_idx & curve_df$corrected_count_is_finite]
    corrected_ratio_valid <- corrected_ratio_valid[is.finite(corrected_ratio_valid)]
    if (length(corrected_ratio_valid) >= 2L) {
      out_metrics$nipter_gc_correction_bin_scale_mean <- mean(corrected_ratio_valid)
      out_metrics$nipter_gc_correction_bin_scale_sd <- stats::sd(corrected_ratio_valid)
      out_metrics$nipter_gc_correction_bin_scale_cv <- .cv_percent(corrected_ratio_valid)
    }

    raw_ratio_valid <- curve_df$raw_ratio_to_fit_eligible_median[auto_fit_idx]
    raw_ratio_valid <- raw_ratio_valid[is.finite(raw_ratio_valid)]
    gc_valid_pre <- gc_auto[valid_for_gc_correction_fit]
    gc_valid_pre <- gc_valid_pre[is.finite(gc_valid_pre)]
    if (length(raw_ratio_valid) == length(gc_valid_pre) && length(gc_valid_pre) >= 2L) {
      out_metrics$nipter_gc_correlation_pre_gc_correction <- suppressWarnings(stats::cor(gc_valid_pre, raw_ratio_valid))
    }
    gc_valid_post <- gc_auto[valid_for_post_gc_summary]
    gc_valid_post <- gc_valid_post[is.finite(gc_valid_post)]
    if (length(corrected_ratio_valid) == length(gc_valid_post) && length(gc_valid_post) >= 2L) {
      out_metrics$nipter_gc_correlation_post_gc_correction <- suppressWarnings(stats::cor(gc_valid_post, corrected_ratio_valid))
    }
  }

  list(
    metrics = out_metrics,
    curve_data = curve_df
  )
}

.write_nipter_gc_curve_data_bgz <- function(curve_data,
                                            out_bgz,
                                            con = NULL,
                                            metadata = list(artifact_type = "nipter_gc_curve_data")) {
  stopifnot(is.data.frame(curve_data))
  stopifnot(all(c("chrom", "start", "end") %in% names(curve_data)))
  stopifnot(is.character(out_bgz), length(out_bgz) == 1L, nzchar(out_bgz))

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)

  .write_tabix_body(curve_data, tmp, metadata = metadata)
  Rduckhts::rduckhts_bgzip(con, tmp, out_bgz, threads = 1L, keep = TRUE, overwrite = TRUE)
  Rduckhts::rduckhts_tabix_index(con, out_bgz, preset = "bed", threads = 1L, comment_char = "#")

  invisible(out_bgz)
}

.write_nipter_gc_curve_plot <- function(curve_data,
                                        sample_name,
                                        out_png,
                                        theme = c("minimal", "light", "bw", "classic"),
                                        base_size = 10,
                                        loess_span = 0.75) {
  .require_ggplot2(".write_nipter_gc_curve_plot()")
  stopifnot(is.data.frame(curve_data))
  stopifnot(is.character(sample_name), length(sample_name) == 1L, nzchar(sample_name))
  stopifnot(is.character(out_png), length(out_png) == 1L, nzchar(out_png))
  stopifnot(is.numeric(loess_span), length(loess_span) == 1L, is.finite(loess_span),
            loess_span > 0, loess_span <= 1)

  theme <- match.arg(theme)
  valid_df <- curve_data[curve_data$valid_for_gc_correction_fit, , drop = FALSE]
  if (!nrow(valid_df)) {
    return(invisible(NULL))
  }

  valid_df$gc_bucket <- round(valid_df$gc, 3)
  curve_summary <- stats::aggregate(
    cbind(raw_ratio_to_fit_eligible_median, corrected_ratio_to_fit_eligible_median) ~ gc_bucket,
    data = valid_df[, c(
      "gc_bucket",
      "raw_ratio_to_fit_eligible_median",
      "corrected_ratio_to_fit_eligible_median"
    ), drop = FALSE],
    FUN = function(x) stats::median(x[is.finite(x)], na.rm = TRUE)
  )

  point_df <- rbind(
    data.frame(
      stage = "Raw",
      gc = valid_df$gc,
      ratio = valid_df$raw_ratio_to_fit_eligible_median
    ),
    data.frame(
      stage = "Corrected",
      gc = valid_df$gc,
      ratio = valid_df$corrected_ratio_to_fit_eligible_median
    )
  )
  point_df <- point_df[is.finite(point_df$gc) & is.finite(point_df$ratio), , drop = FALSE]

  line_df <- rbind(
    data.frame(stage = "Raw", gc = curve_summary$gc_bucket, ratio = curve_summary$raw_ratio_to_fit_eligible_median),
    data.frame(stage = "Corrected", gc = curve_summary$gc_bucket, ratio = curve_summary$corrected_ratio_to_fit_eligible_median)
  )
  line_df <- line_df[is.finite(line_df$gc) & is.finite(line_df$ratio), , drop = FALSE]

  loess_df <- do.call(
    rbind,
    lapply(split(point_df, point_df$stage), function(df) {
      df <- df[is.finite(df$gc) & is.finite(df$ratio), , drop = FALSE]
      if (nrow(df) < 10L || length(unique(df$gc)) < 5L) {
        return(NULL)
      }
      fit <- tryCatch(
        stats::loess(ratio ~ gc, data = df, span = loess_span),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(NULL)
      }
      gc_grid <- seq(min(df$gc), max(df$gc), length.out = 250L)
      fitted_ratio <- suppressWarnings(stats::predict(fit, newdata = data.frame(gc = gc_grid)))
      out <- data.frame(
        stage = unique(df$stage)[1L],
        gc = gc_grid,
        ratio = fitted_ratio
      )
      out[is.finite(out$gc) & is.finite(out$ratio), , drop = FALSE]
    })
  )

  stage_palette <- c(Raw = "#c06c00", Corrected = "#006d77")
  p <- ggplot2::ggplot(point_df, ggplot2::aes(x = gc, y = ratio)) +
    ggplot2::geom_point(
      ggplot2::aes(color = stage),
      alpha = 0.10,
      size = 0.35,
      stroke = 0
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = gc, y = ratio),
      inherit.aes = FALSE,
      color = "#6c757d",
      linewidth = 0.55,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(
      data = loess_df,
      ggplot2::aes(x = gc, y = ratio, color = stage),
      linewidth = 0.95,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~ stage, ncol = 1, scales = "free_y") +
    ggplot2::scale_color_manual(values = stage_palette, guide = "none") +
    ggplot2::labs(
      title = paste0(sample_name, " GC Curve"),
      subtitle = paste0(
        "Autosomal bins used for GC fit: ", sum(curve_data$valid_for_gc_correction_fit),
        " | dashed: GC-bucket median | solid: LOESS"
      ),
      x = "GC fraction",
      y = "Bin ratio to autosomal median"
    ) +
    .nipter_plot_theme(theme, base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = 7.5,
    height = 7.0,
    dpi = 180
  )

  invisible(out_png)
}

.write_nipter_corrected_bin_ratio_genome_plot <- function(curve_data,
                                                          sample_name,
                                                          out_png,
                                                          theme = c("minimal", "light", "bw", "classic"),
                                                          base_size = 10,
                                                          running_median_bins = 51L) {
  .require_ggplot2(".write_nipter_corrected_bin_ratio_genome_plot()")
  stopifnot(is.data.frame(curve_data))
  stopifnot(is.character(sample_name), length(sample_name) == 1L, nzchar(sample_name))
  stopifnot(is.character(out_png), length(out_png) == 1L, nzchar(out_png))

  theme <- match.arg(theme)
  running_median_bins <- as.integer(running_median_bins[[1L]])
  if (!is.finite(running_median_bins) || running_median_bins < 3L) {
    stop("running_median_bins must be an integer >= 3.", call. = FALSE)
  }
  if ((running_median_bins %% 2L) == 0L) {
    running_median_bins <- running_median_bins + 1L
  }

  chrom_levels <- c(as.character(seq_len(22L)), "X", "Y")
  present_levels <- chrom_levels[chrom_levels %in% unique(as.character(curve_data$chrom))]
  df <- curve_data[
    curve_data$chrom %in% present_levels &
      curve_data$in_chromosome_range &
      is.finite(curve_data$corrected_ratio_to_fit_eligible_median) &
      curve_data$corrected_ratio_to_fit_eligible_median > 0,
    c("chrom", "bin", "start", "end", "corrected_ratio_to_fit_eligible_median"),
    drop = FALSE
  ]
  if (!nrow(df)) {
    return(invisible(NULL))
  }

  df$chrom_index <- match(df$chrom, chrom_levels)
  df <- df[order(df$chrom_index, df$bin), , drop = FALSE]
  has_sex <- any(df$chrom %in% c("X", "Y"))
  chrom_sizes <- stats::aggregate(
    bin ~ chrom_index + chrom,
    data = df,
    FUN = max
  )
  chrom_sizes <- chrom_sizes[order(chrom_sizes$chrom_index), , drop = FALSE]
  chrom_offsets <- c(0L, cumsum(head(chrom_sizes$bin, -1L)))
  names(chrom_offsets) <- chrom_sizes$chrom
  df$genome_bin_index <- chrom_offsets[df$chrom] + df$bin
  df$chrom_parity <- ifelse(
    df$chrom %in% c("X", "Y"),
    "sex",
    ifelse((df$chrom_index %% 2L) == 0L, "even", "odd")
  )
  df$log2_corrected_ratio <- log2(df$corrected_ratio_to_fit_eligible_median)

  smooth_k <- min(running_median_bins, nrow(df))
  if ((smooth_k %% 2L) == 0L) {
    smooth_k <- max(3L, smooth_k - 1L)
  }
  smooth_df <- do.call(
    rbind,
    lapply(split(df, df$chrom), function(chr_df) {
      chr_df <- chr_df[order(chr_df$bin), , drop = FALSE]
      chr_k <- min(smooth_k, nrow(chr_df))
      if (chr_k < 3L) {
        return(NULL)
      }
      if ((chr_k %% 2L) == 0L) {
        chr_k <- chr_k - 1L
      }
      if (chr_k < 3L) {
        return(NULL)
      }
      data.frame(
        chrom = chr_df$chrom,
        genome_bin_index = chr_df$genome_bin_index,
        log2_corrected_ratio = stats::runmed(
          chr_df$log2_corrected_ratio,
          k = chr_k,
          endrule = "median"
        )
      )
    })
  )
  if (is.null(smooth_df) || !nrow(smooth_df)) {
    smooth_df <- NULL
  }

  boundary_df <- data.frame(
    boundary = cumsum(chrom_sizes$bin)
  )
  label_df <- data.frame(
    chrom = chrom_sizes$chrom,
    center = chrom_offsets + (chrom_sizes$bin / 2)
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = genome_bin_index,
      y = log2_corrected_ratio
    )
  ) +
    ggplot2::geom_line(
      ggplot2::aes(color = chrom_parity, group = chrom),
      linewidth = 0.20,
      alpha = 0.55,
      lineend = "round"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "#c1121f",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      data = boundary_df,
      ggplot2::aes(xintercept = boundary),
      color = "#dee2e6",
      linewidth = 0.30
    ) +
    ggplot2::scale_color_manual(
      values = c(odd = "#4f5d75", even = "#8d99ae", sex = "#c1121f"),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = label_df$center,
      labels = label_df$chrom,
      expand = c(0.002, 0)
    ) +
    ggplot2::labs(
      title = paste0(sample_name, " Corrected Bin Log2 Ratios"),
      subtitle = paste0(
        if (has_sex) {
          "Post-GC corrected genome-wide bin trace (1-22,X,Y) | out-of-range bins excluded"
        } else {
          "Post-GC corrected autosomal bin trace | out-of-range bins excluded"
        },
        if (!is.null(smooth_df)) paste0(" | running median k=", smooth_k) else ""
      ),
      x = "Chromosome",
      y = "log2(corrected bin ratio to fit-eligible median)"
    ) +
    .nipter_plot_theme(theme, base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold")
    )

  if (!is.null(smooth_df)) {
    p <- p + ggplot2::geom_line(
      data = smooth_df,
      ggplot2::aes(
        x = genome_bin_index,
        y = log2_corrected_ratio,
        group = chrom
      ),
      inherit.aes = FALSE,
      color = "#212529",
      linewidth = 0.45
    )
  }

  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = 11,
    height = 4.8,
    dpi = 180
  )

  invisible(out_png)
}

.coerce_tabix_metadata_column <- function(x) {
  vals <- as.character(x)
  vals[vals %in% c("NA", "")] <- NA_character_
  non_missing <- vals[!is.na(vals)]
  if (!length(non_missing)) {
    return(vals)
  }

  lower_vals <- tolower(non_missing)
  if (all(lower_vals %in% c("true", "false"))) {
    out <- rep(NA, length(vals))
    keep <- !is.na(vals)
    out[keep] <- lower_vals == "true"
    return(out)
  }

  num_vals <- suppressWarnings(as.numeric(non_missing))
  if (all(is.finite(num_vals))) {
    out <- rep(NA_real_, length(vals))
    out[!is.na(vals)] <- num_vals
    return(out)
  }

  vals
}

.tabix_metadata_frame <- function(paths,
                                  sample_names = NULL,
                                  con = NULL) {
  stopifnot(is.character(paths))
  if (!length(paths)) {
    return(data.frame())
  }
  stopifnot(all(file.exists(paths)))

  if (is.null(sample_names)) {
    sample_names <- sub("\\.bed\\.gz$", "", basename(paths), ignore.case = TRUE)
  }
  stopifnot(length(sample_names) == length(paths))

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  rows <- lapply(seq_along(paths), function(i) {
    md <- as.list(tabix_metadata(paths[[i]], con = con))
    c(
      list(
        sample_name = as.character(sample_names[[i]]),
        artifact_path = normalizePath(paths[[i]], winslash = "/", mustWork = TRUE)
      ),
      md
    )
  })

  keys <- unique(unlist(lapply(rows, names), use.names = FALSE))
  mat <- lapply(rows, function(row) {
    out <- stats::setNames(rep(NA_character_, length(keys)), keys)
    for (nm in names(row)) {
      out[[nm]] <- as.character(row[[nm]])
    }
    out
  })
  df <- as.data.frame(do.call(rbind, mat), stringsAsFactors = FALSE, check.names = FALSE)
  convert_cols <- setdiff(names(df), c("sample_name", "artifact_path"))
  for (nm in convert_cols) {
    df[[nm]] <- .coerce_tabix_metadata_column(df[[nm]])
  }
  df
}
