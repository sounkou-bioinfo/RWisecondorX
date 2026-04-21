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

.expand_gc_autosomal_vector <- function(gc_table, n_bins) {
  stopifnot(is.list(gc_table), length(n_bins) == 1L, is.numeric(n_bins), n_bins >= 1L)
  unlist(lapply(as.character(seq_len(22L)), function(chr) {
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
                                  binsize = 50000L,
                                  con = NULL) {
  stopifnot(!is.null(gc_table) || !is.null(fasta))

  gc_tbl <- .resolve_gc_table(gc_table, fasta, binsize, con)
  raw_auto <- autosomal_matrix(sample)
  corrected_auto <- if (is.null(corrected)) raw_auto else autosomal_matrix(corrected)
  n_bins <- ncol(raw_auto)

  gc_auto <- .expand_gc_autosomal_vector(gc_tbl, n_bins = n_bins)
  raw_flat <- as.numeric(t(raw_auto))
  corrected_flat <- as.numeric(t(corrected_auto))
  valid <- !is.na(gc_auto) & gc_auto > 0 & is.finite(raw_flat) & raw_flat > 0

  out_metrics <- list(
    gc_loess_valid_bins = sum(valid),
    nipter_autosomal_bin_cv_pre = NA_real_,
    nipter_autosomal_bin_cv_post = NA_real_,
    nipter_corrected_bin_ratio_mean = NA_real_,
    nipter_corrected_bin_ratio_sd = NA_real_,
    nipter_corrected_bin_ratio_cv = NA_real_,
    nipter_gc_correlation_pre = NA_real_,
    nipter_gc_correlation_post = NA_real_
  )

  chromosomes <- rep(as.character(seq_len(22L)), each = n_bins)
  bin_index <- rep(seq_len(n_bins), times = 22L)
  curve_df <- data.frame(
    chromosome = chromosomes,
    bin = as.integer(bin_index),
    gc = as.numeric(gc_auto),
    raw_count = as.numeric(raw_flat),
    corrected_count = as.numeric(corrected_flat),
    gc_valid = as.logical(valid),
    raw_ratio = NA_real_,
    corrected_ratio = NA_real_
  )

  if (any(valid)) {
    raw_med <- stats::median(raw_flat[valid])
    corrected_med <- stats::median(corrected_flat[valid])

    if (is.finite(raw_med) && raw_med > 0) {
      curve_df$raw_ratio <- raw_flat / raw_med
    }
    if (is.finite(corrected_med) && corrected_med > 0) {
      curve_df$corrected_ratio <- corrected_flat / corrected_med
    }

    out_metrics$nipter_autosomal_bin_cv_pre <- .cv_percent(raw_flat[valid])
    out_metrics$nipter_autosomal_bin_cv_post <- .cv_percent(corrected_flat[valid])

    corrected_ratio_valid <- curve_df$corrected_ratio[valid]
    corrected_ratio_valid <- corrected_ratio_valid[is.finite(corrected_ratio_valid)]
    if (length(corrected_ratio_valid) >= 2L) {
      out_metrics$nipter_corrected_bin_ratio_mean <- mean(corrected_ratio_valid)
      out_metrics$nipter_corrected_bin_ratio_sd <- stats::sd(corrected_ratio_valid)
      out_metrics$nipter_corrected_bin_ratio_cv <- .cv_percent(corrected_ratio_valid)
    }

    raw_ratio_valid <- curve_df$raw_ratio[valid]
    raw_ratio_valid <- raw_ratio_valid[is.finite(raw_ratio_valid)]
    gc_valid <- gc_auto[valid]
    gc_valid <- gc_valid[is.finite(gc_valid)]
    if (length(raw_ratio_valid) == length(gc_valid) && length(gc_valid) >= 2L) {
      out_metrics$nipter_gc_correlation_pre <- suppressWarnings(stats::cor(gc_valid, raw_ratio_valid))
    }
    if (length(corrected_ratio_valid) == length(gc_valid) && length(gc_valid) >= 2L) {
      out_metrics$nipter_gc_correlation_post <- suppressWarnings(stats::cor(gc_valid, corrected_ratio_valid))
    }
  }

  list(
    metrics = out_metrics,
    curve_data = curve_df
  )
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
  valid_df <- curve_data[curve_data$gc_valid, , drop = FALSE]
  if (!nrow(valid_df)) {
    return(invisible(NULL))
  }

  valid_df$gc_bucket <- round(valid_df$gc, 3)
  curve_summary <- stats::aggregate(
    cbind(raw_ratio, corrected_ratio) ~ gc_bucket,
    data = valid_df[, c("gc_bucket", "raw_ratio", "corrected_ratio"), drop = FALSE],
    FUN = function(x) stats::median(x[is.finite(x)], na.rm = TRUE)
  )

  point_df <- rbind(
    data.frame(
      stage = "Raw",
      gc = valid_df$gc,
      ratio = valid_df$raw_ratio
    ),
    data.frame(
      stage = "Corrected",
      gc = valid_df$gc,
      ratio = valid_df$corrected_ratio
    )
  )
  point_df <- point_df[is.finite(point_df$gc) & is.finite(point_df$ratio), , drop = FALSE]

  line_df <- rbind(
    data.frame(stage = "Raw", gc = curve_summary$gc_bucket, ratio = curve_summary$raw_ratio),
    data.frame(stage = "Corrected", gc = curve_summary$gc_bucket, ratio = curve_summary$corrected_ratio)
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
        "Autosomal bins used for GC fit: ", sum(curve_data$gc_valid),
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
