#' Plot chromosome-level NIPTeR control-group QC
#'
#' Draws the chromosome-level coefficient of variation profile from a
#' \code{\link{NIPTControlGroupQC}} object.
#'
#' @param qc A \code{NIPTControlGroupQC} from
#'   \code{\link{nipter_control_group_qc}}.
#' @param theme Theme name. One of \code{"minimal"} (default),
#'   \code{"light"}, \code{"bw"}, or \code{"classic"}.
#' @param base_size Base font size passed to the selected ggplot2 theme.
#' @param method_palette Optional named color vector keyed by method labels.
#' @param shapiro_shapes Optional named shape vector keyed by Shapiro labels
#'   \code{"<0.05"}, \code{">=0.05"}, and \code{"NA"}.
#' @param cv_step Base CV increment used for zero-based y-axis scaling.
#'   Labeled major breaks adapt when the facet range becomes large. Default
#'   \code{0.05}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
nipter_plot_qc_chromosomes <- function(qc,
                                       theme = c("minimal", "light", "bw", "classic"),
                                       base_size = 11,
                                       method_palette = NULL,
                                       shapiro_shapes = NULL,
                                       cv_step = 0.05) {
  .require_ggplot2("nipter_plot_qc_chromosomes()")
  stopifnot(.is_nipt_control_group_qc(qc))
  theme <- match.arg(theme)
  stopifnot(is.numeric(base_size), length(base_size) == 1L, base_size > 0)
  stopifnot(is.numeric(cv_step), length(cv_step) == 1L, cv_step > 0)
  shapiro_shapes <- .merge_named_defaults(
    shapiro_shapes,
    c("<0.05" = 17, ">=0.05" = 16, "NA" = 4)
  )

  df <- .chromosome_plot_data(qc)
  method_palette <- .complete_named_palette(
    .merge_named_defaults(method_palette, .chrom_method_palette()),
    levels = unique(df$series_label)
  )
  df$chromosome <- factor(
    as.character(df$chromosome),
    levels = c(as.character(1:22), "X_XX", "Y_XX", "X_XY", "Y_XY")
  )
  df$method_family <- factor(
    df$method_family,
    levels = c("z_fraction", "ncv", "rbz")
  )
  df$series_label <- factor(df$series_label, levels = unique(df$series_label))
  df$shapiro_flag <- ifelse(
    is.na(df$shapiro_p_value),
    "NA",
    ifelse(df$shapiro_p_value < 0.05, "<0.05", ">=0.05")
  )
  df <- df[!is.na(df$method_family), , drop = FALSE]
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = chromosome, y = cv_fraction, color = series_label, group = series_label)
  ) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(shape = shapiro_flag),
      size = 2.1,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = method_palette,
      drop = FALSE,
      name = "Series"
    ) +
    ggplot2::scale_shape_manual(
      values = shapiro_shapes,
      name = "Shapiro p"
    ) +
    ggplot2::scale_y_continuous(
      limits = function(x) .cv_axis_limits(x, step = cv_step),
      breaks = function(x) .cv_axis_major_breaks(x, step = cv_step),
      minor_breaks = function(x) .cv_axis_minor_breaks(x, step = cv_step),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = "Chromosome CV QC Across Controls",
      x = "Chromosome / sex model",
      y = "CV (%)"
    ) +
    .nipter_plot_theme(theme, base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 70, hjust = 1, vjust = 1),
      legend.position = "bottom"
    )

  p + ggplot2::facet_wrap(
    ~ method_family,
    ncol = 1,
    scales = "free_y",
    labeller = ggplot2::labeller(
      method_family = c(
        z_fraction = "Z fractions",
        ncv = "NCV",
        rbz = "RBZ"
      )
    )
  )
}

.chromosome_plot_data <- function(qc) {
  chr_df <- as.data.frame(qc$chromosome_summary, stringsAsFactors = FALSE)
  if (!nrow(chr_df)) {
    stop("qc$chromosome_summary is empty.", call. = FALSE)
  }
  chr_df <- chr_df[is.finite(chr_df$cv_fraction), , drop = FALSE]
  chr_df$chromosome <- as.character(chr_df$chromosome)
  chr_df$method_family <- if ("method_family" %in% names(chr_df)) {
    as.character(chr_df$method_family)
  } else {
    "z_fraction"
  }
  chr_df$series_label <- if ("method_label" %in% names(chr_df)) {
    as.character(chr_df$method_label)
  } else {
    "combined"
  }
  chr_df[, c("chromosome", "method_family", "series_label",
             "cv_fraction", "shapiro_p_value"), drop = FALSE]
}

.cv_axis_step <- 0.05

.cv_axis_upper <- function(x, step = .cv_axis_step) {
  finite_x <- x[is.finite(x)]
  upper <- if (length(finite_x)) max(finite_x) else 0
  max(step, ceiling(upper / step) * step)
}

.cv_axis_limits <- function(x, step = .cv_axis_step) {
  c(0, .cv_axis_upper(x, step = step))
}

.cv_axis_breaks <- function(x, step = .cv_axis_step) {
  seq(0, .cv_axis_upper(x, step = step), by = step)
}

.cv_axis_major_step <- function(x, step = .cv_axis_step, max_labels = 10L) {
  upper <- .cv_axis_upper(x, step = step)
  if (!is.finite(upper) || upper <= (step * max_labels)) {
    return(step)
  }
  raw_step <- upper / max_labels
  max(step, ceiling(raw_step / step) * step)
}

.cv_axis_major_breaks <- function(x, step = .cv_axis_step, max_labels = 10L) {
  seq(
    0,
    .cv_axis_upper(x, step = step),
    by = .cv_axis_major_step(x, step = step, max_labels = max_labels)
  )
}

.cv_axis_minor_breaks <- function(x, step = .cv_axis_step, max_minor = 200L) {
  upper <- .cv_axis_upper(x, step = step)
  n_minor <- ceiling(upper / step)
  if (!is.finite(n_minor) || n_minor > max_minor) {
    return(NULL)
  }
  seq(0, upper, by = step)
}

#' Plot sample-level NIPTeR control-group QC
#'
#' Visualises the control-group matching summary from a
#' \code{\link{NIPTControlGroupQC}} object.
#'
#' @param qc A \code{NIPTControlGroupQC} from
#'   \code{\link{nipter_control_group_qc}}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
nipter_plot_qc_samples <- function(qc) {
  .require_ggplot2("nipter_plot_qc_samples()")
  stopifnot(.is_nipt_control_group_qc(qc))

  df <- as.data.frame(qc$sample_summary, stringsAsFactors = FALSE)
  has_alert <- if ("n_aberrant_rows" %in% names(df)) {
    df$n_aberrant_rows > 0L
  } else {
    rep(FALSE, nrow(df))
  }
  has_prune_rule <- if ("is_chromosomal_outlier" %in% names(df)) {
    df$is_chromosomal_outlier
  } else {
    rep(FALSE, nrow(df))
  }
  df$status <- "clean"
  df$status[has_alert] <- "retained_alert"
  df$status[df$is_matching_outlier] <- "matching_outlier"
  df$status[has_alert & df$is_matching_outlier] <- "matching_outlier+alert"
  df$status[has_prune_rule] <- "prune_rule_outlier"
  df$status[has_prune_rule & df$is_matching_outlier] <- "prune_rule_outlier+matching"
  label_prune <- if ("is_chromosomal_outlier" %in% names(df)) {
    df$is_chromosomal_outlier
  } else {
    rep(FALSE, nrow(df))
  }
  label_idx <- df$is_matching_outlier |
    has_alert |
    label_prune |
    (!is.na(df$max_abs_z) & df$max_abs_z > 3)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = mean_ssd, y = max_abs_z, color = status)
  ) +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggplot2::geom_hline(yintercept = 3, linetype = "dashed", color = "grey50") +
    ggplot2::geom_text(
      data = df[label_idx, , drop = FALSE],
      ggplot2::aes(label = sample_name),
      nudge_y = 0.15,
      check_overlap = TRUE,
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        clean = "#1f6f8b",
        retained_alert = "#c73e1d",
        matching_outlier = "#8c510a",
        "matching_outlier+alert" = "#5c3a21",
        prune_rule_outlier = "#8b1e3f",
        "prune_rule_outlier+matching" = "#3b0f1b"
      ),
      name = "Sample status"
    ) +
    ggplot2::labs(
      title = "Control Matching QC",
      subtitle = if (identical(qc$settings$outlier_rule, "bidirectional_or_multichromosome")) {
        paste(
          "Only samples present in the supplied control group are plotted.",
          "This rule can retain points above |Z| when they do not meet the bidirectional-or-multi-chromosome drop condition."
        )
      } else {
        paste(
          "Only samples present in the supplied control group are plotted.",
          "Under any_aberrant_score, retained samples should not exceed the pruning z cutoff."
        )
      },
      x = "Mean SSD to Other Controls",
      y = "Max |Z| Across Chromosomes"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Plot autosomal bin-level NIPTeR control-group QC
#'
#' Draws either the scaled-count CV or chi-normalised overdispersion profile
#' across autosomal bins from a \code{\link{NIPTControlGroupQC}} object built
#' with \code{include_bins = TRUE}.
#'
#' @param qc A \code{NIPTControlGroupQC} from
#'   \code{\link{nipter_control_group_qc}} with \code{$bin_summary}.
#' @param metric Either \code{"cv_scaled"} (default) or \code{"chi_z"}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
nipter_plot_qc_bins <- function(qc, metric = c("cv_scaled", "chi_z")) {
  .require_ggplot2("nipter_plot_qc_bins()")
  stopifnot(.is_nipt_control_group_qc(qc))
  metric <- match.arg(metric)

  if (is.null(qc$bin_summary)) {
    stop(
      "qc$bin_summary is NULL. Recompute nipter_control_group_qc(include_bins = TRUE).",
      call. = FALSE
    )
  }

  df_all <- as.data.frame(qc$bin_summary, stringsAsFactors = FALSE)
  out_of_range_hidden <- sum(!df_all$in_chromosome_range, na.rm = TRUE)
  df <- df_all[df_all$in_chromosome_range, , drop = FALSE]
  if (!nrow(df)) {
    stop("No in-range autosomal bins remain to plot.", call. = FALSE)
  }
  chr_levels <- as.character(1:22)
  df$chromosome <- factor(df$chromosome, levels = chr_levels)
  df$chromosome_index <- match(as.character(df$chromosome), chr_levels)
  chr_bin_counts <- tapply(df$bin, as.character(df$chromosome), max)
  chr_bin_counts <- stats::setNames(
    as.integer(chr_bin_counts[chr_levels]),
    chr_levels
  )
  chr_offsets <- cumsum(c(0L, head(chr_bin_counts, -1L)))
  names(chr_offsets) <- chr_levels
  df$global_bin <- df$bin + chr_offsets[as.character(df$chromosome)]
  df$plot_value <- df[[metric]]

  valid_df <- df[df$valid_for_chi & is.finite(df$plot_value), , drop = FALSE]
  invalid_df <- df[!df$valid_for_chi, , drop = FALSE]
  over_df <- valid_df[valid_df$overdispersed, , drop = FALSE]

  y_label <- if (identical(metric, "cv_scaled")) {
    "Scaled Bin CV (%)"
  } else {
    "Chi Normalised Z"
  }
  title <- if (identical(metric, "cv_scaled")) {
    "Autosomal Bin CV Profile"
  } else {
    "Autosomal Bin Chi Profile"
  }
  subtitle <- sprintf(
    "plotted valid=%d invalid=%d (%s) overdispersed=%d; hidden out_of_range=%d",
    sum(df$valid_for_chi, na.rm = TRUE),
    sum(!df$valid_for_chi, na.rm = TRUE),
    .format_invalid_reason_counts(invalid_df$invalid_reason),
    sum(df$overdispersed, na.rm = TRUE),
    out_of_range_hidden
  )

  band_df <- data.frame(
    chromosome = factor(chr_levels, levels = chr_levels),
    chromosome_index = seq_along(chr_levels),
    xmin = chr_offsets + 0.5,
    xmax = chr_offsets + chr_bin_counts + 0.5,
    fill_band = rep(c("odd", "even"), length.out = length(chr_levels)),
    stringsAsFactors = FALSE
  )
  x_breaks <- chr_offsets + (chr_bin_counts / 2)

  valid_range <- range(valid_df$plot_value, na.rm = TRUE)
  if (!all(is.finite(valid_range))) {
    valid_range <- c(0, 1)
  }
  span <- diff(valid_range)
  if (!is.finite(span) || span <= .Machine$double.eps) {
    span <- max(abs(valid_range), na.rm = TRUE)
  }
  if (!is.finite(span) || span <= .Machine$double.eps) {
    span <- 1
  }
  invalid_y <- valid_range[[1L]] - 0.08 * span
  if (nrow(invalid_df)) {
    invalid_df$plot_value <- invalid_y
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = band_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_band),
      inherit.aes = FALSE,
      alpha = 0.55,
      color = NA
    ) +
    ggplot2::scale_fill_manual(values = c(odd = "#f5f7fa", even = "#ebf0f5"), guide = "none") +
    ggplot2::geom_line(
      data = valid_df,
      ggplot2::aes(x = global_bin, y = plot_value, group = chromosome),
      color = "#355c7d",
      linewidth = 0.22,
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = over_df,
      ggplot2::aes(x = global_bin, y = plot_value, color = "overdispersed"),
      size = 0.35,
      alpha = 0.9,
      show.legend = nrow(over_df) > 0L
    ) +
    ggplot2::geom_vline(
      xintercept = band_df$xmax + 0.5,
      color = "white",
      linewidth = 0.35,
      alpha = 0.9
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "grey75",
      linewidth = 0.25
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Chromosome",
      y = y_label
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = chr_levels,
      expand = ggplot2::expansion(mult = c(0.002, 0.002))
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.margin = ggplot2::margin(10, 12, 18, 10)
    )

  if (nrow(invalid_df) > 0L) {
    p <- p + ggplot2::geom_point(
      data = invalid_df,
      ggplot2::aes(x = global_bin, y = plot_value, color = invalid_reason),
      inherit.aes = FALSE,
      shape = 15,
      size = 0.5,
      alpha = 0.95
    ) +
      ggplot2::scale_color_manual(
        values = .bin_flag_palette(),
        name = "Bin flag",
        breaks = unique(c("overdispersed", invalid_df$invalid_reason))
      )
  } else if (nrow(over_df) > 0L) {
    p <- p + ggplot2::scale_color_manual(
      values = .bin_flag_palette(),
      name = "Bin flag",
      breaks = "overdispersed"
    )
  }
  if (identical(metric, "chi_z") && !is.null(qc$settings$chi_cutoff)) {
    p <- p + ggplot2::geom_hline(
      yintercept = as.numeric(qc$settings$chi_cutoff),
      linetype = "dashed",
      color = "#c73e1d",
      linewidth = 0.35
    )
  }

  p
}

#' Plot reference sex-model distributions
#'
#' Draws the per-sex X/Y fraction distributions used by the reference sex
#' models. When the reference frame carries \code{YUniqueRatio}, it is included
#' as an additional facet.
#'
#' @param reference A \code{NIPTReferenceModel} from
#'   \code{\link{nipter_build_reference}}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
nipter_plot_reference_sex_boxplots <- function(reference) {
  .require_ggplot2("nipter_plot_reference_sex_boxplots()")
  frame <- .reference_frame_for_sex_plot(
    reference,
    caller = "nipter_plot_reference_sex_boxplots()"
  )

  rows <- list(
    data.frame(
      Sample_name = frame$Sample_name,
      ConsensusGender = frame$ConsensusGender,
      IsRefSexOutlier = frame$IsRefSexOutlier,
      metric = "-log10(FrChrReads_X)",
      value = .negative_log10_safe(frame$FrChrReads_X),
      is_z_metric = FALSE,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Sample_name = frame$Sample_name,
      ConsensusGender = frame$ConsensusGender,
      IsRefSexOutlier = frame$IsRefSexOutlier,
      metric = "-log10(FrChrReads_Y)",
      value = .negative_log10_safe(frame$FrChrReads_Y),
      is_z_metric = FALSE,
      stringsAsFactors = FALSE
    )
  )
  if ("YUniqueRatio" %in% names(frame)) {
    rows[[length(rows) + 1L]] <- data.frame(
      Sample_name = frame$Sample_name,
      ConsensusGender = frame$ConsensusGender,
      IsRefSexOutlier = frame$IsRefSexOutlier,
      metric = "-log10(YUniqueRatio)",
      value = .negative_log10_safe(frame$YUniqueRatio),
      is_z_metric = FALSE,
      stringsAsFactors = FALSE
    )
  }
  z_metric_map <- c(
    Z_X_XX = "Z_X_XX",
    Z_X_XY = "Z_X_XY",
    Z_Y_XX = "Z_Y_XX",
    Z_Y_XY = "Z_Y_XY"
  )
  for (col in names(z_metric_map)) {
    if (col %in% names(frame)) {
      rows[[length(rows) + 1L]] <- data.frame(
        Sample_name = frame$Sample_name,
        ConsensusGender = frame$ConsensusGender,
        IsRefSexOutlier = frame$IsRefSexOutlier,
        metric = unname(z_metric_map[[col]]),
        value = frame[[col]],
        is_z_metric = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, rows)
  df <- df[is.finite(df$value), , drop = FALSE]
  df$metric <- factor(
    df$metric,
    levels = c(
      "-log10(FrChrReads_X)",
      "-log10(FrChrReads_Y)",
      "-log10(YUniqueRatio)",
      "Z_X_XX", "Z_X_XY", "Z_Y_XX", "Z_Y_XY"
    )
  )
  df$ConsensusGender <- factor(
    df$ConsensusGender,
    levels = c("female", "male", "ambiguous", "unknown")
  )

  z_guides <- unique(df[df$is_z_metric, "metric", drop = FALSE])

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = ConsensusGender, y = value)
  ) +
    ggplot2::geom_hline(
      data = z_guides,
      ggplot2::aes(yintercept = -3),
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.3
    ) +
    ggplot2::geom_hline(
      data = z_guides,
      ggplot2::aes(yintercept = 3),
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.3
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = ConsensusGender),
      outlier.shape = NA,
      alpha = 0.65,
      linewidth = 0.25
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = IsRefSexOutlier),
      width = 0.14,
      alpha = 0.8,
      size = 1.5
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    ggplot2::scale_fill_manual(
      values = .sex_group_palette(),
      name = "Consensus gender",
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#1f6f8b", `TRUE` = "#c73e1d"),
      name = "Sex outlier"
    ) +
    ggplot2::labs(
      title = "Reference Sex-Model Distributions",
      x = "Consensus gender",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
    )
}

#' Plot reference sex-model scatter spaces
#'
#' Draws either the X/Y fraction space or the RR_X/RR_Y ratio space used by the
#' reference sex and gaunosome models.
#'
#' @param reference A \code{NIPTReferenceModel} from
#'   \code{\link{nipter_build_reference}}.
#' @param space Either \code{"fraction"} (default) or \code{"ratio"}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
nipter_plot_reference_sex_scatter <- function(reference,
                                              space = c("fraction", "ratio")) {
  .require_ggplot2("nipter_plot_reference_sex_scatter()")
  space <- match.arg(space)
  frame <- .reference_frame_for_sex_plot(
    reference,
    caller = "nipter_plot_reference_sex_scatter()"
  )
  frame <- frame[frame$ConsensusGender %in% c("female", "male"), , drop = FALSE]
  if (!nrow(frame)) {
    stop("No female/male reference samples available for sex scatter plots.",
         call. = FALSE)
  }

  if (identical(space, "fraction")) {
    frame$plot_x <- .negative_log10_safe(frame$FrChrReads_Y)
    frame$plot_y <- .negative_log10_safe(frame$FrChrReads_X)
    x_label <- "-log10(FrChrReads_Y)"
    y_label <- "-log10(FrChrReads_X)"
    title <- "Reference Sex Fraction Space"
  } else {
    frame <- .reference_ratio_frame(frame)
    frame$plot_x <- frame$RR_Y
    frame$plot_y <- frame$RR_X
    x_label <- "RR_Y"
    y_label <- "RR_X"
    title <- "Reference Sex Ratio Space"
  }

  p <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = plot_x, y = plot_y, color = ConsensusGender, shape = IsRefSexOutlier)
  ) +
    ggplot2::geom_point(size = 2.1, alpha = 0.9) +
    ggplot2::scale_color_manual(
      values = c(female = "#1f6f8b", male = "#c73e1d"),
      name = "Consensus gender"
    ) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 1),
      name = "Sex outlier"
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  ellipse_df <- frame[!frame$IsRefSexOutlier, , drop = FALSE]
  ellipse_counts <- table(ellipse_df$ConsensusGender)
  ellipse_groups <- names(ellipse_counts[ellipse_counts >= 4L])
  ellipse_df <- ellipse_df[ellipse_df$ConsensusGender %in% ellipse_groups, , drop = FALSE]
  if (nrow(ellipse_df)) {
    p <- p + ggplot2::stat_ellipse(
      data = ellipse_df,
      ggplot2::aes(x = plot_x, y = plot_y, color = ConsensusGender),
      inherit.aes = FALSE,
      linewidth = 0.4,
      level = 0.99,
      type = "norm"
    )
  }

  p
}

#' Write NIPTeR reference QC plots
#'
#' Writes PNG plots for the existing NIPTeR QC bundle and reference sex-model
#' spaces.
#'
#' @param qc A \code{NIPTControlGroupQC} from
#'   \code{\link{nipter_control_group_qc}}.
#' @param reference A \code{NIPTReferenceModel} from
#'   \code{\link{nipter_build_reference}}.
#' @param outprefix Output file prefix.
#' @param dpi Output PNG DPI. Default \code{150L}.
#'
#' @return A named character vector of written file paths, invisibly.
#'
#' @export
write_nipter_reference_plots <- function(qc,
                                         reference,
                                         outprefix,
                                         dpi = 150L) {
  .require_ggplot2("write_nipter_reference_plots()")
  stopifnot(.is_nipt_control_group_qc(qc))
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  stopifnot(is.character(outprefix), length(outprefix) == 1L, nzchar(outprefix))
  stopifnot(is.numeric(dpi), length(dpi) == 1L, dpi > 0)

  plots <- list(
    chromosome_cv = nipter_plot_qc_chromosomes(qc),
    sample_qc = nipter_plot_qc_samples(qc),
    sex_boxplots = nipter_plot_reference_sex_boxplots(reference),
    sex_fraction_scatter = nipter_plot_reference_sex_scatter(reference, space = "fraction"),
    sex_ratio_scatter = nipter_plot_reference_sex_scatter(reference, space = "ratio")
  )
  if (!is.null(qc$bin_summary)) {
    plots$bin_cv <- nipter_plot_qc_bins(qc, metric = "cv_scaled")
    plots$bin_chi <- nipter_plot_qc_bins(qc, metric = "chi_z")
  }

  paths <- stats::setNames(character(length(plots)), names(plots))
  for (nm in names(plots)) {
    path <- paste0(outprefix, ".", nm, ".png")
    paths[[nm]] <- path
    ggplot2::ggsave(
      filename = path,
      plot = plots[[nm]],
      width = if (nm %in% c("bin_cv", "bin_chi")) 14 else 10,
      height = if (nm %in% c("bin_cv", "bin_chi")) {
        12
      } else if (nm %in% c("sex_boxplots", "chromosome_cv")) {
        9
      } else {
        7
      },
      dpi = as.integer(dpi)
    )
  }

  invisible(paths)
}

.require_ggplot2 <- function(caller) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      caller,
      " requires the 'ggplot2' package. Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
}

.reference_frame_for_sex_plot <- function(reference, caller) {
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  frame <- as.data.frame(reference$reference_frame, stringsAsFactors = FALSE)
  needed <- c("Sample_name", "ConsensusGender", "IsRefSexOutlier", "FrChrReads_X", "FrChrReads_Y")
  missing <- setdiff(needed, names(frame))
  if (length(missing)) {
    stop(
      caller,
      " requires reference_frame columns: ",
      paste(missing, collapse = ", "),
      ". Rebuild the reference with nipter_build_reference().",
      call. = FALSE
    )
  }
  frame
}

.bin_flag_palette <- function() {
  c(
    overdispersed = "#c73e1d",
    out_of_range = "#7a7a7a",
    nonpositive_expected = "#8b1e3f",
    nonfinite_expected = "#5c4b8a",
    nonfinite_scaled = "#946211",
    invalid_other = "#3d3d3d"
  )
}

.format_invalid_reason_counts <- function(x) {
  x <- x[!is.na(x) & x != "valid"]
  if (!length(x)) {
    return("none")
  }
  tb <- sort(table(x), decreasing = TRUE)
  paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}

.negative_log10_safe <- function(x, floor = 1e-12) {
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x_num))
  ok <- is.finite(x_num)
  out[ok] <- -log10(pmax(x_num[ok], floor))
  out
}

.sex_group_palette <- function() {
  c(
    female = "#1f6f8b",
    male = "#c73e1d",
    ambiguous = "#7f8c8d",
    unknown = "#a6a6a6"
  )
}

.chrom_method_palette <- function() {
  c(
    combined = "#1f6f8b",
    forward = "#4da3c7",
    reverse = "#144b5f",
    ncv = "#3b7a57",
    XX = "#7a3e9d",
    XY = "#b84a62",
    XX_predictor_set_1 = "#7a3e9d",
    XX_predictor_set_2 = "#9467bd",
    XX_predictor_set_3 = "#b08ad9",
    XX_predictor_set_4 = "#5e2a84",
    XY_predictor_set_1 = "#b84a62",
    XY_predictor_set_2 = "#d96f5a",
    XY_predictor_set_3 = "#f49d6e",
    XY_predictor_set_4 = "#8f2d3c",
    predictor_set_1 = "#c73e1d",
    predictor_set_2 = "#d98b2b",
    predictor_set_3 = "#7c4dff",
    predictor_set_4 = "#8b1e3f"
  )
}

.merge_named_defaults <- function(values, defaults) {
  if (is.null(values)) {
    return(defaults)
  }
  stopifnot(!is.null(names(values)))
  out <- defaults
  out[names(values)] <- values
  out
}

.complete_named_palette <- function(values, levels) {
  levels <- unique(as.character(levels))
  missing <- setdiff(levels, names(values))
  if (length(missing)) {
    values[missing] <- stats::setNames(
      grDevices::hcl.colors(length(missing), palette = "Dark 3"),
      missing
    )
  }
  values
}

.nipter_plot_theme <- function(theme = c("minimal", "light", "bw", "classic"),
                               base_size = 11) {
  theme <- match.arg(theme)
  switch(
    theme,
    minimal = ggplot2::theme_minimal(base_size = base_size),
    light = ggplot2::theme_light(base_size = base_size),
    bw = ggplot2::theme_bw(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size)
  )
}
