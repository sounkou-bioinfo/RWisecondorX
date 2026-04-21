#' Build a typed QC report for a NIPTeR control group
#'
#' Summarises the stability of a \code{NIPTControlGroup} at the chromosome,
#' sample, and optional autosomal-bin level. The chromosome summary exposes the
#' mean, standard deviation, coefficient of variation, and Shapiro-Wilk
#' normality statistic for each autosome, across multiple scoring spaces in
#' one long table: post-chi chromosomal fractions, NCV, and the four RBZ
#' predictor sets. The sample summary adds SSD-based matching scores plus
#' strand-aware chromosomal alert counts from the post-chi diagnostic pass.
#' When sex labels are available, the sample summary is enriched with the same
#' per-sample sex-cluster metrics written into the typed reference frame
#' (\code{ConsensusGender}, sex-outlier flag, and \code{Z_X_XX}/\code{Z_X_XY}/
#' \code{Z_Y_XX}/\code{Z_Y_XY}), and the report also includes sex-stratified
#' \code{RR_X}/\code{RR_Y} spread metrics relevant for gaunosome model
#' readiness plus a dedicated long-form sex-model summary for XX/XY
#' sex-chromosome Z, NCV, and regression models. Optional bin-level output
#' exposes the scaled-count CV and chi-correction profile that underlies
#' \code{\link{nipter_chi_correct}}.
#'
#' @param control_group A \code{NIPTControlGroup}.
#' @param sample_sex Optional character vector overriding the sex labels stored
#'   on \code{control_group}. Accepted values are \code{"female"},
#'   \code{"male"}, \code{"ambiguous"}, and \code{"unknown"}.
#' @param reference_model Optional \code{NIPTReferenceModel} built from the
#'   same controls. When supplied, the QC bundle also extracts the explicit
#'   XX/XY sex-chromosome Z, NCV, and regression model summaries from it.
#' @param chi_cutoff Numeric scalar used to flag overdispersed bins in the
#'   optional chi profile. Default \code{3.5}.
#' @param z_cutoff Absolute post-chi z-score threshold used when counting
#'   aberrant chromosome rows in the sample summary. Default \code{3}.
#' @param collapse_strands Logical; when \code{FALSE} (default), separated
#'   strands retain their legacy F/R-specific diagnostic rows and the sample
#'   summary reports both distinct chromosomes and bidirectional aberrations.
#'   When \code{TRUE}, sample-level aberration counts collapse to 22 autosomes.
#' @param max_aberrant_chromosomes Maximum number of distinct aberrant
#'   chromosomes allowed before \code{is_chromosomal_outlier} is set in the
#'   sample summary. Default \code{2L}.
#' @param outlier_rule Character scalar describing how to convert
#'   \code{abberant_scores} into sample-level outlier flags.
#'   \code{"any_aberrant_score"} (default) removes any sample with at least
#'   one post-chi chromosome row beyond \code{z_cutoff}.
#'   \code{"bidirectional_or_multichromosome"} drops a sample only when both
#'   strands of one chromosome are aberrant or when too many distinct
#'   chromosomes are aberrant.
#' @param include_bins Logical; when \code{TRUE}, compute the autosomal
#'   bin-level scaled-count CV and chi profile. Default \code{FALSE}.
#' @param rbz_train_fraction Fraction of controls used to fit autosomal RBZ
#'   QC models. Values below \code{1} use a train/stat split; values greater
#'   than or equal to \code{1} fit and score on all retained controls.
#'   Default \code{1}.
#' @param rbz_seed Integer seed used when \code{rbz_train_fraction < 1}.
#'   Default \code{1995L}.
#' @param rbz_exclude_chromosomes Integer vector excluded from autosomal RBZ
#'   predictor search. Default \code{c(13L, 18L, 21L)}.
#'
#' @return A typed \code{NIPTControlGroupQC} object.
#'
#' @seealso [nipter_diagnose_control_group()], [nipter_match_matrix()],
#'   [nipter_chi_correct()]
#'
#' @export
nipter_control_group_qc <- function(control_group,
                                    sample_sex = NULL,
                                    reference_model = NULL,
                                    chi_cutoff = 3.5,
                                    z_cutoff = 3,
                                    collapse_strands = FALSE,
                                    max_aberrant_chromosomes = 2L,
                                    outlier_rule = c("any_aberrant_score", "bidirectional_or_multichromosome"),
                                    include_bins = FALSE,
                                    rbz_train_fraction = 1,
                                    rbz_seed = 1995L,
                                    rbz_exclude_chromosomes = c(13L, 18L, 21L)) {
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.numeric(chi_cutoff), length(chi_cutoff) == 1L)
  stopifnot(is.numeric(z_cutoff), length(z_cutoff) == 1L, z_cutoff > 0)
  stopifnot(is.logical(collapse_strands), length(collapse_strands) == 1L)
  stopifnot(is.numeric(max_aberrant_chromosomes),
            length(max_aberrant_chromosomes) == 1L,
            max_aberrant_chromosomes >= 0)
  stopifnot(is.numeric(rbz_train_fraction),
            length(rbz_train_fraction) == 1L,
            is.finite(rbz_train_fraction),
            rbz_train_fraction > 0)
  if (!is.null(rbz_seed)) {
    stopifnot(is.numeric(rbz_seed), length(rbz_seed) == 1L, is.finite(rbz_seed))
  }
  if (!is.null(rbz_exclude_chromosomes)) {
    stopifnot(is.numeric(rbz_exclude_chromosomes))
  }
  if (!is.null(reference_model)) {
    stopifnot(.is_nipt_reference_model(reference_model))
  }
  outlier_rule <- match.arg(outlier_rule)
  stopifnot(is.logical(include_bins), length(include_bins) == 1L)

  chi_corrected_control_group <- nipter_chi_correct(
    sample = control_group$samples[[1L]],
    control_group = control_group,
    chi_cutoff = chi_cutoff
  )[["control_group"]]
  diag_post_chi <- nipter_diagnose_control_group(
    chi_corrected_control_group,
    collapse_strands = collapse_strands,
    z_cutoff = z_cutoff
  )
  chromosome_summary <- .qc_chromosome_summary(
    chi_corrected_control_group = chi_corrected_control_group,
    z_cutoff = z_cutoff,
    rbz_train_fraction = rbz_train_fraction,
    rbz_seed = rbz_seed,
    rbz_exclude_chromosomes = rbz_exclude_chromosomes
  )

  ssd_mat <- nipter_match_matrix(control_group)
  sample_names <- unname(as.character(control_names(control_group)))
  mean_ssd <- vapply(seq_len(nrow(ssd_mat)), function(i) {
    mean(ssd_mat[i, setdiff(seq_len(ncol(ssd_mat)), i), drop = TRUE])
  }, numeric(1L))
  ssd_cutoff <- mean(mean_ssd) + 3 * stats::sd(mean_ssd)
  is_matching_outlier <- if (is.finite(ssd_cutoff)) mean_ssd > ssd_cutoff else rep(FALSE, length(mean_ssd))
  sample_summary <- data.frame(
    sample_name = sample_names,
    mean_ssd = mean_ssd,
    max_abs_z = vapply(seq_len(ncol(diag_post_chi$z_scores)), function(i) {
      vals <- abs(diag_post_chi$z_scores[, i])
      vals <- vals[is.finite(vals)]
      if (!length(vals)) NA_real_ else max(vals)
    }, numeric(1L)),
    is_matching_outlier = is_matching_outlier,
    stringsAsFactors = FALSE
  )
  ab_summary <- .qc_sample_aberration_summary(
    sample_names = sample_names,
    diagnostics = diag_post_chi,
    collapse_strands = collapse_strands,
    max_aberrant_chromosomes = max_aberrant_chromosomes,
    outlier_rule = outlier_rule
  )
  sample_summary <- merge(
    sample_summary,
    ab_summary,
    by = "sample_name",
    all.x = TRUE,
    sort = FALSE
  )
  sample_summary <- sample_summary[match(sample_names, sample_summary$sample_name), , drop = FALSE]

  sample_sex <- if (is.null(sample_sex)) {
    .control_sample_sex(control_group)
  } else {
    .normalize_sample_sex(sample_sex, control_group$samples)
  }
  sex_frame <- if (is.null(reference_model)) {
    .qc_annotated_sex_frame(control_group, sample_sex = sample_sex)
  } else {
    as.data.frame(reference_model$reference_frame, stringsAsFactors = FALSE)
  }
  sex_summary <- .qc_sex_summary(sex_frame)
  sex_model_summary <- .qc_sex_model_summary(
    reference_model = reference_model,
    ref = sex_frame
  )
  chromosome_summary <- .append_sex_model_chromosome_summary(
    chromosome_summary,
    sex_model_summary
  )
  if (!is.null(sex_frame)) {
    sex_cols <- sex_frame[, c(
      "Sample_name", "ConsensusGender",
      "RR_X_SexClassMAD", "RR_Y_SexClassMAD",
      "IsRefSexOutlier",
      "Z_X_XX", "Z_X_XY", "Z_Y_XX", "Z_Y_XY"
    ), drop = FALSE]
    names(sex_cols) <- c(
      "sample_name", "consensus_gender",
      "rr_x_sexclass_mad", "rr_y_sexclass_mad",
      "is_reference_sex_outlier",
      "z_x_xx", "z_x_xy", "z_y_xx", "z_y_xy"
    )
    sample_summary <- merge(
      sample_summary,
      sex_cols,
      by = "sample_name",
      all.x = TRUE,
      sort = FALSE
    )
    sample_summary <- sample_summary[match(sample_names, sample_summary$sample_name), , drop = FALSE]
  }

  chi_profile <- .chi_profile(
    control_group,
    chi_cutoff = chi_cutoff,
    include_bin_stats = include_bins
  )
  bin_summary <- if (isTRUE(include_bins)) .chi_profile_bin_summary(chi_profile) else NULL

  n_matching_outliers <- sum(sample_summary$is_matching_outlier %in% TRUE, na.rm = TRUE)
  n_chromosomal_outliers <- sum(sample_summary$is_chromosomal_outlier %in% TRUE, na.rm = TRUE)
  n_samples_with_aberrant_rows <- sum(sample_summary$n_aberrant_rows > 0L, na.rm = TRUE)
  n_reference_sex_outliers <- if ("is_reference_sex_outlier" %in% names(sample_summary)) {
    sum(sample_summary$is_reference_sex_outlier %in% TRUE, na.rm = TRUE)
  } else {
    NA_integer_
  }
  sex_readiness <- if (!is.null(sex_summary) && nrow(sex_summary)) {
    list(
      n_sex_groups_with_ncv_ready = sum(sex_summary$can_build_ncv %in% TRUE, na.rm = TRUE),
      n_sex_groups_with_regression_ready = sum(sex_summary$can_build_regression %in% TRUE, na.rm = TRUE)
    )
  } else {
    list(
      n_sex_groups_with_ncv_ready = NA_integer_,
      n_sex_groups_with_regression_ready = NA_integer_
    )
  }

  .as_nipt_control_group_qc(list(
    sample_names = sample_names,
    chromosome_summary = chromosome_summary,
    sample_summary = sample_summary,
    sex_summary = sex_summary,
    sex_model_summary = sex_model_summary,
    bin_summary = bin_summary,
    settings = list(
      chi_cutoff = chi_cutoff,
      z_cutoff = z_cutoff,
      collapse_strands = collapse_strands,
      max_aberrant_chromosomes = max_aberrant_chromosomes,
      outlier_rule = outlier_rule,
      include_bins = include_bins,
      n_controls = n_controls(control_group),
      n_bins = chi_profile$n_bins,
      binsize = chi_profile$binsize,
      n_valid_chi_bins = sum(chi_profile$valid_bins),
      n_invalid_chi_bins = sum(!chi_profile$valid_bins),
      invalid_chi_reason_counts = as.list(table(chi_profile$invalid_reason[!chi_profile$valid_bins])),
      n_overdispersed_bins = sum(chi_profile$overdispersed),
      chromosome_summary_source = "post_chi_with_model_qc",
      rbz_train_fraction = as.numeric(rbz_train_fraction),
      rbz_seed = if (is.null(rbz_seed)) NA_integer_ else as.integer(rbz_seed),
      rbz_exclude_chromosomes = as.integer(rbz_exclude_chromosomes),
      sample_sex_counts = if (is.null(sample_sex)) NULL else as.list(table(sample_sex)),
      n_matching_outliers = n_matching_outliers,
      n_chromosomal_outliers = n_chromosomal_outliers,
      n_samples_with_aberrant_rows = n_samples_with_aberrant_rows,
      n_reference_sex_outliers = n_reference_sex_outliers,
      n_sex_groups_with_ncv_ready = sex_readiness$n_sex_groups_with_ncv_ready,
      n_sex_groups_with_regression_ready = sex_readiness$n_sex_groups_with_regression_ready
    )
  ))
}

.qc_chromosome_summary <- function(chi_corrected_control_group,
                                   z_cutoff = 3,
                                   rbz_train_fraction = 1,
                                   rbz_seed = 1995L,
                                   rbz_exclude_chromosomes = c(13L, 18L, 21L)) {
  out <- rbind(
    .qc_post_chi_fraction_chromosome_summary(
      chi_corrected_control_group = chi_corrected_control_group,
      z_cutoff = z_cutoff
    ),
    .qc_ncv_chromosome_summary(chi_corrected_control_group),
    .qc_rbz_chromosome_summary(
      chi_corrected_control_group = chi_corrected_control_group,
      train_fraction = rbz_train_fraction,
      seed = rbz_seed,
      exclude_chromosomes = rbz_exclude_chromosomes
    )
  )

  rownames(out) <- NULL
  out
}

.new_qc_chromosome_row <- function(chromosome,
                                   mean_value = NA_real_,
                                   sd_value = NA_real_,
                                   cv_percent = NA_real_,
                                   shapiro_p_value = NA_real_,
                                   method_family,
                                   method_label,
                                   value_space,
                                   cv_type = NA_character_,
                                   model_terms = NA_character_,
                                   model_status = "ok",
                                   model_message = NA_character_,
                                   n_reference_samples = NA_integer_,
                                   n_training_samples = NA_integer_,
                                   n_stat_samples = NA_integer_,
                                   train_fraction_requested = NA_real_,
                                   train_fraction_effective = NA_real_,
                                   split_mode = NA_character_,
                                   split_seed = NA_integer_) {
  data.frame(
    chromosome = as.character(chromosome),
    mean_fraction = as.numeric(mean_value),
    sd_fraction = as.numeric(sd_value),
    cv_fraction = as.numeric(cv_percent),
    shapiro_p_value = as.numeric(shapiro_p_value),
    method_family = as.character(method_family),
    method_label = as.character(method_label),
    value_space = as.character(value_space),
    cv_type = as.character(cv_type),
    model_terms = as.character(model_terms),
    model_status = as.character(model_status),
    model_message = as.character(model_message),
    n_reference_samples = as.integer(n_reference_samples),
    n_training_samples = as.integer(n_training_samples),
    n_stat_samples = as.integer(n_stat_samples),
    train_fraction_requested = as.numeric(train_fraction_requested),
    train_fraction_effective = as.numeric(train_fraction_effective),
    split_mode = as.character(split_mode),
    split_seed = as.integer(split_seed),
    stringsAsFactors = FALSE
  )
}

.append_sex_model_chromosome_summary <- function(chromosome_summary,
                                                 sex_model_summary) {
  if (is.null(sex_model_summary) || !nrow(sex_model_summary)) {
    return(chromosome_summary)
  }

  sex_df <- as.data.frame(sex_model_summary, stringsAsFactors = FALSE)
  sex_df <- sex_df[
    sex_df$model_status == "ok" &
      sex_df$reference_sex %in% c("female", "male") &
      sex_df$focus_chromosome %in% c("X", "Y"),
    ,
    drop = FALSE
  ]
  if (!nrow(sex_df)) {
    return(chromosome_summary)
  }

  sex_group <- ifelse(sex_df$reference_sex == "female", "XX", "XY")
  method_family <- ifelse(
    sex_df$method_family == "regression",
    "rbz",
    sex_df$method_family
  )
  method_label <- ifelse(
    sex_df$method_family == "regression",
    paste0(sex_group, "_", sex_df$method_label),
    sex_group
  )

  rows <- data.frame(
    chromosome = paste0(sex_df$focus_chromosome, "_", sex_group),
    mean_fraction = sex_df$mean_value,
    sd_fraction = sex_df$sd_value,
    cv_fraction = sex_df$cv_percent,
    shapiro_p_value = sex_df$shapiro_p_value,
    method_family = method_family,
    method_label = method_label,
    value_space = sex_df$value_space,
    cv_type = NA_character_,
    model_terms = sex_df$model_terms,
    model_status = sex_df$model_status,
    model_message = sex_df$model_message,
    n_reference_samples = sex_df$n_reference_samples,
    n_training_samples = sex_df$n_training_samples,
    n_stat_samples = sex_df$n_stat_samples,
    train_fraction_requested = NA_real_,
    train_fraction_effective = NA_real_,
    split_mode = "all_samples",
    split_seed = NA_integer_,
    stringsAsFactors = FALSE
  )

  rownames(rows) <- NULL
  rbind(chromosome_summary, rows)
}

.qc_post_chi_fraction_chromosome_summary <- function(chi_corrected_control_group,
                                                     z_cutoff = 3) {
  diag_combined <- nipter_diagnose_control_group(
    chi_corrected_control_group,
    collapse_strands = TRUE,
    z_cutoff = z_cutoff
  )
  combined_stats <- as.data.frame(diag_combined$statistics, stringsAsFactors = FALSE)
  combined_stats$chromosome <- rownames(combined_stats)
  rownames(combined_stats) <- NULL
  combined_rows <- do.call(rbind, lapply(seq_len(nrow(combined_stats)), function(i) {
    mu <- combined_stats$mean[[i]]
    sdv <- combined_stats$SD[[i]]
    cv <- if (is.finite(mu) && abs(mu) > .Machine$double.eps) 100 * sdv / mu else NA_real_
    .new_qc_chromosome_row(
      chromosome = combined_stats$chromosome[[i]],
      mean_value = mu,
      sd_value = sdv,
      cv_percent = cv,
      shapiro_p_value = combined_stats$shapiro_p_value[[i]],
      method_family = "z_fraction",
      method_label = "combined",
      value_space = "fraction",
      n_reference_samples = n_controls(chi_corrected_control_group),
      n_training_samples = n_controls(chi_corrected_control_group),
      n_stat_samples = n_controls(chi_corrected_control_group)
    )
  }))

  if (.strand_type_of(chi_corrected_control_group) != "separated") {
    return(combined_rows)
  }

  diag_strand <- nipter_diagnose_control_group(
    chi_corrected_control_group,
    collapse_strands = FALSE,
    z_cutoff = z_cutoff
  )
  strand_stats <- as.data.frame(diag_strand$statistics, stringsAsFactors = FALSE)
  strand_stats$chromosome <- rownames(strand_stats)
  rownames(strand_stats) <- NULL
  strand_rows <- do.call(rbind, lapply(seq_len(nrow(strand_stats)), function(i) {
    raw_chr <- strand_stats$chromosome[[i]]
    strand <- sub("^.*([FR])$", "\\1", raw_chr)
    method_label <- if (identical(strand, "F")) "forward" else "reverse"
    mu <- strand_stats$mean[[i]]
    sdv <- strand_stats$SD[[i]]
    cv <- if (is.finite(mu) && abs(mu) > .Machine$double.eps) 100 * sdv / mu else NA_real_
    .new_qc_chromosome_row(
      chromosome = sub("[FR]$", "", raw_chr),
      mean_value = mu,
      sd_value = sdv,
      cv_percent = cv,
      shapiro_p_value = strand_stats$shapiro_p_value[[i]],
      method_family = "z_fraction",
      method_label = method_label,
      value_space = "fraction",
      n_reference_samples = n_controls(chi_corrected_control_group),
      n_training_samples = n_controls(chi_corrected_control_group),
      n_stat_samples = n_controls(chi_corrected_control_group)
    )
  }))

  rbind(combined_rows, strand_rows)
}

.qc_ncv_chromosome_summary <- function(chi_corrected_control_group) {
  control_sample <- chi_corrected_control_group$samples[[1L]]
  rows <- lapply(as.character(1:22), function(chrom) {
    res <- tryCatch(
      nipter_ncv_score(
        sample = control_sample,
        control_group = chi_corrected_control_group,
        chromo_focus = as.integer(chrom)
      ),
      error = identity
    )
    if (inherits(res, "error")) {
      return(.new_qc_chromosome_row(
        chromosome = chrom,
        method_family = "ncv",
        method_label = "ncv",
        value_space = "ratio",
        model_status = "error",
        model_message = conditionMessage(res),
        n_reference_samples = n_controls(chi_corrected_control_group),
        n_training_samples = n_controls(chi_corrected_control_group),
        n_stat_samples = NA_integer_
      ))
    }
    stats <- res$control_statistics
    .new_qc_chromosome_row(
      chromosome = chrom,
      mean_value = stats[["mean"]],
      sd_value = stats[["sd"]],
      cv_percent = 100 * stats[["sd"]] / stats[["mean"]],
      shapiro_p_value = stats[["shapiro_p_value"]],
      method_family = "ncv",
      method_label = "ncv",
      value_space = "ratio",
      model_terms = paste(res$denominators, collapse = " "),
      n_reference_samples = n_controls(chi_corrected_control_group),
      n_training_samples = n_controls(chi_corrected_control_group),
      n_stat_samples = length(res$control_z_scores)
    )
  })
  do.call(rbind, rows)
}

.qc_rbz_chromosome_summary <- function(chi_corrected_control_group,
                                       train_fraction = 1,
                                       seed = 1995L,
                                       exclude_chromosomes = c(13L, 18L, 21L)) {
  control_sample <- chi_corrected_control_group$samples[[1L]]
  rows <- vector("list", 22L * 4L)
  k <- 1L
  split_mode <- if (isTRUE(train_fraction >= 1)) "all_samples" else "train_test"

  for (chrom in 1:22) {
    reg <- tryCatch(
      nipter_regression(
        sample = control_sample,
        control_group = chi_corrected_control_group,
        chromo_focus = chrom,
        train_fraction = train_fraction,
        exclude_chromosomes = as.integer(exclude_chromosomes),
        seed = seed
      ),
      error = identity
    )
    if (inherits(reg, "error")) {
      for (i in seq_len(4L)) {
        rows[[k]] <- .new_qc_chromosome_row(
          chromosome = chrom,
          method_family = "rbz",
          method_label = paste0("predictor_set_", i),
          value_space = "ratio",
          model_status = "error",
          model_message = conditionMessage(reg),
          n_reference_samples = n_controls(chi_corrected_control_group),
          n_training_samples = NA_integer_,
          n_stat_samples = NA_integer_,
          train_fraction_requested = as.numeric(train_fraction),
          train_fraction_effective = NA_real_,
          split_mode = split_mode,
          split_seed = if (is.null(seed)) NA_integer_ else as.integer(seed)
        )
        k <- k + 1L
      }
      next
    }

    for (i in seq_len(4L)) {
      if (i > length(reg$models) || is.null(reg$models[[i]])) {
        rows[[k]] <- .new_qc_chromosome_row(
          chromosome = chrom,
          method_family = "rbz",
          method_label = paste0("predictor_set_", i),
          value_space = "ratio",
          model_status = "missing",
          model_message = "Regression builder returned fewer than 4 predictor sets.",
          n_reference_samples = n_controls(chi_corrected_control_group),
          n_training_samples = NA_integer_,
          n_stat_samples = NA_integer_,
          train_fraction_requested = as.numeric(train_fraction),
          train_fraction_effective = NA_real_,
          split_mode = split_mode,
          split_seed = if (is.null(seed)) NA_integer_ else as.integer(seed)
        )
        k <- k + 1L
        next
      }

      mdl <- reg$models[[i]]
      ctrl_ratio <- mdl$control_z_scores * mdl$cv + 1
      rows[[k]] <- .new_qc_chromosome_row(
        chromosome = chrom,
        mean_value = mean(ctrl_ratio),
        sd_value = stats::sd(ctrl_ratio),
        cv_percent = 100 * mdl$cv,
        shapiro_p_value = mdl$shapiro_p_value,
        method_family = "rbz",
        method_label = paste0("predictor_set_", i),
        value_space = "ratio",
        cv_type = mdl$cv_type,
        model_terms = paste(mdl$predictors, collapse = " "),
        n_reference_samples = mdl$n_reference_samples %||% n_controls(chi_corrected_control_group),
        n_training_samples = mdl$n_training_samples %||% NA_integer_,
        n_stat_samples = mdl$n_stat_samples %||% length(mdl$control_z_scores),
        train_fraction_requested = mdl$train_fraction_requested %||% as.numeric(train_fraction),
        train_fraction_effective = mdl$train_fraction_effective %||% NA_real_,
        split_mode = mdl$split_mode %||% split_mode,
        split_seed = mdl$split_seed %||% if (is.null(seed)) NA_integer_ else as.integer(seed)
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}

.qc_sample_aberration_summary <- function(sample_names,
                                         diagnostics,
                                         collapse_strands = FALSE,
                                         max_aberrant_chromosomes = 2L,
                                         outlier_rule = c("any_aberrant_score", "bidirectional_or_multichromosome")) {
  outlier_rule <- match.arg(outlier_rule)
  out <- data.frame(
    sample_name = sample_names,
    n_aberrant_rows = integer(length(sample_names)),
    n_aberrant_chromosomes = integer(length(sample_names)),
    has_bidirectional_aberration = rep(FALSE, length(sample_names)),
    is_chromosomal_outlier = rep(FALSE, length(sample_names)),
    stringsAsFactors = FALSE
  )

  ab <- diagnostics$aberrant_scores
  if (is.null(ab) || !nrow(ab)) {
    return(out)
  }

  flagged <- .flag_control_group_outliers(
    diagnostics,
    collapse_strands = collapse_strands,
    max_aberrant_chromosomes = max_aberrant_chromosomes,
    outlier_rule = outlier_rule
  )
  out$is_chromosomal_outlier[out$sample_name %in% flagged] <- TRUE

  ab$chromosome_core <- if (isTRUE(collapse_strands)) {
    as.character(ab$chromosome)
  } else {
    sub("[FR]$", "", ab$chromosome)
  }

  split_ab <- split(ab, ab$sample_name)
  for (nm in names(split_ab)) {
    idx <- match(nm, out$sample_name)
    if (is.na(idx)) {
      next
    }
    df <- split_ab[[nm]]
    out$n_aberrant_rows[[idx]] <- nrow(df)
    out$n_aberrant_chromosomes[[idx]] <- length(unique(df$chromosome_core))
    if (!isTRUE(collapse_strands)) {
      out$has_bidirectional_aberration[[idx]] <- any(table(df$chromosome_core) > 1L)
    }
  }

  out
}

.qc_annotated_sex_frame <- function(control_group, sample_sex = NULL) {
  if (is.null(sample_sex)) {
    return(NULL)
  }

  ref <- nipter_reference_frame(control_group, sample_sex = sample_sex)
  ref$RR_X <- ref$FrChrReads_X
  ref$RR_Y <- ref$FrChrReads_Y
  consensus_gender <- sample_sex

  if (!all(sample_sex %in% c("female", "male"))) {
    predicted_gender <- tryCatch(
      .predict_reference_sample_sex(
        .control_with_sample_sex(control_group, sample_sex = sample_sex),
        sex_models = list(
          y_fraction = nipter_sex_model(control_group, method = "y_fraction"),
          xy_fraction = nipter_sex_model(control_group, method = "xy_fraction")
        )
      ),
      error = function(e) NULL
    )
    resolved <- .resolve_reference_consensus_gender(
      sample_sex = sample_sex,
      predicted_gender = predicted_gender
    )
    if (!is.null(resolved)) {
      consensus_gender <- resolved
    }
  }

  ref <- .annotate_reference_frame_sex(
    ref,
    consensus_gender = consensus_gender,
    outlier_threshold = 3
  )

  ref
}

.qc_sex_summary <- function(ref) {
  if (is.null(ref)) {
    return(NULL)
  }

  out <- lapply(c("female", "male"), function(sex) {
    idx_all <- ref$ConsensusGender == sex
    idx_use <- idx_all & !ref$IsRefSexOutlier
    rr_x <- ref$RR_X[idx_use]
    rr_y <- ref$RR_Y[idx_use]
    n_use <- sum(idx_use)
    rr_x_sd <- if (n_use >= 2L) stats::sd(rr_x) else NA_real_
    rr_y_sd <- if (n_use >= 2L) stats::sd(rr_y) else NA_real_
    rr_x_mean <- if (n_use) mean(rr_x) else NA_real_
    rr_y_mean <- if (n_use) mean(rr_y) else NA_real_
    data.frame(
      sex_group = sex,
      n_samples = sum(idx_all),
      n_non_outlier_samples = n_use,
      x_model_mean = rr_x_mean,
      x_model_sd = rr_x_sd,
      x_model_cv = if (is.finite(rr_x_mean) && abs(rr_x_mean) > .Machine$double.eps) 100 * rr_x_sd / rr_x_mean else NA_real_,
      y_model_mean = rr_y_mean,
      y_model_sd = rr_y_sd,
      y_model_cv = if (is.finite(rr_y_mean) && abs(rr_y_mean) > .Machine$double.eps) 100 * rr_y_sd / rr_y_mean else NA_real_,
      rr_x_mean = rr_x_mean,
      rr_x_sd = rr_x_sd,
      rr_x_cv = if (is.finite(rr_x_mean) && abs(rr_x_mean) > .Machine$double.eps) 100 * rr_x_sd / rr_x_mean else NA_real_,
      rr_y_mean = rr_y_mean,
      rr_y_sd = rr_y_sd,
      rr_y_cv = if (is.finite(rr_y_mean) && abs(rr_y_mean) > .Machine$double.eps) 100 * rr_y_sd / rr_y_mean else NA_real_,
      rr_x_outliers = sum(idx_all & ref$IsRefSexOutlier),
      rr_y_outliers = sum(idx_all & ref$IsRefSexOutlier),
      can_build_ncv = n_use >= 2L,
      can_build_regression = n_use >= 4L,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

.new_qc_sex_model_row <- function(reference_sex,
                                  focus_chromosome,
                                  method_family,
                                  method_label,
                                  mean_value = NA_real_,
                                  sd_value = NA_real_,
                                  cv_percent = NA_real_,
                                  shapiro_p_value = NA_real_,
                                  value_space = "ratio",
                                  model_terms = NA_character_,
                                  adj_r_squared = NA_real_,
                                  model_status = "ok",
                                  model_message = NA_character_,
                                  n_reference_samples = NA_integer_,
                                  n_training_samples = NA_integer_,
                                  n_stat_samples = NA_integer_) {
  data.frame(
    reference_sex = as.character(reference_sex),
    focus_chromosome = as.character(focus_chromosome),
    method_family = as.character(method_family),
    method_label = as.character(method_label),
    mean_value = as.numeric(mean_value),
    sd_value = as.numeric(sd_value),
    cv_percent = as.numeric(cv_percent),
    shapiro_p_value = as.numeric(shapiro_p_value),
    value_space = as.character(value_space),
    model_terms = as.character(model_terms),
    adj_r_squared = as.numeric(adj_r_squared),
    model_status = as.character(model_status),
    model_message = as.character(model_message),
    n_reference_samples = as.integer(n_reference_samples),
    n_training_samples = as.integer(n_training_samples),
    n_stat_samples = as.integer(n_stat_samples),
    stringsAsFactors = FALSE
  )
}

.qc_sex_model_counts <- function(ref, sex) {
  idx_all <- ref$ConsensusGender == sex
  idx_use <- idx_all & !ref$IsRefSexOutlier
  c(
    n_reference_samples = sum(idx_all),
    n_training_samples = sum(idx_use),
    n_stat_samples = sum(idx_use)
  )
}

.qc_sex_fraction_model_summary <- function(ref) {
  rows <- vector("list", 4L)
  k <- 1L
  for (sex in c("female", "male")) {
    counts <- .qc_sex_model_counts(ref, sex)
    idx_use <- ref$ConsensusGender == sex & !ref$IsRefSexOutlier
    for (focus in c("X", "Y")) {
      values <- ref[[paste0("FrChrReads_", focus)]][idx_use]
      shapiro_p <- if (length(values) >= 3L && length(unique(values)) >= 3L) {
        stats::shapiro.test(values)$p.value
      } else {
        NA_real_
      }
      mean_value <- if (length(values)) mean(values) else NA_real_
      sd_value <- if (length(values) >= 2L) stats::sd(values) else NA_real_
      cv_percent <- if (is.finite(mean_value) && abs(mean_value) > .Machine$double.eps) {
        100 * sd_value / mean_value
      } else {
        NA_real_
      }
      rows[[k]] <- .new_qc_sex_model_row(
        reference_sex = sex,
        focus_chromosome = focus,
        method_family = "z_fraction",
        method_label = "sex_matched",
        mean_value = mean_value,
        sd_value = sd_value,
        cv_percent = cv_percent,
        shapiro_p_value = shapiro_p,
        value_space = "fraction",
        n_reference_samples = counts[["n_reference_samples"]],
        n_training_samples = counts[["n_training_samples"]],
        n_stat_samples = counts[["n_stat_samples"]]
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

.qc_sex_ncv_model_summary <- function(reference_model, ref) {
  rows <- vector("list", 4L)
  k <- 1L
  for (sex in c("female", "male")) {
    counts <- .qc_sex_model_counts(ref, sex)
    for (focus in c("X", "Y")) {
      if (is.null(reference_model) || !("sex_ncv_models" %in% names(reference_model))) {
        rows[[k]] <- .new_qc_sex_model_row(
          reference_sex = sex,
          focus_chromosome = focus,
          method_family = "ncv",
          method_label = "ncv",
          value_space = "ratio",
          model_status = "missing",
          model_message = "reference_model$sex_ncv_models is missing.",
          n_reference_samples = counts[["n_reference_samples"]],
          n_training_samples = counts[["n_training_samples"]],
          n_stat_samples = NA_integer_
        )
        k <- k + 1L
        next
      }
      mdl <- reference_model$sex_ncv_models[[sex]][[focus]]
      stats <- mdl$control_statistics
      rows[[k]] <- .new_qc_sex_model_row(
        reference_sex = sex,
        focus_chromosome = focus,
        method_family = "ncv",
        method_label = "ncv",
        mean_value = stats[["mean"]],
        sd_value = stats[["sd"]],
        cv_percent = stats[["cv"]],
        shapiro_p_value = stats[["shapiro_p_value"]],
        value_space = "ratio",
        model_terms = paste(mdl$denominators, collapse = " "),
        n_reference_samples = counts[["n_reference_samples"]],
        n_training_samples = length(mdl$reference_sample_names),
        n_stat_samples = length(mdl$reference_sample_names)
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

.qc_sex_regression_model_summary <- function(reference_model, ref) {
  rows <- vector("list", 16L)
  k <- 1L
  for (sex in c("female", "male")) {
    counts <- .qc_sex_model_counts(ref, sex)
    for (focus in c("X", "Y")) {
      if (is.null(reference_model) || !("sex_regression_models" %in% names(reference_model))) {
        for (i in seq_len(4L)) {
          rows[[k]] <- .new_qc_sex_model_row(
            reference_sex = sex,
            focus_chromosome = focus,
            method_family = "regression",
            method_label = paste0("predictor_set_", i),
            value_space = "ratio",
            model_status = "missing",
            model_message = "reference_model$sex_regression_models is missing.",
            n_reference_samples = counts[["n_reference_samples"]],
            n_training_samples = counts[["n_training_samples"]],
            n_stat_samples = NA_integer_
          )
          k <- k + 1L
        }
        next
      }

      mdl_set <- reference_model$sex_regression_models[[sex]][[focus]]
      for (i in seq_len(4L)) {
        if (i > length(mdl_set) || is.null(mdl_set[[i]])) {
          rows[[k]] <- .new_qc_sex_model_row(
            reference_sex = sex,
            focus_chromosome = focus,
            method_family = "regression",
            method_label = paste0("predictor_set_", i),
            value_space = "ratio",
            model_status = "missing",
            model_message = "Regression builder returned fewer than 4 predictor sets.",
            n_reference_samples = counts[["n_reference_samples"]],
            n_training_samples = counts[["n_training_samples"]],
            n_stat_samples = NA_integer_
          )
          k <- k + 1L
          next
        }

        mdl <- mdl_set[[i]]
        stats <- mdl$control_statistics
        rows[[k]] <- .new_qc_sex_model_row(
          reference_sex = sex,
          focus_chromosome = focus,
          method_family = "regression",
          method_label = paste0("predictor_set_", i),
          mean_value = stats[["mean_ratio"]],
          sd_value = stats[["sd_ratio"]],
          cv_percent = stats[["cv"]],
          shapiro_p_value = stats[["shapiro_p_value"]],
          value_space = "ratio",
          model_terms = paste(mdl$predictors, collapse = " "),
          adj_r_squared = stats[["adj_r_squared"]],
          n_reference_samples = counts[["n_reference_samples"]],
          n_training_samples = length(mdl$reference_sample_names),
          n_stat_samples = length(mdl$reference_sample_names)
        )
        k <- k + 1L
      }
    }
  }
  do.call(rbind, rows)
}

.qc_sex_model_summary <- function(reference_model = NULL, ref = NULL) {
  if (is.null(ref)) {
    return(NULL)
  }
  do.call(rbind, list(
    .qc_sex_fraction_model_summary(ref),
    .qc_sex_ncv_model_summary(reference_model, ref),
    .qc_sex_regression_model_summary(reference_model, ref)
  ))
}

.chi_profile_bin_summary <- function(profile) {
  if (is.null(profile$scaled_mean) || is.null(profile$scaled_sd) ||
      is.null(profile$scaled_cv)) {
    stop("Chi profile does not contain bin-level scaled statistics.", call. = FALSE)
  }

  data.frame(
    chromosome = rep(as.character(1:22), each = profile$n_bins),
    bin = rep(seq_len(profile$n_bins), times = 22L),
    binsize = profile$binsize,
    chromosome_length_bp = profile$chromosome_length_bp,
    bin_start_bp = profile$bin_start_bp,
    bin_end_bp = profile$bin_end_bp,
    in_chromosome_range = profile$in_chromosome_range,
    valid_for_chi = profile$valid_bins,
    invalid_reason = profile$invalid_reason,
    invalid_reason_detail = profile$invalid_reason_detail,
    expected_finite = profile$expected_finite,
    expected_positive = profile$expected_positive,
    scaled_finite = profile$scaled_finite,
    n_nonzero_scaled = profile$scaled_n_nonzero,
    mean_scaled = profile$scaled_mean,
    sd_scaled = profile$scaled_sd,
    cv_scaled = profile$scaled_cv,
    expected_count = profile$expected,
    chi_square = profile$chi_square,
    chi_z = profile$chi_z,
    overdispersed = profile$overdispersed,
    correction_factor = profile$correction_factor,
    stringsAsFactors = FALSE
  )
}
