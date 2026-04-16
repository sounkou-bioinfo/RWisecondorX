#' Build a typed QC report for a NIPTeR control group
#'
#' Summarises the stability of a \code{NIPTControlGroup} at the chromosome,
#' sample, and optional autosomal-bin level. The chromosome summary exposes the
#' mean, standard deviation, coefficient of variation, and Shapiro-Wilk
#' normality statistic for each autosome. The sample summary adds SSD-based
#' matching scores so aberrant controls can be identified before scoring. When
#' sex labels are available, the report also includes sex-stratified
#' \code{RR_X}/\code{RR_Y} spread metrics relevant for gaunosome model
#' readiness. Optional bin-level output exposes the scaled-count CV and
#' chi-correction profile that underlies \code{\link{nipter_chi_correct}}.
#'
#' @param control_group A \code{NIPTControlGroup}.
#' @param sample_sex Optional character vector overriding the sex labels stored
#'   on \code{control_group}. Accepted values are \code{"female"},
#'   \code{"male"}, \code{"ambiguous"}, and \code{"unknown"}.
#' @param chi_cutoff Numeric scalar used to flag overdispersed bins in the
#'   optional chi profile. Default \code{3.5}.
#' @param include_bins Logical; when \code{TRUE}, compute the autosomal
#'   bin-level scaled-count CV and chi profile. Default \code{FALSE}.
#'
#' @return A typed \code{NIPTControlGroupQC} object.
#'
#' @seealso [nipter_diagnose_control_group()], [nipter_match_matrix()],
#'   [nipter_chi_correct()]
#'
#' @export
nipter_control_group_qc <- function(control_group,
                                    sample_sex = NULL,
                                    chi_cutoff = 3.5,
                                    include_bins = FALSE) {
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.numeric(chi_cutoff), length(chi_cutoff) == 1L)
  stopifnot(is.logical(include_bins), length(include_bins) == 1L)

  diag <- nipter_diagnose_control_group(control_group)
  chr_stats <- as.data.frame(diag$statistics, stringsAsFactors = FALSE)
  chr_stats$chromosome <- rownames(chr_stats)
  rownames(chr_stats) <- NULL
  chr_stats$cv_fraction <- NA_real_
  ok_chr <- is.finite(chr_stats$mean) & abs(chr_stats$mean) > .Machine$double.eps
  chr_stats$cv_fraction[ok_chr] <- 100 * chr_stats$SD[ok_chr] / chr_stats$mean[ok_chr]
  chromosome_summary <- chr_stats[, c(
    "chromosome", "mean", "SD", "cv_fraction", "shapiro_p_value"
  )]
  names(chromosome_summary) <- c(
    "chromosome", "mean_fraction", "sd_fraction", "cv_fraction", "shapiro_p_value"
  )

  ssd_mat <- nipter_match_matrix(control_group)
  sample_names <- unname(as.character(control_names(control_group)))
  mean_ssd <- vapply(seq_len(nrow(ssd_mat)), function(i) {
    mean(ssd_mat[i, setdiff(seq_len(ncol(ssd_mat)), i), drop = TRUE])
  }, numeric(1L))
  max_abs_z <- vapply(seq_len(ncol(diag$z_scores)), function(i) {
    vals <- abs(diag$z_scores[, i])
    vals <- vals[is.finite(vals)]
    if (!length(vals)) NA_real_ else max(vals)
  }, numeric(1L))
  ssd_cutoff <- mean(mean_ssd) + 3 * stats::sd(mean_ssd)
  is_matching_outlier <- if (is.finite(ssd_cutoff)) mean_ssd > ssd_cutoff else rep(FALSE, length(mean_ssd))
  sample_summary <- data.frame(
    sample_name = sample_names,
    mean_ssd = mean_ssd,
    max_abs_z = max_abs_z,
    is_matching_outlier = is_matching_outlier,
    stringsAsFactors = FALSE
  )

  sample_sex <- if (is.null(sample_sex)) {
    .control_sample_sex(control_group)
  } else {
    .normalize_sample_sex(sample_sex, control_group$samples)
  }
  sex_summary <- .qc_sex_summary(control_group, sample_sex = sample_sex)

  chi_profile <- .chi_profile(
    control_group,
    chi_cutoff = chi_cutoff,
    include_bin_stats = include_bins
  )
  bin_summary <- if (isTRUE(include_bins)) .chi_profile_bin_summary(chi_profile) else NULL

  .as_nipt_control_group_qc(list(
    sample_names = sample_names,
    chromosome_summary = chromosome_summary,
    sample_summary = sample_summary,
    sex_summary = sex_summary,
    bin_summary = bin_summary,
    settings = list(
      chi_cutoff = chi_cutoff,
      include_bins = include_bins,
      n_controls = n_controls(control_group),
      n_bins = chi_profile$n_bins,
      n_overdispersed_bins = sum(chi_profile$overdispersed),
      sample_sex_counts = if (is.null(sample_sex)) NULL else as.list(table(sample_sex))
    )
  ))
}

.qc_sex_summary <- function(control_group, sample_sex = NULL) {
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

.chi_profile_bin_summary <- function(profile) {
  if (is.null(profile$scaled_mean) || is.null(profile$scaled_sd) ||
      is.null(profile$scaled_cv)) {
    stop("Chi profile does not contain bin-level scaled statistics.", call. = FALSE)
  }

  data.frame(
    chromosome = rep(as.character(1:22), each = profile$n_bins),
    bin = rep(seq_len(profile$n_bins), times = 22L),
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
