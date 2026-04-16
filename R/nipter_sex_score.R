# Internal helpers for sex-chromosome scoring against a typed reference model.

.sample_reference_metrics <- function(sample) {
  row <- .sample_reference_frame_row(sample)
  c(
    FrChrReads_X = as.numeric(row$FrChrReads_X[[1L]]),
    FrChrReads_Y = as.numeric(row$FrChrReads_Y[[1L]]),
    RR_X = as.numeric(row$FrChrReads_X[[1L]]),
    RR_Y = as.numeric(row$FrChrReads_Y[[1L]])
  )
}

.score_metric_against_reference <- function(value, reference_values, min_controls) {
  ref <- reference_values[is.finite(reference_values)]
  n_ref <- length(ref)
  if (n_ref < min_controls) {
    return(list(z = NA_real_, cv = NA_real_, n = n_ref))
  }

  ref_mean <- mean(ref)
  ref_sd <- stats::sd(ref)
  ref_cv <- if (is.finite(ref_mean) && abs(ref_mean) > .Machine$double.eps) {
    100 * ref_sd / ref_mean
  } else {
    NA_real_
  }

  ref_z <- if (is.finite(ref_sd) && ref_sd > .Machine$double.eps) {
    (value - ref_mean) / ref_sd
  } else {
    NA_real_
  }

  list(z = ref_z, cv = ref_cv, n = n_ref)
}

.score_sex_reference_group <- function(sample_metrics, reference_frame, min_controls) {
  x_stats <- .score_metric_against_reference(
    sample_metrics[["FrChrReads_X"]],
    reference_frame$FrChrReads_X,
    min_controls = min_controls
  )
  y_stats <- .score_metric_against_reference(
    sample_metrics[["FrChrReads_Y"]],
    reference_frame$FrChrReads_Y,
    min_controls = min_controls
  )

  list(
    z_x = x_stats$z,
    z_y = y_stats$z,
    cv_x = x_stats$cv,
    cv_y = y_stats$cv,
    n = min(x_stats$n, y_stats$n)
  )
}

#' Score sex chromosomes against a typed NIPT reference model
#'
#' Computes sex-matched X/Y z-scores from the enriched
#' \code{\link{NIPTReferenceModel}} built by \code{\link{nipter_build_reference}}.
#' The scoring mirrors the production pipeline's gaunosome z-score stage:
#' predicted fetal sex is used to select the relevant non-outlier reference
#' subset, and the sample's X/Y fractions are standardized against that subset.
#'
#' The returned object also includes the alternate XX and XY reference-group
#' scores so downstream code can inspect both hypotheses without rebuilding the
#' reference statistics.
#'
#' @param sample A \code{NIPTeRSample} or typed \code{NIPTSample}.
#' @param reference A \code{NIPTReferenceModel} built by
#'   \code{\link{nipter_build_reference}}.
#' @param y_unique_ratio Optional numeric scalar passed through to
#'   \code{\link{nipter_predict_sex}} when the reference model includes a
#'   \code{"y_unique"} sex model.
#' @param min_controls Minimum number of non-outlier reference samples required
#'   in a sex class to compute a score. Default \code{2L}.
#'
#' @return A list-like \code{NIPTSexScore} S7 object with:
#' \describe{
#'   \item{sample_name}{The sample identifier.}
#'   \item{predicted_sex}{Consensus male/female call from the reference model's
#'     sex classifiers.}
#'   \item{sex_prediction}{The underlying \code{NIPTeRSexPrediction} object.}
#'   \item{sample_metrics}{Named numeric vector containing
#'     \code{FrChrReads_X}, \code{FrChrReads_Y}, \code{RR_X}, and \code{RR_Y}.}
#'   \item{z_scores}{Named numeric vector containing the selected,
#'     XX-reference, and XY-reference X/Y z-scores.}
#'   \item{cv}{Named numeric vector containing the corresponding
#'     coefficients of variation.}
#'   \item{reference_sizes}{Named numeric vector with female, male, and
#'     selected same-sex reference counts after outlier exclusion.}
#'   \item{reference_sample_names}{Character vector of the same-sex control
#'     samples actually used for the selected score.}
#' }
#'
#' @seealso [nipter_build_reference()], [nipter_predict_sex()]
#'
#' @export
nipter_sex_score <- function(sample,
                             reference,
                             y_unique_ratio = NULL,
                             min_controls = 2L) {
  stopifnot(.is_nipt_sample_object(sample))
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  stopifnot(is.numeric(min_controls), length(min_controls) == 1L)
  min_controls <- as.integer(min_controls)
  if (is.na(min_controls) || min_controls < 2L) {
    stop("'min_controls' must be an integer >= 2.", call. = FALSE)
  }

  reference_frame <- as.data.frame(reference$reference_frame,
                                   stringsAsFactors = FALSE)
  required_cols <- c("Sample_name", "ConsensusGender", "FrChrReads_X",
                     "FrChrReads_Y", "IsRefSexOutlier")
  missing_cols <- setdiff(required_cols, names(reference_frame))
  if (length(missing_cols)) {
    stop(
      "reference$reference_frame is missing sex-scoring columns: ",
      paste(missing_cols, collapse = ", "),
      ". Rebuild it with nipter_build_reference().",
      call. = FALSE
    )
  }
  if (!length(reference$sex_models)) {
    stop("reference$sex_models is empty; cannot predict sample sex.",
         call. = FALSE)
  }

  sex_prediction <- nipter_predict_sex(
    sample,
    reference,
    y_unique_ratio = y_unique_ratio
  )
  predicted_sex <- sex_prediction$prediction

  reference_frame <- reference_frame[
    !reference_frame$IsRefSexOutlier &
      reference_frame$ConsensusGender %in% c("female", "male"),
    ,
    drop = FALSE
  ]

  female_ref <- reference_frame[
    reference_frame$ConsensusGender == "female",
    ,
    drop = FALSE
  ]
  male_ref <- reference_frame[
    reference_frame$ConsensusGender == "male",
    ,
    drop = FALSE
  ]
  same_sex_ref <- if (predicted_sex == "male") male_ref else female_ref
  if (nrow(same_sex_ref) < min_controls) {
    stop(
      sprintf(
        "Need at least %d non-outlier %s controls to score the sample; found %d.",
        min_controls,
        predicted_sex,
        nrow(same_sex_ref)
      ),
      call. = FALSE
    )
  }

  sample_metrics <- .sample_reference_metrics(sample)
  score_xx <- .score_sex_reference_group(sample_metrics, female_ref, min_controls)
  score_xy <- .score_sex_reference_group(sample_metrics, male_ref, min_controls)
  selected <- if (predicted_sex == "male") score_xy else score_xx

  .as_nipt_sex_score(list(
    sample_name = .sample_name(sample),
    predicted_sex = predicted_sex,
    sex_prediction = sex_prediction,
    sample_metrics = sample_metrics,
    z_scores = c(
      Z_FrChrReads_X = selected$z_x,
      Z_FrChrReads_Y = selected$z_y,
      Z_FrChrReads_X_XX = score_xx$z_x,
      Z_FrChrReads_Y_XX = score_xx$z_y,
      Z_FrChrReads_X_XY = score_xy$z_x,
      Z_FrChrReads_Y_XY = score_xy$z_y
    ),
    cv = c(
      Z_FrChrReads_CV_X = selected$cv_x,
      Z_FrChrReads_CV_Y = selected$cv_y,
      Z_FrChrReads_CV_X_XX = score_xx$cv_x,
      Z_FrChrReads_CV_Y_XX = score_xx$cv_y,
      Z_FrChrReads_CV_X_XY = score_xy$cv_x,
      Z_FrChrReads_CV_Y_XY = score_xy$cv_y
    ),
    reference_sizes = c(
      female = nrow(female_ref),
      male = nrow(male_ref),
      same_sex = nrow(same_sex_ref)
    ),
    reference_sample_names = same_sex_ref$Sample_name
  ))
}
