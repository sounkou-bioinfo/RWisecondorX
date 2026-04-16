# Internal helpers for aggregate gaunosome scoring.

.gaunosome_summary_frame <- function(predicted_sex,
                                     sex_score,
                                     ncv_scores,
                                     regression_scores,
                                     focus_chromosomes) {
  focus_chromosomes <- unique(match.arg(
    focus_chromosomes,
    c("X", "Y"),
    several.ok = TRUE
  ))

  out <- lapply(focus_chromosomes, function(chrom) {
    ncv <- ncv_scores[[chrom]]
    reg <- regression_scores[[chrom]]
    data.frame(
      chromosome = chrom,
      predicted_sex = predicted_sex,
      z_score = unname(sex_score$z_scores[[paste0("Z_FrChrReads_", chrom)]]),
      cv = unname(sex_score$cv[[paste0("Z_FrChrReads_CV_", chrom)]]),
      ncv_score_female = unname(ncv$sample_scores[["female"]]),
      ncv_score_male = unname(ncv$sample_scores[["male"]]),
      ncv_score_selected = unname(ncv$sample_scores[["selected"]]),
      regression_score_female = unname(reg$aggregate_scores[["female_mean"]]),
      regression_score_male = unname(reg$aggregate_scores[["male_mean"]]),
      regression_score_selected = unname(reg$aggregate_scores[["selected_mean"]]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

#' Build both gaunosome model families on a typed reference
#'
#' Convenience wrapper that attaches both sex-chromosome NCV models and
#' sex-chromosome regression models to a \code{\link{NIPTReferenceModel}}.
#' This is the package-owned precomputation step for downstream X/Y scoring.
#'
#' @param reference A \code{NIPTReferenceModel}.
#' @param candidate_chromosomes Integer vector of autosomes allowed in the
#'   denominator/predictor pool. Defaults to the production-style set
#'   \code{c(1:12, 14:16, 20, 22)}.
#' @param ncv_min_elements Minimum denominator-set size for the NCV search.
#' @param ncv_max_elements Maximum denominator-set size for the NCV search.
#' @param regression_n_models Number of regression models to build per
#'   sex/focus combination.
#' @param regression_n_predictors Maximum predictors per regression model.
#' @param regression_extra_predictors Optional character vector of additional
#'   numeric columns to use when present in \code{reference$reference_frame}.
#' @param focus_chromosomes Character vector; any subset of \code{c("X", "Y")}.
#'
#' @return The input \code{NIPTReferenceModel}, enriched with both
#'   \code{$sex_ncv_models} and \code{$sex_regression_models}.
#'
#' @export
nipter_build_gaunosome_models <- function(reference,
                                          candidate_chromosomes = c(1:12, 14:16, 20, 22),
                                          ncv_min_elements = 6L,
                                          ncv_max_elements = 9L,
                                          regression_n_models = 4L,
                                          regression_n_predictors = 4L,
                                          regression_extra_predictors = "GCPCTAfterFiltering",
                                          focus_chromosomes = c("X", "Y")) {
  reference <- nipter_build_sex_ncv_models(
    reference = reference,
    candidate_chromosomes = candidate_chromosomes,
    min_elements = ncv_min_elements,
    max_elements = ncv_max_elements,
    focus_chromosomes = focus_chromosomes
  )

  nipter_build_sex_regression_models(
    reference = reference,
    candidate_chromosomes = candidate_chromosomes,
    n_models = regression_n_models,
    n_predictors = regression_n_predictors,
    extra_predictors = regression_extra_predictors,
    focus_chromosomes = focus_chromosomes
  )
}

#' Score gaunosomes from a typed reference model
#'
#' Computes the package-level X/Y scoring bundle against a prepared
#' \code{\link{NIPTReferenceModel}}: sex-matched z-scores, sex-matched NCV
#' scores, and sex-matched regression scores. The result is returned as one
#' validated object with a compact per-chromosome summary table.
#'
#' @param sample A \code{NIPTeRSample} or typed \code{NIPTSample}.
#' @param reference A \code{NIPTReferenceModel} prepared with
#'   \code{\link{nipter_build_gaunosome_models}}.
#' @param y_unique_ratio Optional numeric scalar passed through to
#'   \code{\link{nipter_predict_sex}} when the reference includes a
#'   \code{"y_unique"} sex model.
#' @param min_controls Minimum number of non-outlier same-sex controls required
#'   for the z-score component. Default \code{2L}.
#' @param sample_predictors Optional named list of extra predictor values for
#'   the sample, used when the regression models include extra columns such as
#'   \code{GCPCTAfterFiltering}.
#' @param focus_chromosomes Character vector; any subset of \code{c("X", "Y")}.
#'
#' @return A typed \code{NIPTGaunosomeScore} object containing the component
#'   sex, NCV, and regression scores plus a compact summary data frame.
#'
#' @export
nipter_gaunosome_score <- function(sample,
                                   reference,
                                   y_unique_ratio = NULL,
                                   min_controls = 2L,
                                   sample_predictors = NULL,
                                   focus_chromosomes = c("X", "Y")) {
  stopifnot(.is_nipt_sample_object(sample))
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  if (!("sex_ncv_models" %in% names(reference))) {
    stop(
      "reference$sex_ncv_models is missing; build them with nipter_build_gaunosome_models().",
      call. = FALSE
    )
  }
  if (!("sex_regression_models" %in% names(reference))) {
    stop(
      "reference$sex_regression_models is missing; build them with nipter_build_gaunosome_models().",
      call. = FALSE
    )
  }

  focus_chromosomes <- unique(match.arg(
    focus_chromosomes,
    c("X", "Y"),
    several.ok = TRUE
  ))

  sex_score <- nipter_sex_score(
    sample = sample,
    reference = reference,
    y_unique_ratio = y_unique_ratio,
    min_controls = min_controls
  )

  ncv_scores <- stats::setNames(vector("list", length(focus_chromosomes)), focus_chromosomes)
  regression_scores <- stats::setNames(vector("list", length(focus_chromosomes)), focus_chromosomes)
  for (chrom in focus_chromosomes) {
    ncv_scores[[chrom]] <- nipter_ncv_sex_score(
      sample = sample,
      reference = reference,
      focus_chromosome = chrom,
      y_unique_ratio = y_unique_ratio
    )
    regression_scores[[chrom]] <- nipter_regression_sex_score(
      sample = sample,
      reference = reference,
      focus_chromosome = chrom,
      y_unique_ratio = y_unique_ratio,
      sample_predictors = sample_predictors
    )
  }

  .as_nipt_gaunosome_score(list(
    sample_name = .sample_name(sample),
    predicted_sex = sex_score$predicted_sex,
    sex_prediction = sex_score$sex_prediction,
    sex_score = sex_score,
    ncv_scores = ncv_scores,
    regression_scores = regression_scores,
    summary = .gaunosome_summary_frame(
      predicted_sex = sex_score$predicted_sex,
      sex_score = sex_score,
      ncv_scores = ncv_scores,
      regression_scores = regression_scores,
      focus_chromosomes = focus_chromosomes
    )
  ))
}
