.reference_with_model_field <- function(reference, field, value) {
  raw <- as.list(reference)
  raw[[field]] <- value
  .as_nipt_reference_model(raw)
}

.reference_frame_required <- function(reference, cols, caller) {
  frame <- as.data.frame(reference$reference_frame, stringsAsFactors = FALSE)
  missing <- setdiff(cols, names(frame))
  if (length(missing)) {
    stop(
      sprintf(
        "%s requires reference$reference_frame columns: %s. Rebuild the reference with nipter_build_reference().",
        caller,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  frame
}

.reference_sex_groups <- function(reference, caller) {
  frame <- .reference_frame_required(
    reference,
    cols = c("Sample_name", "ConsensusGender", "IsRefSexOutlier"),
    caller = caller
  )
  frame <- frame[
    frame$ConsensusGender %in% c("female", "male") &
      !frame$IsRefSexOutlier,
    ,
    drop = FALSE
  ]
  list(
    female = frame[frame$ConsensusGender == "female", , drop = FALSE],
    male = frame[frame$ConsensusGender == "male", , drop = FALSE]
  )
}

.candidate_ratio_columns <- function(candidate_chromosomes, focus_chromosome,
                                     extra_predictors, frame) {
  focus_chromosome <- match.arg(focus_chromosome, c("X", "Y"))
  other_sex <- setdiff(c("X", "Y"), focus_chromosome)
  cols <- c(paste0("RR_", as.character(candidate_chromosomes)),
            paste0("RR_", other_sex))
  cols <- unique(cols[cols %in% names(frame)])
  if (!is.null(extra_predictors)) {
    extra_present <- extra_predictors[
      extra_predictors %in% names(frame) &
        vapply(frame[extra_predictors[extra_predictors %in% names(frame)]],
               is.numeric, logical(1L))
    ]
    cols <- c(cols, extra_present)
  }
  cols
}

.reference_ratio_frame <- function(reference_frame) {
  frame <- as.data.frame(reference_frame, stringsAsFactors = FALSE)
  chroms <- c(as.character(1:22), "X", "Y")
  for (chrom in chroms) {
    rr_col <- paste0("RR_", chrom)
    fr_col <- paste0("FrChrReads_", chrom)
    if (!(rr_col %in% names(frame)) && fr_col %in% names(frame)) {
      frame[[rr_col]] <- frame[[fr_col]]
    }
  }
  frame
}

.sample_ratio_values <- function(sample, sample_predictors = NULL) {
  row <- .sample_reference_frame_row(sample)
  chroms <- c(as.character(1:22), "X", "Y")
  rr <- stats::setNames(
    vapply(chroms, function(chrom) {
      as.numeric(row[[paste0("FrChrReads_", chrom)]][[1L]])
    }, numeric(1L)),
    paste0("RR_", chroms)
  )
  if (!is.null(sample_predictors)) {
    stopifnot(is.list(sample_predictors))
    extra <- unlist(sample_predictors, use.names = TRUE)
    if (length(extra)) {
      rr <- c(rr, extra)
    }
  }
  rr
}

.build_sex_ncv_model_one <- function(reference_frame,
                                     focus_chromosome,
                                     reference_sex,
                                     candidate_chromosomes,
                                     min_elements,
                                     max_elements) {
  stopifnot(nrow(reference_frame) >= 2L)
  focus_col <- paste0("NChrReads_", focus_chromosome)
  cand_cols <- paste0("NChrReads_", as.character(candidate_chromosomes))
  cand_cols <- cand_cols[cand_cols %in% names(reference_frame)]
  if (length(cand_cols) < min_elements) {
    stop(
      sprintf(
        "Not enough denominator chromosomes remain to build %s %s NCV model.",
        reference_sex, focus_chromosome
      ),
      call. = FALSE
    )
  }

  min_elements <- min(as.integer(min_elements), length(cand_cols))
  max_elements <- min(as.integer(max_elements), length(cand_cols))
  if (min_elements < 1L || max_elements < min_elements) {
    stop("Invalid NCV denominator search bounds.", call. = FALSE)
  }

  counts <- as.matrix(reference_frame[, cand_cols, drop = FALSE])
  focus <- reference_frame[[focus_col]]
  best_cv <- Inf
  best_denoms <- NULL
  best_stats <- NULL

  for (size in seq.int(min_elements, max_elements)) {
    combos <- utils::combn(seq_along(cand_cols), size)
    for (j in seq_len(ncol(combos))) {
      idx <- combos[, j]
      denom_sum <- rowSums(counts[, idx, drop = FALSE])
      ratios <- focus / denom_sum
      mean_ratio <- mean(ratios)
      sd_ratio <- stats::sd(ratios)
      cv_ratio <- sd_ratio / mean_ratio
      if (is.finite(cv_ratio) && cv_ratio < best_cv) {
        shapiro_p <- if (length(ratios) >= 3L && length(unique(ratios)) >= 3L) {
          stats::shapiro.test(ratios)$p.value
        } else {
          NA_real_
        }
        best_cv <- cv_ratio
        best_denoms <- cand_cols[idx]
        best_stats <- c(
          mean = mean_ratio,
          sd = sd_ratio,
          cv = 100 * cv_ratio,
          shapiro_p_value = shapiro_p
        )
      }
    }
  }

  .as_nipt_sex_ncv_model(list(
    focus_chromosome = focus_chromosome,
    reference_sex = reference_sex,
    denominators = best_denoms,
    control_statistics = best_stats,
    reference_sample_names = reference_frame$Sample_name
  ))
}

#' Build sex-chromosome NCV models on a typed reference
#'
#' Searches denominator sets for X and Y separately on the female and male
#' non-outlier subsets of a \code{\link{NIPTReferenceModel}}. The resulting
#' \code{NIPTSexNCVModel} objects are attached to the reference as
#' \code{$sex_ncv_models}.
#'
#' @param reference A \code{NIPTReferenceModel} from
#'   \code{\link{nipter_build_reference}}.
#' @param candidate_chromosomes Integer vector of autosomes allowed in the
#'   denominator pool. Defaults to the production-style set
#'   \code{c(1:12, 14:16, 20, 22)}.
#' @param min_elements Minimum denominator-set size. Default \code{6L}.
#' @param max_elements Maximum denominator-set size. Default \code{9L}.
#' @param focus_chromosomes Character vector; any subset of \code{c("X", "Y")}.
#'
#' @return The input \code{NIPTReferenceModel}, with a typed
#'   \code{$sex_ncv_models} field containing nested \code{female}/\code{male}
#'   and \code{X}/\code{Y} model entries.
#'
#' @export
nipter_build_sex_ncv_models <- function(reference,
                                        candidate_chromosomes = c(1:12, 14:16, 20, 22),
                                        min_elements = 6L,
                                        max_elements = 9L,
                                        focus_chromosomes = c("X", "Y")) {
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  focus_chromosomes <- unique(match.arg(focus_chromosomes, c("X", "Y"), several.ok = TRUE))
  groups <- .reference_sex_groups(reference, "nipter_build_sex_ncv_models()")

  models <- list(female = list(), male = list())
  for (sex in names(models)) {
    for (focus in focus_chromosomes) {
      models[[sex]][[focus]] <- .build_sex_ncv_model_one(
        reference_frame = groups[[sex]],
        focus_chromosome = focus,
        reference_sex = sex,
        candidate_chromosomes = candidate_chromosomes,
        min_elements = min_elements,
        max_elements = max_elements
      )
    }
  }

  .reference_with_model_field(reference, "sex_ncv_models", models)
}

#' Score sex chromosomes with NCV models from a typed reference
#'
#' Computes sex-matched X or Y NCV scores from the \code{$sex_ncv_models}
#' attached to a \code{NIPTReferenceModel}.
#'
#' @param sample A \code{NIPTeRSample} or typed \code{NIPTSample}.
#' @param reference A \code{NIPTReferenceModel} with \code{$sex_ncv_models}
#'   already built by \code{\link{nipter_build_sex_ncv_models}}.
#' @param focus_chromosome Character scalar, \code{"X"} or \code{"Y"}.
#' @param y_unique_ratio Optional scalar passed through to
#'   \code{\link{nipter_predict_sex}} if the reference carries a
#'   \code{"y_unique"} sex model.
#'
#' @return A typed \code{NIPTSexNCVScore}.
#'
#' @export
nipter_ncv_sex_score <- function(sample,
                                 reference,
                                 focus_chromosome = c("X", "Y"),
                                 y_unique_ratio = NULL) {
  stopifnot(inherits(sample, "NIPTeRSample") || S7::S7_inherits(sample, NIPTSample))
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  if (!("sex_ncv_models" %in% names(reference))) {
    stop("reference$sex_ncv_models is missing; build them with nipter_build_sex_ncv_models().",
         call. = FALSE)
  }
  focus_chromosome <- match.arg(focus_chromosome)
  sex_prediction <- nipter_predict_sex(sample, reference, y_unique_ratio = y_unique_ratio)
  sample_row <- .sample_reference_frame_row(sample)
  scores <- c(female = NA_real_, male = NA_real_, selected = NA_real_)
  denoms <- list()
  stats_out <- list()
  ref_sizes <- c(female = NA_real_, male = NA_real_, same_sex = NA_real_)

  for (sex in c("female", "male")) {
    mdl <- reference$sex_ncv_models[[sex]][[focus_chromosome]]
    sample_ratio <- as.numeric(sample_row[[paste0("NChrReads_", focus_chromosome)]][[1L]]) /
      sum(vapply(mdl$denominators, function(col) {
        as.numeric(sample_row[[col]][[1L]])
      }, numeric(1L)))
    scores[[sex]] <- (sample_ratio - mdl$control_statistics[["mean"]]) /
      mdl$control_statistics[["sd"]]
    denoms[[sex]] <- mdl$denominators
    stats_out[[sex]] <- mdl$control_statistics
    ref_sizes[[sex]] <- length(mdl$reference_sample_names)
  }

  selected_sex <- sex_prediction$prediction
  scores[["selected"]] <- scores[[selected_sex]]
  ref_sizes[["same_sex"]] <- ref_sizes[[selected_sex]]

  .as_nipt_sex_ncv_score(list(
    sample_name = .sample_name(sample),
    focus_chromosome = focus_chromosome,
    predicted_sex = selected_sex,
    sex_prediction = sex_prediction,
    sample_scores = scores,
    denominators = denoms,
    control_statistics = stats_out,
    reference_sizes = ref_sizes
  ))
}

.build_sex_regression_models_one <- function(reference_frame,
                                             focus_chromosome,
                                             reference_sex,
                                             candidate_chromosomes,
                                             n_models,
                                             n_predictors,
                                             extra_predictors = "GCPCTAfterFiltering") {
  stopifnot(nrow(reference_frame) >= 4L)
  frame <- .reference_ratio_frame(reference_frame)
  response_col <- paste0("RR_", focus_chromosome)
  predictor_cols <- .candidate_ratio_columns(
    candidate_chromosomes = candidate_chromosomes,
    focus_chromosome = focus_chromosome,
    extra_predictors = extra_predictors,
    frame = frame
  )
  predictor_cols <- setdiff(predictor_cols, response_col)
  if (!length(predictor_cols)) {
    stop(
      sprintf(
        "No regression predictors available to build %s %s models.",
        reference_sex, focus_chromosome
      ),
      call. = FALSE
    )
  }

  used_predictors <- character(0L)
  out <- vector("list", n_models)
  for (i in seq_len(n_models)) {
    remaining <- setdiff(predictor_cols, used_predictors)
    if (!length(remaining)) {
      out <- out[seq_len(i - 1L)]
      break
    }
    n_this <- min(as.integer(n_predictors), length(remaining))
    mat <- t(as.matrix(frame[, c(remaining, response_col), drop = FALSE]))
    focus_row0 <- match(response_col, rownames(mat)) - 1L
    cand_rows0 <- match(remaining, rownames(mat)) - 1L
    sel_idx0 <- nipter_stepwise_cpp(mat, focus_row0, cand_rows0, n_this)
    selected <- remaining[sel_idx0 + 1L]
    if (!length(selected)) {
      out <- out[seq_len(i - 1L)]
      break
    }

    form <- stats::reformulate(selected, response = response_col)
    fit <- stats::lm(form, data = frame)
    pred <- stats::predict(fit, frame)
    ratios <- frame[[response_col]] / pred
    shapiro_p <- if (length(ratios) >= 3L && length(unique(ratios)) >= 3L) {
      stats::shapiro.test(ratios)$p.value
    } else {
      NA_real_
    }

    out[[i]] <- .as_nipt_sex_regression_model(list(
      fit = fit,
      focus_chromosome = focus_chromosome,
      reference_sex = reference_sex,
      response_column = response_col,
      predictors = selected,
      control_statistics = c(
        mean_ratio = mean(ratios),
        sd_ratio = stats::sd(ratios),
        cv = 100 * stats::sd(ratios) / mean(ratios),
        shapiro_p_value = shapiro_p,
        adj_r_squared = summary(fit)$adj.r.squared
      ),
      reference_sample_names = frame$Sample_name
    ))
    used_predictors <- c(used_predictors, selected)
  }

  Filter(Negate(is.null), out)
}

#' Build sex-chromosome regression models on a typed reference
#'
#' Fits multiple sex-matched X/Y ratio models on the non-outlier subsets of a
#' \code{\link{NIPTReferenceModel}}. The returned reference gains a typed
#' \code{$sex_regression_models} field.
#'
#' @param reference A \code{NIPTReferenceModel}.
#' @param candidate_chromosomes Integer vector of autosomes allowed as ratio
#'   predictors. Defaults to the production-style set
#'   \code{c(1:12, 14:16, 20, 22)}.
#' @param n_models Number of models to build per sex/focus combination.
#' @param n_predictors Maximum predictors per model.
#' @param extra_predictors Optional character vector of additional numeric
#'   columns to use when present in \code{reference$reference_frame}. The
#'   default \code{"GCPCTAfterFiltering"} keeps room for explicit QC metadata
#'   without hard-wiring that requirement into the core sample classes.
#' @param focus_chromosomes Character vector; any subset of \code{c("X", "Y")}.
#'
#' @return The input \code{NIPTReferenceModel}, with a typed
#'   \code{$sex_regression_models} field containing nested
#'   \code{female}/\code{male} and \code{X}/\code{Y} model lists.
#'
#' @export
nipter_build_sex_regression_models <- function(reference,
                                               candidate_chromosomes = c(1:12, 14:16, 20, 22),
                                               n_models = 4L,
                                               n_predictors = 4L,
                                               extra_predictors = "GCPCTAfterFiltering",
                                               focus_chromosomes = c("X", "Y")) {
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  focus_chromosomes <- unique(match.arg(focus_chromosomes, c("X", "Y"), several.ok = TRUE))
  groups <- .reference_sex_groups(reference, "nipter_build_sex_regression_models()")

  models <- list(female = list(), male = list())
  for (sex in names(models)) {
    for (focus in focus_chromosomes) {
      models[[sex]][[focus]] <- .build_sex_regression_models_one(
        reference_frame = groups[[sex]],
        focus_chromosome = focus,
        reference_sex = sex,
        candidate_chromosomes = candidate_chromosomes,
        n_models = as.integer(n_models),
        n_predictors = as.integer(n_predictors),
        extra_predictors = extra_predictors
      )
    }
  }

  .reference_with_model_field(reference, "sex_regression_models", models)
}

#' Score sex chromosomes with regression models from a typed reference
#'
#' Applies the sex-matched X/Y regression models attached to a
#' \code{NIPTReferenceModel}. Each per-model score is the standardized ratio
#' of observed to predicted sex-chromosome fraction, mirroring the production
#' pipeline's regression-based z-score idea.
#'
#' @param sample A \code{NIPTeRSample} or typed \code{NIPTSample}.
#' @param reference A \code{NIPTReferenceModel} with \code{$sex_regression_models}
#'   already built by \code{\link{nipter_build_sex_regression_models}}.
#' @param focus_chromosome Character scalar, \code{"X"} or \code{"Y"}.
#' @param y_unique_ratio Optional scalar passed through to
#'   \code{\link{nipter_predict_sex}} if the reference carries a
#'   \code{"y_unique"} sex model.
#' @param sample_predictors Optional named list of extra predictor values for
#'   the sample, used when the fitted regression models include extra columns
#'   such as \code{GCPCTAfterFiltering}.
#'
#' @return A typed \code{NIPTSexRegressionScore}.
#'
#' @export
nipter_regression_sex_score <- function(sample,
                                        reference,
                                        focus_chromosome = c("X", "Y"),
                                        y_unique_ratio = NULL,
                                        sample_predictors = NULL) {
  stopifnot(inherits(sample, "NIPTeRSample") || S7::S7_inherits(sample, NIPTSample))
  if (!.is_nipt_reference_model(reference)) {
    stop("'reference' must be a NIPTReferenceModel.", call. = FALSE)
  }
  if (!("sex_regression_models" %in% names(reference))) {
    stop(
      "reference$sex_regression_models is missing; build them with nipter_build_sex_regression_models().",
      call. = FALSE
    )
  }
  focus_chromosome <- match.arg(focus_chromosome)
  sex_prediction <- nipter_predict_sex(sample, reference, y_unique_ratio = y_unique_ratio)
  sample_values <- .sample_ratio_values(sample, sample_predictors = sample_predictors)
  response_col <- paste0("RR_", focus_chromosome)

  scores <- list(female = numeric(), male = numeric())
  ratios <- list(female = numeric(), male = numeric())
  predictors <- list(female = list(), male = list())
  reference_sizes <- c(female = NA_real_, male = NA_real_, same_sex = NA_real_)

  for (sex in c("female", "male")) {
    mdl_set <- reference$sex_regression_models[[sex]][[focus_chromosome]]
    scores[[sex]] <- vapply(mdl_set, function(mdl) {
      needed <- c(mdl$predictors)
      if (!all(needed %in% names(sample_values))) {
        stop(
          sprintf(
            "Sample is missing regression predictors required by the %s %s model: %s",
            sex, focus_chromosome, paste(setdiff(needed, names(sample_values)), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      sample_df <- as.data.frame(as.list(sample_values[unique(c(mdl$predictors, mdl$response_column))]),
                                 stringsAsFactors = FALSE)
      pred <- stats::predict(mdl$fit, newdata = sample_df)
      ratio <- sample_df[[mdl$response_column]][[1L]] / pred[[1L]]
      ratio
    }, numeric(1L))
    ratios[[sex]] <- scores[[sex]]
    scores[[sex]] <- vapply(seq_along(mdl_set), function(i) {
      mdl <- mdl_set[[i]]
      (ratios[[sex]][[i]] - 1) / mdl$control_statistics[["sd_ratio"]]
    }, numeric(1L))
    names(scores[[sex]]) <- paste0("model_", seq_along(mdl_set))
    names(ratios[[sex]]) <- names(scores[[sex]])
    predictors[[sex]] <- lapply(mdl_set, `[[`, "predictors")
    reference_sizes[[sex]] <- length(mdl_set[[1L]]$reference_sample_names)
  }

  selected_sex <- sex_prediction$prediction
  reference_sizes[["same_sex"]] <- reference_sizes[[selected_sex]]

  .as_nipt_sex_regression_score(list(
    sample_name = .sample_name(sample),
    focus_chromosome = focus_chromosome,
    predicted_sex = selected_sex,
    sex_prediction = sex_prediction,
    scores = scores,
    ratios = ratios,
    predictors = predictors,
    aggregate_scores = c(
      female_mean = mean(scores$female),
      male_mean = mean(scores$male),
      selected_mean = mean(scores[[selected_sex]])
    ),
    reference_sizes = reference_sizes
  ))
}
