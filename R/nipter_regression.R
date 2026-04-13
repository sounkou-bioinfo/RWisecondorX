#' Regression-based Z-score for trisomy prediction
#'
#' Implements NIPTeR's stepwise linear regression approach. Builds multiple
#' regression models that predict the focus chromosome's fraction from
#' other chromosomal fractions, then computes Z-scores based on the
#' prediction residuals.
#'
#' For \code{SeparatedStrands} samples the predictor pool is doubled (each
#' autosomal chromosome contributes both a forward and reverse fraction),
#' and complementary-strand exclusion within each model prevents selecting
#' both strands of the same chromosome as predictors in the same model.
#' Across models, only the exact predictor string is excluded; the
#' complementary strand remains available.
#'
#' @param sample A \code{NIPTeRSample} object.
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param chromo_focus Integer; the target chromosome (1-22).
#' @param n_models Number of regression models to build (default 4).
#' @param n_predictors Maximum number of predictor chromosomes per model
#'   (default 4).
#' @param exclude_chromosomes Integer vector of chromosomes to exclude from
#'   the predictor pool (default \code{c(13, 18, 21)}).
#' @param include_chromosomes Integer vector of chromosomes to force-include.
#' @param train_fraction Fraction of control samples used for training
#'   (default 0.6).
#' @param overdispersion_rate Overdispersion multiplier for theoretical CV
#'   (default 1.15).
#' @param force_practical_cv Logical; always use practical CV regardless of
#'   theoretical CV? Default \code{FALSE}.
#' @param seed Random seed for the train/test split. Default \code{NULL}
#'   (no seed set).
#'
#' @return A list of class \code{"NIPTeRRegression"} with elements:
#' \describe{
#'   \item{models}{A list of \code{n_models} sublists, each containing:
#'     \code{z_score} (numeric), \code{cv} (selected CV), \code{cv_type}
#'     (\code{"practical"} or \code{"theoretical"}), \code{predictors}
#'     (character vector of predictor chromosomes), \code{shapiro_p_value},
#'     \code{control_z_scores} (named numeric vector).}
#'   \item{focus_chromosome}{Character.}
#'   \item{correction_status}{Character vector.}
#'   \item{sample_name}{Character.}
#' }
#'
#' @details
#' The algorithm:
#'
#' 1. Split the control group 60/40 into train and test sets.
#' 2. Compute chromosomal fractions for all samples.
#' 3. For each of \code{n_models} models, greedily select \code{n_predictors}
#'    chromosomes by maximising adjusted R-squared (forward stepwise
#'    selection). Predictors already used in previous models are excluded.
#' 4. For each model, fit \code{lm(focus ~ predictors)} on the training set,
#'    predict on the test set, and compute:
#'    - Practical CV: \code{sd(obs/pred) / mean(obs/pred)} on the test set.
#'    - Theoretical CV: \code{overdispersion_rate / sqrt(total_reads_focus)}
#'      for the test sample.
#'    - The larger CV is used (unless \code{force_practical_cv = TRUE}).
#'    - Sample Z-score: \code{(sample_ratio - 1) / cv}.
#'
#' @seealso [nipter_z_score()], [nipter_ncv_score()]
#'
#' @examples
#' \dontrun{
#' reg21 <- nipter_regression(sample, cg, chromo_focus = 21)
#' reg21$models[[1]]$z_score
#' }
#'
#' @export
nipter_regression <- function(sample,
                              control_group,
                              chromo_focus,
                              n_models             = 4L,
                              n_predictors         = 4L,
                              exclude_chromosomes  = c(13L, 18L, 21L),
                              include_chromosomes  = NULL,
                              train_fraction       = 0.6,
                              overdispersion_rate  = 1.15,
                              force_practical_cv   = FALSE,
                              seed                 = NULL) {
  stopifnot(inherits(sample, "NIPTeRSample"))
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))
  stopifnot(is.numeric(chromo_focus), length(chromo_focus) == 1L,
            chromo_focus >= 1L, chromo_focus <= 22L)
  chromo_focus <- as.integer(chromo_focus)
  n_models     <- as.integer(n_models)
  n_predictors <- as.integer(n_predictors)

  is_ss <- inherits(sample, "SeparatedStrands")

  if (is_ss) {
    .regression_separated_strands(
      sample, control_group, chromo_focus, n_models, n_predictors,
      exclude_chromosomes, include_chromosomes, train_fraction,
      overdispersion_rate, force_practical_cv, seed)
  } else {
    .regression_combined_strands(
      sample, control_group, chromo_focus, n_models, n_predictors,
      exclude_chromosomes, include_chromosomes, train_fraction,
      overdispersion_rate, force_practical_cv, seed)
  }
}


# ---------------------------------------------------------------------------
# CombinedStrands regression (unchanged logic, extracted to helper)
# ---------------------------------------------------------------------------
.regression_combined_strands <- function(sample, control_group, chromo_focus,
                                         n_models, n_predictors,
                                         exclude_chromosomes,
                                         include_chromosomes,
                                         train_fraction, overdispersion_rate,
                                         force_practical_cv, seed) {
  chr_focus_key <- as.character(chromo_focus)

  # Build predictor pool
  control_chroms <- c(1:12, 14:17, 19:20, 22)
  if (!is.null(include_chromosomes)) {
    control_chroms <- sort(unique(c(control_chroms,
                                    as.integer(include_chromosomes))))
  }
  available_predictors <- setdiff(control_chroms,
                                  c(as.integer(exclude_chromosomes),
                                    chromo_focus))
  predictor_keys <- as.character(available_predictors)

  # Compute chromosomal fractions for all control samples
  n_controls <- length(control_group$samples)
  all_frac   <- .control_group_fractions(control_group)
  # all_frac: 22 x n_controls

  # Sample fractions
  sample_frac <- .sample_chr_fractions(sample)

  # Train/test split
  if (!is.null(seed)) set.seed(seed)
  n_train  <- max(2L, round(n_controls * train_fraction))
  train_idx <- sort(sample(seq_len(n_controls), n_train))
  test_idx  <- setdiff(seq_len(n_controls), train_idx)

  if (length(test_idx) < 2L) {
    stop("Control group too small for train/test split (need >= 4 samples).",
         call. = FALSE)
  }

  train_frac <- all_frac[, train_idx, drop = FALSE]
  test_frac  <- all_frac[, test_idx, drop = FALSE]

  # Total reads for focus chromosome of the test sample (for theoretical CV)
  sample_focus_reads <- sum(sample$autosomal_chromosome_reads[[1L]][chr_focus_key, ])

  # Track used predictors across models (no reuse)
  used_predictors <- character(0L)

  models_out <- vector("list", n_models)
  for (m in seq_len(n_models)) {
    remaining <- setdiff(predictor_keys, used_predictors)
    if (length(remaining) == 0L) {
      models_out <- models_out[seq_len(m - 1L)]
      break
    }

    n_pred_this <- min(n_predictors, length(remaining))

    # Forward stepwise selection on the training set
    selected <- character(0L)
    for (step in seq_len(n_pred_this)) {
      best_adj_r2  <- -Inf
      best_pred    <- NULL

      for (candidate in setdiff(remaining, selected)) {
        trial <- c(selected, candidate)
        # Build data frame for lm
        df_train <- data.frame(
          y = train_frac[chr_focus_key, ],
          t(train_frac[trial, , drop = FALSE])
        )
        fit <- stats::lm(y ~ ., data = df_train)
        adj_r2 <- summary(fit)$adj.r.squared
        if (!is.na(adj_r2) && adj_r2 > best_adj_r2) {
          best_adj_r2 <- adj_r2
          best_pred   <- candidate
        }
      }

      if (is.null(best_pred)) break
      selected <- c(selected, best_pred)
    }

    if (length(selected) == 0L) {
      models_out <- models_out[seq_len(m - 1L)]
      break
    }

    used_predictors <- c(used_predictors, selected)

    # Fit final model on training set
    df_train <- data.frame(
      y = train_frac[chr_focus_key, ],
      t(train_frac[selected, , drop = FALSE])
    )
    final_model <- stats::lm(y ~ ., data = df_train)

    # Predict on test set
    df_test <- data.frame(t(test_frac[selected, , drop = FALSE]))
    test_pred <- stats::predict(final_model, newdata = df_test)
    test_obs  <- test_frac[chr_focus_key, ]
    test_ratio <- test_obs / test_pred

    # Predict for the sample
    df_sample <- data.frame(t(sample_frac[selected]))
    colnames(df_sample) <- colnames(df_train)[-1L]
    sample_pred  <- stats::predict(final_model, newdata = df_sample)
    sample_ratio <- sample_frac[chr_focus_key] / sample_pred

    # CVs
    cv_practical   <- stats::sd(test_ratio) / mean(test_ratio)
    cv_theoretical <- overdispersion_rate / sqrt(sample_focus_reads)

    if (force_practical_cv || cv_practical >= cv_theoretical) {
      cv_used  <- cv_practical
      cv_type  <- "practical"
    } else {
      cv_used  <- cv_theoretical
      cv_type  <- "theoretical"
    }

    # Z-score
    z_sample   <- (sample_ratio - 1) / cv_used
    ctrl_z     <- (test_ratio - 1) / cv_used
    ctrl_names <- vapply(control_group$samples[test_idx],
                         function(s) s$sample_name, character(1L))
    names(ctrl_z) <- ctrl_names

    # Shapiro-Wilk
    shap_p <- if (length(ctrl_z) >= 3L && length(unique(ctrl_z)) >= 3L) {
      stats::shapiro.test(ctrl_z)$p.value
    } else {
      NA_real_
    }

    models_out[[m]] <- list(
      z_score          = unname(z_sample),
      cv               = cv_used,
      cv_type          = cv_type,
      predictors       = selected,
      shapiro_p_value  = shap_p,
      control_z_scores = ctrl_z
    )
  }

  structure(
    list(
      models            = models_out,
      focus_chromosome  = chr_focus_key,
      correction_status = sample$correction_status_autosomal,
      sample_name       = sample$sample_name
    ),
    class = "NIPTeRRegression"
  )
}


# ---------------------------------------------------------------------------
# SeparatedStrands regression
# ---------------------------------------------------------------------------
#
# Key differences from CombinedStrands:
# - Fractions matrix is 44 rows ("1F".."22F","1R".."22R") x N samples.
# - Predictor pool: "1F".."22F","1R".."22R" minus trisomy (both strands)
#   minus focus chromosome (both strands).
# - Complementary exclusion WITHIN a model: when "5F" is selected, both
#   "5F" and "5R" are excluded from remaining selections in that model.
# - ACROSS models: only the exact predictor strings from previous models
#   are excluded (the complementary strand remains available).
# - Focus chromosome response: frac[paste0(focus,"F"),] + frac[paste0(focus,"R"),]
#   (summed F+R fractions).
# - Theoretical CV uses Reduce("+", autosomal_chromosome_reads) to get
#   summed reads for the focus chromosome.

.regression_separated_strands <- function(sample, control_group, chromo_focus,
                                          n_models, n_predictors,
                                          exclude_chromosomes,
                                          include_chromosomes,
                                          train_fraction, overdispersion_rate,
                                          force_practical_cv, seed) {
  chr_focus_key <- as.character(chromo_focus)
  focus_F <- paste0(chromo_focus, "F")
  focus_R <- paste0(chromo_focus, "R")

  # Build predictor pool (both strands)
  control_chroms <- c(1:12, 14:17, 19:20, 22)
  if (!is.null(include_chromosomes)) {
    control_chroms <- sort(unique(c(control_chroms,
                                    as.integer(include_chromosomes))))
  }
  available_base <- setdiff(control_chroms,
                            c(as.integer(exclude_chromosomes), chromo_focus))
  # Expand to both strands
  predictor_keys <- c(paste0(available_base, "F"), paste0(available_base, "R"))

  # Trisomy exclusion: both strands of trisomy + focus chromosomes
  trisomy_exclude <- c(
    paste0(c(as.integer(exclude_chromosomes), chromo_focus), "F"),
    paste0(c(as.integer(exclude_chromosomes), chromo_focus), "R")
  )

  # 44-row fractions for all control samples
  n_controls <- length(control_group$samples)
  all_frac   <- .control_group_fractions(control_group)
  # all_frac: 44 x n_controls

  # Sample fractions (44-element vector)
  sample_frac <- .sample_chr_fractions(sample)

  # Focus chromosome fractions: sum F + R (matching
  # retrieve_fractions_of_interest.SeparatedStrands)
  .focus_fracs <- function(frac_mat) {
    frac_mat[focus_F, ] + frac_mat[focus_R, ]
  }
  .focus_frac_vec <- function(frac_vec) {
    frac_vec[focus_F] + frac_vec[focus_R]
  }

  # Train/test split
  if (!is.null(seed)) set.seed(seed)
  n_train   <- max(2L, round(n_controls * train_fraction))
  train_idx <- sort(sample(seq_len(n_controls), n_train))
  test_idx  <- setdiff(seq_len(n_controls), train_idx)

  if (length(test_idx) < 2L) {
    stop("Control group too small for train/test split (need >= 4 samples).",
         call. = FALSE)
  }

  train_frac <- all_frac[, train_idx, drop = FALSE]
  test_frac  <- all_frac[, test_idx, drop = FALSE]

  # Focus fractions for train/test/sample
  train_focus  <- .focus_fracs(train_frac)
  test_focus   <- .focus_fracs(test_frac)
  sample_focus <- .focus_frac_vec(sample_frac)

  # Total reads for focus chromosome (summed F+R) for theoretical CV
  summed_auto <- Reduce("+", sample$autosomal_chromosome_reads)
  rownames(summed_auto) <- as.character(1:22)
  sample_focus_reads <- sum(summed_auto[chr_focus_key, ])

  # Track used predictors across models (exact strings only — no
  # complementary exclusion across models)
  used_predictors <- character(0L)

  models_out <- vector("list", n_models)
  for (m in seq_len(n_models)) {
    # Remove trisomy, used predictors from previous models
    pool <- setdiff(predictor_keys, c(trisomy_exclude, used_predictors))
    if (length(pool) == 0L) {
      models_out <- models_out[seq_len(m - 1L)]
      break
    }

    n_pred_this <- min(n_predictors, length(pool))

    # Forward stepwise selection with complementary exclusion WITHIN model
    selected <- character(0L)
    complementary_excluded <- character(0L)

    for (step in seq_len(n_pred_this)) {
      # Remaining candidates: pool minus selected, minus complementary of
      # already-selected within this model
      candidates <- setdiff(pool, c(selected, complementary_excluded))
      if (length(candidates) == 0L) break

      best_adj_r2 <- -Inf
      best_pred   <- NULL

      for (candidate in candidates) {
        trial <- c(selected, candidate)
        df_train <- data.frame(
          y = train_focus,
          t(train_frac[trial, , drop = FALSE])
        )
        fit <- stats::lm(y ~ ., data = df_train)
        adj_r2 <- summary(fit)$adj.r.squared
        if (!is.na(adj_r2) && adj_r2 > best_adj_r2) {
          best_adj_r2 <- adj_r2
          best_pred   <- candidate
        }
      }

      if (is.null(best_pred)) break
      selected <- c(selected, best_pred)

      # Complementary exclusion: if "5F" selected, exclude both "5F" and "5R"
      comp <- unique(c(best_pred,
                       sub("F$", "R", best_pred),
                       sub("R$", "F", best_pred)))
      complementary_excluded <- unique(c(complementary_excluded, comp))
    }

    if (length(selected) == 0L) {
      models_out <- models_out[seq_len(m - 1L)]
      break
    }

    # Across models: exclude exact predictor strings (not complements)
    used_predictors <- c(used_predictors, selected)

    # Fit final model on training set
    df_train <- data.frame(
      y = train_focus,
      t(train_frac[selected, , drop = FALSE])
    )
    final_model <- stats::lm(y ~ ., data = df_train)

    # Predict on test set
    df_test <- data.frame(t(test_frac[selected, , drop = FALSE]))
    test_pred  <- stats::predict(final_model, newdata = df_test)
    test_ratio <- test_focus / test_pred

    # Predict for the sample
    df_sample <- data.frame(t(sample_frac[selected]))
    colnames(df_sample) <- colnames(df_train)[-1L]
    sample_pred  <- stats::predict(final_model, newdata = df_sample)
    sample_ratio <- sample_focus / sample_pred

    # CVs
    cv_practical   <- stats::sd(test_ratio) / mean(test_ratio)
    cv_theoretical <- overdispersion_rate / sqrt(sample_focus_reads)

    if (force_practical_cv || cv_practical >= cv_theoretical) {
      cv_used <- cv_practical
      cv_type <- "practical"
    } else {
      cv_used <- cv_theoretical
      cv_type <- "theoretical"
    }

    # Z-score
    z_sample   <- (sample_ratio - 1) / cv_used
    ctrl_z     <- (test_ratio - 1) / cv_used
    ctrl_names <- vapply(control_group$samples[test_idx],
                         function(s) s$sample_name, character(1L))
    names(ctrl_z) <- ctrl_names

    # Shapiro-Wilk
    shap_p <- if (length(ctrl_z) >= 3L && length(unique(ctrl_z)) >= 3L) {
      stats::shapiro.test(ctrl_z)$p.value
    } else {
      NA_real_
    }

    models_out[[m]] <- list(
      z_score          = unname(z_sample),
      cv               = cv_used,
      cv_type          = cv_type,
      predictors       = selected,
      shapiro_p_value  = shap_p,
      control_z_scores = ctrl_z
    )
  }

  structure(
    list(
      models            = models_out,
      focus_chromosome  = chr_focus_key,
      correction_status = sample$correction_status_autosomal,
      sample_name       = sample$sample_name
    ),
    class = "NIPTeRRegression"
  )
}
