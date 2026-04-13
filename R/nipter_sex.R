#' Build a sex prediction model from a NIPTeR control group
#'
#' Fits a two-component Gaussian mixture model (GMM) on sex chromosome
#' fractions derived from a \code{NIPTeRControlGroup}. The model
#' distinguishes male and female samples based on Y-chromosome read fraction
#' (univariate) or X+Y chromosome fractions (bivariate).
#'
#' This mirrors the approach used in clinical NIPT pipelines (see
#' \emph{Details}), ported to R using \code{mclust::Mclust()} in place of
#' Python's \code{sklearn.GaussianMixture}.
#'
#' @param control_group A \code{NIPTeRControlGroup} object. Ideally
#'   chi-squared corrected via \code{\link{nipter_chi_correct}}.
#' @param method Character; the feature space for the GMM:
#'   \describe{
#'     \item{\code{"y_fraction"}}{Univariate: Y-chromosome read fraction
#'       relative to total autosomal reads.}
#'     \item{\code{"xy_fraction"}}{Bivariate: X and Y chromosome fractions.}
#'   }
#'
#' @return An object of class \code{"NIPTeRSexModel"} with elements:
#' \describe{
#'   \item{model}{The \code{mclust::Mclust} fitted object.}
#'   \item{method}{Character; the method used.}
#'   \item{male_cluster}{Integer (1 or 2); which cluster is male
#'     (higher Y fraction).}
#'   \item{classifications}{Named character vector of \code{"male"}/
#'     \code{"female"} labels for each control sample.}
#'   \item{fractions}{Matrix of the input features used for fitting
#'     (samples as rows). Columns are \code{"y_fraction"} or
#'     \code{c("x_fraction", "y_fraction")}.}
#' }
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Compute sex chromosome read fractions for every sample in the
#'     control group. Y fraction = \code{sum(Y bins) / sum(autosomal bins)};
#'     X fraction = \code{sum(X bins) / sum(autosomal bins)}.
#'   \item Fit a two-component Gaussian mixture with equal mixing proportions
#'     (\code{mclust::Mclust(data, G = 2,
#'     control = mclust::emControl(equalPro = TRUE))}).
#'   \item Identify the male cluster as the component with the higher median
#'     Y fraction.
#' }
#'
#' This follows the clinical NIPT pipeline pattern of building sex models
#' from control cohorts. The user's pipeline builds three models (Y-unique
#' ratio from samtools, XY fractions, Y fraction) and takes a majority vote;
#' \code{nipter_predict_sex()} implements the consensus when given multiple
#' models.
#'
#' @seealso [nipter_predict_sex()], [nipter_chi_correct()]
#'
#' @examples
#' \dontrun{
#' # Build sex model from chi-corrected control group
#' sex_model <- nipter_sex_model(chi_cg, method = "y_fraction")
#'
#' # Bivariate X+Y model
#' sex_model_xy <- nipter_sex_model(chi_cg, method = "xy_fraction")
#' }
#'
#' @export
nipter_sex_model <- function(control_group,
                             method = c("y_fraction", "xy_fraction")) {
  method <- match.arg(method)
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))

  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop(
      "Package 'mclust' is required for sex prediction. ",
      "Install it with: install.packages('mclust')",
      call. = FALSE
    )
  }

  n_samples <- length(control_group$samples)
  if (n_samples < 4L) {
    stop("At least 4 samples are needed to fit a sex model.", call. = FALSE)
  }

  # Compute sex chromosome fractions
  fracs <- .sex_chr_fractions(control_group)
  # fracs: n_samples x 2 matrix with columns "x_fraction", "y_fraction"

  # Select features for the GMM
  if (method == "y_fraction") {
    fit_data <- fracs[, "y_fraction"]
  } else {
    fit_data <- fracs[, c("x_fraction", "y_fraction"), drop = FALSE]
  }

  # Fit 2-component GMM with equal mixing proportions.
  # NOTE: mclust::Mclust() internally calls mclustBIC() without namespace
  # qualification, so a bare mclust::Mclust() fails when the package is loaded
  # via requireNamespace() but not attached. We evaluate the call inside the
  # mclust namespace to work around this.
  gmm <- .mclust_fit(fit_data, G = 2L, equalPro = TRUE)

  if (is.null(gmm)) {
    stop("mclust::Mclust() failed to converge. ",
         "The control group may lack sufficient sex chromosome signal.",
         call. = FALSE)
  }

  # Identify the male cluster (higher median Y fraction)
  y_fracs <- fracs[, "y_fraction"]
  median_1 <- stats::median(y_fracs[gmm$classification == 1L], na.rm = TRUE)
  median_2 <- stats::median(y_fracs[gmm$classification == 2L], na.rm = TRUE)
  male_cluster <- if (median_1 > median_2) 1L else 2L

  # Build labels
  sample_names <- vapply(control_group$samples,
                         function(s) s$sample_name, character(1L))
  labels <- ifelse(gmm$classification == male_cluster, "male", "female")
  names(labels) <- sample_names

  structure(
    list(
      model          = gmm,
      method         = method,
      male_cluster   = male_cluster,
      classifications = labels,
      fractions      = fracs
    ),
    class = "NIPTeRSexModel"
  )
}


#' Predict fetal sex for a NIPTeR sample
#'
#' Uses one or more \code{NIPTeRSexModel} objects to classify a sample as
#' male or female. When multiple models are supplied, a majority vote
#' determines the consensus call.
#'
#' @param sample A \code{NIPTeRSample} object.
#' @param ... One or more \code{NIPTeRSexModel} objects, or a single list
#'   of models.
#'
#' @return A list of class \code{"NIPTeRSexPrediction"} with elements:
#' \describe{
#'   \item{prediction}{Character; \code{"male"} or \code{"female"} (consensus
#'     when multiple models).}
#'   \item{model_predictions}{Named character vector of per-model predictions
#'     (names are the method names).}
#'   \item{fractions}{Named numeric vector with \code{x_fraction} and
#'     \code{y_fraction} for the sample.}
#'   \item{sample_name}{Character; the sample name.}
#' }
#'
#' @details
#' For each model, the sample's sex chromosome fractions are computed and
#' classified using \code{predict(model$model, newdata = ...)}. The
#' classification is mapped to \code{"male"}/\code{"female"} using the
#' model's \code{male_cluster}.
#'
#' With multiple models, the consensus is the label that appears most
#' frequently (majority vote). In case of a tie, \code{"female"} is
#' returned (conservative default for NIPT, as false-male calls could mask
#' sex chromosome aneuploidies).
#'
#' @seealso [nipter_sex_model()]
#'
#' @examples
#' \dontrun{
#' # Single model
#' pred <- nipter_predict_sex(test_sample, sex_model)
#' pred$prediction
#'
#' # Consensus from two models (majority vote)
#' pred <- nipter_predict_sex(test_sample, model_y, model_xy)
#' pred$prediction
#' }
#'
#' @export
nipter_predict_sex <- function(sample, ...) {
  stopifnot(inherits(sample, "NIPTeRSample"))

  models <- list(...)
  # Allow a single list of models

  if (length(models) == 1L && is.list(models[[1L]]) &&
      !inherits(models[[1L]], "NIPTeRSexModel")) {
    models <- models[[1L]]
  }

  if (length(models) == 0L) {
    stop("At least one NIPTeRSexModel must be provided.", call. = FALSE)
  }

  for (i in seq_along(models)) {
    if (!inherits(models[[i]], "NIPTeRSexModel")) {
      stop("All models must be NIPTeRSexModel objects.", call. = FALSE)
    }
  }

  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop(
      "Package 'mclust' is required for sex prediction. ",
      "Install it with: install.packages('mclust')",
      call. = FALSE
    )
  }

  # Compute sample's sex chromosome fractions
  sample_fracs <- .sample_sex_fractions(sample)

  model_preds <- character(length(models))
  names(model_preds) <- vapply(models, function(m) m$method, character(1L))

  for (i in seq_along(models)) {
    mdl <- models[[i]]

    if (mdl$method == "y_fraction") {
      newdata <- sample_fracs["y_fraction"]
    } else {
      newdata <- matrix(sample_fracs[c("x_fraction", "y_fraction")],
                        nrow = 1L)
      colnames(newdata) <- c("x_fraction", "y_fraction")
    }

    pred <- stats::predict(mdl$model, newdata = newdata)
    cluster <- pred$classification[1L]
    model_preds[i] <- if (cluster == mdl$male_cluster) "male" else "female"
  }

  # Consensus: majority vote (tie goes to "female" — conservative for NIPT)
  n_male   <- sum(model_preds == "male")
  n_female <- sum(model_preds == "female")
  consensus <- if (n_male > n_female) "male" else "female"

  structure(
    list(
      prediction       = consensus,
      model_predictions = model_preds,
      fractions        = sample_fracs,
      sample_name      = sample$sample_name
    ),
    class = "NIPTeRSexPrediction"
  )
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Compute X and Y chromosome fractions for all samples in a control group.
# Returns an n_samples x 2 matrix with columns "x_fraction", "y_fraction".
# Each fraction = sum(sex_chr bins) / sum(autosomal bins).
.sex_chr_fractions <- function(control_group) {
  n <- length(control_group$samples)
  mat <- matrix(NA_real_, nrow = n, ncol = 2L)
  colnames(mat) <- c("x_fraction", "y_fraction")

  for (i in seq_len(n)) {
    mat[i, ] <- .sample_sex_fractions(control_group$samples[[i]])
  }

  rownames(mat) <- vapply(control_group$samples,
                          function(s) s$sample_name, character(1L))
  mat
}


# Compute X and Y chromosome fractions for a single NIPTeRSample.
# Returns named numeric vector c(x_fraction = ..., y_fraction = ...).
# Fractions are relative to total autosomal reads (sum of all autosomal bins).
.sample_sex_fractions <- function(sample) {
  # Autosomal total
  auto <- sample$autosomal_chromosome_reads
  if (inherits(sample, "SeparatedStrands")) {
    auto_total <- sum(Reduce("+", auto))
  } else {
    auto_total <- sum(auto[[1L]])
  }

  if (auto_total == 0) {
    return(c(x_fraction = 0, y_fraction = 0))
  }

  # Sex chromosome totals
  sex <- sample$sex_chromosome_reads
  if (inherits(sample, "SeparatedStrands")) {
    sex_summed <- Reduce("+", sex)
    # rownames: "XF"+"XR" → sum to get X, "YF"+"YR" → sum to get Y
    # After Reduce, rownames come from the first matrix ("XF", "YF"),
    # but values are summed. Rows: row 1 = X total, row 2 = Y total.
    x_total <- sum(sex_summed[1L, ])
    y_total <- sum(sex_summed[2L, ])
  } else {
    x_total <- sum(sex[[1L]]["X", ])
    y_total <- sum(sex[[1L]]["Y", ])
  }

  c(x_fraction = x_total / auto_total,
    y_fraction = y_total / auto_total)
}


# Wrapper around mclust::Mclust() that evaluates the call inside the mclust
# namespace.  mclust::Mclust() internally calls mclustBIC() without its own
# namespace prefix, so the bare mclust::Mclust() path fails when the package
# is loaded via requireNamespace() but not attached to the search path.
# Evaluating in the namespace makes internal symbol lookup succeed, which is
# the CRAN-friendly alternative to library()/require() inside a function.
.mclust_fit <- function(data, G = 2L, equalPro = TRUE) {
  ns <- asNamespace("mclust")
  evalq(
    Mclust(data, G = G, control = emControl(equalPro = equalPro)),
    envir = list2env(list(data = data, G = G, equalPro = equalPro),
                     parent = ns)
  )
}
