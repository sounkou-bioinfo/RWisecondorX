#' Chromosomal Z-score
#'
#' Computes the Z-score for a given chromosome of a test sample against a
#' control group. The Z-score measures how many standard deviations the
#' sample's chromosomal fraction deviates from the control mean.
#'
#' @param sample A \code{NIPTeRSample} object (the test sample).
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param chromo_focus Integer; the chromosome to test (1-22).
#'
#' @return A list of class \code{"NIPTeRZScore"} with elements:
#' \describe{
#'   \item{sample_z_score}{The test sample's Z-score (numeric scalar).}
#'   \item{focus_chromosome}{Character; the tested chromosome.}
#'   \item{control_statistics}{Named numeric vector with \code{mean},
#'     \code{sd}, and \code{shapiro_p_value}.}
#'   \item{control_z_scores}{Named numeric vector of Z-scores for each
#'     control sample.}
#'   \item{correction_status}{Character vector of correction statuses.}
#'   \item{sample_name}{Character; the test sample name.}
#' }
#'
#' @seealso [nipter_ncv_score()], [nipter_as_control_group()],
#'   [nipter_gc_correct()]
#'
#' @examples
#' \dontrun{
#' z21 <- nipter_z_score(sample, cg, chromo_focus = 21)
#' z21$sample_z_score
#' }
#'
#' @export
nipter_z_score <- function(sample, control_group, chromo_focus) {
  stopifnot(inherits(sample, "NIPTeRSample"))
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))
  stopifnot(is.numeric(chromo_focus), length(chromo_focus) == 1L,
            chromo_focus >= 1L, chromo_focus <= 22L)
  chromo_focus <- as.integer(chromo_focus)

  chr_key <- as.character(chromo_focus)

  # Chromosomal fractions
  sample_frac  <- .sample_chr_fractions(sample)
  control_frac <- .control_group_fractions(control_group)
  # control_frac: 22 x n_controls matrix

  # Focus chromosome fractions
  sample_focus  <- sample_frac[chr_key]
  control_focus <- control_frac[chr_key, ]

  # Z-score
  ctrl_mean <- mean(control_focus)
  ctrl_sd   <- stats::sd(control_focus)
  sample_z  <- (sample_focus - ctrl_mean) / ctrl_sd

  # Control Z-scores (standardised)
  control_z <- (control_focus - ctrl_mean) / ctrl_sd
  names(control_z) <- colnames(control_frac)

  # Shapiro-Wilk normality test
  shap_p <- if (length(control_z) >= 3L && length(unique(control_z)) >= 3L) {
    stats::shapiro.test(control_z)$p.value
  } else {
    NA_real_
  }

  structure(
    list(
      sample_z_score     = unname(sample_z),
      focus_chromosome   = chr_key,
      control_statistics = c(mean = ctrl_mean, sd = ctrl_sd,
                             shapiro_p_value = shap_p),
      control_z_scores   = control_z,
      correction_status  = sample$correction_status_autosomal,
      sample_name        = sample$sample_name
    ),
    class = "NIPTeRZScore"
  )
}


#' Normalised Chromosome Value (NCV) score
#'
#' Computes the NCV score for a test sample. The NCV method selects an
#' optimal set of denominator chromosomes that minimises the coefficient of
#' variation (CV) in the control group, then uses the resulting ratio to
#' Z-score the test sample.
#'
#' @param sample A \code{NIPTeRSample} object.
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param chromo_focus Integer; the target chromosome (1-22).
#' @param max_elements Maximum number of denominator chromosomes to try
#'   (default 5). The algorithm searches all combinations of 1 to
#'   \code{max_elements} from the candidate pool.
#' @param exclude_chromosomes Integer vector of chromosomes to exclude from
#'   the denominator pool (default \code{c(13, 18, 21)} — the trisomy
#'   chromosomes).
#' @param include_chromosomes Integer vector of chromosomes to force-include
#'   in the candidate pool. Default \code{NULL}.
#'
#' @return A list of class \code{"NIPTeRNCV"} with elements:
#' \describe{
#'   \item{sample_score}{The NCV score for the test sample (numeric).}
#'   \item{focus_chromosome}{Character; the tested chromosome.}
#'   \item{denominators}{Integer vector of selected denominator chromosomes.}
#'   \item{control_statistics}{Named numeric vector with \code{mean},
#'     \code{sd}, and \code{shapiro_p_value}.}
#'   \item{control_z_scores}{Named numeric vector of control NCV scores.}
#'   \item{best_cv}{The minimum CV achieved.}
#'   \item{correction_status}{Character vector.}
#'   \item{sample_name}{Character.}
#' }
#'
#' @details
#' The candidate denominator pool is: chromosomes 1-12, 14-17, 19-20, 22
#' (the "control chromosomes") minus \code{exclude_chromosomes} minus
#' \code{chromo_focus}, plus any \code{include_chromosomes}.
#'
#' For each combination size from 1 to \code{max_elements}, all
#' \code{choose(n_candidates, size)} combinations are evaluated. The ratio
#' \code{reads[focus] / sum(reads[denominators])} is computed per control
#' sample, and the combination minimising the coefficient of variation
#' (\code{sd/mean}) is selected.
#'
#' @seealso [nipter_z_score()], [nipter_as_control_group()]
#'
#' @examples
#' \dontrun{
#' ncv21 <- nipter_ncv_score(sample, cg, chromo_focus = 21, max_elements = 5)
#' ncv21$sample_score
#' ncv21$denominators
#' }
#'
#' @export
nipter_ncv_score <- function(sample,
                             control_group,
                             chromo_focus,
                             max_elements          = 5L,
                             exclude_chromosomes   = c(13L, 18L, 21L),
                             include_chromosomes   = NULL) {
  stopifnot(inherits(sample, "NIPTeRSample"))
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))
  stopifnot(is.numeric(chromo_focus), length(chromo_focus) == 1L,
            chromo_focus >= 1L, chromo_focus <= 22L)
  chromo_focus <- as.integer(chromo_focus)
  max_elements <- as.integer(max_elements)

  # Build candidate denominator pool
  # Default control chromosomes: 1-12, 14-17, 19-20, 22
  control_chroms <- c(1:12, 14:17, 19:20, 22)
  if (!is.null(include_chromosomes)) {
    control_chroms <- sort(unique(c(control_chroms,
                                    as.integer(include_chromosomes))))
  }
  candidates <- setdiff(control_chroms,
                        c(as.integer(exclude_chromosomes), chromo_focus))

  if (length(candidates) == 0L) {
    stop("No candidate denominator chromosomes remain after exclusions.",
         call. = FALSE)
  }
  max_elements <- min(max_elements, length(candidates))

  # Per-chromosome total reads for each control sample
  n_controls <- length(control_group$samples)
  ctrl_reads <- vapply(control_group$samples, function(s) {
    .sample_chr_reads(s)
  }, numeric(22L))
  # ctrl_reads: 22 x n_controls matrix, rownames "1"-"22"

  # Same for the test sample
  sample_reads <- .sample_chr_reads(sample)

  chr_focus_key <- as.character(chromo_focus)

  # Brute-force search over all combinations
  best_cv   <- Inf
  best_denom <- NULL

  for (n_el in seq_len(max_elements)) {
    combos <- utils::combn(candidates, n_el, simplify = FALSE)
    for (combo in combos) {
      combo_keys <- as.character(combo)
      # Ratio = focus / sum(denominators) for each control
      denom_sums <- colSums(ctrl_reads[combo_keys, , drop = FALSE])
      ratios <- ctrl_reads[chr_focus_key, ] / denom_sums
      cv <- stats::sd(ratios) / mean(ratios)
      if (cv < best_cv) {
        best_cv    <- cv
        best_denom <- combo
      }
    }
  }

  # Compute NCV statistics using the best denominators
  denom_keys   <- as.character(best_denom)
  ctrl_denom   <- colSums(ctrl_reads[denom_keys, , drop = FALSE])
  ctrl_ratios  <- ctrl_reads[chr_focus_key, ] / ctrl_denom
  ctrl_mean    <- mean(ctrl_ratios)
  ctrl_sd      <- stats::sd(ctrl_ratios)

  # Control NCV Z-scores
  ctrl_z <- (ctrl_ratios - ctrl_mean) / ctrl_sd
  names(ctrl_z) <- vapply(control_group$samples,
                          function(s) s$sample_name, character(1L))

  # Sample NCV score
  sample_denom <- sum(sample_reads[denom_keys])
  sample_ratio <- sample_reads[chr_focus_key] / sample_denom
  sample_score <- (sample_ratio - ctrl_mean) / ctrl_sd

  # Shapiro-Wilk
  shap_p <- if (length(ctrl_z) >= 3L && length(unique(ctrl_z)) >= 3L) {
    stats::shapiro.test(ctrl_z)$p.value
  } else {
    NA_real_
  }

  structure(
    list(
      sample_score       = unname(sample_score),
      focus_chromosome   = chr_focus_key,
      denominators       = as.integer(best_denom),
      control_statistics = c(mean = ctrl_mean, sd = ctrl_sd,
                             shapiro_p_value = shap_p),
      control_z_scores   = ctrl_z,
      best_cv            = best_cv,
      correction_status  = sample$correction_status_autosomal,
      sample_name        = sample$sample_name
    ),
    class = "NIPTeRNCV"
  )
}
