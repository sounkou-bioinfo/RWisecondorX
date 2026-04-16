#' Chi-squared correction for overdispersed bins
#'
#' Identifies bins in the control group with excess variance (overdispersion)
#' using a chi-squared test and downweights them. This reduces the influence
#' of noisy bins on downstream Z-scores and NCV calculations.
#'
#' The correction is applied simultaneously to both the test sample and all
#' control group samples, maintaining consistency.
#'
#' @param sample A \code{NIPTeRSample} object (the test sample).
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param chi_cutoff Normalised chi-squared threshold. Bins with
#'   \code{(chi - df) / sqrt(2*df) > chi_cutoff} are corrected. Default
#'   \code{3.5}.
#' @param include_sex Logical; correct sex chromosomes as well? Default
#'   \code{FALSE}. Sex chromosome bins use the same chi-squared weights
#'   derived from autosomes.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{sample}{The corrected \code{NIPTeRSample}.}
#'   \item{control_group}{The corrected \code{NIPTeRControlGroup}.}
#' }
#'
#' @details
#' The algorithm follows NIPTeR's chi-squared correction:
#'
#' 1. For each control sample, scale autosomal bin counts so that total
#'    reads match the overall mean across all control samples.
#' 2. Compute the expected count per bin (mean of scaled counts).
#' 3. Compute chi-squared per bin:
#'    \eqn{\chi^2 = \sum_i (expected - scaled_i)^2 / expected}.
#' 4. Normalise: \eqn{z = (\chi^2 - df) / \sqrt{2 \cdot df}} where
#'    \eqn{df = n_{controls} - 1}.
#' 5. For bins where \eqn{z > \text{chi\_cutoff}}: divide reads by
#'    \eqn{\chi^2 / df}.
#'
#' @seealso [nipter_gc_correct()], [nipter_z_score()],
#'   [nipter_as_control_group()]
#'
#' @examples
#' \dontrun{
#' result <- nipter_chi_correct(sample, control_group)
#' sample_corrected <- result$sample
#' cg_corrected     <- result$control_group
#' }
#'
#' @export
nipter_chi_correct <- function(sample,
                               control_group,
                               chi_cutoff  = 3.5,
                               include_sex = FALSE) {
  stopifnot(.is_nipt_sample_object(sample))
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.numeric(chi_cutoff), length(chi_cutoff) == 1L)
  stopifnot(is.logical(include_sex), length(include_sex) == 1L)

  # Strand-type compatibility guard
  sample_st <- .strand_type_of(sample)
  cg_st     <- .strand_type_of(control_group)
  if (!identical(sample_st, cg_st)) {
    stop(sprintf(
      "Strand type mismatch: sample is '%s' but control_group is '%s'.",
      sample_st, cg_st
    ), call. = FALSE)
  }

  if (isTRUE(include_sex)) {
    stop(
      "include_sex = TRUE deviates from NIPTeR upstream (chi_correct always uses ",
      "include_XY = FALSE in the production pipeline). ",
      "For GC correction of sex chromosomes use nipter_gc_correct(include_sex = TRUE).",
      call. = FALSE
    )
  }

  profile <- .chi_profile(control_group, chi_cutoff = chi_cutoff)

  # Apply correction to the test sample
  sample <- .chi_correct_sample(sample, profile$correction_factor, profile$n_bins)

  # Apply correction to all control samples
  control_group <- .control_group_with_samples(
    control_group,
    lapply(control_group$samples,
           .chi_correct_sample,
           correction_factor = profile$correction_factor,
           n_bins = profile$n_bins)
  )

  list(sample = sample, control_group = control_group)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.chi_profile <- function(control_group,
                         chi_cutoff = 3.5,
                         include_bin_stats = FALSE) {
  n_controls <- length(control_group$samples)
  n_bins <- ncol(autosomal_matrix(control_group$samples[[1L]]))

  # Flatten each control sample's autosomal reads into a single vector.
  # autosomal_matrix() handles both S7 and S3 objects and sums strands for
  # SeparatedStrands, matching upstream NIPTeR's chi_correct behaviour.
  ctrl_flat <- lapply(control_group$samples, function(s) {
    as.numeric(t(autosomal_matrix(s)))
  })

  total_bins <- 22L * n_bins
  df <- n_controls - 1L
  ctrl_totals <- vapply(ctrl_flat, sum, numeric(1L))
  overall_mean <- mean(ctrl_totals)
  scaling <- overall_mean / ctrl_totals

  scaled <- lapply(seq_len(n_controls), function(i) {
    ctrl_flat[[i]] * scaling[[i]]
  })

  expected <- Reduce(`+`, scaled) / n_controls
  chi_sum <- numeric(total_bins)
  for (i in seq_len(n_controls)) {
    chi_sum <- chi_sum + (expected - scaled[[i]])^2 / pmax(expected, .Machine$double.eps)
  }

  chi_norm <- (chi_sum - df) / sqrt(2 * df)
  overdispersed <- chi_norm > chi_cutoff
  correction_factor <- rep(1.0, total_bins)
  correction_factor[overdispersed] <- chi_sum[overdispersed] / df

  out <- list(
    n_controls = n_controls,
    n_bins = n_bins,
    df = df,
    chi_cutoff = chi_cutoff,
    expected = expected,
    chi_square = chi_sum,
    chi_z = chi_norm,
    overdispersed = overdispersed,
    correction_factor = correction_factor
  )

  if (isTRUE(include_bin_stats)) {
    scaled_mat <- do.call(cbind, scaled)
    sd_scaled <- apply(scaled_mat, 1L, stats::sd)
    mean_scaled <- expected
    cv_scaled <- rep(NA_real_, length(mean_scaled))
    ok <- is.finite(mean_scaled) & abs(mean_scaled) > .Machine$double.eps
    cv_scaled[ok] <- 100 * sd_scaled[ok] / mean_scaled[ok]
    out$scaled_mean <- mean_scaled
    out$scaled_sd <- sd_scaled
    out$scaled_cv <- cv_scaled
  }

  out
}

# Apply chi-squared correction factors to a single NIPTeRSample.
# For SeparatedStrands, applies the same correction factors (derived from
# summed F+R) independently to each strand matrix via lapply, matching
# upstream NIPTeR's correctsamples pattern.
.chi_correct_sample <- function(sample, correction_factor, n_bins) {
  auto_list <- .sample_autosomal_reads(sample)

  corrected_auto <- lapply(auto_list, function(mat) {
    corrected <- mat
    for (i in seq_len(22L)) {
      offset <- (i - 1L) * n_bins
      idx <- seq(offset + 1L, offset + n_bins)
      corrected[i, ] <- mat[i, ] / correction_factor[idx]
    }
    corrected
  })

  sample <- .sample_with_reads(sample, autosomal = corrected_auto)
  sample <- .sample_append_correction_step(
    sample,
    "autosomal",
    .nipt_chi_correction_step()
  )

  sample
}
