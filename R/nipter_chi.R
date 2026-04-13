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
  stopifnot(inherits(sample, "NIPTeRSample"))
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))
  stopifnot(is.numeric(chi_cutoff), length(chi_cutoff) == 1L)
  stopifnot(is.logical(include_sex), length(include_sex) == 1L)

  n_controls <- length(control_group$samples)
  df         <- n_controls - 1L
  n_bins     <- ncol(sample$autosomal_chromosome_reads[[1L]])

  # Flatten each control sample's autosomal reads into a single vector.
  # For SeparatedStrands, sum F+R via Reduce first (matching upstream NIPTeR's
  # chi_correct which does Reduce("+", x[[autosomal_chromosome_reads]])).
  ctrl_flat <- lapply(control_group$samples, function(s) {
    auto <- s$autosomal_chromosome_reads
    if (inherits(s, "SeparatedStrands")) {
      summed <- Reduce("+", auto)
      rownames(summed) <- as.character(1:22)
      return(as.numeric(t(summed)))
    }
    as.numeric(t(auto[[1L]]))
  })

  total_bins <- 22L * n_bins

  # Overall mean of total reads across control samples
  ctrl_totals <- vapply(ctrl_flat, sum, numeric(1L))
  overall_mean <- mean(ctrl_totals)

  # Scale each control to equalise totals
  scaled <- lapply(seq_len(n_controls), function(i) {
    ctrl_flat[[i]] * (overall_mean / ctrl_totals[i])
  })

  # Expected per bin: mean of scaled values
  expected <- Reduce(`+`, scaled) / n_controls

  # Chi-squared per bin: sum over controls of (expected - scaled_i)^2 / expected
  chi_sum <- numeric(total_bins)
  for (i in seq_len(n_controls)) {
    chi_sum <- chi_sum + (expected - scaled[[i]])^2 / pmax(expected, .Machine$double.eps)
  }

  # Normalise chi scores to approximate standard normal
  chi_norm <- (chi_sum - df) / sqrt(2 * df)

  # Identify overdispersed bins
  overdispersed <- chi_norm > chi_cutoff
  correction_factor <- rep(1.0, total_bins)
  correction_factor[overdispersed] <- chi_sum[overdispersed] / df

  # Apply correction to the test sample
  sample <- .chi_correct_sample(sample, correction_factor, n_bins,
                                include_sex)

  # Apply correction to all control samples
  control_group$samples <- lapply(control_group$samples,
                                  .chi_correct_sample,
                                  correction_factor = correction_factor,
                                  n_bins = n_bins,
                                  include_sex = include_sex)

  # Update control group correction statuses
  control_group$correction_status_autosomal <- unique(unlist(lapply(
    control_group$samples, `[[`, "correction_status_autosomal"
  )))
  if (include_sex) {
    control_group$correction_status_sex <- unique(unlist(lapply(
      control_group$samples, `[[`, "correction_status_sex"
    )))
  }

  list(sample = sample, control_group = control_group)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Apply chi-squared correction factors to a single NIPTeRSample.
# For SeparatedStrands, applies the same correction factors (derived from
# summed F+R) independently to each strand matrix via lapply, matching
# upstream NIPTeR's correctsamples pattern.
.chi_correct_sample <- function(sample, correction_factor, n_bins,
                                include_sex) {
  auto_list <- sample$autosomal_chromosome_reads

  corrected_auto <- lapply(auto_list, function(mat) {
    corrected <- mat
    for (i in seq_len(22L)) {
      offset <- (i - 1L) * n_bins
      idx <- seq(offset + 1L, offset + n_bins)
      corrected[i, ] <- mat[i, ] / correction_factor[idx]
    }
    corrected
  })

  sample$autosomal_chromosome_reads <- corrected_auto
  sample$correction_status_autosomal <- .update_correction_status(
    sample$correction_status_autosomal, "Chi square corrected"
  )

  if (include_sex) {
    corrected_sex <- lapply(sample$sex_chromosome_reads, function(sex_mat) {
      corrected <- sex_mat
      # Sex chromosomes use the same autosomal chi weights — the chi
      # distribution is derived from autosomes. We apply the mean correction
      # factor (across autosomal bins at comparable positions) to sex bins.
      for (b in seq_len(n_bins)) {
        auto_idx <- seq(b, by = n_bins, length.out = 22L)
        mean_cf  <- mean(correction_factor[auto_idx])
        corrected[, b] <- sex_mat[, b] / mean_cf
      }
      corrected
    })

    sample$sex_chromosome_reads <- corrected_sex
    sample$correction_status_sex <- .update_correction_status(
      sample$correction_status_sex, "Chi square corrected"
    )
  }

  sample
}
