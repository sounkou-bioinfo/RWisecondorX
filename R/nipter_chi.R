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
  bin_layout <- .chi_bin_layout(control_group, n_bins)

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
  expected_finite <- is.finite(expected)
  expected_positive <- expected > 0
  scaled_finite <- rep(TRUE, total_bins)
  for (i in seq_len(n_controls)) {
    scaled_finite <- scaled_finite & is.finite(scaled[[i]])
  }
  valid_bins <- expected_finite & expected_positive & scaled_finite
  if (!all(is.na(bin_layout$in_chromosome_range))) {
    valid_bins <- valid_bins & bin_layout$in_chromosome_range
  }

  invalid_reason_detail <- rep("valid", total_bins)
  invalid_reason <- rep("valid", total_bins)
  invalid_idx <- which(!valid_bins)
  if (length(invalid_idx)) {
    detail_parts <- vector("list", length = total_bins)
    if (!all(is.na(bin_layout$in_chromosome_range))) {
      idx <- invalid_idx[!bin_layout$in_chromosome_range[invalid_idx]]
      if (length(idx)) {
        for (j in idx) detail_parts[[j]] <- c(detail_parts[[j]], "out_of_range")
      }
    }
    idx <- invalid_idx[!scaled_finite[invalid_idx]]
    if (length(idx)) {
      for (j in idx) detail_parts[[j]] <- c(detail_parts[[j]], "nonfinite_scaled")
    }
    idx <- invalid_idx[!expected_finite[invalid_idx]]
    if (length(idx)) {
      for (j in idx) detail_parts[[j]] <- c(detail_parts[[j]], "nonfinite_expected")
    }
    idx <- invalid_idx[expected_finite[invalid_idx] & !expected_positive[invalid_idx]]
    if (length(idx)) {
      for (j in idx) detail_parts[[j]] <- c(detail_parts[[j]], "nonpositive_expected")
    }

    invalid_reason_detail[invalid_idx] <- vapply(invalid_idx, function(i) {
      parts <- unique(detail_parts[[i]])
      if (!length(parts)) "invalid_other" else paste(parts, collapse = ";")
    }, character(1L))
    invalid_reason[invalid_idx] <- vapply(invalid_idx, function(i) {
      parts <- unique(detail_parts[[i]])
      if (!length(parts)) {
        return("invalid_other")
      }
      if ("out_of_range" %in% parts) {
        return("out_of_range")
      }
      if ("nonfinite_scaled" %in% parts) {
        return("nonfinite_scaled")
      }
      if ("nonfinite_expected" %in% parts) {
        return("nonfinite_expected")
      }
      if ("nonpositive_expected" %in% parts) {
        return("nonpositive_expected")
      }
      "invalid_other"
    }, character(1L))
  }

  chi_sum <- rep(NA_real_, total_bins)
  chi_norm <- rep(NA_real_, total_bins)
  if (any(valid_bins)) {
    chi_sum[valid_bins] <- 0
    denom <- expected[valid_bins]
    for (i in seq_len(n_controls)) {
      chi_sum[valid_bins] <- chi_sum[valid_bins] +
        (expected[valid_bins] - scaled[[i]][valid_bins])^2 / denom
    }
    chi_norm[valid_bins] <- (chi_sum[valid_bins] - df) / sqrt(2 * df)
  }
  overdispersed <- rep(FALSE, total_bins)
  overdispersed[valid_bins] <- chi_norm[valid_bins] > chi_cutoff
  correction_factor <- rep(1.0, total_bins)
  correction_factor[overdispersed] <- chi_sum[overdispersed] / df

  out <- list(
    n_controls = n_controls,
    n_bins = n_bins,
    df = df,
    chi_cutoff = chi_cutoff,
    binsize = bin_layout$binsize[[1L]],
    chromosome_length_bp = bin_layout$chromosome_length_bp,
    bin_start_bp = bin_layout$bin_start_bp,
    bin_end_bp = bin_layout$bin_end_bp,
    in_chromosome_range = bin_layout$in_chromosome_range,
    expected = expected,
    expected_finite = expected_finite,
    expected_positive = expected_positive,
    scaled_finite = scaled_finite,
    chi_square = chi_sum,
    chi_z = chi_norm,
    valid_bins = valid_bins,
    invalid_reason = invalid_reason,
    invalid_reason_detail = invalid_reason_detail,
    overdispersed = overdispersed,
    correction_factor = correction_factor
  )

  if (isTRUE(include_bin_stats)) {
    scaled_mat <- do.call(cbind, scaled)
    sd_scaled <- apply(scaled_mat, 1L, stats::sd)
    mean_scaled <- expected
    cv_scaled <- rep(NA_real_, length(mean_scaled))
    n_nonzero_scaled <- rowSums(scaled_mat != 0)
    ok <- valid_bins & is.finite(mean_scaled) & abs(mean_scaled) > .Machine$double.eps
    cv_scaled[ok] <- 100 * sd_scaled[ok] / mean_scaled[ok]
    out$scaled_mean <- mean_scaled
    out$scaled_sd <- sd_scaled
    out$scaled_cv <- cv_scaled
    out$scaled_n_nonzero <- as.integer(n_nonzero_scaled)
  }

  out
}

.chi_bin_layout <- function(control_group, n_bins) {
  template <- control_group$samples[[1L]]
  chr_names <- as.character(1:22)
  binsize <- .sample_binsize(template)
  if (!is.finite(binsize) || length(binsize) != 1L || is.na(binsize) || binsize < 1L) {
    binsize <- NA_integer_
  } else {
    binsize <- as.integer(binsize)
  }

  chrom_lengths <- .sample_chrom_lengths(template)
  chr_lengths <- rep(NA_integer_, length(chr_names))
  names(chr_lengths) <- chr_names
  if (!is.null(chrom_lengths)) {
    matched <- match(chr_names, names(chrom_lengths))
    keep <- !is.na(matched)
    chr_lengths[keep] <- as.integer(chrom_lengths[matched[keep]])
  }

  chromosome <- rep(chr_names, each = n_bins)
  bin <- rep(seq_len(n_bins), times = length(chr_names))
  chromosome_length_bp <- rep(chr_lengths, each = n_bins)
  if (is.na(binsize)) {
    bin_start_bp <- rep(NA_integer_, length(bin))
    bin_end_bp <- rep(NA_integer_, length(bin))
    in_chromosome_range <- rep(NA, length(bin))
  } else {
    bin_offsets <- as.integer((seq_len(n_bins) - 1L) * binsize)
    bin_start_bp <- rep(bin_offsets, times = length(chr_names))
    in_chromosome_range <- !is.na(chromosome_length_bp) & bin_start_bp < chromosome_length_bp
    bin_end_bp <- rep(NA_integer_, length(bin))
    keep <- which(in_chromosome_range)
    if (length(keep)) {
      bin_end_bp[keep] <- pmin(bin_start_bp[keep] + binsize, chromosome_length_bp[keep])
    }
  }

  data.frame(
    chromosome = chromosome,
    bin = as.integer(bin),
    binsize = rep(binsize, length(bin)),
    chromosome_length_bp = chromosome_length_bp,
    bin_start_bp = bin_start_bp,
    bin_end_bp = bin_end_bp,
    in_chromosome_range = in_chromosome_range,
    stringsAsFactors = FALSE
  )
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
