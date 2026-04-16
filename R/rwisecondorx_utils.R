# WisecondorX native R implementation — shared utilities
#
# Ports overall_tools.py and parts of newref_tools.py from the upstream
# WisecondorX Python package. Original authors: Lennart Raman, Roy Straver,
# Wim Audenaert. Ported to R with credit under GPL-3.


#' Scale a binned sample to a new bin size
#'
#' Rescales per-chromosome bin-count vectors by aggregating bins.
#' The new bin size must be a positive integer multiple of the original.
#' This is the R equivalent of `overall_tools.scale_sample()` in the upstream
#' WisecondorX Python package.
#'
#' @param sample Named list of integer vectors keyed by chromosome
#'   (`"1"`--`"24"`), as returned by [bam_convert()].
#' @param from_size Current bin size in base pairs.
#' @param to_size Target bin size in base pairs. Must be an integer multiple
#'   of `from_size`.
#'
#' @return Named list structured like `sample` but with bins aggregated
#'   to the new size.
#'
#' @examples
#' \dontrun{
#' bins_5k  <- bam_convert("sample.bam", binsize = 5000L)
#' bins_100k <- scale_sample(bins_5k, from_size = 5000, to_size = 100000)
#' }
#'
#' @export
scale_sample <- function(sample, from_size, to_size) {
  from_size <- as.numeric(from_size)
  to_size   <- as.numeric(to_size)
  stopifnot(is.list(sample))
  stopifnot(from_size > 0, to_size > 0)

  if (from_size == to_size) return(sample)

  if (to_size < from_size || (to_size %% from_size) != 0) {
    stop(sprintf("Impossible binsize scaling: %d to %d", as.integer(from_size),
                 as.integer(to_size)), call. = FALSE)
  }

  scale <- as.integer(to_size / from_size)
  result <- vector("list", length(sample))
  names(result) <- names(sample)

  for (chr_name in names(sample)) {
    chr_data <- sample[[chr_name]]
    if (is.null(chr_data)) next
    n_src   <- length(chr_data)
    new_len <- as.integer(ceiling(n_src / scale))
    # Pad to an exact multiple of scale so matrix reshape works cleanly,
    # then sum each column (= group of `scale` original bins).
    padded_len <- new_len * scale
    padded <- integer(padded_len)
    padded[seq_len(n_src)] <- chr_data
    result[[chr_name]] <- as.integer(colSums(matrix(padded, nrow = scale)))
  }

  if (.is_wcx_sample(sample)) .as_wcx_sample(result) else result
}


#' Correct gonosomal read counts for male samples
#'
#' Doubles chrX (key `"23"`) and chrY (key `"24"`) counts for male samples
#' to level with autosomal diploid coverage. This is a no-op for female
#' samples. Mirrors `overall_tools.gender_correct()` in upstream WisecondorX.
#'
#' @param sample Named list of integer/numeric vectors keyed by chromosome.
#' @param gender Character; `"M"` for male, `"F"` for female.
#'
#' @return Modified `sample` with doubled chrX/Y if male.
#'
#' @keywords internal
.gender_correct <- function(sample, gender) {
  if (gender == "M") {
    if (!is.null(sample[["23"]])) sample[["23"]] <- sample[["23"]] * 2L
    if (!is.null(sample[["24"]])) sample[["24"]] <- sample[["24"]] * 2L
  }
  sample
}


#' Predict gender from Y-chromosome read fraction
#'
#' Classifies a sample as male or female based on the fraction of reads
#' mapping to chrY, using a cutoff from the reference GMM model.
#' Mirrors `predict_tools.predict_gender()` in upstream WisecondorX.
#'
#' @param sample Named list of integer/numeric vectors keyed by chromosome.
#' @param trained_cutoff Numeric; Y-fraction threshold from reference.
#'
#' @return `"M"` or `"F"`.
#'
#' @keywords internal
.predict_gender <- function(sample, trained_cutoff) {
  total <- sum(vapply(sample, function(x) {
    if (is.null(x)) 0 else sum(as.numeric(x))
  }, numeric(1)))
  y_sum <- if (is.null(sample[["24"]])) 0 else sum(as.numeric(sample[["24"]]))
  y_frac <- y_sum / total
  if (y_frac > trained_cutoff) "M" else "F"
}


#' Compute the global bin mask from a set of reference samples
#'
#' Identifies bins with sufficient coverage across all reference samples.
#' Bins where the summed normalized coverage is less than 5% of the median
#' are masked out. Mirrors `newref_tools.get_mask()` in upstream WisecondorX.
#'
#' @param samples List of sample objects (each a named list of integer vectors
#'   keyed by chromosome `"1"`--`"24"`).
#' @param ref_bins_per_chr Optional integer vector of length 24 giving the
#'   minimum number of bins per chromosome. When supplied (e.g. from the global
#'   mask), each chromosome is zero-padded to at least this many bins so that
#'   the returned mask has the same length as the global mask. This avoids
#'   length mismatches when combining gender-specific masks with `&`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{mask}{Logical vector of length `sum(bins_per_chr)`; `TRUE` for
#'       bins with sufficient coverage.}
#'     \item{bins_per_chr}{Integer vector of length 24 giving the number of
#'       bins per chromosome.}
#'   }
#'
#' @keywords internal
.get_mask <- function(samples, ref_bins_per_chr = NULL) {
  n_samples <- length(samples)
  bins_per_chr <- integer(24L)
  by_chr <- vector("list", 24L)

  for (chr in 1:24) {
    chr_key <- as.character(chr)
    max_len <- max(vapply(samples, function(s) {
      x <- s[[chr_key]]
      if (is.null(x)) 0L else length(x)
    }, integer(1)))
    # Ensure at least ref_bins_per_chr bins so the mask length matches the
    # global mask when combining gender-specific masks with &.
    if (!is.null(ref_bins_per_chr)) {
      max_len <- max(max_len, ref_bins_per_chr[chr])
    }
    bins_per_chr[chr] <- max_len

    mat <- matrix(0, nrow = max_len, ncol = n_samples)
    for (j in seq_len(n_samples)) {
      x <- samples[[j]][[chr_key]]
      if (!is.null(x)) mat[seq_along(x), j] <- as.numeric(x)
    }
    by_chr[[chr]] <- mat
  }

  all_data <- do.call(rbind, by_chr)

  # Normalize each sample column by total reads
  col_sums <- colSums(all_data)
  col_sums[col_sums == 0] <- 1  # avoid division by zero
  all_data <- sweep(all_data, 2, col_sums, "/")

  sum_per_bin <- rowSums(all_data)
  median_cov <- stats::median(sum_per_bin[sum_per_bin > 0])
  mask <- sum_per_bin > (0.05 * median_cov)

  list(mask = mask, bins_per_chr = bins_per_chr)
}


#' Train a gender model using GMM on Y-fractions
#'
#' Fits a 2-component Gaussian mixture model to the Y-chromosome read fractions
#' of a set of reference samples. Uses `mclust::Mclust()` as the R-native
#' replacement for `sklearn.GaussianMixture`. The local minimum of the fitted
#' density is used as the male/female cutoff.
#'
#' Mirrors `newref_tools.train_gender_model()` in upstream WisecondorX.
#'
#' @param samples List of sample objects.
#' @param yfrac Optional numeric; if given, used as the manual cutoff instead
#'   of the GMM-derived one.
#'
#' @return A list with:
#'   \describe{
#'     \item{genders}{Character vector of `"M"` / `"F"` per sample.}
#'     \item{cutoff}{Numeric; the Y-fraction cutoff used.}
#'   }
#'
#' @keywords internal
.train_gender_model <- function(samples, yfrac = NULL) {
  y_fractions <- vapply(samples, function(s) {
    total <- sum(vapply(s, function(x) {
      if (is.null(x)) 0 else sum(as.numeric(x))
    }, numeric(1)))
    y_sum <- if (is.null(s[["24"]])) 0 else sum(as.numeric(s[["24"]]))
    y_sum / total
  }, numeric(1))

  if (!is.null(yfrac)) {
    cutoff <- yfrac
  } else {
    fit <- .mclust_gender_fit(y_fractions)

    if (is.null(fit)) {
      # Mclust returned NULL (rare edge case). Use the gap between the
      # smallest nonzero fraction and zero as the cutoff.
      nonzero <- y_fractions[y_fractions > 0]
      cutoff <- if (length(nonzero) > 0L) min(nonzero) / 2 else 0.001
      warning("GMM gender model did not converge; using fallback cutoff ", cutoff,
              ". Verify gender assignments before trusting gonosomal CNV calls.",
              call. = FALSE)
    } else {
      # Find local minimum of the density on a fine grid [0, 0.02]
      gmm_x <- seq(0, 0.02, length.out = 5000)
      gmm_y <- .gmm_density(fit, gmm_x)

      local_mins <- which(diff(sign(diff(gmm_y))) == 2) + 1L
      if (length(local_mins) == 0L) {
        cutoff <- mean(fit$parameters$mean)
      } else {
        cutoff <- gmm_x[local_mins[1]]
      }
    }
  }

  genders <- ifelse(y_fractions > cutoff, "M", "F")
  list(genders = genders, cutoff = cutoff)
}


#' Fit a 2-component GMM using mclust (gender model)
#'
#' @param y_fractions Numeric vector of Y-chromosome fractions.
#' @return mclust model object.
#' @keywords internal
.mclust_gender_fit <- function(y_fractions) {
  # Use the same namespace trick as .mclust_fit() in nipter_sex.R
  env <- list2env(
    list(data = matrix(y_fractions, ncol = 1)),
    parent = asNamespace("mclust")
  )
  evalq(Mclust(data, G = 2, modelNames = "V"), envir = env)
}


#' Evaluate the density of a mclust GMM
#'
#' @param fit mclust model object.
#' @param x Numeric vector of evaluation points.
#' @return Numeric vector of density values.
#' @keywords internal
.gmm_density <- function(fit, x) {
  params <- fit$parameters
  n_comp <- length(params$pro)
  dens <- numeric(length(x))
  for (k in seq_len(n_comp)) {
    dens <- dens + params$pro[k] * stats::dnorm(
      x, mean = params$mean[k],
      sd = sqrt(as.numeric(params$variance$sigmasq[k]))
    )
  }
  dens
}


#' Stack samples into a matrix for a subset of chromosomes and apply a mask
#'
#' @param samples List of sample objects.
#' @param chr_range Integer vector of chromosome indices (1-based).
#' @param mask Logical vector covering all bins in `chr_range`.
#' @return Numeric matrix of shape `(n_masked_bins, n_samples)`.
#' @keywords internal
.normalize_and_mask <- function(samples, chr_range, mask) {
  n_samples <- length(samples)
  by_chr <- vector("list", length(chr_range))

  for (ii in seq_along(chr_range)) {
    chr <- chr_range[ii]
    chr_key <- as.character(chr)
    max_len <- max(vapply(samples, function(s) {
      x <- s[[chr_key]]
      if (is.null(x)) 0L else length(x)
    }, integer(1)))

    mat <- matrix(0, nrow = max_len, ncol = n_samples)
    for (j in seq_len(n_samples)) {
      x <- samples[[j]][[chr_key]]
      if (!is.null(x)) mat[seq_along(x), j] <- as.numeric(x)
    }
    by_chr[[ii]] <- mat
  }

  all_data <- do.call(rbind, by_chr)

  # Normalize each sample by total reads
  col_sums <- colSums(all_data)
  col_sums[col_sums == 0] <- 1
  all_data <- sweep(all_data, 2, col_sums, "/")

  # Apply mask
  all_data[mask, , drop = FALSE]
}


#' Train PCA and return ratio-corrected data
#'
#' Fits a PCA model with `n_comp` components to the masked reference data,
#' then corrects the data by dividing by the PCA reconstruction (ratio
#' correction, not subtractive). Mirrors `newref_tools.train_pca()`.
#'
#' @param ref_data Numeric matrix `(n_masked_bins, n_samples)`.
#' @param n_comp Number of PCA components (default 5).
#'
#' @return A list with:
#'   \describe{
#'     \item{corrected}{Numeric matrix `(n_masked_bins, n_samples)` — ratio-
#'       corrected data.}
#'     \item{components}{Numeric matrix `(n_comp, n_masked_bins)` — PCA
#'       rotation (loadings).}
#'     \item{center}{Numeric vector of length `n_masked_bins` — PCA column
#'       means.}
#'   }
#'
#' @keywords internal
.train_pca <- function(ref_data, n_comp = 5L) {
  # PCA is fitted on transposed data: (n_samples, n_bins)
  t_data <- t(ref_data)

  pca <- stats::prcomp(t_data, center = TRUE, scale. = FALSE, rank. = n_comp)
  # Reconstruct: transform then inverse-transform
  transformed <- pca$x[, seq_len(n_comp), drop = FALSE]
  reconstructed <- transformed %*% t(pca$rotation[, seq_len(n_comp), drop = FALSE])
  reconstructed <- sweep(reconstructed, 2, pca$center, "+")

  # Ratio correction: data / reconstruction
  corrected <- t_data / reconstructed
  # Handle NaN/Inf from division
  corrected[!is.finite(corrected)] <- 0

  list(
    corrected  = t(corrected),  # back to (n_bins, n_samples)
    components = t(pca$rotation[, seq_len(n_comp), drop = FALSE]),  # (n_comp, n_bins)
    center     = pca$center  # length n_bins
  )
}


#' Project a single sample through stored PCA and apply ratio correction
#'
#' Mirrors `predict_tools.project_pc()` in upstream WisecondorX.
#'
#' @param sample_data Numeric vector of length `n_masked_bins` (coverage-
#'   normalized and masked).
#' @param pca_components Matrix `(n_comp, n_masked_bins)`.
#' @param pca_mean Numeric vector of length `n_masked_bins`.
#'
#' @return Numeric vector of length `n_masked_bins` — ratio-corrected.
#'
#' @keywords internal
.project_pc <- function(sample_data, pca_components, pca_mean) {
  # Transform: sample into PCA space
  centered <- sample_data - pca_mean
  transformed <- centered %*% t(pca_components)  # (1, n_comp)

  # Reconstruct
  reconstructed <- drop(transformed %*% pca_components) + pca_mean

  result <- sample_data / reconstructed
  result[!is.finite(result)] <- 0
  result
}
