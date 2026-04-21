# WisecondorX native R implementation — prediction (predict stage)
#
# Ports predict_control.py, predict_tools.py, and the orchestration from
# main.py tool_test() from the upstream WisecondorX Python package.
# Original authors: Lennart Raman, Roy Straver, Wim Audenaert.
# Ported to R with credit under GPL-3.
#
# This file implements the complete `predict` pipeline:
# 1. Load sample and reference.
# 2. Predict gender from Y-fraction.
# 3. Normalize autosomes (coverage + PCA + within-sample).
# 4. Normalize gonosomes (gender-specific sub-reference).
# 5. Combine autosomal and gonosomal results.
# 6. Post-process (inflate results, log-transform, blacklist).
# 7. Circular binary segmentation (CBS) via DNAcopy.
# 8. Aberration calling.
# 9. Output generation (BED files, statistics).


#' Predict copy number aberrations using WisecondorX
#'
#' Native R implementation of the WisecondorX `predict` pipeline. Takes a
#' binned sample and a reference (from [rwisecondorx_newref()]) and detects
#' copy number aberrations using within-sample normalization, PCA correction,
#' and circular binary segmentation (CBS).
#'
#' CBS is performed using the DNAcopy Bioconductor package (or optionally
#' ParDNAcopy for parallel segmentation). Both must be available in
#' `Suggests`. Unlike [wisecondorx_predict()], this function's `parallel`
#' argument controls native R CBS execution through `ParDNAcopy::parSegment()`;
#' the upstream Python CLI wrapper delegates segmentation to upstream
#' WisecondorX and does not expose this native switch.
#'
#' @param sample Named list of integer vectors keyed by chromosome
#'   (`"1"`--`"24"`), as returned by [bam_convert()].
#' @param reference A `WisecondorXReference` object from [rwisecondorx_newref()],
#'   or a list with equivalent structure.
#' @param sample_binsize Integer; bin size of the input sample. Used for
#'   rescaling to the reference bin size.
#' @param outprefix Character; path prefix for output files. BED and statistics
#'   files will be written as `<outprefix>_bins.bed`, `<outprefix>_segments.bed`,
#'   `<outprefix>_aberrations.bed`, `<outprefix>_statistics.txt`.
#'   If `NULL` (default), no files are written and results are returned as a list.
#' @param minrefbins Integer; minimum reference bins per target. Bins with
#'   fewer references are zeroed out. Default `150L`.
#' @param maskrepeats Integer; number of iterative distance-masking cycles.
#'   Default `5L`.
#' @param alpha Numeric; CBS breakpoint p-value threshold. Default `1e-4`.
#' @param zscore Numeric; Z-score cutoff for aberration calling. Default `5`.
#' @param optimal_cutoff_sd_multiplier Numeric; number of population standard
#'   deviations added to the mean distance when iteratively deriving the
#'   optimal within-reference cutoff. Default `3`.
#' @param within_sample_mask_iterations Integer; number of iterative
#'   within-sample masking passes. Default `3L`.
#' @param within_sample_mask_quantile Numeric scalar in `(0, 1)`; bins with
#'   \code{|z| >= qnorm(within_sample_mask_quantile)} are excluded from the
#'   next within-sample normalization pass. Default `0.99`.
#' @param cbs_split_min_gap_bp Integer; NA stretches spanning more than this
#'   many base pairs trigger CBS segment splitting. Default `2000000L`.
#' @param segment_zscore_cap Numeric; absolute cap applied to segment z-scores
#'   after they are derived from the null-ratio distributions. Default `1000`.
#' @param beta Optional numeric; if given, ratio-based cutoff is used instead
#'   of Z-score. Should approximate purity (0, 1]. Default `NULL`.
#' @param blacklist Optional character; path to a headerless BED file of
#'   regions to mask.
#' @param gender Optional character; force gender (`"F"` or `"M"`).
#' @param seed Optional integer; RNG seed for CBS reproducibility.
#' @param parallel Logical; use `ParDNAcopy::parSegment()` for CBS when
#'   `TRUE` (default). Requires the `ParDNAcopy` package. Set
#'   `parallel = FALSE` to use `DNAcopy::segment()` explicitly.
#' @param cpus Integer; number of threads for parallel CBS (`parSegment`) and
#'   any other OpenMP-accelerated steps. Default `4L`.
#'
#' @return A list with class `"WisecondorXPrediction"` containing:
#'   \describe{
#'     \item{results_r}{List of per-chromosome log2-ratio vectors.}
#'     \item{results_z}{List of per-chromosome Z-score vectors.}
#'     \item{results_w}{List of per-chromosome weight vectors.}
#'     \item{results_c}{Data frame of CBS segments (`chr`, `start`, `end`,
#'       `zscore`, `ratio`).}
#'     \item{aberrations}{Data frame of called aberrations.}
#'     \item{statistics}{Data frame of per-chromosome statistics.}
#'     \item{gender}{Predicted (or forced) gender.}
#'     \item{algorithm_params}{Named list of native prediction parameters used
#'       for this result.}
#'     \item{n_reads}{Total read count.}
#'     \item{binsize}{Reference bin size.}
#'   }
#'
#' @seealso [rwisecondorx_newref()], [scale_sample()]
#'
#' @export
rwisecondorx_predict <- function(sample,
                                reference,
                                sample_binsize  = NULL,
                                outprefix       = NULL,
                                minrefbins      = 150L,
                                maskrepeats     = 5L,
                                alpha           = 1e-4,
                                zscore          = 5,
                                optimal_cutoff_sd_multiplier = 3,
                                within_sample_mask_iterations = 3L,
                                within_sample_mask_quantile = 0.99,
                                cbs_split_min_gap_bp = 2000000L,
                                segment_zscore_cap = 1000,
                                beta            = NULL,
                                blacklist       = NULL,
                                gender          = NULL,
                                seed            = NULL,
                                parallel        = TRUE,
                                cpus            = 4L) {
  # ---------- validation ----------
  sample <- .as_wcx_sample(sample)
  reference <- .as_wcx_reference(reference)
  stopifnot(is.list(sample), is.list(reference))
  minrefbins  <- as.integer(minrefbins)
  maskrepeats <- as.integer(maskrepeats)
  within_sample_mask_iterations <- as.integer(within_sample_mask_iterations)
  cbs_split_min_gap_bp <- as.integer(cbs_split_min_gap_bp)
  stopifnot(alpha > 0, alpha <= 1)
  stopifnot(zscore > 0)
  stopifnot(
    is.numeric(optimal_cutoff_sd_multiplier),
    length(optimal_cutoff_sd_multiplier) == 1L,
    is.finite(optimal_cutoff_sd_multiplier),
    optimal_cutoff_sd_multiplier > 0
  )
  stopifnot(within_sample_mask_iterations >= 1L)
  stopifnot(within_sample_mask_quantile > 0, within_sample_mask_quantile < 1)
  stopifnot(cbs_split_min_gap_bp >= 1L)
  stopifnot(
    is.numeric(segment_zscore_cap),
    length(segment_zscore_cap) == 1L,
    is.finite(segment_zscore_cap),
    segment_zscore_cap > 0
  )
  if (!is.null(beta)) stopifnot(beta > 0, beta <= 1)
  if (!is.null(gender)) stopifnot(gender %in% c("F", "M"))
  if (!is.null(blacklist)) stopifnot(file.exists(blacklist))
  .validate_rwisecondorx_predict_threshold(reference, minrefbins)

  binsize <- as.integer(reference$binsize)

  # ---------- Step 1: rescale sample ----------
  n_reads <- sum(vapply(sample, function(x) {
    if (is.null(x)) 0 else sum(as.numeric(x))
  }, numeric(1)))

  if (!is.null(sample_binsize)) {
    sample <- scale_sample(sample, from_size = sample_binsize, to_size = binsize)
  }

  # ---------- Step 2: predict gender ----------
  pred_gender <- .predict_gender(sample, reference$trained_cutoff)
  is_nipt <- isTRUE(reference$is_nipt)

  if (!is_nipt) {
    if (!is.null(gender)) pred_gender <- gender
    sample <- .gender_correct(sample, pred_gender)
    ref_gender <- pred_gender
  } else {
    if (!is.null(gender)) pred_gender <- gender
    ref_gender <- "F"
  }

  # Handle missing gonosomal references
  if (!is_nipt) {
    if (!isTRUE(reference$has_male) && ref_gender == "M") {
      stop(
        "Sample was classified as male but the reference lacks a male gonosomal partition. ",
        "Refuse to substitute the female partition.",
        call. = FALSE
      )
    } else if (!isTRUE(reference$has_female) && ref_gender == "F") {
      stop(
        "Sample was classified as female but the reference lacks a female gonosomal partition. ",
        "Refuse to substitute the male partition.",
        call. = FALSE
      )
    }
  }

  # ---------- Step 3: normalize autosomes ----------
  aut <- .normalize(
    sample,
    reference,
    "A",
    maskrepeats,
    optimal_cutoff_sd_multiplier = optimal_cutoff_sd_multiplier,
    within_sample_mask_iterations = within_sample_mask_iterations,
    within_sample_mask_quantile = within_sample_mask_quantile
  )

  # ---------- Step 4: normalize gonosomes ----------
  gon <- .normalize(
    sample,
    reference,
    ref_gender,
    maskrepeats,
    optimal_cutoff_sd_multiplier = optimal_cutoff_sd_multiplier,
    within_sample_mask_iterations = within_sample_mask_iterations,
    within_sample_mask_quantile = within_sample_mask_quantile
  )

  # .normalize() with ref_gender="F"/"M" already returns ONLY the gonosomal
  # bins (it sets ct = masked_bins_per_chr_cum.{F/M}[22] internally and starts
  # iteration from there). No further slicing is needed — use results directly.
  # This matches upstream main.py lines 221-223: results_r_2 from normalize()
  # is appended directly with np.append(results_r, results_r_2).
  gon_results_r   <- gon$results_r
  gon_results_z   <- gon$results_z
  gon_results_w   <- gon$results_w
  gon_ref_sizes   <- gon$ref_sizes

  # ---------- Step 4b: remap autosomal results to gonosomal mask space ----------
  aut_mapped <- .remap_autosomal_to_partition(aut, reference, ref_gender)
  aut_r_remapped <- aut_mapped$r
  aut_z_remapped <- aut_mapped$z
  aut_w_remapped <- aut_mapped$w
  aut_rs_remapped <- aut_mapped$ref_sizes

  # Null ratios: autosomal from "A", gonosomal from gender-specific (sliced).
  # Upstream slices gonosomal null_ratios at len(null_ratios_aut), which
  # implicitly assumes A_aut == F_aut. When they differ, use the gonosomal
  # partition's own autosomal cum count for a correct split.
  null_ratios_gon_key <- .ref_key("null_ratios.", ref_gender)
  null_ratios_gon_full <- reference[[null_ratios_gon_key]]
  null_ratios_gon <- null_ratios_gon_full[-seq_len(aut_mapped$n_aut_gon), , drop = FALSE]
  null_ratios_aut <- aut_mapped$null_ratios_aut

  # Align null ratio column counts (partitions may have different sample counts)
  n_null_min <- min(ncol(null_ratios_aut), ncol(null_ratios_gon))
  null_ratios_aut <- null_ratios_aut[, seq_len(n_null_min), drop = FALSE]
  null_ratios_gon <- null_ratios_gon[, seq_len(n_null_min), drop = FALSE]

  # ---------- Step 5: combine ----------
  results_r <- c(aut_r_remapped, gon_results_r)
  results_z <- c(aut_z_remapped, gon_results_z) - aut$m_z
  results_w <- c(aut_w_remapped * mean(gon_results_w, na.rm = TRUE),
                 gon_results_w * mean(aut_w_remapped, na.rm = TRUE))
  results_w <- results_w / mean(results_w, na.rm = TRUE)

  if (any(!is.finite(results_w))) {
    stop(
      "RWisecondorX normalization produced non-finite bin weights. ",
      "Refuse to replace them with uniform weights.",
      call. = FALSE
    )
  }

  ref_sizes <- c(aut_rs_remapped, gon_ref_sizes)
  null_ratios <- rbind(null_ratios_aut, null_ratios_gon)

  # ---------- Step 6: post-process ----------
  # Build rem_input structure
  gon_suffix <- paste0(".", ref_gender)
  gon_mask_key <- .ref_key("mask", gon_suffix)
  gon_bpc_key  <- .ref_key("bins_per_chr", gon_suffix)
  gon_mbpc_key <- .ref_key("masked_bins_per_chr", gon_suffix)

  rem_mask             <- reference[[gon_mask_key]]
  rem_bins_per_chr     <- reference[[gon_bpc_key]]
  rem_masked_bins_per_chr <- reference[[gon_mbpc_key]]

  # Zero out bins with insufficient references
  insuff_mask <- ref_sizes < minrefbins
  results_r[insuff_mask] <- 0
  results_z[insuff_mask] <- 0
  results_w[insuff_mask] <- 0

  # Inflate results (un-mask)
  results_r_inflated <- .inflate_results(results_r, rem_mask)
  results_z_inflated <- .inflate_results(results_z, rem_mask)
  results_w_inflated <- .inflate_results(results_w, rem_mask)

  # Handle null_ratios inflation — each column independently
  nr_inflated <- matrix(0, nrow = length(rem_mask), ncol = ncol(null_ratios))
  for (col_i in seq_len(ncol(null_ratios))) {
    nr_inflated[, col_i] <- .inflate_results(null_ratios[, col_i], rem_mask)
  }

  # Split into per-chromosome lists
  results_r_chr <- .split_by_chr(results_r_inflated, rem_bins_per_chr)
  results_z_chr <- .split_by_chr(results_z_inflated, rem_bins_per_chr)
  results_w_chr <- .split_by_chr(results_w_inflated, rem_bins_per_chr)
  results_nr_chr <- .split_nr_by_chr(nr_inflated, rem_bins_per_chr)

  # Log-transform ratios
  m_lr <- aut$m_lr
  lt <- .log_trans(results_r_chr, results_z_chr, results_w_chr, m_lr)
  results_r_chr <- lt$results_r
  results_z_chr <- lt$results_z
  results_w_chr <- lt$results_w

  # Apply blacklist
  if (!is.null(blacklist)) {
    bl <- .apply_blacklist(blacklist, results_r_chr, results_z_chr, results_w_chr, binsize)
    results_r_chr <- bl$results_r
    results_z_chr <- bl$results_z
    results_w_chr <- bl$results_w
  }

  # ---------- Step 7: CBS ----------
  cbs_result <- .exec_cbs(
    results_r_chr,
    results_w_chr,
    ref_gender,
    alpha,
    binsize,
    seed,
    parallel,
    cpus = cpus,
    split_min_gap_bp = cbs_split_min_gap_bp
  )

  # Compute segment Z-scores
  segment_zscores <- .get_segment_zscores(
    cbs_result,
    results_nr_chr,
    results_r_chr,
    results_w_chr,
    zscore_cap = segment_zscore_cap
  )

  # Build results_c
  results_c <- data.frame(
    chr     = cbs_result$chr,
    start   = cbs_result$start,
    end     = cbs_result$end,
    zscore  = segment_zscores,
    ratio   = cbs_result$ratio
  )

  # ---------- Step 8: aberrations ----------
  aberrations <- .call_aberrations(results_c, zscore, beta, ref_gender)

  # ---------- Step 9: statistics ----------
  statistics <- .compute_statistics(results_r_chr, results_w_chr, results_c,
                                   results_nr_chr, rem_bins_per_chr, binsize,
                                   ref_gender, pred_gender, n_reads)

  # ---------- output ----------
  prediction <- list(
    results_r    = results_r_chr,
    results_z    = results_z_chr,
    results_w    = results_w_chr,
    results_c    = results_c,
    aberrations  = aberrations,
    statistics   = statistics,
    gender       = pred_gender,
    ref_gender   = ref_gender,
    n_reads      = n_reads,
    binsize      = binsize,
    bins_per_chr = rem_bins_per_chr,
    algorithm_params = list(
      minrefbins = minrefbins,
      maskrepeats = maskrepeats,
      alpha = alpha,
      zscore = zscore,
      optimal_cutoff_sd_multiplier = optimal_cutoff_sd_multiplier,
      within_sample_mask_iterations = within_sample_mask_iterations,
      within_sample_mask_quantile = within_sample_mask_quantile,
      cbs_split_min_gap_bp = cbs_split_min_gap_bp,
      segment_zscore_cap = segment_zscore_cap,
      beta = beta,
      gender = if (is.null(gender)) NULL else as.character(gender),
      blacklist = if (is.null(blacklist)) NULL else normalizePath(blacklist, mustWork = TRUE),
      seed = if (is.null(seed)) NULL else as.integer(seed),
      parallel = isTRUE(parallel),
      cpus = as.integer(cpus)
    )
  )
  prediction <- .as_wcx_prediction(prediction)

  if (!is.null(outprefix)) {
    .write_prediction_output(prediction, outprefix)
  }

  prediction
}

.WCX_AUTOSOMAL_CHR_COUNT <- 22L

.remap_autosomal_to_partition <- function(aut, reference, ref_gender) {
  gon_suffix <- paste0(".", ref_gender)
  gon_mbpc_cum <- reference[[.ref_key("masked_bins_per_chr_cum", gon_suffix)]]
  n_aut_gon <- gon_mbpc_cum[.WCX_AUTOSOMAL_CHR_COUNT]
  n_aut_a <- reference$masked_bins_per_chr_cum[.WCX_AUTOSOMAL_CHR_COUNT]
  null_ratios_aut_full <- reference$null_ratios

  if (n_aut_a == n_aut_gon) {
    return(list(
      r = aut$results_r,
      z = aut$results_z,
      w = aut$results_w,
      ref_sizes = aut$ref_sizes,
      null_ratios_aut = null_ratios_aut_full,
      n_aut_gon = n_aut_gon
    ))
  }

  n_aut_unmask <- sum(reference$bins_per_chr[seq_len(.WCX_AUTOSOMAL_CHR_COUNT)])
  mask_a <- reference$mask[seq_len(n_aut_unmask)]
  mask_gon <- reference[[.ref_key("mask", gon_suffix)]][seq_len(n_aut_unmask)]

  a_idx <- which(mask_a)
  g_idx <- which(mask_gon)
  a_rank <- integer(n_aut_unmask)
  a_rank[a_idx] <- seq_along(a_idx)
  mapped_ranks <- a_rank[g_idx]
  keep <- mapped_ranks > 0L

  remap_aut <- function(aut_vec) {
    out <- numeric(length(g_idx))
    out[keep] <- aut_vec[mapped_ranks[keep]]
    out
  }

  null_ratios_aut <- matrix(
    0,
    nrow = length(g_idx),
    ncol = ncol(null_ratios_aut_full)
  )
  if (any(keep)) {
    null_ratios_aut[keep, ] <- null_ratios_aut_full[mapped_ranks[keep], , drop = FALSE]
  }

  list(
    r = remap_aut(aut$results_r),
    z = remap_aut(aut$results_z),
    w = remap_aut(aut$results_w),
    ref_sizes = remap_aut(aut$ref_sizes),
    null_ratios_aut = null_ratios_aut,
    n_aut_gon = n_aut_gon
  )
}

.validate_rwisecondorx_predict_threshold <- function(reference, minrefbins) {
  stopifnot(is.list(reference))
  stopifnot(is.numeric(minrefbins), length(minrefbins) == 1L, minrefbins >= 1L)

  index_keys <- grep("^indexes(\\.|$)", names(reference), value = TRUE)
  if (!length(index_keys)) {
    return(invisible(NULL))
  }

  refsize <- max(vapply(index_keys, function(key) {
    idx <- reference[[key]]
    dims <- dim(idx)
    if (is.null(dims) || length(dims) < 2L) {
      return(as.integer(length(idx)))
    }
    as.integer(dims[[2L]])
  }, integer(1L)))

  if (!is.finite(refsize) || minrefbins <= refsize) {
    return(invisible(NULL))
  }

  stop(
    sprintf(
      "Invalid RWisecondorX predict settings: minrefbins=%d exceeds the reference refsize=%d. ",
      as.integer(minrefbins),
      as.integer(refsize)
    ),
    "No target bin can retain enough reference bins under this combination. ",
    "Lower 'minrefbins' or rebuild the reference with a larger 'refsize'.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# Internal normalization helpers
# ---------------------------------------------------------------------------

#' Run the full normalization chain for one gender partition
#'
#' @param sample Named list of bin counts.
#' @param ref Reference list.
#' @param ref_gender `"A"`, `"F"`, or `"M"`.
#' @param maskrepeats Number of optimal cutoff iterations.
#' @return List with results_r, results_z, results_w, ref_sizes, m_lr, m_z.
#' @keywords internal
.normalize <- function(sample,
                       ref,
                       ref_gender,
                       maskrepeats,
                       optimal_cutoff_sd_multiplier = 3,
                       within_sample_mask_iterations = 3L,
                       within_sample_mask_quantile = 0.99) {
  ap <- if (ref_gender == "A") "" else paste0(".", ref_gender)
  cp <- if (ref_gender == "A") 0L else .WCX_AUTOSOMAL_CHR_COUNT
  ct <- if (ref_gender == "A") 0L else ref[[.ref_key("masked_bins_per_chr_cum", ap)]][cp]

  # Coverage normalization + masking
  sample_masked <- .coverage_normalize_and_mask(sample, ref, ap)

  # PCA projection
  pca_comp <- ref[[paste0("pca_components", ap)]]
  pca_mean <- ref[[paste0("pca_mean", ap)]]
  sample_corrected <- .project_pc(sample_masked, pca_comp, pca_mean)

  # Weights
  results_w <- .get_weights(ref, ap)
  results_w <- results_w[(ct + 1L):length(results_w)]

  # Optimal distance cutoff
  optimal_cutoff <- .get_optimal_cutoff(
    ref,
    maskrepeats,
    sd_multiplier = optimal_cutoff_sd_multiplier
  )

  # Within-sample normalization (3 iterations)
  norm_result <- .normalize_repeat(
    sample_corrected,
    ref,
    optimal_cutoff,
    ct,
    cp,
    ap,
    iterations = within_sample_mask_iterations,
    mask_quantile = within_sample_mask_quantile
  )

  list(
    results_r  = norm_result$results_r,
    results_z  = norm_result$results_z,
    results_w  = results_w,
    ref_sizes  = norm_result$ref_sizes,
    m_lr       = norm_result$m_lr,
    m_z        = norm_result$m_z
  )
}


#' Coverage normalization and masking
#' @keywords internal
.coverage_normalize_and_mask <- function(sample, ref, ap) {
  bins_per_chr <- ref[[.ref_key("bins_per_chr", ap)]]
  mask         <- ref[[.ref_key("mask", ap)]]
  n_chrs <- length(bins_per_chr)

  by_chr <- vector("list", n_chrs)
  for (chr in seq_len(n_chrs)) {
    chr_key <- as.character(chr)
    target_len <- bins_per_chr[chr]
    this_chr <- numeric(target_len)
    x <- sample[[chr_key]]
    if (!is.null(x)) {
      min_len <- min(target_len, length(x))
      this_chr[seq_len(min_len)] <- as.numeric(x[seq_len(min_len)])
    }
    by_chr[[chr]] <- this_chr
  }

  all_data <- unlist(by_chr, use.names = FALSE)
  total <- sum(all_data)
  if (total > 0) all_data <- all_data / total

  all_data[mask]
}


#' Get weights from reference distances
#' @keywords internal
.get_weights <- function(ref, ap) {
  distances <- ref[[.ref_key("distances", ap)]]
  # weights = 1 / mean(sqrt(distances_per_bin))
  inv_w <- apply(distances, 1, function(x) mean(sqrt(x)))
  1 / inv_w
}


#' Iterative optimal distance cutoff
#' @keywords internal
.get_optimal_cutoff <- function(ref, repeats, sd_multiplier = 3) {
  distances <- ref$distances
  cutoff <- Inf
  for (i in seq_len(repeats)) {
    vals <- distances[distances < cutoff]
    if (length(vals) == 0L) break
    avg <- mean(vals)
    sdev <- .wcx_population_sd(vals)
    cutoff <- avg + as.numeric(sd_multiplier) * sdev
  }
  cutoff
}


#' Within-sample normalization with iterative aberration masking (3 iterations)
#' @keywords internal
.normalize_repeat <- function(test_data,
                              ref,
                              optimal_cutoff,
                              ct,
                              cp,
                              ap,
                              iterations = 3L,
                              mask_quantile = 0.99) {
  test_copy <- test_data  # mutable copy
  results_z <- NULL
  results_r <- NULL
  ref_sizes <- NULL

  z_mask_cutoff <- stats::qnorm(mask_quantile)
  for (iter in seq_len(as.integer(iterations))) {
    norm_once <- .normalize_once(test_data, test_copy, ref, optimal_cutoff, ct, cp, ap)
    results_z <- norm_once$results_z
    results_r <- norm_once$results_r
    ref_sizes <- norm_once$ref_sizes

    # Mask aberrant bins for the next within-sample normalization pass.
    # NaN/Inf z-scores from zero-sd bins: Inf gets masked,
    # NaN stays unmasked (matches numpy behavior: np.abs(nan) >= x is False).
    aberrant <- abs(results_z) >= z_mask_cutoff
    aberrant[is.na(aberrant)] <- FALSE
    # ct offsets into the full array; aberrant is already ct-relative
    idx_start <- ct + 1L
    idx_end   <- ct + length(results_z)
    test_copy[idx_start:idx_end][aberrant] <- -1
  }

  m_lr <- stats::median(log2(results_r), na.rm = TRUE)
  m_z  <- stats::median(results_z, na.rm = TRUE)

  list(results_z = results_z, results_r = results_r,
       ref_sizes = ref_sizes, m_lr = m_lr, m_z = m_z)
}


#' Single within-sample normalization pass
#' @keywords internal
.normalize_once <- function(test_data, test_copy, ref, optimal_cutoff, ct, cp, ap) {
  masked_bins_per_chr     <- ref[[.ref_key("masked_bins_per_chr", ap)]]
  masked_bins_per_chr_cum <- ref[[.ref_key("masked_bins_per_chr_cum", ap)]]
  indexes   <- ref[[.ref_key("indexes", ap)]]
  distances <- ref[[.ref_key("distances", ap)]]

  total_bins <- masked_bins_per_chr_cum[length(masked_bins_per_chr_cum)]
  n_out <- total_bins - ct

  results_z <- numeric(n_out)
  results_r <- numeric(n_out)
  ref_sizes <- numeric(n_out)

  i  <- ct + 1L   # 1-based position in the full masked array
  i2 <- 1L        # 1-based position in the output

  for (chr_idx in (cp + 1L):length(masked_bins_per_chr)) {
    n_chr <- masked_bins_per_chr[chr_idx]
    if (n_chr == 0L) next

    chr_cum <- masked_bins_per_chr_cum[chr_idx]
    chr_start <- chr_cum - n_chr + 1L

    # Exclude same-chromosome data from test_copy
    before <- if (chr_start > 1L) test_copy[seq_len(chr_start - 1L)] else numeric(0)
    after  <- if (chr_cum < length(test_copy)) test_copy[(chr_cum + 1L):length(test_copy)] else numeric(0)
    chr_data <- c(before, after)

    for (local_bin in seq_len(n_chr)) {
      bin_idx <- chr_start + local_bin - 1L

      # Get reference indices, filtering by distance
      idx_vec  <- indexes[bin_idx, ]
      dist_vec <- distances[bin_idx, ]
      close_enough <- dist_vec < optimal_cutoff
      ref_idx <- idx_vec[close_enough]

      # Convert global indexes to chr_data-local indexes.
      # chr_data = c(test_copy[1..(chr_start-1)], test_copy[(chr_cum+1)..end])
      # Global index g maps to: g if g < chr_start, g - n_chr if g > chr_cum.
      # Indexes within [chr_start, chr_cum] should never appear (same-chr bins
      # are excluded during KNN reference building).
      bad <- ref_idx >= chr_start & ref_idx <= chr_cum
      if (any(bad)) {
        stop(sprintf(
          "KNN reference contains same-chromosome indexes for chr %d (bin %d): %s",
          chr_idx, bin_idx, paste(ref_idx[bad], collapse = ", ")
        ))
      }
      ref_idx_local <- ref_idx - n_chr * (ref_idx >= chr_start)

      # Get reference values from within the sample (from excluded-chr data)
      ref_vals <- chr_data[ref_idx_local]
      ref_vals <- ref_vals[ref_vals >= 0]  # exclude masked (-1) bins

      ref_mean <- if (length(ref_vals)) mean(ref_vals) else NaN
      ref_sd <- .wcx_population_sd(ref_vals)
      ref_median <- if (length(ref_vals)) stats::median(ref_vals) else NaN

      results_z[i2] <- (test_data[i] - ref_mean) / ref_sd
      results_r[i2] <- test_data[i] / ref_median
      ref_sizes[i2] <- length(ref_vals)

      i  <- i + 1L
      i2 <- i2 + 1L
    }
  }

  list(results_z = results_z, results_r = results_r, ref_sizes = ref_sizes)
}


# ---------------------------------------------------------------------------
# Post-processing helpers
# ---------------------------------------------------------------------------

#' Inflate (un-mask) a compressed result vector
#' @keywords internal
.inflate_results <- function(results, mask) {
  out <- numeric(length(mask))
  out[mask] <- results
  out
}


#' Split a full-genome vector into per-chromosome lists
#' @keywords internal
.split_by_chr <- function(data, bins_per_chr) {
  result <- vector("list", length(bins_per_chr))
  chr_names <- vapply(seq_along(bins_per_chr), .idx_to_chr_name, character(1L))
  names(result) <- chr_names
  pos <- 1L
  for (i in seq_along(bins_per_chr)) {
    n <- bins_per_chr[i]
    result[[i]] <- data[pos:(pos + n - 1L)]
    pos <- pos + n
  }
  result
}


#' Split a null-ratio matrix into per-chromosome list of matrices
#' @keywords internal
.split_nr_by_chr <- function(nr_matrix, bins_per_chr) {
  result <- vector("list", length(bins_per_chr))
  chr_names <- vapply(seq_along(bins_per_chr), .idx_to_chr_name, character(1L))
  names(result) <- chr_names
  pos <- 1L
  for (i in seq_along(bins_per_chr)) {
    n <- bins_per_chr[i]
    result[[i]] <- nr_matrix[pos:(pos + n - 1L), , drop = FALSE]
    pos <- pos + n
  }
  result
}


#' Log-transform ratios and clean non-finite values
#'
#' Modifies lists in place. Matches upstream `predict_tools.log_trans()`.
#' @keywords internal
.log_trans <- function(results_r, results_z, results_w, m_lr) {
  for (chr in seq_along(results_r)) {
    results_r[[chr]] <- log2(results_r[[chr]])
    # Clean non-finite values
    bad <- !is.finite(results_r[[chr]])
    results_r[[chr]][bad] <- 0
    results_z[[chr]][bad] <- 0
    results_w[[chr]][bad] <- 0
    # Center by autosomal median log-ratio
    nonzero <- results_r[[chr]] != 0
    results_r[[chr]][nonzero] <- results_r[[chr]][nonzero] - m_lr
  }
  list(results_r = results_r, results_z = results_z, results_w = results_w)
}


#' Apply a blacklist BED to zero out regions
#' @keywords internal
.apply_blacklist <- function(blacklist_path, results_r, results_z, results_w, binsize) {
  lines <- readLines(blacklist_path)
  for (line in lines) {
    parts <- strsplit(trimws(line), "\t")[[1]]
    if (length(parts) < 3L) next
    chr_name <- .normalize_chr_name(parts[1])
    chr_idx <- suppressWarnings(as.integer(chr_name))
    if (is.na(chr_idx) || chr_idx < 1L || chr_idx > length(results_r)) next

    s_bin <- as.integer(as.numeric(parts[2]) / binsize) + 1L  # 1-based R index
    e_bin <- as.integer(as.numeric(parts[3]) / binsize) + 1L  # 1-based; upstream int(e/b)+1 is exclusive, R : is inclusive

    n <- length(results_r[[chr_idx]])
    s_bin <- max(1L, s_bin)
    e_bin <- min(n, e_bin)
    if (s_bin <= e_bin) {
      results_r[[chr_idx]][s_bin:e_bin] <- 0
      results_z[[chr_idx]][s_bin:e_bin] <- 0
      results_w[[chr_idx]][s_bin:e_bin] <- 0
    }
  }
  list(results_r = results_r, results_z = results_z, results_w = results_w)
}


#' Call aberrations from CBS segments
#' @keywords internal
.call_aberrations <- function(results_c, zscore_cutoff, beta, ref_gender) {
  if (nrow(results_c) == 0L) {
    return(data.frame(chr = integer(0), start = integer(0), end = integer(0),
                      ratio = numeric(0), zscore = numeric(0),
                      type = character(0)))
  }

  aberrations <- vector("list", nrow(results_c))
  n_ab <- 0L

  for (i in seq_len(nrow(results_c))) {
    seg <- results_c[i, ]
    chr_name <- seg$chr
    # Determine ploidy for beta mode
    ploidy <- 2L
    if (chr_name %in% c(23L, 24L) && ref_gender == "M") ploidy <- 1L

    ab_type <- NULL

    if (!is.null(beta)) {
      loss_cutoff <- log2((ploidy - beta / 2) / ploidy)
      gain_cutoff <- log2((ploidy + beta / 2) / ploidy)
      if (is.finite(seg$ratio) && seg$ratio > gain_cutoff) {
        ab_type <- "gain"
      } else if (is.finite(seg$ratio) && seg$ratio < loss_cutoff) {
        ab_type <- "loss"
      }
    } else {
      if (is.finite(seg$zscore) && seg$zscore > zscore_cutoff) {
        ab_type <- "gain"
      } else if (is.finite(seg$zscore) && seg$zscore < -zscore_cutoff) {
        ab_type <- "loss"
      }
    }

    if (!is.null(ab_type)) {
      n_ab <- n_ab + 1L
      aberrations[[n_ab]] <- data.frame(
        chr    = seg$chr,
        start  = seg$start,
        end    = seg$end,
        ratio  = seg$ratio,
        zscore = seg$zscore,
        type   = ab_type
      )
    }
  }

  if (n_ab == 0L) {
    return(data.frame(chr = integer(0), start = integer(0), end = integer(0),
                      ratio = numeric(0), zscore = numeric(0),
                      type = character(0)))
  }

  do.call(rbind, aberrations[seq_len(n_ab)])
}
