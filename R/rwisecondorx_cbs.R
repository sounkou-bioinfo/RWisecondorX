# WisecondorX native R implementation — CBS segmentation
#
# Ports the CBS logic from CBS.R and predict_tools.exec_cbs() from upstream
# WisecondorX. This is a direct R-native wrapper around DNAcopy::segment()
# (or optionally ParDNAcopy::parSegment()), replacing the Python->Rscript
# subprocess call used in the upstream implementation.
#
# Original WisecondorX authors: Lennart Raman, Roy Straver, Wim Audenaert.
# DNAcopy: Adam B. Olshen, E. S. Venkatraman, et al.
# ParDNAcopy: Alex Krasnitz, Guoli Sun.
#
# CBS Coordinate Convention
# -------------------------
# .exec_cbs() returns segments with:
#   start = DNAcopy_loc.start - 1  (0-based inclusive)
#   end   = DNAcopy_loc.end        (0-based exclusive; equivalently 1-based inclusive)
#
# This forms a half-open [start, end) interval in 0-based terms, matching
# Python's array[start:end] slicing semantics.
#
# R consumers must use (start + 1):end for 1-based indexing.
# Segment length in bins is end - start.
#
# Upstream bugs intentionally replicated for conformance:
# - CPA overcounts by 1 bin per segment (overall_tools.py:146 uses
#   segment[2] - segment[1] + 1 instead of segment[2] - segment[1]).
# - Whole-chromosome Z-scores drop the last bin (predict_output.py:210
#   sets end = bins_per_chr - 1 instead of bins_per_chr).
# Both are documented in @note sections on the affected functions.


#' Execute CBS on WisecondorX results
#'
#' Runs circular binary segmentation on log2-ratios using weights, then
#' splits segments spanning large NA gaps and recalculates segmental ratios.
#' Exact port of upstream `CBS.R` logic.
#'
#' @param results_r List of per-chromosome log2-ratio vectors.
#' @param results_w List of per-chromosome weight vectors.
#' @param ref_gender `"F"` or `"M"` — determines whether chrY is included.
#' @param alpha CBS breakpoint p-value threshold.
#' @param binsize Reference bin size in bp.
#' @param seed Optional RNG seed.
#' @param parallel Logical; use ParDNAcopy when \code{TRUE} (default). This
#'   requires the \code{ParDNAcopy} package to be installed. Set
#'   \code{parallel = FALSE} to use \code{DNAcopy::segment()} explicitly.
#' @param cpus Integer; number of threads passed to \code{parSegment()}. Only
#'   used when \code{parallel = TRUE} and ParDNAcopy is available.
#'
#' @return Data frame with columns `chr` (integer, 1-based), `start` (integer,
#'   0-based bin index), `end` (integer, exclusive bin index), `ratio` (numeric,
#'   weighted mean segment ratio).
#'
#' @keywords internal
.exec_cbs <- function(results_r, results_w, ref_gender, alpha, binsize, seed,
                      parallel = TRUE, cpus = 1L,
                      split_min_gap_bp = 2000000L) {
  # Determine which chromosomes to include
  chrs <- if (ref_gender == "M") seq_len(24L) else seq_len(23L)

  # Pre-allocate then fill: avoids O(n²) copies from repeated c().
  chr_lens <- vapply(chrs, function(chr) {
    if (chr > length(results_r)) return(0L)
    as.integer(length(results_r[[chr]]))
  }, integer(1L))
  total_n <- sum(chr_lens)

  ratio_vec  <- numeric(total_n)
  weight_vec <- numeric(total_n)
  chr_vec    <- integer(total_n)
  pos_vec    <- integer(total_n)

  fill_pos <- 1L
  for (ii in seq_along(chrs)) {
    chr <- chrs[ii]
    n <- chr_lens[ii]
    if (n == 0L) next
    idx <- fill_pos:(fill_pos + n - 1L)
    ratio_vec[idx]  <- results_r[[chr]]
    weight_vec[idx] <- results_w[[chr]]
    chr_vec[idx]    <- chr
    pos_vec[idx]    <- seq_len(n)
    fill_pos <- fill_pos + n
  }

  # Upstream convention: 0 ratios become NA (blacklisted bins)
  ratio_vec[ratio_vec == 0] <- NA_real_
  # Weights cannot be 0 or NA for DNAcopy
  weight_vec[weight_vec == 0 | !is.finite(weight_vec)] <- 1e-99

  # Remove chromosomes that are entirely NA
  keep <- !is.na(ratio_vec)
  if (sum(keep) == 0L) {
    return(data.frame(chr = integer(0), start = integer(0),
                      end = integer(0), ratio = numeric(0)))
  }

  ratio_cbs  <- ratio_vec[keep]
  weight_cbs <- weight_vec[keep]
  chr_cbs    <- chr_vec[keep]
  pos_cbs    <- pos_vec[keep]

  # Set seed temporarily if provided, but restore the caller RNG state.
  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit(
      {
        if (has_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    set.seed(as.integer(seed))
  }

  # Create CNA object
  cna_obj <- DNAcopy::CNA(ratio_cbs, chr_cbs, pos_cbs,
                           data.type = "logratio", sampleid = "X")

  # Run CBS — ParDNAcopy when parallel = TRUE and available, else DNAcopy::segment().
  use_parallel <- isTRUE(parallel) && requireNamespace("ParDNAcopy", quietly = TRUE)
  if (isTRUE(parallel) && !use_parallel) {
    stop(
      "parallel = TRUE requires the 'ParDNAcopy' package. ",
      "Refuse to fall back to serial DNAcopy::segment(); either install ParDNAcopy or call with parallel = FALSE.",
      call. = FALSE
    )
  }

  if (use_parallel) {
    seg_result <- ParDNAcopy::parSegment(cna_obj,
                                         distrib = "Rparallel",
                                         njobs = as.integer(cpus),
                                         alpha = alpha,
                                         verbose = 0,
                                         weights = weight_cbs)
  } else {
    seg_result <- DNAcopy::segment(
      cna_obj,
      alpha = alpha,
      verbose = 0,
      weights = weight_cbs
    )
  }

  seg_df <- seg_result$output
  seg_df <- seg_df[, c("chrom", "loc.start", "loc.end", "seg.mean"),
                   drop = FALSE]
  colnames(seg_df) <- c("chr", "s", "e", "r")

  # Build the full position-indexed data frame for segment splitting
  full_df <- data.frame(chromosome = chr_cbs, x = pos_cbs,
                        y = ratio_cbs, w = weight_cbs)

  # Split segments spanning large NA gaps
  na_threshold <- as.integer(ceiling(as.numeric(split_min_gap_bp) / as.numeric(binsize)))
  new_segs <- .split_na_segments(seg_df, full_df, na_threshold)

  if (nrow(new_segs) == 0L) {
    return(data.frame(chr = integer(0), start = integer(0),
                      end = integer(0), ratio = numeric(0)))
  }

  # Recalculate segmental ratios using weighted means
  for (row_i in seq_len(nrow(new_segs))) {
    chr_sub <- full_df[full_df$chromosome == new_segs$chr[row_i], ]
    s <- new_segs$s[row_i]
    e <- new_segs$e[row_i]
    # Positions in chr_sub$x are 1-based sequential indices
    sel <- chr_sub$x >= s & chr_sub$x <= e
    if (any(sel)) {
      new_segs$r[row_i] <- stats::weighted.mean(chr_sub$y[sel], chr_sub$w[sel],
                                                  na.rm = TRUE)
    }
  }

  # Convert to 0-based start (upstream subtracts 1 for Python compatibility)
  data.frame(
    chr   = as.integer(new_segs$chr),
    start = as.integer(new_segs$s) - 1L,
    end   = as.integer(new_segs$e),
    ratio = as.numeric(new_segs$r)
  )
}


#' Split segments that span large NA gaps
#'
#' Exact port of the upstream CBS.R segment-splitting logic.
#'
#' @param seg_df Data frame with chr, s, e, r.
#' @param full_df Data frame with chromosome, x, y, w.
#' @param na_threshold Consecutive NA stretch that triggers a split.
#' @return Modified segment data frame.
#' @keywords internal
.split_na_segments <- function(seg_df, full_df, na_threshold) {
  new_rows <- vector("list", nrow(seg_df) * 2L)
  n_new <- 0L

  for (row_i in seq_len(nrow(seg_df))) {
    start_i <- seg_df$s[row_i]
    end_i   <- seg_df$e[row_i]
    chr_i   <- seg_df$chr[row_i]

    sub_frame <- full_df[full_df$chromosome == chr_i, ]
    if (nrow(sub_frame) == 0L) next

    # Get segment values in position order
    sel <- sub_frame$x >= start_i & sub_frame$x <= end_i
    segment_y <- sub_frame$y[sel]

    if (length(segment_y) < 2L) {
      n_new <- n_new + 1L
      new_rows[[n_new]] <- data.frame(chr = chr_i, s = start_i,
                                       e = end_i, r = seg_df$r[row_i])
      next
    }

    diff_na <- diff(is.na(segment_y))

    # Start positions of NA stretches (diff goes from 0/FALSE to 1/TRUE)
    start_pos <- which(diff_na == 1) + start_i - 1L
    # End positions of NA stretches (diff goes from TRUE to FALSE)
    end_pos   <- which(diff_na == -1) + start_i - 1L

    # Filter by threshold
    if (length(start_pos) > 0L && length(end_pos) > 0L) {
      selection <- (end_pos - start_pos) > na_threshold
      start_pos <- start_pos[selection]
      end_pos   <- end_pos[selection]
    } else {
      start_pos <- integer(0)
      end_pos   <- integer(0)
    }

    if (length(start_pos) != length(end_pos)) {
      warning(
        sprintf(
          "CBS NA-gap split mismatch on chr %s segment [%d, %d]; keeping the original segment unsplit.",
          chr_i, start_i, end_i
        ),
        call. = FALSE
      )
      n_new <- n_new + 1L
      new_rows[[n_new]] <- data.frame(chr = chr_i, s = start_i,
                                      e = end_i, r = seg_df$r[row_i])
      next
    }
    stopifnot(length(start_pos) == length(end_pos))

    # Compute inverse segments (the non-NA stretches)
    inv_start <- c(start_i, end_pos)
    inv_end   <- c(start_pos, end_i)

    # Segments should be at least two in length
    valid <- (inv_end - inv_start) > 0
    if (!any(valid)) next
    inv_start <- inv_start[valid]
    inv_end   <- inv_end[valid]

    for (j in seq_along(inv_start)) {
      n_new <- n_new + 1L
      new_rows[[n_new]] <- data.frame(chr = chr_i, s = inv_start[j],
                                       e = inv_end[j], r = seg_df$r[row_i])
    }
  }

  if (n_new == 0L) {
    return(data.frame(chr = integer(0), s = integer(0), e = integer(0),
                      r = numeric(0)))
  }

  do.call(rbind, new_rows[seq_len(n_new)])
}


#' Compute between-sample Z-scores for CBS segments
#'
#' For each segment, computes a Z-score by comparing the observed segment
#' ratio to a null distribution derived from the stored null ratios.
#' Mirrors `overall_tools.get_z_score()`.
#'
#' @param cbs_result Data frame with chr, start, end, ratio.
#' @param results_nr List of per-chromosome null-ratio matrices.
#' @param results_r List of per-chromosome log2-ratio vectors.
#' @param results_w List of per-chromosome weight vectors.
#' @param zscore_cap Numeric; absolute cap applied to finite segment z-scores.
#'
#' @return Numeric vector of Z-scores, one per segment.
#' @keywords internal
.get_segment_zscores <- function(cbs_result,
                                 results_nr,
                                 results_r,
                                 results_w,
                                 zscore_cap = 1000) {
  if (nrow(cbs_result) == 0L) return(numeric(0))

  zscores <- numeric(nrow(cbs_result))

  for (i in seq_len(nrow(cbs_result))) {
    chr <- cbs_result$chr[i]
    s   <- cbs_result$start[i] + 1L  # convert 0-based back to 1-based
    e   <- cbs_result$end[i]

    if (chr > length(results_nr)) {
      zscores[i] <- NaN
      next
    }

    # Get segment data
    seg_nr <- results_nr[[chr]][s:e, , drop = FALSE]
    seg_rr <- results_r[[chr]][s:e]
    seg_w  <- results_w[[chr]][s:e]

    # Filter out blacklisted bins (where ratio is 0)
    valid <- seg_rr != 0
    if (sum(valid) == 0L) {
      zscores[i] <- NaN
      next
    }

    seg_nr <- seg_nr[valid, , drop = FALSE]
    seg_w  <- seg_w[valid]

    # Clean non-finite values in null ratios
    seg_nr[!is.finite(seg_nr)] <- NA_real_

    # For each null sample, compute weighted average across segment bins
    n_null <- ncol(seg_nr)
    null_segments <- numeric(n_null)
    for (j in seq_len(n_null)) {
      vals <- seg_nr[, j]
      if (all(is.na(vals))) {
        null_segments[j] <- NA_real_
      } else {
        null_segments[j] <- stats::weighted.mean(vals, seg_w, na.rm = TRUE)
      }
    }

    finite_null <- null_segments[is.finite(null_segments)]
    if (length(finite_null) < 2L) {
      zscores[i] <- NaN
      next
    }

    null_mean <- mean(finite_null)
    # Upstream uses np.ma.std() which defaults to ddof=0 (population SD).
    # R's sd() uses ddof=1 (sample SD). Use population SD for conformance.
    n_fn <- length(finite_null)
    null_sd <- sqrt(sum((finite_null - null_mean)^2) / n_fn)

    if (is.na(null_mean) || is.na(null_sd) || null_sd == 0) {
      zscores[i] <- NaN
    } else {
      z <- (cbs_result$ratio[i] - null_mean) / null_sd
      zscores[i] <- max(min(z, zscore_cap), -zscore_cap)
    }
  }

  zscores
}


#' Compute per-chromosome statistics
#'
#' Mirrors `predict_output._generate_chr_statistics_file()`.
#'
#' @note **Upstream bug replicated for conformance (last-bin drop).**
#' Upstream Python (`predict_output.py:210`) sets `end = bins_per_chr - 1`
#' when constructing whole-chromosome segments. Combined with `get_z_score`'s
#' `array[s:e]` slicing (Python half-open), this drops the last bin of every
#' chromosome from the whole-chromosome Z-score calculation. We replicate this
#' exactly (`end = bins_per_chr - 1L`) for conformance with Python WisecondorX
#' output.
#'
#' @keywords internal
.compute_statistics <- function(results_r, results_w, results_c, results_nr,
                                bins_per_chr, binsize, ref_gender, gender,
                                n_reads) {
  n_chr <- length(results_r)

  # Per-chromosome ratio means and medians
  chr_ratio_means <- rep(NaN, n_chr)
  chr_ratio_medians <- rep(NaN, n_chr)
  for (chr in seq_len(n_chr)) {
    r <- results_r[[chr]]
    w <- results_w[[chr]]
    if (length(r) > 0 && any(w != 0)) {
      chr_ratio_means[chr] <- stats::weighted.mean(r, w, na.rm = TRUE)
    }
    nonzero <- r[r != 0]
    if (length(nonzero) > 0) {
      chr_ratio_medians[chr] <- stats::median(nonzero, na.rm = TRUE)
    }
  }

  # Whole-chromosome Z-scores (treat each chr as one segment)
  chr_segments <- data.frame(
    chr   = seq_len(n_chr),
    start = rep(0L, n_chr),
    # NOTE: upstream Python uses `bins_per_chr - 1` here, which combined with
    # get_z_score's `array[s:e]` slicing drops the last bin per chromosome.
    # This is a known upstream bug (predict_output.py:210) but we replicate it
    # for exact conformance with Python WisecondorX output.
    end   = bins_per_chr - 1L,
    ratio = chr_ratio_means
  )
  chr_zscores <- .get_segment_zscores(chr_segments, results_nr, results_r, results_w)

  # MSV (median segment variance)
  msv <- .get_msv(results_c, results_r)

  # CPA (copy number profile abnormality)
  cpa <- .get_cpa(results_c, binsize)

  # Chr names
  chr_names <- character(n_chr)
  for (i in seq_len(n_chr)) {
    chr_names[i] <- if (i == 23L) "X" else if (i == 24L) "Y" else as.character(i)
  }

  data.frame(
    chr           = chr_names,
    ratio_mean    = chr_ratio_means,
    ratio_median  = chr_ratio_medians,
    zscore        = chr_zscores,
    msv           = rep(msv, n_chr),
    cpa           = rep(cpa, n_chr)
  )
}


#' Median segment variance (MSV)
#' @keywords internal
.get_msv <- function(results_c, results_r) {
  if (nrow(results_c) == 0L) return(NA_real_)
  vars <- numeric(nrow(results_c))
  n_valid <- 0L
  for (i in seq_len(nrow(results_c))) {
    chr <- results_c$chr[i]
    s   <- results_c$start[i] + 1L
    e   <- results_c$end[i]
    if (chr > length(results_r)) next
    seg_r <- results_r[[chr]][s:e]
    seg_r <- seg_r[seg_r != 0]
    if (length(seg_r) > 0L) {
      n_valid <- n_valid + 1L
      seg_mean <- mean(seg_r)
      vars[n_valid] <- mean((seg_r - seg_mean)^2)
    }
  }
  if (n_valid == 0L) return(NA_real_)
  stats::median(vars[seq_len(n_valid)])
}


#' Copy number profile abnormality (CPA)
#'
#' @note **Upstream bug replicated for conformance (CPA +1 overcount).**
#' Upstream Python (`overall_tools.py:146`) computes segment length as
#' `segment[2] - segment[1] + 1`. Given that CBS returns half-open
#' `[start, end)` intervals, the correct length is `end - start`, so the
#' `+ 1` overcounts by one bin per segment. We replicate this exactly
#' (`end - start + 1L`) for conformance with Python WisecondorX output.
#'
#' @keywords internal
.get_cpa <- function(results_c, binsize) {
  if (nrow(results_c) == 0L) return(NA_real_)
  x <- 0
  for (i in seq_len(nrow(results_c))) {
    # NOTE: upstream Python uses `segment[2] - segment[1] + 1` (overall_tools.py:146).
    # Given the half-open CBS convention, this overcounts by 1 bin per segment.
    # We replicate this for exact conformance with Python WisecondorX output.
    x <- x + (results_c$end[i] - results_c$start[i] + 1L) * binsize * abs(results_c$ratio[i])
  }
  round(x / nrow(results_c) * 1e-8, 5)
}
