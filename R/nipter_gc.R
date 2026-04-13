#' GC-correct a NIPTeR sample or control group
#'
#' Adjusts bin counts for GC-content bias using either LOESS regression or
#' bin-weight normalisation. GC content per bin is computed on-the-fly from the
#' reference FASTA via \code{\link[Rduckhts]{rduckhts_fasta_nuc}} rather than
#' from bundled precomputed tables (as the original NIPTeR does with
#' \code{sysdata.rda}).
#'
#' @param object A \code{NIPTeRSample} or \code{NIPTeRControlGroup}.
#' @param fasta Path to an indexed reference FASTA file (.fa/.fasta with
#'   .fai index).
#' @param method GC correction method: \code{"loess"} (default) or
#'   \code{"bin"} (bin-weight).
#' @param span LOESS smoothing parameter (only used when
#'   \code{method = "loess"}). Default \code{0.75}.
#' @param include_sex Logical; correct sex chromosomes (X, Y) as well?
#'   Default \code{FALSE}.
#' @param binsize Bin size used when binning the sample (default 50000L).
#'   Must match the binsize of the sample.
#' @param con Optional open DBI connection with duckhts loaded. If \code{NULL},
#'   a temporary connection is created.
#'
#' @return A corrected copy of \code{object} with the same class. Correction
#'   status is updated from \code{"Uncorrected"} to \code{"GC corrected"}.
#'
#' @details
#' **LOESS method** (default): Fits a LOESS curve of read counts vs GC
#' percentage across all autosomal bins with known GC and non-zero reads.
#' Each bin is then scaled by \code{median(counts) / fitted(loess)}, so that
#' all bins are normalised to the genome-wide median. This is the NIPTeR
#' default method.
#'
#' **Bin-weight method**: Groups bins by GC percentage (0.1\% resolution),
#' computes the mean read count per GC bucket, then scales each bin by
#' \code{global_mean / bucket_mean}. Faster than LOESS but less smooth.
#'
#' Sex chromosome correction (when \code{include_sex = TRUE}) uses a
#' nearest-neighbour lookup against the autosomal LOESS curve (LOESS method)
#' or the same GC bucket weights (bin-weight method).
#'
#' @seealso [nipter_bin_bam()], [nipter_chi_correct()],
#'   [nipter_as_control_group()]
#'
#' @examples
#' \dontrun{
#' sample <- nipter_bin_bam("sample.bam")
#' corrected <- nipter_gc_correct(sample, fasta = "hg38.fa")
#'
#' # Correct an entire control group
#' cg <- nipter_as_control_group(samples)
#' cg_corrected <- nipter_gc_correct(cg, fasta = "hg38.fa")
#' }
#'
#' @export
nipter_gc_correct <- function(object,
                              fasta,
                              method      = c("loess", "bin"),
                              span        = 0.75,
                              include_sex = FALSE,
                              binsize     = 50000L,
                              con         = NULL) {
  method <- match.arg(method)
  stopifnot(is.character(fasta), length(fasta) == 1L, nzchar(fasta))
  stopifnot(file.exists(fasta))
  stopifnot(is.numeric(span), length(span) == 1L, span > 0, span <= 1)
  stopifnot(is.logical(include_sex), length(include_sex) == 1L)
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  binsize <- as.integer(binsize)

  if (inherits(object, "NIPTeRControlGroup")) {
    object$samples <- lapply(object$samples, nipter_gc_correct,
                             fasta = fasta, method = method, span = span,
                             include_sex = include_sex, binsize = binsize,
                             con = con)
    object$correction_status_autosomal <- unique(unlist(lapply(
      object$samples, `[[`, "correction_status_autosomal"
    )))
    if (include_sex) {
      object$correction_status_sex <- unique(unlist(lapply(
        object$samples, `[[`, "correction_status_sex"
      )))
    }
    return(object)
  }

  stopifnot(inherits(object, "NIPTeRSample"))

  # Fetch GC content for all bins via rduckhts_fasta_nuc
  gc_table <- .get_gc_table(fasta, binsize, con)

  if (method == "loess") {
    .gc_correct_loess(object, gc_table, span, include_sex)
  } else {
    .gc_correct_bin(object, gc_table, include_sex)
  }
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Fetch GC percentages for 50kb bins tiling all chromosomes.
# Returns a named list: keys "1"-"22","X","Y", each a numeric vector of
# pct_gc values (one per bin), with NA for bins with no sequence or all-N.
.get_gc_table <- function(fasta, binsize, con) {
  own_con <- is.null(con)
  if (own_con) {
    if (!requireNamespace("Rduckhts", quietly = TRUE)) {
      stop("Rduckhts is required for GC correction.", call. = FALSE)
    }
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  nuc <- Rduckhts::rduckhts_fasta_nuc(con, fasta, bin_width = binsize)

  # Normalise chrom names: strip chr prefix, accept 1-22,X,Y
  nuc$chrom <- sub("^[Cc][Hh][Rr]", "", nuc$chrom)

  valid_chroms <- c(as.character(1:22), "X", "Y")
  nuc <- nuc[nuc$chrom %in% valid_chroms, , drop = FALSE]

  # Build per-chromosome GC vectors, ordered by start position
  gc_list <- list()
  for (chr in valid_chroms) {
    sub_nuc <- nuc[nuc$chrom == chr, , drop = FALSE]
    if (nrow(sub_nuc) == 0L) {
      gc_list[[chr]] <- numeric(0L)
      next
    }
    sub_nuc <- sub_nuc[order(sub_nuc$start), , drop = FALSE]
    gc <- as.numeric(sub_nuc$pct_gc)
    # Mark bins with all Ns or no sequence as NA
    seq_len_col <- as.numeric(sub_nuc$seq_len)
    num_n <- as.numeric(sub_nuc$num_n)
    gc[seq_len_col == 0L | num_n == seq_len_col] <- NA_real_
    gc_list[[chr]] <- gc
  }

  gc_list
}


# LOESS GC correction for a single NIPTeRSample.
.gc_correct_loess <- function(sample, gc_table, span, include_sex) {
  auto_list <- sample$autosomal_chromosome_reads
  is_ss     <- inherits(sample, "SeparatedStrands")

  # For fitting: always use the summed (F+R) matrix.
  # CombinedStrands: auto_list has 1 element.
  # SeparatedStrands: Reduce("+", auto_list) sums fwd + rev.
  if (is_ss) {
    summed_auto <- Reduce("+", auto_list)
    rownames(summed_auto) <- as.character(1:22)
  } else {
    summed_auto <- auto_list[[1L]]
  }
  n_bins <- ncol(summed_auto)

  # Build combined GC vector (all autosomal bins concatenated)
  gc_auto <- unlist(lapply(as.character(1:22), function(chr) {
    gc <- gc_table[[chr]]
    if (length(gc) < n_bins) {
      gc <- c(gc, rep(NA_real_, n_bins - length(gc)))
    } else if (length(gc) > n_bins) {
      gc <- gc[seq_len(n_bins)]
    }
    gc
  }))

  # Flatten summed autosomal reads row by row (chr1 bins, chr2 bins, ...)
  reads_flat <- as.numeric(t(summed_auto))  # 22*n_bins vector

  # Valid bins: known GC and non-zero reads
  valid <- !is.na(gc_auto) & reads_flat > 0
  if (sum(valid) < 10L) {
    warning("Too few valid bins for LOESS GC correction; returning uncorrected.",
            call. = FALSE)
    return(sample)
  }

  median_reads <- stats::median(reads_flat[valid])

  fit <- stats::loess(reads_flat[valid] ~ gc_auto[valid], span = span)
  fitted_vals <- stats::predict(fit)

  # Correction factor: median / fitted (so bins normalise to the median)
  correction <- rep(1.0, length(reads_flat))
  correction[valid] <- median_reads / fitted_vals

  # Apply correction to each autosomal matrix in the list
  corrected_auto <- lapply(auto_list, function(mat) {
    corrected <- mat
    for (i in seq_len(22L)) {
      offset <- (i - 1L) * n_bins
      idx <- seq(offset + 1L, offset + n_bins)
      corrected[i, ] <- mat[i, ] * correction[idx]
    }
    corrected
  })

  sample$autosomal_chromosome_reads <- corrected_auto
  sample$correction_status_autosomal <- .update_correction_status(
    sample$correction_status_autosomal, "GC corrected"
  )

  # Sex chromosome correction (nearest-neighbour lookup)
  if (include_sex) {
    sex_list <- sample$sex_chromosome_reads

    # For each sex chromosome bin, find the nearest GC in the LOESS curve
    gc_fitted_vals <- gc_auto[valid]
    fitted_reads   <- fitted_vals

    corrected_sex <- lapply(sex_list, function(sex_mat) {
      corrected <- sex_mat
      for (chr_label in c("X", "Y")) {
        gc_sex <- gc_table[[chr_label]]
        if (length(gc_sex) < n_bins) {
          gc_sex <- c(gc_sex, rep(NA_real_, n_bins - length(gc_sex)))
        } else if (length(gc_sex) > n_bins) {
          gc_sex <- gc_sex[seq_len(n_bins)]
        }
        # Find the row — for SeparatedStrands rows are "XF"/"YF" or "XR"/"YR"
        row_idx <- grep(paste0("^", chr_label), rownames(sex_mat))
        if (length(row_idx) == 0L) next

        for (b in seq_len(n_bins)) {
          if (is.na(gc_sex[b]) || sex_mat[row_idx, b] == 0) next
          dists <- abs(gc_fitted_vals - gc_sex[b])
          nearest <- which.min(dists)
          cf <- median_reads / fitted_reads[nearest]
          corrected[row_idx, b] <- sex_mat[row_idx, b] * cf
        }
      }
      corrected
    })

    sample$sex_chromosome_reads <- corrected_sex
    sample$correction_status_sex <- .update_correction_status(
      sample$correction_status_sex, "GC corrected"
    )
  }

  sample
}


# Bin-weight GC correction for a single NIPTeRSample.
.gc_correct_bin <- function(sample, gc_table, include_sex) {
  auto_list <- sample$autosomal_chromosome_reads
  is_ss     <- inherits(sample, "SeparatedStrands")

  # Sum F+R for SeparatedStrands; use single matrix for CombinedStrands
  if (is_ss) {
    summed_auto <- Reduce("+", auto_list)
    rownames(summed_auto) <- as.character(1:22)
  } else {
    summed_auto <- auto_list[[1L]]
  }
  n_bins <- ncol(summed_auto)

  # Flatten summed autosomal reads and GC values
  reads_flat <- as.numeric(t(summed_auto))

  gc_auto <- unlist(lapply(as.character(1:22), function(chr) {
    gc <- gc_table[[chr]]
    if (length(gc) < n_bins) {
      gc <- c(gc, rep(NA_real_, n_bins - length(gc)))
    } else if (length(gc) > n_bins) {
      gc <- gc[seq_len(n_bins)]
    }
    gc
  }))

  # Zero out bins with unknown GC
  reads_flat[is.na(gc_auto)] <- 0

  # Bucket bins by GC at 0.1% resolution (multiply GC by 1000, round)
  gc_bucket <- round(gc_auto * 1000)

  # Compute mean reads per non-empty bin in each GC bucket
  buckets <- unique(gc_bucket[!is.na(gc_bucket)])
  bucket_mean <- stats::setNames(numeric(length(buckets)), buckets)

  for (b in buckets) {
    idx <- which(gc_bucket == b)
    nonzero <- reads_flat[idx] > 0
    if (any(nonzero)) {
      bucket_mean[as.character(b)] <- sum(reads_flat[idx]) / sum(nonzero)
    }
  }

  # Global mean: total reads / total non-empty bins
  total_reads   <- sum(reads_flat)
  total_nonzero <- sum(reads_flat > 0)
  if (total_nonzero == 0L) {
    warning("No non-zero autosomal bins; returning uncorrected.", call. = FALSE)
    return(sample)
  }
  global_mean <- total_reads / total_nonzero

  # Weight per bucket: global_mean / bucket_mean
  weights <- rep(1.0, length(reads_flat))
  for (b in buckets) {
    bm <- bucket_mean[as.character(b)]
    if (bm > 0) {
      idx <- which(gc_bucket == b)
      weights[idx] <- global_mean / bm
    }
  }

  # Apply weights to each autosomal matrix in the list
  corrected_auto <- lapply(auto_list, function(mat) {
    corrected <- mat
    for (i in seq_len(22L)) {
      offset <- (i - 1L) * n_bins
      idx <- seq(offset + 1L, offset + n_bins)
      corrected[i, ] <- mat[i, ] * weights[idx]
    }
    corrected
  })

  sample$autosomal_chromosome_reads <- corrected_auto
  sample$correction_status_autosomal <- .update_correction_status(
    sample$correction_status_autosomal, "GC corrected"
  )

  # Sex chromosome correction
  if (include_sex) {
    corrected_sex <- lapply(sample$sex_chromosome_reads, function(sex_mat) {
      corrected <- sex_mat
      for (chr_label in c("X", "Y")) {
        gc_sex <- gc_table[[chr_label]]
        if (length(gc_sex) < n_bins) {
          gc_sex <- c(gc_sex, rep(NA_real_, n_bins - length(gc_sex)))
        } else if (length(gc_sex) > n_bins) {
          gc_sex <- gc_sex[seq_len(n_bins)]
        }
        row_idx <- grep(paste0("^", chr_label), rownames(sex_mat))
        if (length(row_idx) == 0L) next
        gc_sex_bucket <- round(gc_sex * 1000)

        for (b_pos in seq_len(n_bins)) {
          if (is.na(gc_sex_bucket[b_pos]) || sex_mat[row_idx, b_pos] == 0) next
          bm <- bucket_mean[as.character(gc_sex_bucket[b_pos])]
          if (!is.na(bm) && bm > 0) {
            corrected[row_idx, b_pos] <- sex_mat[row_idx, b_pos] * (global_mean / bm)
          }
        }
      }
      corrected
    })

    sample$sex_chromosome_reads <- corrected_sex
    sample$correction_status_sex <- .update_correction_status(
      sample$correction_status_sex, "GC corrected"
    )
  }

  sample
}


# Update correction status vector: add new status, remove "Uncorrected".
.update_correction_status <- function(current, new_status) {
  updated <- unique(c(current, new_status))
  updated <- updated[updated != "Uncorrected"]
  updated
}
