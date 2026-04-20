#' Pre-compute and save per-bin GC content to a TSV.bgz file
#'
#' Runs \code{rduckhts_fasta_nuc()} once and writes the GC percentage table to
#' a bgzipped, tabix-indexed TSV file. Pass the resulting path to
#' \code{nipter_gc_correct(gc_table = ...)} to avoid recomputing GC content for
#' every sample in a large cohort.
#'
#' @param fasta Path to an indexed reference FASTA file (.fa/.fasta with .fai).
#' @param binsize Bin size in base pairs (default \code{50000L}).
#' @param out Path for the output file. The tabix index is written alongside as
#'   \code{<out>.tbi}.
#' @param con Optional open DBI connection with duckhts loaded.
#'
#' @return \code{out} invisibly.
#'
#' @details
#' The output is a 5-column, tab-delimited TSV.bgz:
#' \code{chrom}, \code{start}, \code{end}, \code{pct_gc}, \code{seq_len}.
#' Coordinates are 0-based half-open intervals (BED convention). Chromosomes
#' use no \code{chr} prefix (\code{1}--\code{22}, \code{X}, \code{Y}).
#' Bins where all bases are N are written with \code{pct_gc = NA}.
#'
#' @seealso [nipter_gc_correct()]
#'
#' @examples
#' \dontrun{
#' nipter_gc_precompute("hg38.fa", binsize = 50000L, out = "hg38_gc_50k.tsv.bgz")
#' cg <- nipter_gc_correct(cg, gc_table = "hg38_gc_50k.tsv.bgz")
#' }
#'
#' @export
nipter_gc_precompute <- function(fasta, binsize = 50000L, out, con = NULL) {
  stopifnot(is.character(fasta), length(fasta) == 1L, nzchar(fasta),
            file.exists(fasta))
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  stopifnot(is.character(out), length(out) == 1L, nzchar(out))
  binsize <- as.integer(binsize)

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  gc_list <- .get_gc_table(fasta, binsize, con)
  valid_chroms <- c(as.character(1:22), "X", "Y")

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)

  rows <- vector("list", length(valid_chroms))
  for (i in seq_along(valid_chroms)) {
    chr <- valid_chroms[i]
    gc  <- gc_list[[chr]]
    n   <- length(gc)
    if (n == 0L) next
    starts <- (seq_len(n) - 1L) * binsize
    rows[[i]] <- data.frame(
      chrom   = chr,
      start   = starts,
      end     = starts + binsize,
      pct_gc  = gc,
      seq_len = binsize,
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, Filter(Negate(is.null), rows))
  utils::write.table(df, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE, na = "NA")

  Rduckhts::rduckhts_bgzip(con, tmp, out, overwrite = TRUE)
  Rduckhts::rduckhts_tabix_index(con, out, preset = "bed", threads = 1L)

  message("GC table written to ", out, " (", nrow(df), " bins)")
  invisible(out)
}


#' GC-correct a NIPTeR sample or control group
#'
#' Adjusts bin counts for GC-content bias using either LOESS regression or
#' bin-weight normalisation.
#'
#' GC content can be supplied in three ways via the \code{gc_table} parameter:
#' \describe{
#'   \item{Pre-computed file}{Path to a TSV.bgz produced by
#'     \code{\link{nipter_gc_precompute}}. Fastest for large cohorts — compute
#'     once, reuse for every sample.}
#'   \item{In-memory list}{Named list of numeric vectors (one per chromosome)
#'     as returned by \code{.get_gc_table()}. Useful when chaining corrections
#'     within a session.}
#'   \item{FASTA path via \code{fasta}}{Compute GC on-the-fly for every call.
#'     Convenient for single-sample use; slow for many samples.}
#' }
#'
#' @param object A \code{NIPTeRSample} or \code{NIPTeRControlGroup}.
#' @param fasta Path to an indexed reference FASTA file (.fa/.fasta with
#'   .fai index). Ignored when \code{gc_table} is supplied.
#' @param method GC correction method: \code{"loess"} (default) or
#'   \code{"bin"} (bin-weight).
#' @param span LOESS smoothing parameter (only used when
#'   \code{method = "loess"}). Default \code{0.75}.
#' @param include_sex Logical; correct sex chromosomes (X, Y) as well?
#'   Default \code{FALSE}.
#' @param binsize Bin size used when binning the sample (default 50000L).
#'   Must match the binsize of the sample. Ignored when \code{gc_table} is a
#'   list (bin size is already encoded in the vector lengths).
#' @param gc_table Pre-computed GC table. Either a path to a TSV.bgz file (from
#'   \code{\link{nipter_gc_precompute}}) or the in-memory named list returned by
#'   a previous \code{.get_gc_table()} call. When \code{NULL} (default), GC
#'   content is computed from \code{fasta}.
#' @param con Optional open DBI connection with duckhts loaded. If \code{NULL},
#'   a temporary connection is created.
#'
#' @return A corrected copy of \code{object} with the same class. Correction
#'   status is updated from \code{"Uncorrected"} to \code{"GC Corrected"}.
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
#' @seealso [nipter_gc_precompute()], [nipter_bin_bam()], [nipter_chi_correct()],
#'   [nipter_as_control_group()]
#'
#' @examples
#' \dontrun{
#' # One-shot: compute GC and correct in one call
#' corrected <- nipter_gc_correct(sample, fasta = "hg38.fa")
#'
#' # Recommended for cohorts: precompute once, reuse
#' nipter_gc_precompute("hg38.fa", binsize = 50000L, out = "hg38_gc_50k.tsv.bgz")
#' cg <- nipter_gc_correct(cg, gc_table = "hg38_gc_50k.tsv.bgz")
#' test_sample <- nipter_gc_correct(test_sample, gc_table = "hg38_gc_50k.tsv.bgz")
#' }
#'
#' @export
nipter_gc_correct <- function(object,
                              fasta       = NULL,
                              method      = c("loess", "bin"),
                              span        = 0.75,
                              include_sex = FALSE,
                              binsize     = 50000L,
                              gc_table    = NULL,
                              con         = NULL) {
  method <- match.arg(method)
  stopifnot(is.numeric(span), length(span) == 1L, span > 0, span <= 1)
  stopifnot(is.logical(include_sex), length(include_sex) == 1L)
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  binsize <- as.integer(binsize)

  if (is.null(gc_table)) {
    # Fallback: compute from FASTA
    stopifnot(!is.null(fasta), is.character(fasta), length(fasta) == 1L,
              nzchar(fasta), file.exists(fasta))
  }

  if (.is_nipt_control_group_object(object)) {
    # Resolve gc_table once for all samples; convert path → list so we only
    # hit disk / run rduckhts_fasta_nuc once.
    resolved_gc <- .resolve_gc_table(gc_table, fasta, binsize, con)
    object <- .control_group_with_samples(
      object,
      lapply(object$samples, nipter_gc_correct,
             method = method, span = span,
             include_sex = include_sex, binsize = binsize,
             gc_table = resolved_gc, con = con)
    )
    return(object)
  }

  stopifnot(.is_nipt_sample_object(object))

  gc_tbl <- .resolve_gc_table(gc_table, fasta, binsize, con)

  if (method == "loess") {
    .gc_correct_loess(object, gc_tbl, span, include_sex)
  } else {
    .gc_correct_bin(object, gc_tbl, include_sex)
  }
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Resolve the gc_table argument to an in-memory named list.
#
# gc_table can be:
#   NULL          → compute from FASTA (requires fasta + binsize + con)
#   character     → path to a TSV.bgz produced by nipter_gc_precompute()
#   list          → already resolved; return as-is
#
.resolve_gc_table <- function(gc_table, fasta, binsize, con) {
  if (is.list(gc_table)) return(gc_table)      # already resolved

  if (is.character(gc_table)) {
    # Load from pre-computed TSV.bgz
    stopifnot(length(gc_table) == 1L, nzchar(gc_table), file.exists(gc_table))
    return(.load_gc_table(gc_table, con))
  }

  # gc_table is NULL → compute on-the-fly from FASTA
  .get_gc_table(fasta, binsize, con)
}


# Load a pre-computed GC table from a TSV.bgz file.
# Format: chrom  start  end  pct_gc  seq_len  (no header, NA written as "NA")
.load_gc_table <- function(path, con) {
  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  path_sql <- gsub("'", "''", normalizePath(path, mustWork = TRUE))
  sql <- sprintf("
    SELECT
      column0                       AS chrom,
      CAST(column1 AS INTEGER)      AS start_pos,
      CAST(column2 AS INTEGER)      AS end_pos,
      TRY_CAST(column3 AS DOUBLE)   AS pct_gc
    FROM read_tabix('%s')
  ", path_sql)

  rows <- DBI::dbGetQuery(con, sql)
  rows$chrom <- .normalize_chr_name(rows$chrom, xy_to_numeric = FALSE)

  valid_chroms <- c(as.character(1:22), "X", "Y")
  rows <- rows[rows$chrom %in% valid_chroms, , drop = FALSE]

  gc_list <- list()
  for (chr in valid_chroms) {
    sub <- rows[rows$chrom == chr, , drop = FALSE]
    if (nrow(sub) == 0L) {
      gc_list[[chr]] <- numeric(0L)
      next
    }
    sub <- sub[order(sub$start_pos), , drop = FALSE]
    gc_list[[chr]] <- sub$pct_gc   # NA already preserved via TRY_CAST
  }

  gc_list
}


# Fetch GC percentages for bins tiling all chromosomes directly from FASTA.
# Returns a named list: keys "1"-"22","X","Y", each a numeric vector of
# pct_gc values (one per bin), with NA for bins with no sequence or all-N.
.get_gc_table <- function(fasta, binsize, con) {
  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  nuc <- Rduckhts::rduckhts_fasta_nuc(con, fasta, bin_width = binsize)

  # Normalise chrom names: strip chr prefix, accept 1-22,X,Y
  nuc$chrom <- .normalize_chr_name(nuc$chrom, xy_to_numeric = FALSE)

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


.clamp_nonnegative_read_matrices <- function(mats) {
  lapply(mats, function(mat) {
    out <- mat
    bad <- !is.finite(out) | out < 0
    if (any(bad)) {
      out[bad] <- 0
    }
    out
  })
}


# LOESS GC correction for a single NIPTeRSample.
.gc_correct_loess <- function(sample, gc_table, span, include_sex) {
  auto_list <- .sample_autosomal_reads(sample)

  # For fitting: always use the summed (F+R) matrix via autosomal_matrix().
  summed_auto <- autosomal_matrix(sample)
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

  # Match upstream NIPTeR LOESS semantics: only strictly positive GC bins
  # with non-zero reads contribute to the fit. Our GC table uses NA for
  # all-N/empty bins, but true 0.0 GC bins should still be excluded here.
  valid <- !is.na(gc_auto) & gc_auto > 0 & reads_flat > 0
  if (sum(valid) < 10L) {
    stop(
      "Too few valid autosomal bins for LOESS GC correction in sample '",
      .sample_name(sample),
      "' (valid bins = ",
      sum(valid),
      "). Refuse to return an uncorrected sample.",
      call. = FALSE
    )
  }

  median_reads <- stats::median(reads_flat[valid])

  fit <- stats::loess(reads_flat[valid] ~ gc_auto[valid], span = span)
  fitted_vals <- stats::predict(fit)

  # Correction factor: median / fitted (so bins normalise to the median).
  # Keep the upstream fit-selection semantics (gc > 0 and count > 0), but do
  # not allow non-finite or exact-zero fitted values to propagate.
  correction <- rep(1.0, length(reads_flat))
  safe_fitted <- fitted_vals
  bad_fitted <- !is.finite(safe_fitted) | safe_fitted == 0
  safe_fitted[bad_fitted] <- median_reads
  correction[valid] <- median_reads / safe_fitted

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
  corrected_auto <- .clamp_nonnegative_read_matrices(corrected_auto)

  sample <- .sample_with_reads(sample, autosomal = corrected_auto)
  sample <- .sample_append_correction_step(sample, "autosomal", .nipt_gc_correction_step())

  # Sex chromosome correction (nearest-neighbour lookup, vectorized per chromosome)
  if (include_sex) {
    sex_list <- .sample_sex_reads(sample)

    # Pre-sort the autosomal GC/fitted-values for binary-search nearest-neighbour.
    gc_fitted_vals <- gc_auto[valid]
    fitted_reads   <- fitted_vals
    sort_idx       <- order(gc_fitted_vals)
    gc_sorted      <- gc_fitted_vals[sort_idx]
    fit_sorted     <- fitted_reads[sort_idx]

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

        # Vectorised nearest-neighbour: binary search via findInterval, then
        # check adjacent positions.  Replaces the O(n_bins * n_valid) loop.
        cf_vec <- vapply(gc_sex, function(g) {
          if (is.na(g) || g <= 0) return(1.0)
          pos <- findInterval(g, gc_sorted)
          candidates <- c(pos, pos + 1L)
          candidates <- candidates[candidates >= 1L & candidates <= length(gc_sorted)]
          if (length(candidates) == 0L) return(1.0)
          best <- candidates[which.min(abs(gc_sorted[candidates] - g))]
          median_reads / fit_sorted[best]
        }, numeric(1L))

        corrected[row_idx, ] <- sex_mat[row_idx, ] * cf_vec
      }
      corrected
    })
    corrected_sex <- .clamp_nonnegative_read_matrices(corrected_sex)

    sample <- .sample_with_reads(sample, sex = corrected_sex)
    sample <- .sample_append_correction_step(sample, "sex", .nipt_gc_correction_step())
  }

  sample
}


# Bin-weight GC correction for a single NIPTeRSample.
.gc_correct_bin <- function(sample, gc_table, include_sex) {
  auto_list <- .sample_autosomal_reads(sample)

  # Sum F+R via autosomal_matrix() for both S3 and S7 objects
  summed_auto <- autosomal_matrix(sample)
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
    stop(
      "No non-zero autosomal bins remain for GC correction in sample '",
      .sample_name(sample),
      "'. Refuse to return an uncorrected sample.",
      call. = FALSE
    )
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
  corrected_auto <- .clamp_nonnegative_read_matrices(corrected_auto)

  sample <- .sample_with_reads(sample, autosomal = corrected_auto)
  sample <- .sample_append_correction_step(sample, "autosomal", .nipt_gc_correction_step())

  # Sex chromosome correction (vectorised per chromosome)
  if (include_sex) {
    corrected_sex <- lapply(.sample_sex_reads(sample), function(sex_mat) {
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

        # Vectorised: look up bucket weights for all bins at once
        bm_vals <- bucket_mean[as.character(gc_sex_bucket)]
        use <- !is.na(gc_sex_bucket) & !is.na(bm_vals) & bm_vals > 0
        cf_vec <- rep(1.0, n_bins)
        cf_vec[use] <- global_mean / bm_vals[use]
        corrected[row_idx, ] <- sex_mat[row_idx, ] * cf_vec
      }
      corrected
    })
    corrected_sex <- .clamp_nonnegative_read_matrices(corrected_sex)

    sample <- .sample_with_reads(sample, sex = corrected_sex)
    sample <- .sample_append_correction_step(sample, "sex", .nipt_gc_correction_step())
  }

  sample
}
