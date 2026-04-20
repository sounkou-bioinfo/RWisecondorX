# WisecondorX native R implementation — reference building (newref stage)
#
# Ports newref_control.py and newref_tools.py from the upstream WisecondorX
# Python package. Original authors: Lennart Raman, Roy Straver, Wim Audenaert.
# Ported to R with credit under GPL-3.
#
# This file implements the complete `newref` pipeline:
# 1. Load and rescale binned samples to a common bin size.
# 2. Train a gender model (2-component GMM on Y-fractions).
# 3. Compute the global bin mask.
# 4. For each gender partition (Autosomal, Female gonosomal, Male gonosomal):
#    a. Normalize and mask.
#    b. Train PCA (5 components, ratio correction).
#    c. PCA distance filtering (remove anomalous bins).
#    d. Find within-sample reference bins (K-nearest neighbour).
#    e. Compute null ratios for between-sample Z-scoring.
# 5. Merge autosomal and gonosomal sub-references.
#
# The result is a list (the "reference object") that can be passed to
# `rwisecondorx_predict()` and optionally serialized to disk via `saveRDS()`.


#' Build a WisecondorX reference from binned samples
#'
#' Native R implementation of the WisecondorX `newref` pipeline. Takes a
#' list of binned samples (as returned by [bam_convert()] or loaded from NPZ
#' via reticulate) and builds a PCA-based reference suitable for
#' [rwisecondorx_predict()].
#'
#' The pipeline trains a gender model (2-component GMM on Y-fractions),
#' optionally applies gender correction for non-NIPT workflows, computes a
#' global bin mask, then builds three sub-references: autosomes (A), female
#' gonosomes (F), and male gonosomes (M). Each sub-reference includes PCA
#' components, within-sample reference bin indices and distances, and null
#' ratios for between-sample Z-scoring.
#'
#' This is a faithful port of the upstream Python `wisecondorx newref`,
#' crediting the original WisecondorX authors.
#'
#' @param samples List of sample objects, each a named list of integer vectors
#'   keyed by chromosome (`"1"`--`"24"`), as returned by [bam_convert()].
#'   At least 10 samples are required. Mutually exclusive with `bed_dir`.
#' @param binsize Integer; the target bin size in base pairs. All samples are
#'   rescaled to this size. Default `100000L`. Inferred from BED files when
#'   `bed_dir` is supplied and `sample_binsizes` is `NULL`.
#' @param sample_binsizes Optional integer vector of per-sample bin sizes. If
#'   `NULL` (default), all samples are assumed to already be at `binsize`.
#' @param nipt Logical; if `TRUE`, NIPT mode (no gender correction, no male
#'   gonosomal reference). Default `FALSE`.
#' @param refsize Integer; number of reference bin locations per target bin.
#'   Default `300L`.
#' @param yfrac Optional numeric; manual Y-fraction cutoff for gender
#'   classification. If `NULL` (default), the cutoff is derived from a GMM.
#' @param cpus Integer; number of threads for reference bin finding.
#'   Default `4L`.
#' @param bed_dir Optional character; path to a directory of 4-column bgzipped
#'   BED files (as written by [bam_convert_bed()]). All files matching
#'   `bed_pattern` are loaded in a single DuckDB pass via
#'   `rduckhts_tabix_multi()`. Mutually exclusive with `samples`.
#' @param bed_pattern Glob pattern for matching BED files inside `bed_dir`.
#'   Default `"*.bed.gz"`.
#' @param con Optional existing DuckDB connection. Used only when `bed_dir` is
#'   supplied. A temporary connection is created (and closed) if `NULL`.
#'
#' @return A list-like `WisecondorXReference` S7 object containing autosomal
#'   and (optionally) gonosomal sub-references. See Details for the full
#'   structure.
#'
#' @details
#' The returned reference object contains:
#' \describe{
#'   \item{binsize}{Integer; the reference bin size.}
#'   \item{is_nipt}{Logical; whether NIPT mode was used.}
#'   \item{trained_cutoff}{Numeric; Y-fraction gender cutoff.}
#'   \item{has_female}{Logical; whether a female gonosomal reference exists.}
#'   \item{has_male}{Logical; whether a male gonosomal reference exists.}
#'   \item{mask, bins_per_chr, masked_bins_per_chr, masked_bins_per_chr_cum,
#'     pca_components, pca_mean, indexes, distances, null_ratios}{Autosomal
#'     reference components.}
#'   \item{mask.F, ..., null_ratios.F}{Female gonosomal reference (if present).}
#'   \item{mask.M, ..., null_ratios.M}{Male gonosomal reference (if present).}
#' }
#'
#' @seealso [rwisecondorx_predict()], [scale_sample()], [bam_convert()]
#'
#' @export
rwisecondorx_newref <- function(samples         = NULL,
                               binsize         = 100000L,
                               sample_binsizes = NULL,
                               nipt            = FALSE,
                               refsize         = 300L,
                               yfrac           = NULL,
                               cpus            = 4L,
                               bed_dir         = NULL,
                               bed_pattern     = "*.bed.gz",
                               con             = NULL) {
  # ---------- bed_dir loading (alternative to samples list) ----------
  if (!is.null(bed_dir)) {
    if (!is.null(samples))
      stop("Provide 'samples' OR 'bed_dir', not both.")
    files <- sort(Sys.glob(file.path(bed_dir, bed_pattern)))
    if (!length(files))
      stop("No files matching '", bed_pattern, "' found in bed_dir: ", bed_dir)
    own_con <- is.null(con)
    if (own_con) {
      drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
      con <- DBI::dbConnect(drv)
      Rduckhts::rduckhts_load(con)
    }
    on.exit(if (own_con) DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    tbl <- "rwx_newref_beds"
    Rduckhts::rduckhts_tabix_multi(con, tbl, files, overwrite = TRUE)
    rows <- DBI::dbGetQuery(con, sprintf(
      "SELECT column0 AS chrom,
              CAST(column1 AS INTEGER) AS start_pos,
              CAST(column2 AS INTEGER) AS end_pos,
              CAST(column3 AS INTEGER) AS count,
              filename
       FROM %s ORDER BY filename, chrom, start_pos", tbl))
    rows$sample_name <- sub("\\.bed(\\.gz)?$", "", basename(rows$filename))
    sample_names <- unique(rows$sample_name)
    # Infer binsize from coordinate spacing if not supplied
    if (is.null(sample_binsizes)) {
      first_rows <- rows[rows$sample_name == sample_names[1L], ]
      if (nrow(first_rows) >= 2L) {
        inferred <- first_rows$end_pos[1L] - first_rows$start_pos[1L]
        if (inferred > 0L) binsize <- as.integer(inferred)
      }
    }
    samples <- lapply(sample_names, function(nm) {
      .bed_rows_to_wcx_sample(rows[rows$sample_name == nm, ])
    })
    message("Loaded ", length(samples), " samples from bed_dir.")
  }

  # ---------- input validation ----------
  stopifnot(is.list(samples), length(samples) >= 10L)
  binsize <- as.integer(binsize)
  refsize <- as.integer(refsize)
  cpus    <- as.integer(cpus)
  stopifnot(binsize > 0L, refsize > 0L, cpus >= 1L)
  if (!is.null(yfrac)) {
    stopifnot(is.numeric(yfrac), length(yfrac) == 1L, yfrac > 0, yfrac <= 1)
  }

  samples <- lapply(samples, .as_wcx_sample)

  # ---------- Step 1: rescale all samples ----------
  if (!is.null(sample_binsizes)) {
    stopifnot(length(sample_binsizes) == length(samples))
    for (i in seq_along(samples)) {
      samples[[i]] <- scale_sample(samples[[i]],
                                   from_size = sample_binsizes[i],
                                   to_size   = binsize)
    }
  }

  # ---------- Step 2: train gender model ----------
  gender_model <- .train_gender_model(samples, yfrac = yfrac)
  genders        <- gender_model$genders
  trained_cutoff <- gender_model$cutoff

  n_female <- sum(genders == "F")
  n_male   <- sum(genders == "M")

  if (n_female < 5L && isTRUE(nipt)) {
    stop(
      "NIPT-mode RWisecondorX reference building requires at least 5 female feti samples. ",
      "Refuse to silently disable NIPT mode.",
      call. = FALSE
    )
  }

  # ---------- Step 3: gender correction (non-NIPT) ----------
  if (!nipt) {
    for (i in seq_along(samples)) {
      samples[[i]] <- .gender_correct(samples[[i]], genders[i])
    }
  }

  # ---------- Step 4: compute global mask ----------
  mask_info    <- .get_mask(samples)
  total_mask   <- mask_info$mask
  bins_per_chr <- mask_info$bins_per_chr

  # AND with gender-specific masks
  if (n_female > 4L) {
    female_samples <- samples[genders == "F"]
    female_mask <- .get_mask(female_samples, ref_bins_per_chr = bins_per_chr)$mask
    total_mask <- total_mask & female_mask
  }
  if (n_male > 4L && !nipt) {
    male_samples <- samples[genders == "M"]
    male_mask <- .get_mask(male_samples, ref_bins_per_chr = bins_per_chr)$mask
    total_mask <- total_mask & male_mask
  }

  # ---------- Step 5: build sub-references ----------
  result <- list(
    binsize        = binsize,
    is_nipt        = nipt,
    trained_cutoff = trained_cutoff,
    has_female     = FALSE,
    has_male       = FALSE
  )

  # 5a: Autosomal reference (all samples)
  message("Building autosomal reference...")
  auto_ref <- .build_sub_reference(samples, "A", total_mask, bins_per_chr,
                                   refsize = refsize, cpus = cpus)
  for (nm in names(auto_ref)) result[[nm]] <- auto_ref[[nm]]

  # 5b: Female gonosomal reference
  if (n_female > 4L) {
    message("Building female gonosomal reference...")
    female_samples <- samples[genders == "F"]
    gon_f_ref <- .build_sub_reference(female_samples, "F", total_mask, bins_per_chr,
                                      refsize = refsize, cpus = cpus)
    result$has_female <- TRUE
    for (nm in names(gon_f_ref)) result[[paste0(nm, ".F")]] <- gon_f_ref[[nm]]
  }

  # 5c: Male gonosomal reference (non-NIPT only)
  if (!nipt && n_male > 4L) {
    message("Building male gonosomal reference...")
    male_samples <- samples[genders == "M"]
    gon_m_ref <- .build_sub_reference(male_samples, "M", total_mask, bins_per_chr,
                                      refsize = refsize, cpus = cpus)
    result$has_male <- TRUE
    for (nm in names(gon_m_ref)) result[[paste0(nm, ".M")]] <- gon_m_ref[[nm]]
  }

  .as_wcx_reference(result)
}


#' Build a sub-reference for one gender partition
#'
#' @param samples List of sample objects.
#' @param gender `"A"` (autosomes), `"F"` (female gonosomes), `"M"` (male gonosomes).
#' @param total_mask Logical mask over all 24 chromosomes.
#' @param bins_per_chr Integer vector of length 24.
#' @param refsize Number of reference bins per target.
#' @param cpus Number of threads.
#' @return Named list with mask, bins_per_chr, masked_bins_per_chr, etc.
#' @keywords internal
.build_sub_reference <- function(samples, gender, total_mask, bins_per_chr,
                                 refsize, cpus) {
  last_chr <- switch(gender, "A" = 22L, "F" = 23L, "M" = 24L)
  chr_range <- seq_len(last_chr)

  # Trim mask and bins_per_chr to the relevant chromosomes
  sub_bins <- bins_per_chr[chr_range]
  total_bins <- sum(sub_bins)
  sub_mask <- total_mask[seq_len(total_bins)]

  # Normalize and mask
  masked_data <- .normalize_and_mask(samples, chr_range, sub_mask)

  # Train PCA
  pca_result <- .train_pca(masked_data)
  corrected <- pca_result$corrected

  # PCA distance filtering: remove anomalous bins by comparing each bin's
  # cross-sample profile to the per-bin median profile, matching upstream.
  med_prof <- apply(corrected, 1, stats::median)
  dist_to_med <- rowSums((corrected - med_prof)^2)
  mad_val <- stats::median(abs(dist_to_med - stats::median(dist_to_med)))
  cutoff_pca <- max(stats::median(dist_to_med) + 10 * mad_val, 5.0)
  bad_bins <- dist_to_med > cutoff_pca

  if (any(bad_bins)) {
    n_removed <- sum(bad_bins)
    message(sprintf("  Removing %d anomalous bins (PCA distance cutoff=%.4f)",
                    n_removed, cutoff_pca))
    # Update mask: find which masked positions are bad
    masked_indices <- which(sub_mask)
    sub_mask[masked_indices[bad_bins]] <- FALSE

    # Re-train PCA on cleaned set
    masked_data <- .normalize_and_mask(samples, chr_range, sub_mask)
    pca_result <- .train_pca(masked_data)
    corrected <- pca_result$corrected
  }

  # Compute masked_bins_per_chr
  masked_bins_per_chr <- integer(length(sub_bins))
  cumpos <- 0L
  for (i in seq_along(sub_bins)) {
    masked_bins_per_chr[i] <- sum(sub_mask[(cumpos + 1L):(cumpos + sub_bins[i])])
    cumpos <- cumpos + sub_bins[i]
  }
  masked_bins_per_chr_cum <- cumsum(masked_bins_per_chr)

  # Find within-sample reference bins
  message("  Finding reference bins...")
  ref_result <- .get_reference(corrected, masked_bins_per_chr, masked_bins_per_chr_cum,
                               refsize, gender, cpus)

  list(
    mask                   = sub_mask,
    bins_per_chr           = sub_bins,
    masked_bins_per_chr     = masked_bins_per_chr,
    masked_bins_per_chr_cum = masked_bins_per_chr_cum,
    pca_components         = pca_result$components,
    pca_mean               = pca_result$center,
    indexes                = ref_result$indexes,
    distances              = ref_result$distances,
    null_ratios            = ref_result$null_ratios
  )
}


#' Find within-sample reference bins (KNN)
#'
#' For each target bin, finds the `refsize` most similar bins from other
#' chromosomes, measured by Euclidean distance in PCA-corrected space.
#' Also computes null ratios for between-sample Z-scoring.
#'
#' Mirrors `newref_tools.get_reference()` and `get_ref_for_bins()`.
#' Uses Rcpp + OpenMP for performance-critical distance computation and
#' null ratio calculation.
#'
#' @param pca_corrected Numeric matrix `(n_masked_bins, n_samples)`.
#' @param masked_bins_per_chr Integer vector.
#' @param masked_bins_per_chr_cum Integer vector (cumulative).
#' @param refsize Number of reference bins per target.
#' @param gender `"A"`, `"F"`, or `"M"`.
#' @param cpus Number of threads for parallel computation.
#'
#' @return List with `indexes` (integer matrix), `distances` (numeric matrix),
#'   `null_ratios` (numeric matrix).
#'
#' @keywords internal
.get_reference <- function(pca_corrected, masked_bins_per_chr,
                           masked_bins_per_chr_cum, refsize, gender, cpus) {
  n_bins <- nrow(pca_corrected)
  n_samples <- ncol(pca_corrected)

  # KNN reference bin finding via Rcpp (OpenMP-parallelized)
  knn_result <- knn_reference_cpp(
    pca_corrected,
    as.integer(masked_bins_per_chr),
    as.integer(masked_bins_per_chr_cum),
    as.integer(refsize),
    gender,
    as.integer(cpus)
  )

  indexes   <- knn_result$indexes
  distances <- knn_result$distances

  # Null ratios via Rcpp (OpenMP-parallelized)
  message("  Computing null ratios...")
  n_null <- min(n_samples, 100L)
  null_sample_idx <- sample(n_samples, n_null)
  null_ratios <- null_ratios_cpp(pca_corrected, indexes, null_sample_idx,
                                 as.integer(cpus))

  list(indexes = indexes, distances = distances, null_ratios = null_ratios)
}


#' Convert BED rows (from rduckhts_tabix_multi) to WisecondorX sample format
#'
#' @param rows data.frame with columns chrom, start_pos, end_pos, count.
#' @return Named list of integer vectors keyed by chromosome ("1"–"24").
#' @keywords internal
.bed_rows_to_wcx_sample <- function(rows) {
  sample <- stats::setNames(vector("list", 24L), as.character(1:24))

  # Normalize all unique chromosome names once
  unique_chroms <- unique(rows$chrom)
  normed <- stats::setNames(.normalize_chr_name(unique_chroms), unique_chroms)

  for (chr_raw in unique_chroms) {
    chr_key <- normed[chr_raw]
    if (!(chr_key %in% as.character(1:24))) next
    r <- rows[rows$chrom == chr_raw, ]
    r <- r[order(r$start_pos), ]
    sample[[chr_key]] <- as.integer(r$count)
  }

  # Fill missing chromosomes with zero-length vectors
  for (k in as.character(1:24)) {
    if (is.null(sample[[k]])) sample[[k]] <- integer(0L)
  }

  .as_wcx_sample(sample)
}
