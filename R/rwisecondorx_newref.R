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
#' @param gender_model_names Character scalar; mclust covariance model name used
#'   for the Y-fraction Gaussian mixture during gender-model fitting. Default
#'   `"V"`, matching the current upstream WisecondorX implementation.
#' @param female_partition_min_samples Integer; minimum female samples required
#'   before applying female-specific masking and building the female gonosomal
#'   partition. In NIPT mode this is also the minimum female count required to
#'   proceed. Default `5L`.
#' @param male_partition_min_samples Integer; minimum male samples required
#'   before applying male-specific masking and building the male gonosomal
#'   partition in non-NIPT mode. Default `5L`.
#' @param mask_min_median_coverage_fraction Numeric scalar; bins with summed
#'   normalized coverage below this fraction of the global median are masked
#'   before model fitting. Default `0.05`.
#' @param gender_model_grid_min Numeric scalar; lower bound of the density grid
#'   used to locate the Y-fraction cutoff when `yfrac` is not provided.
#'   Default `0`.
#' @param gender_model_grid_max Numeric scalar; upper bound of that density
#'   grid. Default `0.02`.
#' @param gender_model_grid_length Integer; number of grid points used to
#'   locate the Y-fraction cutoff. Default `5000L`.
#' @param pca_components Integer; PCA components retained during between-sample
#'   normalization. Default `5L`.
#' @param pca_distance_min_class_bins Integer; minimum masked bins required
#'   before applying class-specific PCA-distance pruning. Default `10L`.
#' @param pca_distance_autosome_mad_multiplier Numeric; MAD multiplier used for
#'   autosomal PCA-distance pruning. Default `20`.
#' @param pca_distance_autosome_floor Numeric; minimum autosomal PCA-distance
#'   cutoff. Default `10`.
#' @param pca_distance_chrX_mad_multiplier Numeric; MAD multiplier used for
#'   chrX PCA-distance pruning. Default `20`.
#' @param pca_distance_chrX_floor Numeric; minimum chrX PCA-distance cutoff.
#'   Default `10`.
#' @param pca_distance_chrY_mad_multiplier Numeric; MAD multiplier used for
#'   chrY PCA-distance pruning. Default `50`.
#' @param pca_distance_chrY_floor Numeric; minimum chrY PCA-distance cutoff.
#'   Default `15`.
#' @param null_ratio_max_samples Integer; maximum number of cohort samples used
#'   when sampling columns for null-ratio estimation. Default `100L`, matching
#'   the current upstream implementation.
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
#'   \item{algorithm_params}{Named list of native reference-building
#'     parameters used to create the object.}
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
                               gender_model_names = "V",
                               female_partition_min_samples = 5L,
                               male_partition_min_samples = 5L,
                               mask_min_median_coverage_fraction = 0.05,
                               gender_model_grid_min = 0,
                               gender_model_grid_max = 0.02,
                               gender_model_grid_length = 5000L,
                               pca_components = 5L,
                               pca_distance_min_class_bins = 10L,
                               pca_distance_autosome_mad_multiplier = 20,
                               pca_distance_autosome_floor = 10,
                               pca_distance_chrX_mad_multiplier = 20,
                               pca_distance_chrX_floor = 10,
                               pca_distance_chrY_mad_multiplier = 50,
                               pca_distance_chrY_floor = 15,
                               null_ratio_max_samples = 100L,
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
    # dbGetQuery() fully materializes the result into R before the on.exit
    # disconnect fires, so downstream sample construction does not retain any
    # lazy dependency on the DuckDB connection.
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
  null_ratio_max_samples <- as.integer(null_ratio_max_samples)
  female_partition_min_samples <- as.integer(female_partition_min_samples)
  male_partition_min_samples <- as.integer(male_partition_min_samples)
  gender_model_grid_length <- as.integer(gender_model_grid_length)
  pca_components <- as.integer(pca_components)
  pca_distance_min_class_bins <- as.integer(pca_distance_min_class_bins)
  stopifnot(is.character(gender_model_names), length(gender_model_names) == 1L,
            nzchar(gender_model_names))
  stopifnot(binsize > 0L, refsize > 0L, cpus >= 1L)
  if (!is.null(yfrac)) {
    stopifnot(is.numeric(yfrac), length(yfrac) == 1L, yfrac > 0, yfrac <= 1)
  }
  stopifnot(
    is.numeric(mask_min_median_coverage_fraction),
    length(mask_min_median_coverage_fraction) == 1L,
    is.finite(mask_min_median_coverage_fraction),
    mask_min_median_coverage_fraction >= 0
  )
  stopifnot(
    is.numeric(gender_model_grid_min),
    is.numeric(gender_model_grid_max),
    length(gender_model_grid_min) == 1L,
    length(gender_model_grid_max) == 1L,
    is.finite(gender_model_grid_min),
    is.finite(gender_model_grid_max),
    gender_model_grid_max > gender_model_grid_min
  )
  stopifnot(
    female_partition_min_samples >= 1L,
    male_partition_min_samples >= 1L,
    gender_model_grid_length >= 2L,
    pca_components >= 1L,
    pca_distance_min_class_bins >= 1L
  )
  stopifnot(null_ratio_max_samples >= 1L)

  samples <- lapply(samples, .as_wcx_sample)

  pca_distance_params <- list(
    min_class_bins = pca_distance_min_class_bins,
    autosome_mad_multiplier = pca_distance_autosome_mad_multiplier,
    autosome_floor = pca_distance_autosome_floor,
    chrX_mad_multiplier = pca_distance_chrX_mad_multiplier,
    chrX_floor = pca_distance_chrX_floor,
    chrY_mad_multiplier = pca_distance_chrY_mad_multiplier,
    chrY_floor = pca_distance_chrY_floor
  )

  # ---------- Step 1: rescale all samples ----------
  if (!is.null(sample_binsizes)) {
    stopifnot(length(sample_binsizes) == length(samples))
    stopifnot(is.numeric(sample_binsizes), all(is.finite(sample_binsizes)),
              all(sample_binsizes > 0))
    for (i in seq_along(samples)) {
      samples[[i]] <- scale_sample(samples[[i]],
                                   from_size = sample_binsizes[i],
                                   to_size   = binsize)
    }
  }

  # ---------- Step 2: train gender model ----------
  gender_model <- .train_gender_model(
    samples,
    yfrac = yfrac,
    gender_model_names = gender_model_names,
    density_grid_min = gender_model_grid_min,
    density_grid_max = gender_model_grid_max,
    density_grid_length = gender_model_grid_length
  )
  genders        <- gender_model$genders
  trained_cutoff <- gender_model$cutoff

  n_female <- sum(genders == "F")
  n_male   <- sum(genders == "M")

  if (n_female < female_partition_min_samples && isTRUE(nipt)) {
    stop(
      sprintf(
        "NIPT-mode RWisecondorX reference building requires at least %d female fetal samples. ",
        female_partition_min_samples
      ),
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
  mask_info <- .get_mask(
    samples,
    min_median_coverage_fraction = mask_min_median_coverage_fraction
  )
  total_mask   <- mask_info$mask
  bins_per_chr <- mask_info$bins_per_chr

  # AND with gender-specific masks
  if (n_female >= female_partition_min_samples) {
    female_samples <- samples[genders == "F"]
    female_mask <- .get_mask(
      female_samples,
      ref_bins_per_chr = bins_per_chr,
      min_median_coverage_fraction = mask_min_median_coverage_fraction
    )$mask
    total_mask <- total_mask & female_mask
  }
  if (n_male >= male_partition_min_samples && !nipt) {
    male_samples <- samples[genders == "M"]
    male_mask <- .get_mask(
      male_samples,
      ref_bins_per_chr = bins_per_chr,
      min_median_coverage_fraction = mask_min_median_coverage_fraction
    )$mask
    total_mask <- total_mask & male_mask
  }

  # ---------- Step 5: build sub-references ----------
  result <- list(
    binsize        = binsize,
    is_nipt        = nipt,
    trained_cutoff = trained_cutoff,
    has_female     = FALSE,
    has_male       = FALSE,
    algorithm_params = list(
      refsize = refsize,
      yfrac = yfrac,
      gender_model_names = gender_model_names,
      female_partition_min_samples = female_partition_min_samples,
      male_partition_min_samples = male_partition_min_samples,
      mask_min_median_coverage_fraction = mask_min_median_coverage_fraction,
      gender_model_grid_min = gender_model_grid_min,
      gender_model_grid_max = gender_model_grid_max,
      gender_model_grid_length = gender_model_grid_length,
      pca_components = pca_components,
      pca_distance_min_class_bins = pca_distance_min_class_bins,
      pca_distance_autosome_mad_multiplier = pca_distance_autosome_mad_multiplier,
      pca_distance_autosome_floor = pca_distance_autosome_floor,
      pca_distance_chrX_mad_multiplier = pca_distance_chrX_mad_multiplier,
      pca_distance_chrX_floor = pca_distance_chrX_floor,
      pca_distance_chrY_mad_multiplier = pca_distance_chrY_mad_multiplier,
      pca_distance_chrY_floor = pca_distance_chrY_floor,
      null_ratio_max_samples = null_ratio_max_samples,
      cpus = cpus
    )
  )

  # 5a: Autosomal reference (all samples)
  message("Building autosomal reference...")
  auto_ref <- .build_sub_reference(
    samples,
    "A",
    total_mask,
    bins_per_chr,
    refsize = refsize,
    cpus = cpus,
    pca_components = pca_components,
    pca_distance_params = pca_distance_params,
    null_ratio_max_samples = null_ratio_max_samples
  )
  result <- modifyList(result, auto_ref)

  # Keep autosomal masks aligned across branches: once the autosomal
  # reference has applied its own PCA-distance pruning, reuse that finalized
  # autosomal prefix when building sex-specific branches.
  autosomal_span <- as.integer(sum(bins_per_chr[seq_len(22L)]))
  aligned_total_mask <- total_mask
  aligned_total_mask[seq_len(autosomal_span)] <- auto_ref$mask

  # 5b: Female gonosomal reference
  if (n_female >= female_partition_min_samples) {
    message("Building female gonosomal reference...")
    female_samples <- samples[genders == "F"]
    gon_f_ref <- .build_sub_reference(
      female_samples,
      "F",
      aligned_total_mask,
      bins_per_chr,
      refsize = refsize,
      cpus = cpus,
      pca_components = pca_components,
      pca_distance_params = pca_distance_params,
      null_ratio_max_samples = null_ratio_max_samples
    )
    result$has_female <- TRUE
    result <- modifyList(
      result,
      stats::setNames(gon_f_ref, paste0(names(gon_f_ref), ".F"))
    )
  }

  # 5c: Male gonosomal reference (non-NIPT only)
  if (!nipt && n_male >= male_partition_min_samples) {
    message("Building male gonosomal reference...")
    male_samples <- samples[genders == "M"]
    gon_m_ref <- .build_sub_reference(
      male_samples,
      "M",
      aligned_total_mask,
      bins_per_chr,
      refsize = refsize,
      cpus = cpus,
      pca_components = pca_components,
      pca_distance_params = pca_distance_params,
      null_ratio_max_samples = null_ratio_max_samples
    )
    result$has_male <- TRUE
    result <- modifyList(
      result,
      stats::setNames(gon_m_ref, paste0(names(gon_m_ref), ".M"))
    )
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
                                 refsize, cpus,
                                 pca_components,
                                 pca_distance_params,
                                 null_ratio_max_samples) {
  last_chr <- switch(gender, "A" = 22L, "F" = 23L, "M" = 24L)
  chr_range <- seq_len(last_chr)

  # Trim mask and bins_per_chr to the relevant chromosomes
  sub_bins <- bins_per_chr[chr_range]
  total_bins <- sum(sub_bins)
  sub_mask <- total_mask[seq_len(total_bins)]

  # Normalize and mask
  masked_data <- .normalize_and_mask(samples, chr_range, sub_mask)

  # Train PCA
  pca_result <- .train_pca(masked_data, n_comp = pca_components)
  corrected <- pca_result$corrected

  masked_indices <- which(sub_mask)
  masked_chr_ids <- rep(seq_along(sub_bins), times = sub_bins)[sub_mask]
  bad_bins_mask <- rep(FALSE, length(masked_indices))

  class_settings <- switch(
    gender,
    "A" = list(list(name = "autosomes", mask = masked_chr_ids <= 22L,
                    mad_multiplier = pca_distance_params$autosome_mad_multiplier,
                    floor = pca_distance_params$autosome_floor)),
    "F" = list(list(name = "chrX", mask = masked_chr_ids == 23L,
                    mad_multiplier = pca_distance_params$chrX_mad_multiplier,
                    floor = pca_distance_params$chrX_floor)),
    "M" = list(
      list(name = "chrX", mask = masked_chr_ids == 23L,
           mad_multiplier = pca_distance_params$chrX_mad_multiplier,
           floor = pca_distance_params$chrX_floor),
      list(name = "chrY", mask = masked_chr_ids == 24L,
           mad_multiplier = pca_distance_params$chrY_mad_multiplier,
           floor = pca_distance_params$chrY_floor)
    ),
    stop("Unsupported sub-reference gender: ", gender, call. = FALSE)
  )

  for (cfg in class_settings) {
    class_mask <- cfg$mask
    if (sum(class_mask) < pca_distance_params$min_class_bins) {
      next
    }
    class_data <- corrected[class_mask, , drop = FALSE]
    med_prof <- apply(class_data, 2L, stats::median)
    dist_to_med <- rowSums(sweep(class_data, 2L, med_prof, "-")^2)
    cutoff <- .wcx_robust_cutoff(dist_to_med, cfg$mad_multiplier, cfg$floor)
    class_bad_bins <- dist_to_med > cutoff
    bad_bins_mask[class_mask] <- class_bad_bins
    if (any(class_bad_bins)) {
      message(sprintf(
        "  Removing %d anomalous %s bins (PCA distance cutoff=%.4f)",
        sum(class_bad_bins), cfg$name, cutoff
      ))
    }
  }

  if (any(bad_bins_mask)) {
    sub_mask[masked_indices[bad_bins_mask]] <- FALSE

    masked_data <- .normalize_and_mask(samples, chr_range, sub_mask)
    pca_result <- .train_pca(masked_data, n_comp = pca_components)
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
                               refsize, gender, cpus,
                               null_ratio_max_samples = null_ratio_max_samples)

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

.wcx_robust_cutoff <- function(distances, mad_multiplier, floor) {
  med <- stats::median(distances)
  mad <- stats::median(abs(distances - med))
  max(med + mad_multiplier * mad, floor)
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
                           masked_bins_per_chr_cum, refsize, gender, cpus,
                           null_ratio_max_samples = 100L) {
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
  n_null <- min(n_samples, as.integer(null_ratio_max_samples))
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
  sample <- stats::setNames(vector("list", 24L), as.character(seq_len(24L)))

  # Normalize all unique chromosome names once
  unique_chroms <- unique(rows$chrom)
  normed <- stats::setNames(.normalize_chr_name(unique_chroms), unique_chroms)

  for (chr_raw in unique_chroms) {
    chr_key <- normed[chr_raw]
    if (!(chr_key %in% as.character(seq_len(24L)))) next
    r <- rows[rows$chrom == chr_raw, ]
    r <- r[order(r$start_pos), ]
    sample[[chr_key]] <- as.integer(r$count)
  }

  # Fill missing chromosomes with zero-length vectors
  for (k in as.character(seq_len(24L))) {
    if (is.null(sample[[k]])) sample[[k]] <- integer(0L)
  }

  .as_wcx_sample(sample)
}
