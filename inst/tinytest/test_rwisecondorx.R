library(tinytest)
library(RWisecondorX)

# ---------------------------------------------------------------------------
# Tests for the native WisecondorX R implementation:
#   rwisecondorx_utils.R   — scale_sample, gender_correct, predict_gender,
#                            get_mask, train_pca, project_pc
#   rwisecondorx_newref.R  — rwisecondorx_newref, rwisecondorx_ref_qc
#   rwisecondorx_predict.R — rwisecondorx_predict, .log_trans, .inflate_results
#   rwisecondorx_cbs.R     — .exec_cbs, .split_na_segments
#   rwisecondorx_output.R  — write_wisecondorx_output
#
# Tests use synthetic bin-count data (no BAM fixture or real genome required).
# ---------------------------------------------------------------------------

# Expose internal functions for unit testing
for (.fn in c(".gender_correct", ".predict_gender", ".get_mask",
              ".train_pca", ".project_pc", ".log_trans",
              ".inflate_results", ".normalize_and_mask",
              ".exec_cbs", ".split_na_segments",
              ".get_segment_zscores", ".call_aberrations",
              ".idx_to_chr_name", ".train_gender_model")) {
  tryCatch(
    assign(.fn, getFromNamespace(.fn, "RWisecondorX")),
    error = function(e) NULL
  )
}
rm(.fn)

# ---------------------------------------------------------------------------
# Helpers: build synthetic samples
# ---------------------------------------------------------------------------

# Create a fake sample: named list of integer vectors keyed "1" through "24".
# Each chromosome gets `n_bins` bins with counts drawn from Poisson(lambda).
.make_sample <- function(n_bins = 50L, lambda = 100, seed = NULL,
                         y_extra = 0) {
  if (!is.null(seed)) set.seed(seed)
  sample <- vector("list", 24L)
  names(sample) <- as.character(1:24)
  for (chr in 1:24) {
    lam <- lambda
    if (chr == 24L) lam <- lambda * (0.001 + y_extra)
    sample[[as.character(chr)]] <- as.integer(rpois(n_bins, lam))
  }
  sample
}

# Create a "female" sample (low Y counts)
.make_female_sample <- function(n_bins = 50L, lambda = 100, seed = NULL) {
  .make_sample(n_bins = n_bins, lambda = lambda, seed = seed, y_extra = 0)
}

# Create a "male" sample (higher Y counts)
.make_male_sample <- function(n_bins = 50L, lambda = 100, seed = NULL) {
  .make_sample(n_bins = n_bins, lambda = lambda, seed = seed, y_extra = 0.5)
}

# ===========================================================================
# Tests for scale_sample()
# ===========================================================================

expect_silent(
  s <- .make_sample(n_bins = 100L, lambda = 50, seed = 42)
)
s <- getFromNamespace(".as_wcx_sample", "RWisecondorX")(s)

expect_true(is.list(s),
            info = "synthetic sample stays list-like after wrapping")
expect_true(S7::S7_inherits(s, WisecondorXSample),
            info = "sample wrapper uses the WisecondorXSample S7 class")

# Identity scaling
s_same <- scale_sample(s, from_size = 5000, to_size = 5000)
expect_identical(s_same, s, info = "identity scaling returns original")

# 2× scaling: 100 bins of 5kb → 50 bins of 10kb
s2 <- scale_sample(s, from_size = 5000, to_size = 10000)
expect_equal(length(s2[["1"]]), 50L,
             info = "2× scaling halves number of bins")
expect_equal(sum(s2[["1"]]), sum(s[["1"]]),
             info = "2× scaling preserves total counts")

# Check first scaled bin = sum of first two original bins
expect_equal(s2[["1"]][1], s[["1"]][1] + s[["1"]][2],
             info = "first scaled bin = sum of first two source bins")

# 5× scaling: 100 bins of 5kb → 20 bins of 25kb
s5 <- scale_sample(s, from_size = 5000, to_size = 25000)
expect_equal(length(s5[["1"]]), 20L,
             info = "5× scaling: 100 bins → 20 bins")
expect_equal(sum(s5[["1"]]), sum(s[["1"]]),
             info = "5× scaling preserves total counts")

# Impossible scaling
expect_error(scale_sample(s, from_size = 5000, to_size = 3000),
             pattern = "Impossible",
             info = "scale_sample rejects non-multiple target")

expect_error(scale_sample(s, from_size = 5000, to_size = 7000),
             pattern = "Impossible",
             info = "scale_sample rejects non-divisible target")

# ===========================================================================
# Tests for .gender_correct()
# ===========================================================================

s_gc <- .make_sample(n_bins = 10L, lambda = 50, seed = 1)
original_x <- s_gc[["23"]]
original_y <- s_gc[["24"]]

s_male <- .gender_correct(s_gc, "M")
expect_equal(s_male[["23"]], original_x * 2L,
             info = "male gender correction doubles chrX")
expect_equal(s_male[["24"]], original_y * 2L,
             info = "male gender correction doubles chrY")

s_female <- .gender_correct(s_gc, "F")
expect_equal(s_female[["23"]], original_x,
             info = "female gender correction is a no-op for chrX")
expect_equal(s_female[["24"]], original_y,
             info = "female gender correction is a no-op for chrY")

# ===========================================================================
# Tests for .predict_gender()
# ===========================================================================

# Sample with nearly zero Y
s_low_y <- .make_female_sample(n_bins = 50L, seed = 10)
gender_low <- .predict_gender(s_low_y, trained_cutoff = 0.001)
expect_equal(gender_low, "F",
             info = "low Y-fraction sample classified as female")

# Sample with substantial Y
s_high_y <- .make_male_sample(n_bins = 50L, seed = 10)
gender_high <- .predict_gender(s_high_y, trained_cutoff = 0.001)
expect_equal(gender_high, "M",
             info = "high Y-fraction sample classified as male")

# ===========================================================================
# Tests for .get_mask()
# ===========================================================================

# All samples identical — all bins should be masked in (above 5% median)
samples_uniform <- replicate(12, .make_sample(n_bins = 20L, lambda = 100, seed = NULL),
                             simplify = FALSE)
# Hack: make all samples identical for determinism
set.seed(99)
one_sample <- .make_sample(n_bins = 20L, lambda = 100, seed = 99)
samples_same <- replicate(12, one_sample, simplify = FALSE)

mask_info <- .get_mask(samples_same)
expect_true(is.logical(mask_info$mask),
            info = "mask is logical vector")
expect_equal(length(mask_info$bins_per_chr), 24L,
             info = "bins_per_chr has 24 elements")
expect_equal(mask_info$bins_per_chr[1], 20L,
             info = "bins_per_chr matches input bin count")

# With uniform coverage, most bins should be kept (but chr24 has very low coverage)
n_kept <- sum(mask_info$mask)
n_total <- sum(mask_info$bins_per_chr)
expect_true(n_kept > n_total * 0.8,
            info = "uniform-coverage mask keeps most bins")

# ===========================================================================
# Tests for .normalize_and_mask()
# ===========================================================================

set.seed(42)
norm_samples <- replicate(10, .make_sample(n_bins = 20L, lambda = 100), simplify = FALSE)
mask_result <- .get_mask(norm_samples)
nm_data <- .normalize_and_mask(norm_samples, chr_range = 1:22,
                               mask = mask_result$mask[seq_len(sum(mask_result$bins_per_chr[1:22]))])
expect_true(is.matrix(nm_data),
            info = ".normalize_and_mask returns a matrix")
expect_equal(ncol(nm_data), 10L,
             info = ".normalize_and_mask has one column per sample")

# ===========================================================================
# Tests for .train_pca() and .project_pc()
# ===========================================================================

# Create a small matrix for PCA
set.seed(123)
n_bins <- 100L
n_samp <- 15L
pca_data <- matrix(rnorm(n_bins * n_samp, mean = 1e-4, sd = 1e-5),
                   nrow = n_bins, ncol = n_samp)
pca_data <- abs(pca_data)  # ensure positive

pca_result <- .train_pca(pca_data, n_comp = 3L)
expect_true(is.matrix(pca_result$corrected),
            info = ".train_pca returns corrected matrix")
expect_equal(dim(pca_result$corrected), c(n_bins, n_samp),
             info = ".train_pca corrected has same dims as input")
expect_equal(nrow(pca_result$components), 3L,
             info = ".train_pca components has n_comp rows")
expect_equal(ncol(pca_result$components), n_bins,
             info = ".train_pca components has n_bins cols")
expect_equal(length(pca_result$center), n_bins,
             info = ".train_pca center has length n_bins")

# Project a single sample through the PCA
test_vec <- pca_data[, 1]
projected <- .project_pc(test_vec, pca_result$components, pca_result$center)
expect_equal(length(projected), n_bins,
             info = ".project_pc returns vector of same length")
expect_true(all(is.finite(projected)),
            info = ".project_pc returns finite values")

# ===========================================================================
# Tests for .log_trans()
# ===========================================================================

results_r <- list(c(1.5, 2.0, 0, 0.5), c(1.0, 0, 3.0))
results_z <- list(c(1.0, 2.0, 3.0, 4.0), c(5.0, 6.0, 7.0))
results_w <- list(c(1.0, 1.0, 1.0, 1.0), c(1.0, 1.0, 1.0))

lt <- .log_trans(results_r, results_z, results_w, m_lr = 0.1)

# Zeros in ratio should produce zeros everywhere
expect_equal(lt$results_r[[1]][3], 0,
             info = ".log_trans zeroes ratio where input is 0")
expect_equal(lt$results_z[[1]][3], 0,
             info = ".log_trans zeroes z-score where ratio was 0")
expect_equal(lt$results_w[[1]][3], 0,
             info = ".log_trans zeroes weight where ratio was 0")

# Non-zero values should be log2-transformed and centered
expect_equal(lt$results_r[[1]][1], log2(1.5) - 0.1, tolerance = 1e-10,
             info = ".log_trans applies log2 and subtracts m_lr")

# Second chromosome
expect_equal(lt$results_z[[2]][2], 0,
             info = ".log_trans zeroes z-score at ratio=0 on chr2")

# ===========================================================================
# Tests for .inflate_results()
# ===========================================================================

mask <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
vals <- c(10, 20, 30)
inflated <- .inflate_results(vals, mask)
expect_equal(inflated, c(10, 0, 20, 30, 0),
             info = ".inflate_results un-masks correctly")

# ===========================================================================
# Tests for .idx_to_chr_name()
# ===========================================================================

expect_equal(.idx_to_chr_name(1L), "1",
             info = "chr 1 → '1'")
expect_equal(.idx_to_chr_name(22L), "22",
             info = "chr 22 → '22'")
expect_equal(.idx_to_chr_name(23L), "X",
             info = "chr 23 → 'X'")
expect_equal(.idx_to_chr_name(24L), "Y",
             info = "chr 24 → 'Y'")

# ===========================================================================
# Tests for .call_aberrations()
# ===========================================================================

results_c <- data.frame(
  chr = c(1L, 2L, 3L),
  start = c(0L, 0L, 0L),
  end = c(10L, 10L, 10L),
  zscore = c(6.0, -6.0, 2.0),
  ratio = c(0.3, -0.3, 0.01),
  stringsAsFactors = FALSE
)

# Z-score mode (default)
ab_z <- .call_aberrations(results_c, zscore_cutoff = 5, beta = NULL, ref_gender = "F")
expect_equal(nrow(ab_z), 2L,
             info = "two aberrations called with z-cutoff=5")
expect_equal(ab_z$type, c("gain", "loss"),
             info = "gain and loss correctly identified")

# Beta mode
ab_b <- .call_aberrations(results_c, zscore_cutoff = 5, beta = 0.5, ref_gender = "F")
expect_true(nrow(ab_b) >= 1L,
            info = "beta mode calls at least one aberration")

# No aberrations
results_c_clean <- data.frame(
  chr = 1L, start = 0L, end = 10L, zscore = 1.0, ratio = 0.01,
  stringsAsFactors = FALSE
)
ab_clean <- .call_aberrations(results_c_clean, zscore_cutoff = 5, beta = NULL, ref_gender = "F")
expect_equal(nrow(ab_clean), 0L,
             info = "no aberrations when z-score below cutoff")

# ===========================================================================
# Tests for CBS (.exec_cbs) — requires DNAcopy
# ===========================================================================

has_dnacopy <- TRUE

if (has_dnacopy) {
  # Simple CBS test: flat signal with one spike
  set.seed(77)
  n <- 200L
  ratio_flat <- rnorm(n, mean = 0, sd = 0.05)
  # Insert a "gain" in the middle
  ratio_flat[80:120] <- ratio_flat[80:120] + 0.5

  results_r_cbs <- list(ratio_flat)
  results_w_cbs <- list(rep(1, n))

  cbs_out <- .exec_cbs(results_r_cbs, results_w_cbs, ref_gender = "F",
                        alpha = 0.01, binsize = 100000L, seed = 42,
                        parallel = FALSE)

  expect_true(is.data.frame(cbs_out),
              info = ".exec_cbs returns a data frame")
  expect_true(nrow(cbs_out) >= 2L,
              info = ".exec_cbs detects at least 2 segments (gain + baseline)")
  expect_true(all(c("chr", "start", "end", "ratio") %in% names(cbs_out)),
              info = ".exec_cbs returns correct columns")

  # The segment covering the gain region should have a positive ratio
  gain_seg <- cbs_out[cbs_out$ratio > 0.2, ]
  expect_true(nrow(gain_seg) >= 1L,
              info = ".exec_cbs detects the gain segment")
} else {
  expect_message(TRUE, pattern = NA,
                 info = "DNAcopy not available; CBS tests skipped")
}

# ===========================================================================
# Tests for full pipeline: rwisecondorx_newref() + rwisecondorx_predict()
# ===========================================================================

# Build 15 synthetic "female" samples + 5 "male" samples (minimum 10 required)
set.seed(2024)
n_bins <- 30L
n_female <- 12L
n_male <- 8L
all_samples <- vector("list", n_female + n_male)

for (i in seq_len(n_female)) {
  all_samples[[i]] <- .make_female_sample(n_bins = n_bins, lambda = 200)
}
for (i in seq_len(n_male)) {
  all_samples[[n_female + i]] <- .make_male_sample(n_bins = n_bins, lambda = 200)
}

# rwisecondorx_newref() requires mclust for gender model
has_mclust <- TRUE

if (has_mclust && has_dnacopy) {
  # --- Test newref ---
  # Suppress messages during test
  ref <- suppressMessages(
    rwisecondorx_newref(
      samples     = all_samples,
      binsize     = 100000L,
      nipt        = FALSE,
      refsize     = 10L,  # small for testing speed
      cpus        = 1L
    )
  )

  expect_true(is.list(ref),
              info = "rwisecondorx_newref returns a list-like reference object")
  expect_true(S7::S7_inherits(ref, WisecondorXReference),
              info = "rwisecondorx_newref returns the typed WisecondorXReference S7 class")
  expect_equal(ref$binsize, 100000L,
               info = "reference binsize matches input")
  expect_true(ref$trained_cutoff > 0 && ref$trained_cutoff < 1,
              info = "gender cutoff is in (0, 1)")
  expect_identical(sum(ref$masked_bins_per_chr), nrow(ref$indexes),
                   info = "reference rows match masked autosomal bins")
  expect_identical(dim(ref$indexes), dim(ref$distances),
                   info = "reference indexes and distances share the same shape")
  expect_identical(nrow(ref$null_ratios), nrow(ref$indexes),
                   info = "null ratios have one row per target bin")
  expect_identical(ncol(ref$pca_components), length(ref$pca_mean),
                   info = "PCA components and center describe the same feature space")

  qc <- rwisecondorx_ref_qc(ref, min_ref_bins = 5L)
  expect_true(qc$overall_verdict %in% c("PASS", "WARN", "FAIL"),
              info = "QC report has a valid overall verdict")
  expect_identical(sort(names(qc$metrics)), c("F", "M"),
                   info = "non-NIPT QC reports female and male branches")
  expect_true(all(vapply(qc$metrics, function(branch) {
    is.finite(branch$mean_of_means) &&
      is.finite(branch$std_of_means) &&
      branch$n_bins > 0L
  }, logical(1L))),
  info = "QC metrics contain finite branch summaries")

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    qc_json <- tempfile(fileext = ".json")
    qc_json_result <- rwisecondorx_ref_qc(
      ref,
      min_ref_bins = 5L,
      output_json = qc_json
    )
    qc_parsed <- jsonlite::read_json(qc_json, simplifyVector = TRUE)
    expect_identical(qc_parsed$overall_verdict, qc_json_result$overall_verdict,
                     info = "JSON output preserves the QC verdict")
    expect_identical(names(qc_parsed$metrics), names(qc_json_result$metrics),
                     info = "JSON output preserves QC branch labels")
    unlink(qc_json)
  }

  # --- Test predict with a female sample ---
  test_sample_f <- .make_female_sample(n_bins = n_bins, lambda = 200, seed = 999)

  pred <- suppressMessages(suppressWarnings(
    rwisecondorx_predict(
      sample    = test_sample_f,
      reference = ref,
      zscore    = 5,
      minrefbins = 5L,
      alpha     = 0.01,
      seed      = 42,
      parallel  = FALSE
    )
  ))

  expect_true(is.list(pred),
              info = "rwisecondorx_predict returns a list-like prediction object")
  expect_true(S7::S7_inherits(pred, WisecondorXPrediction),
              info = "rwisecondorx_predict returns the typed WisecondorXPrediction S7 class")
  expect_true(pred$gender %in% c("F", "M"),
              info = "predicted gender is F or M")
  expect_true(pred$n_reads > 0,
              info = "n_reads is positive")
  expect_identical(length(pred$results_r), length(pred$bins_per_chr),
                   info = "prediction returns one ratio vector per chromosome")
  expect_identical(vapply(pred$results_r, length, integer(1L)),
                   vapply(pred$results_z, length, integer(1L)),
                   info = "ratio and z-score vectors stay aligned per chromosome")
  expect_true(all(pred$results_c$start <= pred$results_c$end),
              info = "CBS segments have valid start/end ordering")
  expect_true(all(pred$statistics$chr %in% c(as.character(1:22), "X", "Y")),
              info = "statistics rows are chromosome-labelled")
  expect_true(nrow(pred$statistics) >= 22L,
              info = "statistics include chromosome-level summaries")
  expect_true(nrow(pred$aberrations) == 0L ||
              all(pred$aberrations$type %in% c("gain", "loss")),
              info = "aberrations use gain/loss labels only")

  # --- Test output writing ---
  tmpdir <- tempfile("rwx_test_")
  dir.create(tmpdir)
  outprefix <- file.path(tmpdir, "test_output")

  out_path <- write_wisecondorx_output(pred, outprefix)
  expect_equal(out_path, outprefix,
               info = "write_wisecondorx_output returns outprefix")

  bins_file <- paste0(outprefix, "_bins.bed")
  segs_file <- paste0(outprefix, "_segments.bed")
  aber_file <- paste0(outprefix, "_aberrations.bed")
  stat_file <- paste0(outprefix, "_statistics.txt")

  bins_lines <- readLines(bins_file)
  seg_lines <- readLines(segs_file)
  aber_lines <- readLines(aber_file)
  stat_lines <- readLines(stat_file)
  expect_identical(bins_lines[1L], "chr\tstart\tend\tid\tratio\tzscore",
                   info = "bins BED header matches the public format")
  expect_identical(seg_lines[1L], "chr\tstart\tend\tratio\tzscore",
                   info = "segments BED header matches the public format")
  expect_identical(aber_lines[1L], "chr\tstart\tend\tratio\tzscore\ttype",
                   info = "aberrations BED header matches the public format")
  expect_identical(length(bins_lines), sum(vapply(pred$results_r, length, integer(1L))) + 1L,
                   info = "bins BED writes one data row per bin")
  expect_identical(length(seg_lines), nrow(pred$results_c) + 1L,
                   info = "segments BED writes one data row per segment")
  expect_identical(length(aber_lines), nrow(pred$aberrations) + 1L,
                   info = "aberrations BED writes one data row per call")
  expect_true(any(grepl(paste0("^Number of reads: ", pred$n_reads, "$"), stat_lines)),
              info = "statistics file records the prediction read count")
  expect_true(any(grepl(paste0(": ", pred$gender, "$"), stat_lines)),
              info = "statistics file records the predicted gender")

  unlink(tmpdir, recursive = TRUE)

  # --- Test NIPT mode ---
  ref_nipt <- suppressMessages(
    rwisecondorx_newref(
      samples = all_samples,
      binsize = 100000L,
      nipt    = TRUE,
      refsize = 10L,
      cpus    = 1L
    )
  )

  expect_true(isTRUE(ref_nipt$is_nipt),
              info = "NIPT mode reference has is_nipt=TRUE")

  qc_nipt <- rwisecondorx_ref_qc(ref_nipt, min_ref_bins = 5L)
  expect_identical(names(qc_nipt$metrics), "F",
                   info = "NIPT QC follows upstream female-only branch selection")

} else {
  if (!has_mclust) {
    message("mclust not available; full pipeline tests skipped")
  }
  if (!has_dnacopy) {
    message("DNAcopy not available; full pipeline tests skipped")
  }
}

# ===========================================================================
# Tests for scale_sample edge cases
# ===========================================================================

# Non-divisible bin count: 7 bins / 3 = 3 (ceiling)
s_odd <- list("1" = as.integer(1:7))
s_odd_scaled <- scale_sample(s_odd, from_size = 1000, to_size = 3000)
expect_equal(length(s_odd_scaled[["1"]]), 3L,
             info = "non-divisible bin count uses ceiling")
expect_equal(s_odd_scaled[["1"]][1], 1L + 2L + 3L,
             info = "first scaled bin sums correctly")
expect_equal(s_odd_scaled[["1"]][3], 7L,
             info = "last partial bin sums correctly")

# Empty chromosome
s_empty <- list("1" = integer(0), "2" = as.integer(1:10))
s_empty_scaled <- scale_sample(s_empty, from_size = 1000, to_size = 2000)
expect_equal(length(s_empty_scaled[["1"]]), 0L,
             info = "empty chromosome stays empty after scaling")
expect_equal(length(s_empty_scaled[["2"]]), 5L,
             info = "non-empty chromosome scaled correctly")
