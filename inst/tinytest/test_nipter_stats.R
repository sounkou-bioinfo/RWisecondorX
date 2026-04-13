library(tinytest)

# ---------------------------------------------------------------------------
# Tests for the NIPTeR statistical layer:
#   nipter_control.R — control group construction & diagnostics
#   nipter_score.R   — Z-score and NCV
#   nipter_chi.R     — chi-squared correction
#   nipter_regression.R — regression Z-score
#
# These tests use synthetic NIPTeRSample objects (no BAM fixture required).
# ---------------------------------------------------------------------------

.source_candidate <- function(path) {
  candidates <- c(path, file.path("..", "..", path))
  candidates[file.exists(candidates)][1L]
}

src_files <- c("R/convert.R", "R/nipter_bin.R", "R/nipter_control.R",
               "R/nipter_gc.R", "R/nipter_chi.R", "R/nipter_score.R",
               "R/nipter_regression.R")
src_paths <- vapply(src_files, .source_candidate, character(1L))

if (all(!is.na(src_paths))) {
  for (p in src_paths) source(p)
} else if (requireNamespace("RWisecondorX", quietly = TRUE)) {
  for (fn in c("nipter_as_control_group", "nipter_diagnose_control_group",
               "nipter_match_control_group", "nipter_z_score",
               "nipter_ncv_score", "nipter_chi_correct",
               "nipter_regression")) {
    assign(fn, getExportedValue("RWisecondorX", fn))
  }
} else {
  stop("Unable to locate source files or installed RWisecondorX.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Helpers: build synthetic NIPTeRSample objects
# ---------------------------------------------------------------------------

# Create a fake NIPTeRSample with specified per-chromosome total reads.
# chr_totals: named integer vector with keys "1"-"22" (autosomal).
# n_bins: number of bins (reads spread evenly then jittered).
.make_sample <- function(chr_totals, name, n_bins = 100L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  auto_mat <- matrix(0L, nrow = 22L, ncol = n_bins)
  rownames(auto_mat) <- as.character(1:22)
  colnames(auto_mat) <- as.character(seq_len(n_bins))

  for (chr in names(chr_totals)) {
    total <- chr_totals[chr]
    if (total == 0L) next
    # Spread reads roughly evenly across bins with some noise
    base   <- total %/% n_bins
    remainder <- total - base * n_bins
    counts <- rep(base, n_bins)
    if (remainder > 0L) {
      bump <- sample(seq_len(n_bins), remainder)
      counts[bump] <- counts[bump] + 1L
    }
    auto_mat[chr, ] <- as.integer(counts)
  }

  sex_mat <- matrix(0L, nrow = 2L, ncol = n_bins)
  rownames(sex_mat) <- c("X", "Y")
  colnames(sex_mat) <- as.character(seq_len(n_bins))

  structure(
    list(
      autosomal_chromosome_reads  = list(auto_mat),
      sex_chromosome_reads        = list(sex_mat),
      correction_status_autosomal = "Uncorrected",
      correction_status_sex       = "Uncorrected",
      sample_name                 = name
    ),
    class = c("NIPTeRSample", "CombinedStrands")
  )
}

# Generate a "normal" sample: reads roughly proportional to chromosome length.
# Each chromosome gets ~1000 * scale reads per length unit.
.normal_chr_totals <- function(scale = 1, trisomy_chr = NULL,
                               trisomy_frac = 0.05) {
  # Approximate relative chromosome sizes (hg38, Mb, rounded)
  chr_sizes <- c(248, 242, 198, 190, 182, 171, 159, 145, 138, 134, 135, 133,
                 114, 107, 102, 90, 83, 80, 59, 64, 47, 51)
  totals <- as.integer(round(chr_sizes * 10 * scale))
  names(totals) <- as.character(1:22)

  # Add trisomy signal if requested
  if (!is.null(trisomy_chr)) {
    chr_key <- as.character(trisomy_chr)
    extra <- as.integer(round(totals[chr_key] * trisomy_frac))
    totals[chr_key] <- totals[chr_key] + extra
  }

  totals
}

# Build a set of n normal samples as a control group
.make_control_set <- function(n = 10L, scale = 1, seed = 42L) {
  set.seed(seed)
  samples <- vector("list", n)
  for (i in seq_len(n)) {
    # Add small random noise to each sample's totals
    noise <- runif(22, 0.95, 1.05)
    chr_totals <- .normal_chr_totals(scale = scale)
    chr_totals <- as.integer(round(chr_totals * noise))
    names(chr_totals) <- as.character(1:22)
    samples[[i]] <- .make_sample(chr_totals, paste0("ctrl_", i), seed = seed + i)
  }
  samples
}

# ===================================================================
# CONTROL GROUP
# ===================================================================

ctrl_samples <- .make_control_set(10L)

# --- nipter_as_control_group ---

cg <- nipter_as_control_group(ctrl_samples)

expect_true(inherits(cg, "NIPTeRControlGroup"),
            info = "control group has correct class")
expect_true(inherits(cg, "CombinedStrands"),
            info = "control group inherits strand type")
expect_identical(length(cg$samples), 10L,
                 info = "control group has 10 samples")
expect_identical(cg$description, "General control group",
                 info = "default description")
expect_true(is.character(cg$correction_status_autosomal),
            info = "correction_status_autosomal is character")

# Validation: must have at least 2 samples
expect_error(nipter_as_control_group(ctrl_samples[1]),
             info = "rejects single-sample control group")

# Validation: rejects non-NIPTeRSample
expect_error(nipter_as_control_group(list(1, 2)),
             info = "rejects non-NIPTeRSample objects")

# Duplicate removal
dup_samples <- c(ctrl_samples, ctrl_samples[1:2])
cg_dedup <- nipter_as_control_group(dup_samples)
expect_identical(length(cg_dedup$samples), 10L,
                 info = "duplicate samples removed by name")

# --- nipter_diagnose_control_group ---

diag <- nipter_diagnose_control_group(cg)

expect_true(is.matrix(diag$z_scores),   info = "z_scores is a matrix")
expect_identical(nrow(diag$z_scores), 22L,
                 info = "z_scores has 22 rows (chromosomes)")
expect_identical(ncol(diag$z_scores), 10L,
                 info = "z_scores has 10 columns (samples)")

expect_true(is.matrix(diag$statistics), info = "statistics is a matrix")
expect_identical(colnames(diag$statistics),
                 c("mean", "SD", "shapiro_p_value"),
                 info = "statistics has expected columns")
expect_identical(nrow(diag$statistics), 22L,
                 info = "statistics has 22 rows")

# With well-behaved synthetic data, no aberrant scores
# (the noise is small relative to the signal)
# This can be NULL or a data.frame — just check it doesn't error
expect_true(is.null(diag$aberrant_scores) ||
            is.data.frame(diag$aberrant_scores),
            info = "aberrant_scores is NULL or data.frame")

# --- nipter_match_control_group ---

test_sample <- .make_sample(.normal_chr_totals(scale = 1), "test_subject",
                            seed = 999L)

matched <- nipter_match_control_group(test_sample, cg, n = 5L,
                                      mode = "subset")
expect_true(inherits(matched, "NIPTeRControlGroup"),
            info = "matched control group has correct class")
expect_identical(length(matched$samples), 5L,
                 info = "matched control group has n=5 samples")

scores <- nipter_match_control_group(test_sample, cg, n = 10L,
                                     mode = "report")
expect_true(is.numeric(scores),        info = "report mode returns numeric")
expect_identical(length(scores), 10L,  info = "report returns all controls")
# Scores should be sorted ascending
expect_true(all(diff(scores) >= 0),    info = "scores sorted ascending")


# ===================================================================
# Z-SCORE
# ===================================================================

z21 <- nipter_z_score(test_sample, cg, chromo_focus = 21L)

expect_true(inherits(z21, "NIPTeRZScore"), info = "z_score result class")
expect_true(is.numeric(z21$sample_z_score),
            info = "sample_z_score is numeric")
expect_identical(z21$focus_chromosome, "21",
                 info = "focus_chromosome is '21'")
expect_identical(names(z21$control_statistics),
                 c("mean", "sd", "shapiro_p_value"),
                 info = "control_statistics has expected names")
expect_identical(length(z21$control_z_scores), 10L,
                 info = "control_z_scores has n_controls entries")
expect_identical(z21$sample_name, "test_subject",
                 info = "sample_name preserved")

# A normal sample should have Z ~ 0 (within reason for small control group)
expect_true(abs(z21$sample_z_score) < 5,
            info = "normal sample has moderate Z-score")

# Test with a trisomy sample — should have elevated Z
trisomy_sample <- .make_sample(
  .normal_chr_totals(scale = 1, trisomy_chr = 21, trisomy_frac = 0.10),
  "trisomy_21", seed = 777L
)

z21_tri <- nipter_z_score(trisomy_sample, cg, chromo_focus = 21L)
expect_true(z21_tri$sample_z_score > z21$sample_z_score,
            info = "trisomy sample has higher Z-score than normal")
expect_true(z21_tri$sample_z_score > 2,
            info = "trisomy sample Z-score > 2")

# ===================================================================
# NCV SCORE
# ===================================================================

ncv21 <- nipter_ncv_score(test_sample, cg, chromo_focus = 21L,
                          max_elements = 3L)

expect_true(inherits(ncv21, "NIPTeRNCV"), info = "NCV result class")
expect_true(is.numeric(ncv21$sample_score),
            info = "sample_score is numeric")
expect_identical(ncv21$focus_chromosome, "21",
                 info = "NCV focus_chromosome")
expect_true(is.integer(ncv21$denominators),
            info = "denominators is integer vector")
expect_true(length(ncv21$denominators) >= 1L &&
            length(ncv21$denominators) <= 3L,
            info = "denominators between 1 and max_elements")
expect_true(is.numeric(ncv21$best_cv) && ncv21$best_cv >= 0,
            info = "best_cv is non-negative numeric")
expect_identical(length(ncv21$control_z_scores), 10L,
                 info = "NCV control_z_scores has n_controls entries")

# Trisomy sample NCV should be elevated
ncv21_tri <- nipter_ncv_score(trisomy_sample, cg, chromo_focus = 21L,
                              max_elements = 3L)
expect_true(ncv21_tri$sample_score > ncv21$sample_score,
            info = "trisomy NCV score > normal NCV score")


# ===================================================================
# CHI-SQUARED CORRECTION
# ===================================================================

chi_result <- nipter_chi_correct(test_sample, cg, chi_cutoff = 3.5)

expect_true(is.list(chi_result),          info = "chi_correct returns list")
expect_true(inherits(chi_result$sample, "NIPTeRSample"),
            info = "corrected sample is NIPTeRSample")
expect_true(inherits(chi_result$control_group, "NIPTeRControlGroup"),
            info = "corrected control_group is NIPTeRControlGroup")

# Correction status should include "Chi square corrected"
expect_true("Chi square corrected" %in%
            chi_result$sample$correction_status_autosomal,
            info = "sample correction status updated")
expect_true(!"Uncorrected" %in%
            chi_result$sample$correction_status_autosomal,
            info = "'Uncorrected' removed after chi correction")

# All control samples should also be corrected
for (i in seq_along(chi_result$control_group$samples)) {
  s <- chi_result$control_group$samples[[i]]
  expect_true("Chi square corrected" %in% s$correction_status_autosomal,
              info = paste("control", i, "correction status updated"))
}

# Corrected read counts should still be non-negative
corrected_auto <- chi_result$sample$autosomal_chromosome_reads[[1L]]
expect_true(all(corrected_auto >= 0),
            info = "chi-corrected autosomal counts are non-negative")

# Total reads should be similar (correction redistributes, not removes)
orig_total <- sum(test_sample$autosomal_chromosome_reads[[1L]])
corr_total <- sum(corrected_auto)
expect_true(abs(corr_total - orig_total) / orig_total < 0.5,
            info = "chi correction doesn't drastically change total reads")


# ===================================================================
# REGRESSION Z-SCORE
# ===================================================================

# Need a larger control group for train/test split
ctrl_large <- .make_control_set(20L, seed = 100L)
cg_large   <- nipter_as_control_group(ctrl_large)

reg21 <- nipter_regression(test_sample, cg_large, chromo_focus = 21L,
                           n_models = 2L, n_predictors = 3L,
                           seed = 42L)

expect_true(inherits(reg21, "NIPTeRRegression"),
            info = "regression result class")
expect_identical(reg21$focus_chromosome, "21",
                 info = "regression focus_chromosome")
expect_identical(reg21$sample_name, "test_subject",
                 info = "regression sample_name")

expect_true(is.list(reg21$models),        info = "models is a list")
expect_true(length(reg21$models) >= 1L,   info = "at least 1 model")
expect_true(length(reg21$models) <= 2L,   info = "at most n_models models")

m1 <- reg21$models[[1L]]
expect_true(is.numeric(m1$z_score),       info = "model z_score is numeric")
expect_true(is.numeric(m1$cv),            info = "model cv is numeric")
expect_true(m1$cv_type %in% c("practical", "theoretical"),
            info = "model cv_type valid")
expect_true(is.character(m1$predictors),  info = "model predictors is character")
expect_true(length(m1$predictors) >= 1L && length(m1$predictors) <= 3L,
            info = "model has 1-3 predictors")
expect_true(is.numeric(m1$control_z_scores),
            info = "model control_z_scores is numeric")

# No predictor reuse between models
if (length(reg21$models) >= 2L) {
  m2 <- reg21$models[[2L]]
  expect_true(length(intersect(m1$predictors, m2$predictors)) == 0L,
              info = "no predictor reuse between models")
}

# Trisomy regression Z should be elevated
reg21_tri <- nipter_regression(trisomy_sample, cg_large, chromo_focus = 21L,
                               n_models = 2L, n_predictors = 3L,
                               seed = 42L)
expect_true(reg21_tri$models[[1]]$z_score > reg21$models[[1]]$z_score,
            info = "trisomy regression Z > normal regression Z")

# Control group too small should error
expect_error(
  nipter_regression(test_sample,
                    nipter_as_control_group(ctrl_samples[1:3]),
                    chromo_focus = 21L, seed = 42L),
  info = "regression rejects tiny control group"
)
