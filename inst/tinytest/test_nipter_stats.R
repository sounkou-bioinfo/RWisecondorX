library(tinytest)
library(RWisecondorX)

.helper_nipter <- system.file("tinytest", "helper_nipter.R", package = "RWisecondorX")
if (nzchar(.helper_nipter)) {
  sys.source(.helper_nipter, envir = environment())
} else {
  sys.source("inst/tinytest/helper_nipter.R", envir = environment())
}

# ---------------------------------------------------------------------------
# Tests for the NIPTeR statistical layer:
#   nipter_control.R — control group construction & diagnostics
#   nipter_score.R   — Z-score and NCV
#   nipter_chi.R     — chi-squared correction
#   nipter_regression.R — regression Z-score
#
# These tests use synthetic NIPTeRSample objects (no BAM fixture required).
# Simulation helpers come from helper_nipter.R (auto-sourced by tinytest).
# ---------------------------------------------------------------------------

# ===================================================================
# CONTROL GROUP
# ===================================================================

ctrl_samples <- .sim_nipter_control_set(10L)
.sample_names <- function(samples) {
  vapply(samples, `[[`, character(1L), "sample_name")
}

# --- nipter_as_control_group ---

cg <- nipter_as_control_group(ctrl_samples)

expect_identical(length(cg$samples), 10L,
                 info = "control group has 10 samples")
expect_identical(.sample_names(cg$samples), .sample_names(ctrl_samples),
                 info = "control group preserves sample ordering and names")
expect_identical(cg$description, "General control group",
                 info = "default description")
expect_identical(cg$correction_status_autosomal, "Uncorrected",
                 info = "new control groups start uncorrected")

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

expect_identical(dim(diag$z_scores), c(22L, 10L),
                 info = "diagnose returns one Z-score per chromosome and sample")
expect_identical(colnames(diag$statistics),
                 c("mean", "SD", "shapiro_p_value"),
                 info = "statistics has expected columns")
expect_identical(nrow(diag$statistics), 22L,
                 info = "statistics has 22 rows")
expect_true(all(is.finite(diag$statistics[, c("mean", "SD")])),
            info = "diagnose returns finite per-chromosome mean/SD values")

# Aberrant-score output is either empty or only contains |Z| > 3 events.
if (is.null(diag$aberrant_scores)) {
  expect_true(TRUE, info = "aberrant_scores can be NULL for stable synthetic controls")
} else {
  expect_true(all(abs(diag$aberrant_scores$z_score) > 3),
              info = "aberrant_scores only contains true |Z| > 3 events")
}

# --- nipter_match_control_group ---

test_sample <- .sim_nipter_sample(.sim_chr_totals(scale = 1), "test_subject",
                                  seed = 999L)

matched <- nipter_match_control_group(test_sample, cg, n = 5L,
                                      mode = "subset")
expect_identical(length(matched$samples), 5L,
                 info = "matched control group has n=5 samples")

scores <- nipter_match_control_group(test_sample, cg, n = 10L,
                                     mode = "report")
expect_identical(length(scores), 10L,  info = "report returns all controls")
# Scores should be sorted ascending
expect_true(all(diff(scores) >= 0),    info = "scores sorted ascending")
expect_identical(sort(.sample_names(matched$samples)),
                 sort(names(scores)[seq_len(5L)]),
                 info = "subset mode returns the same top controls as report mode")


# ===================================================================
# Z-SCORE
# ===================================================================

z21 <- nipter_z_score(test_sample, cg, chromo_focus = 21L)

expect_identical(z21$focus_chromosome, "21",
                 info = "focus_chromosome is '21'")
expect_identical(names(z21$control_statistics),
                 c("mean", "sd", "shapiro_p_value"),
                 info = "control_statistics has expected names")
expect_identical(length(z21$control_z_scores), 10L,
                 info = "control_z_scores has n_controls entries")
expect_identical(z21$sample_name, "test_subject",
                 info = "sample_name preserved")

# A normal sample should stay reasonably close to the control centre.
expect_true(abs(z21$sample_z_score) < 3,
            info = "normal sample stays within 3 SD of the control centre")

# Test with a trisomy sample — should have elevated Z
trisomy_sample <- .sim_nipter_sample(
  .sim_chr_totals(scale = 1, trisomy_chr = 21, trisomy_frac = 0.10),
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

expect_identical(ncv21$focus_chromosome, "21",
                 info = "NCV focus_chromosome")
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

# Correction status should include "Chi square corrected"
expect_true("Chi square corrected" %in%
            chi_result$sample$correction_status_autosomal,
            info = "sample correction status updated")
expect_true(!"Uncorrected" %in%
            chi_result$sample$correction_status_autosomal,
            info = "'Uncorrected' removed after chi correction")

expect_true(all(vapply(chi_result$control_group$samples, function(s) {
  "Chi square corrected" %in% s$correction_status_autosomal
}, logical(1L))),
info = "chi correction updates every control sample")

# Corrected read counts should still be non-negative
corrected_auto <- chi_result$sample$autosomal_chromosome_reads[[1L]]
expect_true(all(corrected_auto >= 0),
            info = "chi-corrected autosomal counts are non-negative")

# Total reads should stay close (correction reweights noisy bins, not halve coverage)
orig_total <- sum(test_sample$autosomal_chromosome_reads[[1L]])
corr_total <- sum(corrected_auto)
expect_true(abs(corr_total - orig_total) / orig_total < 0.05,
            info = "chi correction keeps total autosomal reads within 5%")

mixed_ctrl <- ctrl_samples[1:2]
mixed_ctrl[[1]] <- RWisecondorX:::.sample_append_correction_step(
  mixed_ctrl[[1]], "autosomal", RWisecondorX:::.nipt_gc_correction_step()
)
mixed_ctrl[[2]] <- RWisecondorX:::.sample_append_correction_step(
  mixed_ctrl[[2]], "autosomal", RWisecondorX:::.nipt_gc_correction_step()
)
mixed_ctrl[[2]] <- RWisecondorX:::.sample_append_correction_step(
  mixed_ctrl[[2]], "autosomal", RWisecondorX:::.nipt_chi_correction_step()
)
mixed_cg <- nipter_as_control_group(mixed_ctrl)
expect_true(all(c("GC Corrected", "Chi square corrected") %in%
                  mixed_cg$correction_status_autosomal),
            info = "control-group correction status returns the union of typed correction steps")


# ===================================================================
# REGRESSION Z-SCORE
# ===================================================================

# Need a larger control group for train/test split
ctrl_large <- .sim_nipter_control_set(20L, seed = 100L)
cg_large   <- nipter_as_control_group(ctrl_large)

reg21 <- nipter_regression(test_sample, cg_large, chromo_focus = 21L,
                           n_models = 2L, n_predictors = 3L,
                           seed = 42L)

expect_identical(reg21$focus_chromosome, "21",
                 info = "regression focus_chromosome")
expect_identical(reg21$sample_name, "test_subject",
                 info = "regression sample_name")

expect_true(length(reg21$models) >= 1L,   info = "at least 1 model")
expect_true(length(reg21$models) <= 2L,   info = "at most n_models models")

m1 <- reg21$models[[1L]]
expect_true(m1$cv_type %in% c("practical", "theoretical"),
            info = "model cv_type valid")
expect_true(length(m1$predictors) >= 1L && length(m1$predictors) <= 3L,
            info = "model has 1-3 predictors")
expect_true(all(is.finite(c(m1$z_score, m1$cv, m1$control_z_scores))),
            info = "model metrics are finite")

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


# ===================================================================
# SEPARATED STRANDS
# ===================================================================

# Helpers from helper_nipter.R (auto-sourced):
#   .sim_nipter_ss_sample(), .sim_nipter_ss_control_set(), .sim_chr_totals()

# --- SeparatedStrands object structure ---

ss_sample <- .sim_nipter_ss_sample(.sim_chr_totals(), "ss_test", seed = 500L)
expect_identical(length(ss_sample$autosomal_chromosome_reads), 2L,
                 info = "SS auto has 2 matrices (fwd, rev)")
expect_identical(nrow(ss_sample$autosomal_chromosome_reads[[1L]]), 22L,
                 info = "SS fwd auto has 22 rows")
expect_identical(nrow(ss_sample$autosomal_chromosome_reads[[2L]]), 22L,
                 info = "SS rev auto has 22 rows")
expect_true(all(grepl("F$", rownames(ss_sample$autosomal_chromosome_reads[[1L]]))),
            info = "SS fwd rows end in F")
expect_true(all(grepl("R$", rownames(ss_sample$autosomal_chromosome_reads[[2L]]))),
            info = "SS rev rows end in R")

# --- SeparatedStrands control group ---

ss_ctrl <- .sim_nipter_ss_control_set(10L)
ss_cg <- nipter_as_control_group(ss_ctrl)

expect_identical(length(ss_cg$samples), 10L,
                 info = "SS control group has 10 samples")
expect_identical(.sample_names(ss_cg$samples), .sample_names(ss_ctrl),
                 info = "SS control group preserves sample ordering and names")

# Mixed strand types should error
expect_error(
  nipter_as_control_group(c(ctrl_samples[1:3], ss_ctrl[1:3])),
  info = "mixed CombinedStrands/SeparatedStrands rejected"
)

# --- SeparatedStrands Z-score ---

ss_z21 <- nipter_z_score(ss_sample, ss_cg, chromo_focus = 21L)
expect_true(abs(ss_z21$sample_z_score) < 5,
            info = "SS normal sample has moderate Z-score")

# Trisomy SS sample should have elevated Z
ss_trisomy <- .sim_nipter_ss_sample(
  .sim_chr_totals(scale = 1, trisomy_chr = 21, trisomy_frac = 0.10),
  "ss_trisomy_21", seed = 777L
)
ss_z21_tri <- nipter_z_score(ss_trisomy, ss_cg, chromo_focus = 21L)
expect_true(ss_z21_tri$sample_z_score > ss_z21$sample_z_score,
            info = "SS trisomy Z > SS normal Z")

# --- SeparatedStrands NCV ---

ss_ncv21 <- nipter_ncv_score(ss_sample, ss_cg, chromo_focus = 21L,
                             max_elements = 3L)
expect_true(length(ss_ncv21$denominators) >= 1L,
            info = "SS NCV chooses at least one denominator chromosome")

ss_ncv21_tri <- nipter_ncv_score(ss_trisomy, ss_cg, chromo_focus = 21L,
                                 max_elements = 3L)
expect_true(ss_ncv21_tri$sample_score > ss_ncv21$sample_score,
            info = "SS trisomy NCV > SS normal NCV")

# --- SeparatedStrands chi-squared correction ---

ss_chi <- nipter_chi_correct(ss_sample, ss_cg, chi_cutoff = 3.5)
expect_identical(length(ss_chi$sample$autosomal_chromosome_reads), 2L,
                 info = "SS chi sample still has 2 auto matrices")
expect_true("Chi square corrected" %in%
            ss_chi$sample$correction_status_autosomal,
            info = "SS chi correction status updated")

# Both strand matrices should have non-negative counts
expect_true(all(ss_chi$sample$autosomal_chromosome_reads[[1L]] >= 0),
            info = "SS chi fwd auto counts non-negative")
expect_true(all(ss_chi$sample$autosomal_chromosome_reads[[2L]] >= 0),
            info = "SS chi rev auto counts non-negative")

# --- SeparatedStrands regression ---

ss_ctrl_large <- .sim_nipter_ss_control_set(20L, seed = 200L)
ss_cg_large   <- nipter_as_control_group(ss_ctrl_large)

ss_reg21 <- nipter_regression(ss_sample, ss_cg_large, chromo_focus = 21L,
                              n_models = 2L, n_predictors = 3L,
                              seed = 42L)

expect_identical(ss_reg21$focus_chromosome, "21",
                 info = "SS regression focus_chromosome")
expect_true(length(ss_reg21$models) >= 1L,
             info = "SS regression has at least 1 model")

ss_m1 <- ss_reg21$models[[1L]]
expect_true(all(is.finite(c(ss_m1$z_score, ss_m1$cv, ss_m1$control_z_scores))),
            info = "SS regression metrics are finite")

# Predictors should be strand-specific (end in F or R)
expect_true(all(grepl("[FR]$", ss_m1$predictors)),
            info = "SS regression predictors are strand-specific")

# Complementary exclusion within model: no two predictors should be the
# same chromosome number (one F and one R)
ss_pred_nums <- sub("[FR]$", "", ss_m1$predictors)
expect_true(length(ss_pred_nums) == length(unique(ss_pred_nums)),
            info = "SS regression complementary exclusion within model")

# No predictor reuse between models
if (length(ss_reg21$models) >= 2L) {
  ss_m2 <- ss_reg21$models[[2L]]
  expect_true(length(intersect(ss_m1$predictors, ss_m2$predictors)) == 0L,
              info = "SS regression no predictor reuse across models")
  # But complementary strands CAN appear across models
  # (e.g., "5F" in model 1, "5R" in model 2 is allowed)
}

# Trisomy SS regression should be elevated
ss_reg21_tri <- nipter_regression(ss_trisomy, ss_cg_large,
                                  chromo_focus = 21L,
                                  n_models = 2L, n_predictors = 3L,
                                  seed = 42L)
expect_true(ss_reg21_tri$models[[1]]$z_score > ss_reg21$models[[1]]$z_score,
            info = "SS trisomy regression Z > SS normal regression Z")

# --- SeparatedStrands diagnose ---

ss_diag <- nipter_diagnose_control_group(ss_cg)
expect_identical(dim(ss_diag$z_scores), c(22L, 10L),
                 info = "SS diagnose collapses to one Z-score matrix")

ss_diag_stranded <- nipter_diagnose_control_group(ss_cg, collapse_strands = FALSE)
expect_identical(dim(ss_diag_stranded$z_scores), c(44L, 10L),
                 info = "SS diagnose can expose strand-resolved diagnostics")
expect_true(all(grepl("^[0-9]+[FR]$", rownames(ss_diag_stranded$z_scores))),
            info = "SS strand-resolved diagnostics use upstream-style row names")

# --- SeparatedStrands match ---

ss_matched <- nipter_match_control_group(ss_sample, ss_cg, n = 5L,
                                         mode = "subset")
expect_identical(length(ss_matched$samples), 5L,
                 info = "SS matched has 5 samples")

# --- Control-group dropping / iterative pruning ---

ss_drop <- nipter_drop_control_group_samples(ss_cg, "ss_ctrl_01")
expect_identical(length(ss_drop$samples), 9L,
                 info = "dropping one control leaves 9 samples")
expect_true(!("ss_ctrl_01" %in% control_names(ss_drop)),
            info = "requested sample is removed from control group")

ss_qc <- data.frame(
  sample_name = control_names(ss_cg),
  TotalUniqueReads = c(6e6, 6.1e6, 3e6, 6.2e6, 6.0e6, 5.9e6, 6.05e6, 6.15e6, 5.95e6, 6.1e6),
  GCPCTAfterFiltering = c(40.0, 40.2, 39.9, 40.1, 39.8, 40.3, 40.0, 40.2, 39.7, 45.0),
  stringsAsFactors = FALSE
)
ss_qc_filtered <- nipter_filter_control_group_qc(
  ss_cg,
  sample_qc = ss_qc,
  min_total_unique_reads = 4.5e6,
  gc_mad_cutoff = 3
)
expect_identical(sort(ss_qc_filtered$excluded_samples),
                 c("ss_ctrl_03", "ss_ctrl_10"),
                 info = "QC filtering removes low-depth and GC-outlier controls")
expect_identical(length(ss_qc_filtered$control_group$samples), 8L,
                 info = "QC filtering returns the retained control subset")

ss_qc_ref <- nipter_build_reference(
  ss_cg,
  sample_qc = ss_qc,
  min_total_unique_reads = 4.5e6,
  gc_mad_cutoff = 3
)
expect_identical(length(ss_qc_ref$control_group$samples), 8L,
                 info = "nipter_build_reference applies QC filtering before fitting")
expect_true(all(!c("ss_ctrl_03", "ss_ctrl_10") %in%
                  control_names(ss_qc_ref$control_group)),
            info = "QC-filtered controls are absent from the built reference")

ss_prune_samples <- .sim_nipter_ss_control_set(10L, seed = 123L)
ss_prune_samples[[1]]@auto_fwd[c("1F", "2F", "3F"), ] <- ss_prune_samples[[1]]@auto_fwd[c("1F", "2F", "3F"), ] * 100
ss_prune_samples[[1]]@auto_rev[c("1R", "2R", "3R"), ] <- ss_prune_samples[[1]]@auto_rev[c("1R", "2R", "3R"), ] * 100
ss_prune_cg <- nipter_as_control_group(ss_prune_samples)
ss_pruned <- nipter_prune_control_group_outliers(
  ss_prune_cg,
  collapse_strands = FALSE,
  z_cutoff = 2.5,
  min_controls = 9L
)
expect_true("ss_ctrl_01" %in% ss_pruned$dropped_samples,
            info = "iterative pruning drops a sample with both-strand aberrant chromosome fractions")
expect_identical(length(ss_pruned$control_group$samples), 9L,
                 info = "iterative pruning removes the flagged outlier sample")
expect_true("Chi square corrected" %in%
              ss_pruned$chi_corrected_control_group$correction_status_autosomal,
            info = "iterative pruning returns the terminal chi-corrected control group")
