library(tinytest)
library(RWisecondorX)

.helper_nipter <- system.file("tinytest", "helper_nipter.R", package = "RWisecondorX")
if (nzchar(.helper_nipter)) {
  sys.source(.helper_nipter, envir = environment())
} else {
  sys.source("inst/tinytest/helper_nipter.R", envir = environment())
}

ctrl_samples <- .sim_nipter_control_set(6L, seed = 21L)
sex_labels <- c(
  ctrl_01 = "female",
  ctrl_02 = "female",
  ctrl_03 = "female",
  ctrl_04 = "male",
  ctrl_05 = "male",
  ctrl_06 = "male"
)

cg <- nipter_as_control_group(
  ctrl_samples,
  description = "QC controls",
  sample_sex = sex_labels,
  sex_source = "synthetic_truth"
)

qc <- nipter_control_group_qc(cg, include_bins = TRUE)

expect_true(S7::S7_inherits(qc, NIPTControlGroupQC),
            info = "control-group QC returns the typed S7 QC object")
expect_identical(qc$sample_names, control_names(cg),
                 info = "QC report preserves control sample names")
expect_true(all(c("chromosome", "mean_fraction", "sd_fraction",
                  "cv_fraction", "shapiro_p_value",
                  "method_family", "method_label", "model_status",
                  "n_reference_samples", "n_training_samples", "n_stat_samples") %in%
                  names(qc$chromosome_summary)),
            info = "chromosome summary exposes the expected CV statistics")
expect_identical(qc$settings$outlier_rule, "any_aberrant_score",
                 info = "control-group QC defaults to the upstream any-aberrant-score outlier rule")
expect_true(all(c("sample_name", "mean_ssd", "max_abs_z",
                  "is_matching_outlier",
                  "n_aberrant_rows", "n_aberrant_chromosomes",
                  "has_bidirectional_aberration", "is_chromosomal_outlier",
                  "consensus_gender", "is_reference_sex_outlier",
                  "z_x_xx", "z_x_xy", "z_y_xx", "z_y_xy") %in%
                  names(qc$sample_summary)),
            info = "sample summary exposes matching and z-score metrics")
expect_true(all(c("sex_group", "n_samples", "rr_x_cv", "rr_y_cv",
                  "x_model_cv", "y_model_cv",
                  "can_build_ncv", "can_build_regression") %in%
                  names(qc$sex_summary)),
            info = "sex summary exposes gaunosome readiness metrics")
expect_true(all(c("reference_sex", "focus_chromosome", "method_family",
                  "method_label", "mean_value", "sd_value", "cv_percent",
                  "n_reference_samples", "n_training_samples", "n_stat_samples",
                  "model_status") %in% names(qc$sex_model_summary)),
            info = "sex-model summary exposes XX/XY z, NCV, and regression metrics")
expect_true("n_non_outlier_samples" %in% names(qc$sex_summary),
            info = "sex summary reports non-outlier subset sizes")
expect_true(all(c("chromosome", "bin", "n_nonzero_scaled", "mean_scaled", "sd_scaled",
                  "valid_for_chi", "invalid_reason", "invalid_reason_detail",
                  "in_chromosome_range", "expected_finite", "expected_positive", "scaled_finite",
                  "cv_scaled", "chi_z", "overdispersed",
                  "correction_factor") %in%
                  names(qc$bin_summary)),
            info = "bin summary exposes chi-profile CV and correction metrics")
expect_true(is.logical(qc$bin_summary$valid_for_chi),
            info = "bin summary exposes a logical chi-validity mask")
expect_true(all(qc$bin_summary$invalid_reason == "valid"),
            info = "clean synthetic controls label every bin as chi-valid")

collapsed_fracs <- RWisecondorX:::.control_group_fractions_collapsed(
  nipter_chi_correct(cg$samples[[1L]], cg, chi_cutoff = 3.5)[["control_group"]]
)
manual_cv_21 <- 100 * stats::sd(collapsed_fracs["21", ]) /
  mean(collapsed_fracs["21", ])
qc_cv_21 <- qc$chromosome_summary$cv_fraction[
  qc$chromosome_summary$chromosome == "21" &
    qc$chromosome_summary$method_family == "z_fraction" &
    qc$chromosome_summary$method_label == "combined"
]
expect_equal(qc_cv_21, manual_cv_21,
             info = "chromosome summary CV matches the manual chr21 fraction CV")

expect_true(all(c("z_fraction", "ncv", "rbz") %in% qc$chromosome_summary$method_family),
            info = "chromosome summary includes post-chi fraction, NCV, and RBZ rows")
expect_true(all(c("combined", "ncv", paste0("predictor_set_", 1:4)) %in%
                  qc$chromosome_summary$method_label),
            info = "chromosome summary includes combined, NCV, and all four RBZ predictor-set labels")
rbz_chr21 <- qc$chromosome_summary[
  qc$chromosome_summary$chromosome == "21" &
    qc$chromosome_summary$method_family == "rbz" &
    qc$chromosome_summary$method_label == "predictor_set_1",
  ,
  drop = FALSE
]
expect_identical(rbz_chr21$n_reference_samples, 6L,
                 info = "RBZ chromosome summary keeps the full control-group size as the reference count")
expect_identical(rbz_chr21$n_training_samples, 4L,
                 info = "RBZ chromosome summary records the regression training-set size separately")
expect_identical(rbz_chr21$n_stat_samples, 2L,
                 info = "RBZ chromosome summary records the held-out statistic-set size separately")

mm <- nipter_match_matrix(cg)
manual_mean_ssd_1 <- mean(mm[1L, -1L])
expect_equal(qc$sample_summary$mean_ssd[[1L]], manual_mean_ssd_1,
             info = "sample summary mean_ssd matches the matching-matrix row mean")
expect_identical(qc$sample_summary$sample_name, qc$sample_names,
                 info = "sample summary rows align with sample_names")
expect_true(!any(qc$sample_summary$is_chromosomal_outlier),
            info = "clean synthetic controls do not trigger the chromosomal prune rule")
expect_true(all(qc$sample_summary$n_aberrant_rows == 0L),
            info = "clean synthetic controls have no post-chi aberrant chromosome rows")
expect_true(all(qc$sample_summary$consensus_gender %in% c("female", "male")),
            info = "sample summary carries consensus gender labels when sample sex is supplied")

expect_identical(sort(qc$sex_summary$sex_group), c("female", "male"),
                 info = "sex summary reports both female and male groups")
expect_true(all(qc$sex_summary$n_samples == 3L),
            info = "sex summary reports the expected per-sex sample counts")
expect_true(all(qc$sex_summary$can_build_ncv),
            info = "three controls per sex are sufficient for NCV model building")
expect_true(!any(qc$sex_summary$can_build_regression),
            info = "three controls per sex are insufficient for regression model building")
expect_true(all(qc$sex_summary$n_non_outlier_samples == 3L),
            info = "clean cohort keeps all sex-stratified controls after outlier filtering")
expect_true(all(c("z_fraction", "ncv", "regression") %in% qc$sex_model_summary$method_family),
            info = "sex-model summary includes fraction, NCV, and regression families")
expect_true(all(c("female", "male") %in% qc$sex_model_summary$reference_sex),
            info = "sex-model summary includes both sex-stratified reference groups")
expect_true(all(c("X", "Y") %in% qc$sex_model_summary$focus_chromosome),
            info = "sex-model summary includes both sex chromosomes")
expect_true(all(qc$sex_model_summary$n_reference_samples == 3L),
            info = "sex-model summary reports the raw same-sex reference counts")
expect_true(all(qc$sex_model_summary$n_training_samples == 3L |
                  is.na(qc$sex_model_summary$n_training_samples)),
            info = "sex-model summary records the filtered same-sex model size")
expect_true(all(qc$sex_model_summary$model_status[
  qc$sex_model_summary$method_family == "z_fraction"
] == "ok"),
info = "sex-matched fraction rows are available from the reference frame")
expect_true(all(qc$sex_model_summary$model_status[
  qc$sex_model_summary$method_family == "ncv"
] == "missing"),
info = "sex-NCV rows stay explicit when no reference-side NCV models are supplied")
expect_true(all(qc$sex_model_summary$model_status[
  qc$sex_model_summary$method_family == "regression"
] == "missing"),
info = "sex-regression rows stay explicit when the reference omits regression models")
expect_true(all(qc$sex_model_summary$n_stat_samples[
  qc$sex_model_summary$method_family == "z_fraction"
] == 3L),
info = "sex-model summary keeps all three filtered controls for sex-matched z statistics")

profile_n_bins <- ncol(ctrl_samples[[1L]]$autosomal_chromosome_reads[[1L]])
expect_identical(nrow(qc$bin_summary), 22L * profile_n_bins,
                 info = "bin summary has one row per autosomal chromosome bin")
expect_identical(qc$bin_summary$chromosome[[1L]], "1",
                 info = "bin summary starts with chromosome 1")
expect_identical(qc$bin_summary$bin[[1L]], 1L,
                 info = "bin summary starts with bin 1")

ctrl_flat <- lapply(cg$samples, function(s) as.numeric(t(autosomal_matrix(s))))
ctrl_totals <- vapply(ctrl_flat, sum, numeric(1L))
overall_mean <- mean(ctrl_totals)
scaling <- overall_mean / ctrl_totals
scaled_bin1 <- vapply(seq_along(ctrl_flat), function(i) {
  ctrl_flat[[i]][1L] * scaling[[i]]
}, numeric(1L))
manual_mean_bin1 <- mean(scaled_bin1)
manual_sd_bin1 <- stats::sd(scaled_bin1)
manual_cv_bin1 <- 100 * manual_sd_bin1 / manual_mean_bin1
manual_n_nonzero_bin1 <- sum(scaled_bin1 != 0)
manual_chi_bin1 <- sum((manual_mean_bin1 - scaled_bin1)^2 / manual_mean_bin1)
manual_chiz_bin1 <- (manual_chi_bin1 - (length(ctrl_flat) - 1L)) /
  sqrt(2 * (length(ctrl_flat) - 1L))

expect_equal(qc$bin_summary$mean_scaled[[1L]], manual_mean_bin1,
             tolerance = 1e-10,
             info = "bin summary mean_scaled matches the manual first-bin value")
expect_equal(qc$bin_summary$sd_scaled[[1L]], manual_sd_bin1,
             tolerance = 1e-10,
             info = "bin summary sd_scaled matches the manual first-bin value")
expect_equal(qc$bin_summary$cv_scaled[[1L]], manual_cv_bin1,
             tolerance = 1e-10,
             info = "bin summary cv_scaled matches the manual first-bin value")
expect_identical(qc$bin_summary$n_nonzero_scaled[[1L]], manual_n_nonzero_bin1,
                 info = "bin summary n_nonzero_scaled matches the manual first-bin value")
expect_equal(qc$bin_summary$chi_z[[1L]], manual_chiz_bin1,
             tolerance = 1e-10,
             info = "bin summary chi_z matches the manual first-bin value")
expect_true(all(qc$bin_summary$valid_for_chi),
            info = "clean synthetic controls have chi-valid bins throughout")

qc_no_bins <- nipter_control_group_qc(cg, include_bins = FALSE)
expect_true(is.null(qc_no_bins$bin_summary),
            info = "bin_summary is omitted when include_bins = FALSE")

padded_samples <- .sim_nipter_control_set(4L, seed = 88L, n_bins = 40L)
for (i in seq_along(padded_samples)) {
  padded_samples[[i]]@chrom_lengths["21"] <- 10L * padded_samples[[i]]@binsize
}
padded_cg <- nipter_as_control_group(
  padded_samples,
  description = "Padded controls"
)
padded_qc <- nipter_control_group_qc(padded_cg, include_bins = TRUE)
padded_21_tail <- padded_qc$bin_summary[
  padded_qc$bin_summary$chromosome == "21" &
    padded_qc$bin_summary$bin > 10L,
  ,
  drop = FALSE
]
expect_true(nrow(padded_21_tail) > 0L,
            info = "shortened chr21 produces a padded tail in the bin summary")
expect_true(all(!padded_21_tail$valid_for_chi),
            info = "padded out-of-range chr21 bins are excluded from chi")
expect_true(all(padded_21_tail$invalid_reason == "out_of_range"),
            info = "padded chr21 bins are classified as out_of_range")

outlier_samples <- .sim_nipter_control_set(9L, seed = 33L)
outlier_labels <- c(
  ctrl_01 = "female",
  ctrl_02 = "female",
  ctrl_03 = "female",
  ctrl_04 = "female",
  ctrl_05 = "male",
  ctrl_06 = "male",
  ctrl_07 = "male",
  ctrl_08 = "male",
  ctrl_09 = "male"
)

for (nm in names(outlier_labels)) {
  idx <- match(nm, vapply(outlier_samples, function(s) s@sample_name, character(1L)))
  sx <- outlier_samples[[idx]]@sex_matrix_
  if (outlier_labels[[nm]] == "female") {
    sx["X", ] <- 25L
    sx["Y", ] <- 0L
  } else {
    sx["X", ] <- 12L
    sx["Y", ] <- 4L
  }
  outlier_samples[[idx]]@sex_matrix_ <- sx
}

# One extreme female should be dropped by the same outlier logic used by the
# gaunosome reference builders, leaving only three usable female controls.
sx <- outlier_samples[[1L]]@sex_matrix_
sx["X", ] <- 2L
sx["Y", ] <- 20L
outlier_samples[[1L]]@sex_matrix_ <- sx

outlier_cg <- nipter_as_control_group(
  outlier_samples,
  description = "QC outlier controls",
  sample_sex = outlier_labels,
  sex_source = "synthetic_truth"
)

qc_outlier <- nipter_control_group_qc(outlier_cg)
female_row <- qc_outlier$sex_summary[qc_outlier$sex_summary$sex_group == "female", , drop = FALSE]
male_row <- qc_outlier$sex_summary[qc_outlier$sex_summary$sex_group == "male", , drop = FALSE]

expect_identical(female_row$n_samples, 4L,
                 info = "outlier cohort reports raw female control count")
expect_identical(female_row$n_non_outlier_samples, 3L,
                 info = "outlier cohort reports filtered female control count")
expect_true(!female_row$can_build_regression,
            info = "regression readiness follows filtered female subset size")
expect_true(male_row$can_build_regression,
            info = "regression readiness remains true for unaffected male subset")
expect_true(qc_outlier$sample_summary$is_reference_sex_outlier[
  qc_outlier$sample_summary$sample_name == "ctrl_01"
],
info = "extreme female control is flagged as a reference sex outlier")
female_model_rows <- qc_outlier$sex_model_summary[
  qc_outlier$sex_model_summary$reference_sex == "female",
  ,
  drop = FALSE
]
male_model_rows <- qc_outlier$sex_model_summary[
  qc_outlier$sex_model_summary$reference_sex == "male",
  ,
  drop = FALSE
]
expect_true(all(female_model_rows$n_reference_samples == 4L),
            info = "sex-model summary preserves the raw female reference count before outlier dropping")
expect_true(all(female_model_rows$n_training_samples == 3L | is.na(female_model_rows$n_training_samples)),
            info = "sex-model summary reports the filtered female model size after outlier dropping")
expect_true(all(male_model_rows$n_reference_samples == 5L),
            info = "sex-model summary preserves the raw male reference count")
expect_true(all(male_model_rows$n_training_samples == 4L | is.na(male_model_rows$n_training_samples)),
            info = "sex-model summary reports the filtered male model size after sex-outlier removal")
