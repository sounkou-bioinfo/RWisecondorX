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
                  "cv_fraction", "shapiro_p_value") %in%
                  names(qc$chromosome_summary)),
            info = "chromosome summary exposes the expected CV statistics")
expect_true(all(c("sample_name", "mean_ssd", "max_abs_z",
                  "is_matching_outlier") %in%
                  names(qc$sample_summary)),
            info = "sample summary exposes matching and z-score metrics")
expect_true(all(c("sex_group", "n_samples", "rr_x_cv", "rr_y_cv",
                  "can_build_ncv", "can_build_regression") %in%
                  names(qc$sex_summary)),
            info = "sex summary exposes gaunosome readiness metrics")
expect_true(all(c("chromosome", "bin", "mean_scaled", "sd_scaled",
                  "cv_scaled", "chi_z", "overdispersed",
                  "correction_factor") %in%
                  names(qc$bin_summary)),
            info = "bin summary exposes chi-profile CV and correction metrics")

collapsed_fracs <- RWisecondorX:::.control_group_fractions_collapsed(cg)
manual_cv_21 <- 100 * stats::sd(collapsed_fracs["21", ]) /
  mean(collapsed_fracs["21", ])
qc_cv_21 <- qc$chromosome_summary$cv_fraction[
  qc$chromosome_summary$chromosome == "21"
]
expect_equal(qc_cv_21, manual_cv_21,
             info = "chromosome summary CV matches the manual chr21 fraction CV")

mm <- nipter_match_matrix(cg)
manual_mean_ssd_1 <- mean(mm[1L, -1L])
expect_equal(qc$sample_summary$mean_ssd[[1L]], manual_mean_ssd_1,
             info = "sample summary mean_ssd matches the matching-matrix row mean")
expect_identical(qc$sample_summary$sample_name, qc$sample_names,
                 info = "sample summary rows align with sample_names")

expect_identical(sort(qc$sex_summary$sex_group), c("female", "male"),
                 info = "sex summary reports both female and male groups")
expect_true(all(qc$sex_summary$n_samples == 3L),
            info = "sex summary reports the expected per-sex sample counts")
expect_true(all(qc$sex_summary$can_build_ncv),
            info = "three controls per sex are sufficient for NCV model building")
expect_true(!any(qc$sex_summary$can_build_regression),
            info = "three controls per sex are insufficient for regression model building")

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
manual_chi_bin1 <- sum((manual_mean_bin1 - scaled_bin1)^2 /
                         max(manual_mean_bin1, .Machine$double.eps))
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
expect_equal(qc$bin_summary$chi_z[[1L]], manual_chiz_bin1,
             tolerance = 1e-10,
             info = "bin summary chi_z matches the manual first-bin value")

qc_no_bins <- nipter_control_group_qc(cg, include_bins = FALSE)
expect_true(is.null(qc_no_bins$bin_summary),
            info = "bin_summary is omitted when include_bins = FALSE")
