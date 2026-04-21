library(tinytest)
library(RWisecondorX)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  exit_file("ggplot2 not available; skipping NIPTeR plotting tests")
}

.helper_nipter <- system.file("tinytest", "helper_nipter.R", package = "RWisecondorX")
if (nzchar(.helper_nipter)) {
  sys.source(.helper_nipter, envir = environment())
} else {
  sys.source("inst/tinytest/helper_nipter.R", envir = environment())
}

plot_samples <- .sim_nipter_control_set(10L, seed = 77L)
plot_sex <- c(
  ctrl_01 = "female",
  ctrl_02 = "female",
  ctrl_03 = "female",
  ctrl_04 = "female",
  ctrl_05 = "female",
  ctrl_06 = "male",
  ctrl_07 = "male",
  ctrl_08 = "male",
  ctrl_09 = "male",
  ctrl_10 = "male"
)

for (nm in names(plot_sex)) {
  idx <- match(nm, vapply(plot_samples, function(s) s@sample_name, character(1L)))
  sx <- plot_samples[[idx]]@sex_matrix_
  if (plot_sex[[nm]] == "female") {
    sx["X", ] <- 25L
    sx["Y", ] <- 1L
  } else {
    sx["X", ] <- 12L
    sx["Y", ] <- 5L
  }
  plot_samples[[idx]]@sex_matrix_ <- sx
}

plot_cg <- nipter_as_control_group(
  plot_samples,
  description = "Plot controls",
  sample_sex = plot_sex,
  sex_source = "synthetic_truth"
)

plot_ref <- nipter_build_reference(
  plot_cg,
  y_unique_ratios = c(
    ctrl_01 = 1e-05,
    ctrl_02 = 1.2e-05,
    ctrl_03 = 1.1e-05,
    ctrl_04 = 0.9e-05,
    ctrl_05 = 1.05e-05,
    ctrl_06 = 1.5e-04,
    ctrl_07 = 1.4e-04,
    ctrl_08 = 1.6e-04,
    ctrl_09 = 1.55e-04,
    ctrl_10 = 1.45e-04
  )
)
plot_ref <- nipter_build_gaunosome_models(
  plot_ref,
  ncv_min_elements = 1L,
  ncv_max_elements = 2L,
  regression_n_predictors = 2L,
  regression_extra_predictors = character()
)
plot_qc <- nipter_control_group_qc(
  plot_cg,
  include_bins = TRUE,
  reference_model = plot_ref
)

p_chr <- nipter_plot_qc_chromosomes(plot_qc)
expect_true(inherits(p_chr, "ggplot"),
            info = "chromosome QC plot returns a ggplot object")
expect_false(all(is.na(p_chr$data$chromosome)),
             info = "chromosome QC plot keeps chromosome labels after strand-aware ordering")
expect_true(all(c("z_fraction", "ncv", "rbz") %in% p_chr$data$method_family),
            info = "chromosome QC plot keeps the three method families after merging in sex-model rows")
expect_true(all(c("X_XX", "Y_XX", "X_XY", "Y_XY") %in% as.character(p_chr$data$chromosome)),
            info = "chromosome QC plot includes the sex-model XX/XY X/Y rows")
chr_y_scale <- p_chr$scales$get_scales("y")
expect_equal(chr_y_scale$breaks(c(0, 0.13)), c(0, 0.05, 0.10, 0.15),
             info = "chromosome QC plot uses 0.05 y-axis increments from zero")
expect_equal(chr_y_scale$limits(c(0.02, 0.13)), c(0, 0.15),
             info = "chromosome QC plot starts each y-axis at zero")
big_chr_breaks <- RWisecondorX:::.cv_axis_major_breaks(c(0, 5), step = 0.05)
expect_true(length(big_chr_breaks) < length(RWisecondorX:::.cv_axis_breaks(c(0, 5), step = 0.05)),
            info = "chromosome QC plot uses coarser labeled breaks when the CV range is large")

p_samples <- nipter_plot_qc_samples(plot_qc)
expect_true(inherits(p_samples, "ggplot"),
            info = "sample QC plot returns a ggplot object")

p_bins_cv <- nipter_plot_qc_bins(plot_qc, metric = "cv_scaled")
expect_true(inherits(p_bins_cv, "ggplot"),
            info = "bin CV plot returns a ggplot object")

p_bins_chi <- nipter_plot_qc_bins(plot_qc, metric = "chi_z")
expect_true(inherits(p_bins_chi, "ggplot"),
            info = "bin chi plot returns a ggplot object")

p_sex_box <- nipter_plot_reference_sex_boxplots(plot_ref)
expect_true(inherits(p_sex_box, "ggplot"),
            info = "sex boxplot panel returns a ggplot object")
expect_true(all(c("Z_X_XX", "Z_X_XY", "Z_Y_XX", "Z_Y_XY") %in% levels(p_sex_box$data$metric)),
            info = "sex boxplot panel includes the sex-cluster Z metrics")
expect_true(all(c("-log10(FrChrReads_X)", "-log10(FrChrReads_Y)") %in% levels(p_sex_box$data$metric)),
            info = "sex boxplot panel uses -log10 scales for fraction metrics")

p_sex_scatter_fraction <- nipter_plot_reference_sex_scatter(plot_ref, space = "fraction")
expect_true(inherits(p_sex_scatter_fraction, "ggplot"),
            info = "sex fraction scatter returns a ggplot object")
expect_identical(p_sex_scatter_fraction$labels$x, "-log10(FrChrReads_Y)",
                 info = "sex fraction scatter keeps Y on the x-axis")
expect_identical(p_sex_scatter_fraction$labels$y, "-log10(FrChrReads_X)",
                 info = "sex fraction scatter keeps X on the y-axis")

p_sex_scatter_ratio <- nipter_plot_reference_sex_scatter(plot_ref, space = "ratio")
expect_true(inherits(p_sex_scatter_ratio, "ggplot"),
            info = "sex ratio scatter returns a ggplot object")
expect_identical(p_sex_scatter_ratio$labels$x, "RR_Y",
                 info = "sex ratio scatter keeps Y on the x-axis")
expect_identical(p_sex_scatter_ratio$labels$y, "RR_X",
                 info = "sex ratio scatter keeps X on the y-axis")

plot_prefix <- file.path(tempdir(), "nipter_qc_plot_test")
plot_paths <- write_nipter_reference_plots(plot_qc, plot_ref, plot_prefix, dpi = 72L)
expect_true(all(file.exists(unname(plot_paths))),
            info = "write_nipter_reference_plots writes all expected PNG files")
expect_true(all(grepl("\\.png$", unname(plot_paths))),
            info = "written reference plot outputs are PNG files")
