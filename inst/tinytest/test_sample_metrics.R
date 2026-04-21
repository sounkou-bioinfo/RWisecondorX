library(tinytest)
library(RWisecondorX)

base_row <- RWisecondorX:::.sample_metrics_row(
  sample_name = "s1",
  bam = "/tmp/s1.bam"
)
expect_true(
  "too_few_reads_for_metrics" %in% names(base_row),
  info = "sample metrics rows expose the low-read guard column"
)
expect_true(
  is.na(base_row$too_few_reads_for_metrics),
  info = "low-read guard stays NA when native read counts are unavailable"
)

tiny_row <- RWisecondorX:::.sample_metrics_row(
  sample_name = "tiny",
  bam = "/tmp/tiny.bam",
  native = list(count_pre = 2, count_post = 1)
)
expect_true(
  isTRUE(tiny_row$too_few_reads_for_metrics),
  info = "low-read guard marks degenerate BAM-like samples"
)

large_row <- RWisecondorX:::.sample_metrics_row(
  sample_name = "large",
  bam = "/tmp/large.bam",
  native = list(count_pre = 6e6, count_post = 5.5e6)
)
expect_true(
  identical(large_row$too_few_reads_for_metrics, FALSE),
  info = "low-read guard stays FALSE for ordinary NIPT-scale samples"
)

qc_row <- RWisecondorX:::.sample_metrics_row(
  sample_name = "qc",
  bam = "/tmp/qc.bam",
  native = list(count_pre = 6e6, count_post = 5.5e6),
  nipter_qc = list(
    gc_loess_valid_bins = 1200,
    nipter_autosomal_bin_cv_pre = 10.5,
    nipter_autosomal_bin_cv_post = 8.2,
    nipter_corrected_bin_ratio_mean = 1.01,
    nipter_corrected_bin_ratio_sd = 0.09,
    nipter_corrected_bin_ratio_cv = 8.91,
    nipter_gc_correlation_pre = 0.42,
    nipter_gc_correlation_post = 0.03,
    gc_curve_plot = "/tmp/qc.gc_curve.png"
  )
)
expect_true(
  all(c(
    "gc_loess_valid_bins",
    "nipter_autosomal_bin_cv_pre",
    "nipter_autosomal_bin_cv_post",
    "nipter_corrected_bin_ratio_mean",
    "nipter_corrected_bin_ratio_sd",
    "nipter_corrected_bin_ratio_cv",
    "nipter_gc_correlation_pre",
    "nipter_gc_correlation_post",
    "nipter_gc_curve_png"
  ) %in% names(qc_row)),
  info = "sample metrics rows expose the NIPTeR preprocessing QC columns"
)
expect_identical(
  qc_row$nipter_gc_curve_png,
  "/tmp/qc.gc_curve.png",
  info = "sample metrics rows preserve the per-sample GC curve path"
)

ok_row <- RWisecondorX:::.sample_metrics_with_nipter_status(
  base_row,
  written = TRUE
)
expect_identical(
  ok_row$nipter_bed_status,
  "ok",
  info = "sample metrics rows mark successful NIPTeR BED writes explicitly"
)
expect_true(
  isTRUE(ok_row$nipter_bed_written),
  info = "successful NIPTeR BED writes record a TRUE written flag"
)
expect_true(
  is.na(ok_row$nipter_bed_error),
  info = "successful NIPTeR BED writes clear any error text"
)

failed_row <- RWisecondorX:::.sample_metrics_with_nipter_status(
  base_row,
  written = FALSE,
  error = "Too few valid autosomal bins"
)
expect_identical(
  failed_row$nipter_bed_status,
  "failed",
  info = "sample metrics rows mark failed NIPTeR BED writes explicitly"
)
expect_true(
  identical(failed_row$nipter_bed_written, FALSE),
  info = "failed NIPTeR BED writes record a FALSE written flag"
)
expect_identical(
  failed_row$nipter_bed_error,
  "Too few valid autosomal bins",
  info = "failed NIPTeR BED writes preserve the failure reason"
)

legacy_row <- RWisecondorX:::.sample_metrics_with_nipter_status(base_row)
expect_true(
  all(c("nipter_bed_status", "nipter_bed_written", "nipter_bed_error") %in% names(legacy_row)),
  info = "legacy sample metrics rows can be normalized to include explicit NIPTeR BED status columns"
)
