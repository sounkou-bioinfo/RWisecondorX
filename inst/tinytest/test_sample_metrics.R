library(tinytest)
library(RWisecondorX)

.helper_nipter <- system.file("tinytest", "helper_nipter.R", package = "RWisecondorX")
if (nzchar(.helper_nipter)) {
  sys.source(.helper_nipter, envir = environment())
} else {
  sys.source("inst/tinytest/helper_nipter.R", envir = environment())
}

base_row <- RWisecondorX:::.sample_qc_row(
  sample_name = "s1",
  bam = "/tmp/s1.bam"
)
expect_true(
  "read_counts_input_total" %in% names(base_row),
  info = "sample metrics rows expose the canonical read-count columns"
)
expect_true(
  is.na(base_row$read_counts_input_total),
  info = "read-count columns stay NA when read-count metrics are unavailable"
)

counts_row <- RWisecondorX:::.sample_qc_row(
  sample_name = "large",
  bam = "/tmp/large.bam",
  read_counts = list(
    input_total = 7e6,
    mapped_total = 6.2e6,
    unmapped_total = 8e5,
    count_pre = 6e6,
    count_post = 5.5e6,
    retained_fraction = 5.5 / 6
  )
)
expect_true(
  isTRUE(all(c(
    "read_counts_input_total",
    "read_counts_mapped_total",
    "read_counts_unmapped_total",
    "read_counts_binned_pre_sum",
    "read_counts_binned_post_sum",
    "read_counts_binned_retained_fraction"
  ) %in% names(counts_row))),
  info = "sample metrics rows expose the semantic read-count fields"
)

qc_row <- RWisecondorX:::.sample_qc_row(
  sample_name = "qc",
  bam = "/tmp/qc.bam",
  read_counts = list(count_pre = 6e6, count_post = 5.5e6),
  coverage = list(
    binsize = 50000,
    pre_mapq = 0,
    pre_exclude_flags = 1796,
    pre_summary_txt = "/tmp/qc.coverage_pre.summary.txt",
    pre_length_bases = 3.1e9,
    pre_bases = 6.9e8,
    pre_covered_fraction = 0.22,
    pre_mean_depth = 0.22,
    pre_min_depth = 0,
    pre_max_depth = 153,
    post_mapq = 40,
    post_exclude_flags = 1796,
    post_summary_txt = "/tmp/qc.coverage_post.summary.txt",
    post_length_bases = 3.1e9,
    post_bases = 5.3e8,
    post_covered_fraction = 0.17,
    post_mean_depth = 0.17,
    post_min_depth = 0,
    post_max_depth = 65,
    source = "native_mosdepth_dense"
  ),
  nipter_qc = list(
    gc_loess_valid_bins = 1200,
    nipter_autosomal_bin_cv_pre_gc_correction = 10.5,
    nipter_autosomal_bin_cv_post_gc_correction = 8.2,
    nipter_gc_correction_bin_scale_mean = 1.01,
    nipter_gc_correction_bin_scale_sd = 0.09,
    nipter_gc_correction_bin_scale_cv = 8.91,
    nipter_gc_correlation_pre_gc_correction = 0.42,
    nipter_gc_correlation_post_gc_correction = 0.03,
    gc_curve_plot = "/tmp/qc.gc_curve.png"
  )
)
expect_true(
  all(c(
    "sample_processing_status",
    "sample_failure_stage",
    "sample_failure_step",
    "sample_failure_message",
    "coverage_binsize",
    "coverage_pre_summary_txt",
    "coverage_post_summary_txt",
    "gc_loess_valid_bins",
    "nipter_autosomal_bin_cv_pre_gc_correction",
    "nipter_autosomal_bin_cv_post_gc_correction",
    "nipter_gc_correction_bin_scale_mean",
    "nipter_gc_correction_bin_scale_sd",
    "nipter_gc_correction_bin_scale_cv",
    "nipter_gc_correlation_pre_gc_correction",
    "nipter_gc_correlation_post_gc_correction",
    "nipter_gc_curve_png",
    "nipter_corrected_bin_ratio_genome_plot_png"
  ) %in% names(qc_row)),
  info = "sample metrics rows expose the NIPTeR preprocessing QC columns"
)
expect_identical(
  qc_row$nipter_gc_curve_png,
  "/tmp/qc.gc_curve.png",
  info = "sample metrics rows preserve the per-sample GC curve path"
)

stats_template <- as.list(setNames(
  rep(0, length(RWisecondorX:::.samtools_stats_summary_fields)),
  RWisecondorX:::.samtools_stats_summary_fields
))
stats_template$raw_total_sequences <- 7e6
stats_template$filtered_sequences <- 0
stats_template$sequences <- 7e6
stats_template$reads_mapped <- 6.2e6
stats_template$reads_unmapped <- 8e5
stats_template$total_length <- 105e8
stats_template$bases_mapped <- 93e8
stats_template$bases_mapped_cigar <- 92e8
stats_template$mismatches_from_nm <- 1.2e6
stats_template$error_rate <- 1.2e6 / 92e8
stats_template$duplicated_read_fraction <- 0.08
stats_template$duplicated_base_fraction <- 0.08

bam_stats <- RWisecondorX:::.samtools_stats_metric_record(
  pre = data.frame(
    c(
      list(
        min_mapq = 0,
        require_flags = 0,
        exclude_flags = 0,
        report_filtered_stream_bookkeeping = FALSE
      ),
      stats_template
    ),
    check.names = FALSE
  ),
  post = data.frame(
    c(
      list(
        min_mapq = 40,
        require_flags = 0,
        exclude_flags = 1024,
        report_filtered_stream_bookkeeping = TRUE
      ),
      stats_template
    ),
    check.names = FALSE
  ),
  source = "native_htslib_samtools_stats"
)

complete_row <- RWisecondorX:::.sample_qc_row(
  sample_name = "complete",
  bam = "/tmp/complete.bam",
  seqff = list(seqff_pre = 0.10, enet_pre = 0.09, wrsc_pre = 0.11,
               seqff_post = 0.09, enet_post = 0.08, wrsc_post = 0.10,
               source = "computed"),
  y_unique = list(ratio_pre = 0.001, reads_pre = 10, total_pre = 10000,
                  ratio_post = 0.0008, reads_post = 8, total_post = 10000,
                  regions_file = "/tmp/y_regions.bed", source = "computed"),
  read_counts = list(
    input_total = 7e6,
    mapped_total = 6.2e6,
    unmapped_total = 8e5,
    count_pre = 6e6,
    count_post = 5.5e6,
    count_fwd = 2.8e6,
    count_rev = 2.7e6,
    n_nonzero_bins_post = 1200,
    retained_fraction = 5.5 / 6,
    gc_pre = 40.1,
    gc_post = 40.0,
    mean_mapq_post = 58.2,
    source = "samtools_idxstats+bam_bin_counts(gc,mq)"
  ),
  bam_stats = bam_stats,
  coverage = list(
    binsize = 50000,
    pre_mapq = 0,
    pre_exclude_flags = 1796,
    pre_summary_txt = "/tmp/complete.coverage_pre.summary.txt",
    pre_length_bases = 3101804739,
    pre_bases = 692370794,
    pre_covered_fraction = 0.223216,
    pre_mean_depth = 0.22,
    pre_min_depth = 0,
    pre_max_depth = 153,
    post_mapq = 40,
    post_exclude_flags = 1796,
    post_summary_txt = "/tmp/complete.coverage_post.summary.txt",
    post_length_bases = 3101804739,
    post_bases = 527658506,
    post_covered_fraction = 0.170112,
    post_mean_depth = 0.17,
    post_min_depth = 0,
    post_max_depth = 65,
    source = "native_mosdepth_dense"
  ),
  nipter_qc = list(
    gc_loess_valid_bins = 1200,
    gc_loess_total_bins = 1300,
    gc_loess_invalid_bins = 100,
    gc_loess_valid_bin_fraction_pct = 92.3,
    gc_loess_invalid_out_of_chromosome_range = 12,
    gc_loess_invalid_gc_missing_or_nonfinite = 5,
    gc_loess_invalid_gc_nonpositive = 0,
    gc_loess_invalid_raw_nonfinite = 0,
    gc_loess_invalid_raw_nonpositive = 95,
    gc_loess_invalid_corrected_nonfinite = 0,
    gc_curve_has_valid_bins = TRUE,
    gc_curve_has_loess_support = TRUE,
    nipter_autosomal_bin_cv_pre_gc_correction = 10.5,
    nipter_autosomal_bin_cv_post_gc_correction = 8.2,
    nipter_gc_correction_bin_scale_mean = 1.01,
    nipter_gc_correction_bin_scale_sd = 0.09,
    nipter_gc_correction_bin_scale_cv = 8.91,
    nipter_gc_correlation_pre_gc_correction = 0.42,
    nipter_gc_correlation_post_gc_correction = 0.03,
    gc_curve_plot = "/tmp/qc.gc_curve.png",
    gc_curve_data_bgz = "/tmp/qc.gc_curve.tsv.bgz",
    corrected_bin_ratio_genome_plot = "/tmp/qc.corrected_ratio.png"
  ),
  filters_pre = list(mapq = 0, require_flags = 0, exclude_flags = 0),
  filters_post = list(mapq = 40, require_flags = 0, exclude_flags = 1024)
)

expect_true(
  all(c(
    "samtools_stats_source",
    "samtools_stats_pre_bookkeeping_mode",
    "samtools_stats_pre_raw_total_sequences",
    "samtools_stats_post_bookkeeping_mode",
    "samtools_stats_post_error_rate"
  ) %in% names(complete_row)),
  info = "sample metrics rows expose the pre/post samtools-stats-compatible fields"
)
expect_identical(
  complete_row$samtools_stats_post_bookkeeping_mode,
  "filtered_stream",
  info = "sample metrics rows preserve the post-filter bookkeeping mode"
)

qc_gc_tbl <- stats::setNames(
  rep(list(seq(0.35, 0.55, length.out = 20L)), 24L),
  c(as.character(1:22), "X", "Y")
)
qc_sample <- .sim_nipter_sample(.sim_chr_totals(scale = 2), "qc_sex_plot", n_bins = 20L)
qc_sex <- sex_matrix(qc_sample)
qc_sex["X", ] <- 5
qc_sex["Y", ] <- 2
qc_sample <- RWisecondorX:::.sample_with_reads(qc_sample, sex = list(qc_sex))
qc_corrected <- RWisecondorX:::.sample_with_reads(
  qc_sample,
  autosomal = list(autosomal_matrix(qc_sample) * 1.05),
  sex = list(sex_matrix(qc_sample) * 1.10)
)
qc_preprocess <- RWisecondorX:::.nipter_preprocess_qc(
  sample = qc_sample,
  corrected = qc_corrected,
  gc_table = qc_gc_tbl,
  include_sex = TRUE,
  binsize = 50000L
)
expect_true(
  all(c("X", "Y") %in% qc_preprocess$curve_data$chrom),
  info = "preprocess QC includes X/Y rows when sex GC correction is enabled"
)
expect_true(
  all(!qc_preprocess$curve_data$valid_for_gc_correction_fit[qc_preprocess$curve_data$chrom %in% c("X", "Y")]),
  info = "sex rows remain excluded from the autosomal GC-fit set"
)
expect_true(
  all(is.finite(
    qc_preprocess$curve_data$corrected_ratio_to_fit_eligible_median[
      qc_preprocess$curve_data$chrom %in% c("X", "Y") &
        qc_preprocess$curve_data$valid_for_corrected_ratio_genome_plot
    ]
  )),
  info = "sex rows receive corrected-ratio values for the genome-order plot"
)

processed_ok_row <- RWisecondorX:::.sample_qc_row_with_processing_status(
  complete_row,
  status = "ok"
)
expect_identical(
  processed_ok_row$sample_processing_status,
  "ok",
  info = "sample metrics rows can be annotated with an explicit successful processing state"
)

processed_failed_row <- RWisecondorX:::.sample_qc_row_with_processing_status(
  base_row,
  status = "failed",
  failure_stage = "seqff",
  failure_step = "compute",
  failure_message = "SeqFF failed"
)
expect_identical(
  processed_failed_row$sample_failure_stage,
  "seqff",
  info = "failed processing rows preserve the failing stage"
)
expect_identical(
  processed_failed_row$sample_failure_message,
  "SeqFF failed",
  info = "failed processing rows preserve the failure message"
)

ok_row <- RWisecondorX:::.sample_qc_row_with_nipter_status(
  processed_ok_row,
  status = "not_run"
)
expect_identical(
  ok_row$nipter_bed_status,
  "not_run",
  info = "successful non-NIPTeR sample rows can explicitly mark the NIPTeR artifact as not run"
)

nipter_ok_row <- RWisecondorX:::.sample_qc_row_with_nipter_status(
  processed_ok_row,
  written = TRUE,
  status = "ok"
)
expect_identical(
  nipter_ok_row$nipter_bed_status,
  "ok",
  info = "sample metrics rows mark successful NIPTeR BED writes explicitly"
)
expect_true(
  isTRUE(nipter_ok_row$nipter_bed_written),
  info = "successful NIPTeR BED writes record a TRUE written flag"
)
expect_true(
  is.na(nipter_ok_row$nipter_bed_error),
  info = "successful NIPTeR BED writes clear any error text"
)

failed_row <- RWisecondorX:::.sample_qc_row_with_nipter_status(
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

legacy_row <- RWisecondorX:::.sample_qc_row_with_nipter_status(base_row)
expect_true(
  all(c("nipter_bed_status", "nipter_bed_written", "nipter_bed_error") %in% names(legacy_row)),
  info = "legacy sample metrics rows can be normalized to include explicit NIPTeR BED status columns"
)

expect_error(
  RWisecondorX:::.sample_qc_row_with_nipter_status(
    RWisecondorX:::.sample_qc_row_with_processing_status(base_row, status = "ok"),
    written = TRUE,
    status = "ok"
  ),
  info = "successful sample rows without mandatory preprocessing metrics are rejected"
)
