# Helpers for per-sample SeqFF / Y-unique metrics used by the preprocessing
# and BED-writing CLIs. These stay internal; the public surface remains the
# existing BED writers plus generic `metadata =`.

.null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

.is_finite_scalar <- function(x) {
  if (is.null(x) || !length(x)) {
    return(FALSE)
  }
  val <- suppressWarnings(as.numeric(x[[1L]]))
  length(val) == 1L && is.finite(val)
}

.sample_qc_value_present <- function(x) {
  if (is.null(x) || !length(x)) {
    return(FALSE)
  }
  value <- x[[1L]]
  if (is.logical(value)) {
    return(!is.na(value))
  }
  if (is.numeric(value) || is.integer(value)) {
    return(length(value) == 1L && is.finite(value))
  }
  if (isTRUE(is.na(value))) {
    return(FALSE)
  }
  txt <- trimws(as.character(value))
  nzchar(txt) && !identical(txt, "NA")
}

.samtools_stats_summary_fields <- c(
  "raw_total_sequences",
  "filtered_sequences",
  "sequences",
  "first_fragments",
  "last_fragments",
  "reads_mapped",
  "reads_mapped_and_paired",
  "reads_unmapped",
  "reads_properly_paired",
  "reads_paired",
  "reads_duplicated",
  "reads_mq0",
  "reads_qc_failed",
  "non_primary_alignments",
  "supplementary_alignments",
  "pairs_on_different_chromosomes",
  "total_length",
  "total_first_fragment_length",
  "total_last_fragment_length",
  "bases_mapped",
  "bases_mapped_cigar",
  "bases_duplicated",
  "mismatches_from_nm",
  "error_rate",
  "duplicated_read_fraction",
  "duplicated_base_fraction"
)

.samtools_stats_summary_labels <- c(
  raw_total_sequences = "raw total sequences",
  filtered_sequences = "filtered sequences",
  sequences = "sequences",
  first_fragments = "1st fragments",
  last_fragments = "last fragments",
  reads_mapped = "reads mapped",
  reads_mapped_and_paired = "reads mapped and paired",
  reads_unmapped = "reads unmapped",
  reads_properly_paired = "reads properly paired",
  reads_paired = "reads paired",
  reads_duplicated = "reads duplicated",
  reads_mq0 = "reads MQ0",
  reads_qc_failed = "reads QC failed",
  non_primary_alignments = "non-primary alignments",
  supplementary_alignments = "supplementary alignments",
  pairs_on_different_chromosomes = "pairs on different chromosomes",
  total_length = "total length",
  total_first_fragment_length = "total first fragment length",
  total_last_fragment_length = "total last fragment length",
  bases_mapped = "bases mapped",
  bases_mapped_cigar = "bases mapped (cigar)",
  bases_duplicated = "bases duplicated",
  mismatches_from_nm = "mismatches from NM",
  error_rate = "error rate",
  duplicated_read_fraction = "duplicated read fraction",
  duplicated_base_fraction = "duplicated base fraction"
)

.sample_qc_column_dictionary <- local({
  samtools_stats_columns <- c(
    "samtools_stats_source",
    "samtools_stats_pre_bookkeeping_mode",
    paste0("samtools_stats_pre_", .samtools_stats_summary_fields),
    "samtools_stats_post_bookkeeping_mode",
    paste0("samtools_stats_post_", .samtools_stats_summary_fields)
  )

  column_name <- c(
    "sample_name", "bam",
    "sample_processing_status", "sample_failure_stage", "sample_failure_step", "sample_failure_message",
    "read_counts_input_total", "read_counts_mapped_total", "read_counts_unmapped_total",
    "read_counts_binned_pre_sum", "read_counts_binned_post_sum",
    "read_counts_binned_fwd_sum", "read_counts_binned_rev_sum",
    "read_counts_binned_nonzero_bins_post",
    "read_counts_binned_retained_fraction",
    "gc_read_perc_pre", "gc_read_perc_post", "mean_mapq_post",
    "coverage_binsize",
    "coverage_pre_mapq", "coverage_pre_exclude_flags",
    "coverage_pre_summary_txt", "coverage_pre_length_bases", "coverage_pre_bases",
    "coverage_pre_covered_fraction", "coverage_pre_mean_depth",
    "coverage_pre_min_depth", "coverage_pre_max_depth",
    "coverage_post_mapq", "coverage_post_exclude_flags",
    "coverage_post_summary_txt", "coverage_post_length_bases", "coverage_post_bases",
    "coverage_post_covered_fraction", "coverage_post_mean_depth",
    "coverage_post_min_depth", "coverage_post_max_depth",
    "coverage_source",
    "nipter_gc_correction_applied",
    "nipter_gc_correction_method",
    "nipter_gc_correction_loess_span",
    "nipter_gc_correction_include_sex",
    "nipter_gc_correction_binsize",
    "nipter_gc_correction_table_bgz",
    "nipter_gc_correction_fasta",
    "gc_loess_valid_bins", "gc_loess_total_bins", "gc_loess_invalid_bins",
    "gc_loess_valid_bin_fraction_pct",
    "gc_loess_invalid_out_of_chromosome_range",
    "gc_loess_invalid_gc_missing_or_nonfinite",
    "gc_loess_invalid_gc_nonpositive",
    "gc_loess_invalid_raw_nonfinite",
    "gc_loess_invalid_raw_nonpositive",
    "gc_loess_invalid_corrected_nonfinite",
    "gc_curve_has_valid_bins", "gc_curve_has_loess_support",
    "nipter_autosomal_bin_cv_pre_gc_correction",
    "nipter_autosomal_bin_cv_post_gc_correction",
    "nipter_gc_correction_bin_scale_mean",
    "nipter_gc_correction_bin_scale_sd",
    "nipter_gc_correction_bin_scale_cv",
    "nipter_gc_correlation_pre_gc_correction",
    "nipter_gc_correlation_post_gc_correction",
    "nipter_gc_curve_png", "nipter_gc_curve_data_bgz",
    "nipter_corrected_bin_ratio_genome_plot_png",
    "seqff_pre", "enet_pre", "wrsc_pre", "fetal_fraction_pre",
    "seqff_post", "enet_post", "wrsc_post", "fetal_fraction_post",
    "y_unique_ratio_pre", "y_unique_reads_pre", "total_nuclear_reads_pre",
    "y_unique_ratio_post", "y_unique_reads_post", "total_nuclear_reads_post",
    "y_unique_regions_file",
    samtools_stats_columns,
    "metrics_pre_mapq", "metrics_pre_require_flags", "metrics_pre_exclude_flags",
    "metrics_post_mapq", "metrics_post_require_flags", "metrics_post_exclude_flags",
    "seqff_source", "y_unique_source", "read_counts_source",
    "nipter_bed_status", "nipter_bed_written", "nipter_bed_error"
  )

  layer <- vapply(column_name, function(x) {
    if (x %in% c("sample_name", "bam")) {
      "identity"
    } else if (x %in% c(
      "sample_processing_status",
      "sample_failure_stage",
      "sample_failure_step",
      "sample_failure_message"
    )) {
      "sample_processing"
    } else if (x %in% c(
      "read_counts_input_total", "read_counts_mapped_total", "read_counts_unmapped_total",
      "read_counts_binned_pre_sum", "read_counts_binned_post_sum",
      "read_counts_binned_fwd_sum", "read_counts_binned_rev_sum",
      "read_counts_binned_nonzero_bins_post",
      "read_counts_binned_retained_fraction", "gc_read_perc_pre", "gc_read_perc_post",
      "mean_mapq_post", "read_counts_source"
    )) {
      "read_counts"
    } else if (x %in% c(
      "coverage_binsize",
      "coverage_pre_mapq", "coverage_pre_exclude_flags",
      "coverage_pre_summary_txt", "coverage_pre_length_bases", "coverage_pre_bases",
      "coverage_pre_covered_fraction", "coverage_pre_mean_depth",
      "coverage_pre_min_depth", "coverage_pre_max_depth",
      "coverage_post_mapq", "coverage_post_exclude_flags",
      "coverage_post_summary_txt", "coverage_post_length_bases", "coverage_post_bases",
      "coverage_post_covered_fraction", "coverage_post_mean_depth",
      "coverage_post_min_depth", "coverage_post_max_depth",
      "coverage_source"
    )) {
      "coverage"
    } else if (x %in% c(
      "nipter_gc_correction_applied",
      "nipter_gc_correction_method",
      "nipter_gc_correction_loess_span",
      "nipter_gc_correction_include_sex",
      "nipter_gc_correction_binsize",
      "nipter_gc_correction_table_bgz",
      "nipter_gc_correction_fasta"
    )) {
      "nipter_gc_correction"
    } else if (x %in% c(
      "seqff_pre", "enet_pre", "wrsc_pre", "fetal_fraction_pre",
      "seqff_post", "enet_post", "wrsc_post", "fetal_fraction_post",
      "seqff_source"
    )) {
      "fetal_fraction"
    } else if (x %in% c(
      "y_unique_ratio_pre", "y_unique_reads_pre", "total_nuclear_reads_pre",
      "y_unique_ratio_post", "y_unique_reads_post", "total_nuclear_reads_post",
      "y_unique_regions_file", "y_unique_source"
    )) {
      "y_unique"
    } else if (x %in% samtools_stats_columns) {
      "samtools_stats"
    } else if (x %in% c(
      "gc_loess_valid_bins", "gc_loess_total_bins", "gc_loess_invalid_bins",
      "gc_loess_valid_bin_fraction_pct",
      "gc_loess_invalid_out_of_chromosome_range",
      "gc_loess_invalid_gc_missing_or_nonfinite",
      "gc_loess_invalid_gc_nonpositive",
      "gc_loess_invalid_raw_nonfinite",
      "gc_loess_invalid_raw_nonpositive",
      "gc_loess_invalid_corrected_nonfinite",
      "gc_curve_has_valid_bins", "gc_curve_has_loess_support",
      "nipter_autosomal_bin_cv_pre_gc_correction",
      "nipter_autosomal_bin_cv_post_gc_correction",
      "nipter_gc_correction_bin_scale_mean",
      "nipter_gc_correction_bin_scale_sd",
      "nipter_gc_correction_bin_scale_cv",
      "nipter_gc_correlation_pre_gc_correction",
      "nipter_gc_correlation_post_gc_correction",
      "nipter_gc_curve_png", "nipter_gc_curve_data_bgz",
      "nipter_corrected_bin_ratio_genome_plot_png"
    )) {
      "nipter_preprocess_qc"
    } else if (startsWith(x, "metrics_")) {
      "metric_filter_contract"
    } else if (startsWith(x, "nipter_bed_")) {
      "nipter_write_status"
    } else {
      "other"
    }
  }, character(1L))

  required_on_success <- !(column_name %in% c(
    "sample_failure_stage",
    "sample_failure_step",
    "sample_failure_message",
    "nipter_gc_correction_loess_span",
    "nipter_gc_correction_binsize",
    "nipter_gc_correction_table_bgz",
    "nipter_gc_correction_fasta",
    "nipter_gc_curve_png",
    "nipter_corrected_bin_ratio_genome_plot_png",
    "nipter_bed_error"
  ))

  description <- vapply(column_name, function(x) {
    if (identical(x, "samtools_stats_source")) {
      return("Provenance for samtools-stats-compatible alignment summary metrics.")
    }
    if (identical(x, "samtools_stats_pre_bookkeeping_mode")) {
      return("Bookkeeping mode used for the pre-filter samtools-stats-compatible summary, for example input_file or filtered_stream.")
    }
    if (identical(x, "samtools_stats_post_bookkeeping_mode")) {
      return("Bookkeeping mode used for the post-filter samtools-stats-compatible summary, for example input_file or filtered_stream.")
    }
    if (startsWith(x, "samtools_stats_pre_")) {
      field <- sub("^samtools_stats_pre_", "", x)
      return(paste0(
        "Pre-filter samtools-stats-compatible summary field '",
        .samtools_stats_summary_labels[[field]],
        "'."
      ))
    }
    if (startsWith(x, "samtools_stats_post_")) {
      field <- sub("^samtools_stats_post_", "", x)
      return(paste0(
        "Post-filter samtools-stats-compatible summary field '",
        .samtools_stats_summary_labels[[field]],
        "'."
      ))
    }
    switch(
      x,
      sample_name = "Canonical sample identifier used across preprocessing and prediction outputs.",
      bam = "Input BAM or CRAM path used to derive this QC row.",
      sample_processing_status = "Final sample-level preprocessing status, for example ok or failed.",
      sample_failure_stage = "Processing stage that produced the terminal sample failure when sample_processing_status is failed.",
      sample_failure_step = "More specific processing step within sample_failure_stage that produced the terminal sample failure.",
      sample_failure_message = "Recorded terminal failure message for the sample when sample_processing_status is failed.",
      read_counts_input_total = "Total reads reported by samtools idxstats-compatible alignment accounting.",
      read_counts_mapped_total = "Total mapped reads reported by samtools idxstats-compatible alignment accounting.",
      read_counts_unmapped_total = "Total unmapped reads reported by samtools idxstats-compatible alignment accounting.",
      read_counts_binned_pre_sum = "Sum of pre-filter per-bin counts from bam_bin_counts().",
      read_counts_binned_post_sum = "Sum of post-filter per-bin counts retained after filtering and deduplication logic.",
      read_counts_binned_fwd_sum = "Sum of forward-strand post-filter per-bin counts.",
      read_counts_binned_rev_sum = "Sum of reverse-strand post-filter per-bin counts.",
      read_counts_binned_nonzero_bins_post = "Number of bins with non-zero post-filter counts.",
      read_counts_binned_retained_fraction = "Post-filter binned count sum divided by the pre-filter binned count sum.",
      gc_read_perc_pre = "Weighted GC percentage from pre-filter native bin counts.",
      gc_read_perc_post = "Weighted GC percentage from post-filter native bin counts.",
      mean_mapq_post = "Weighted mean MAPQ across post-filter native bin counts.",
      coverage_binsize = "Bin size in base pairs used by the cohort coverage summaries.",
      coverage_pre_mapq = "MAPQ threshold used for the pre-filter coverage summary.",
      coverage_pre_exclude_flags = "Excluded SAM flag mask used for the pre-filter coverage summary.",
      coverage_pre_summary_txt = "Path to the pre-filter mosdepth-compatible summary file.",
      coverage_pre_length_bases = "Reference length in bases summarized by the pre-filter coverage pass.",
      coverage_pre_bases = "Summed covered bases reported by the pre-filter coverage pass.",
      coverage_pre_covered_fraction = "Pre-filter covered-base fraction reported as coverage_pre_bases divided by coverage_pre_length_bases.",
      coverage_pre_mean_depth = "Mean depth reported by the pre-filter coverage pass.",
      coverage_pre_min_depth = "Minimum per-region depth reported by the pre-filter coverage pass.",
      coverage_pre_max_depth = "Maximum per-region depth reported by the pre-filter coverage pass.",
      coverage_post_mapq = "MAPQ threshold used for the post-filter coverage summary.",
      coverage_post_exclude_flags = "Excluded SAM flag mask used for the post-filter coverage summary.",
      coverage_post_summary_txt = "Path to the post-filter mosdepth-compatible summary file.",
      coverage_post_length_bases = "Reference length in bases summarized by the post-filter coverage pass.",
      coverage_post_bases = "Summed covered bases reported by the post-filter coverage pass.",
      coverage_post_covered_fraction = "Post-filter covered-base fraction reported as coverage_post_bases divided by coverage_post_length_bases.",
      coverage_post_mean_depth = "Mean depth reported by the post-filter coverage pass.",
      coverage_post_min_depth = "Minimum per-region depth reported by the post-filter coverage pass.",
      coverage_post_max_depth = "Maximum per-region depth reported by the post-filter coverage pass.",
      coverage_source = "Provenance for current coverage summary metrics.",
      nipter_gc_correction_applied = "Whether NIPTeR GC correction was applied before BED export and downstream QC summarization.",
      nipter_gc_correction_method = "GC correction method used for the NIPTeR sample, or 'none' when no correction was applied.",
      nipter_gc_correction_loess_span = "LOESS span used for NIPTeR GC correction when the method is loess.",
      nipter_gc_correction_include_sex = "Whether NIPTeR GC correction was also applied to X and Y bins.",
      nipter_gc_correction_binsize = "Bin size used for the NIPTeR GC correction pass.",
      nipter_gc_correction_table_bgz = "Path to the precomputed GC table used for NIPTeR GC correction, when applicable.",
      nipter_gc_correction_fasta = "Reference FASTA path used directly for NIPTeR GC correction or to derive the referenced GC table.",
      gc_loess_valid_bins = "Number of autosomal bins retained as valid for GC LOESS fitting.",
      gc_loess_total_bins = "Total autosomal bins evaluated for GC LOESS fitting.",
      gc_loess_invalid_bins = "Number of autosomal bins excluded from GC LOESS fitting.",
      gc_loess_valid_bin_fraction_pct = "Percentage of autosomal bins retained as valid for GC LOESS fitting.",
      gc_loess_invalid_out_of_chromosome_range = "Invalid-bin count due to autosomal dense-grid bins extending past the end of a chromosome.",
      gc_loess_invalid_gc_missing_or_nonfinite = "Invalid-bin count due to missing or non-finite GC percentages.",
      gc_loess_invalid_gc_nonpositive = "Invalid-bin count due to GC percentages less than or equal to zero.",
      gc_loess_invalid_raw_nonfinite = "Invalid-bin count due to non-finite raw counts.",
      gc_loess_invalid_raw_nonpositive = "Invalid-bin count due to raw counts less than or equal to zero.",
      gc_loess_invalid_corrected_nonfinite = "Count of fit-eligible autosomal bins with non-finite corrected counts in post-correction summaries.",
      gc_curve_has_valid_bins = "Whether the sample had at least one valid autosomal bin for GC-curve reporting.",
      gc_curve_has_loess_support = "Whether the sample had enough valid bins and GC spread to support LOESS fitting.",
      nipter_autosomal_bin_cv_pre_gc_correction = "Autosomal bin-level coefficient of variation before GC correction.",
      nipter_autosomal_bin_cv_post_gc_correction = "Autosomal bin-level coefficient of variation after GC correction.",
      nipter_gc_correction_bin_scale_mean = "Mean autosomal bin-level GC-correction scale factor computed as corrected_count/raw_count after GC correction.",
      nipter_gc_correction_bin_scale_sd = "Standard deviation of the autosomal bin-level GC-correction scale factor computed as corrected_count/raw_count after GC correction.",
      nipter_gc_correction_bin_scale_cv = "Coefficient of variation of the autosomal bin-level GC-correction scale factor computed as corrected_count/raw_count after GC correction.",
      nipter_gc_correlation_pre_gc_correction = "GC-to-count correlation on autosomal bins before GC correction.",
      nipter_gc_correlation_post_gc_correction = "GC-to-count correlation on autosomal bins after GC correction.",
      nipter_gc_curve_png = "Optional path to the per-sample GC-curve PNG.",
      nipter_gc_curve_data_bgz = "Path to the per-sample GC-curve interval table written as TSV.bgz with a tabix index.",
      nipter_corrected_bin_ratio_genome_plot_png = "Optional path to the per-sample genome-order plot of post-GC corrected bin ratios excluding out-of-range bins.",
      seqff_pre = "SeqFF estimate from the nonfiltered auxiliary metrics pass.",
      enet_pre = "Elastic-net SeqFF component from the nonfiltered auxiliary metrics pass.",
      wrsc_pre = "WRSC SeqFF component from the nonfiltered auxiliary metrics pass.",
      fetal_fraction_pre = "Current nonfiltered fetal-fraction summary used by the QC layer.",
      seqff_post = "SeqFF estimate from the filtered auxiliary metrics pass.",
      enet_post = "Elastic-net SeqFF component from the filtered auxiliary metrics pass.",
      wrsc_post = "WRSC SeqFF component from the filtered auxiliary metrics pass.",
      fetal_fraction_post = "Current filtered fetal-fraction summary used by the QC layer.",
      y_unique_ratio_pre = "Y-unique ratio from the nonfiltered auxiliary metrics pass.",
      y_unique_reads_pre = "Y-unique read count from the nonfiltered auxiliary metrics pass.",
      total_nuclear_reads_pre = "Total nuclear read denominator used for y_unique_ratio_pre.",
      y_unique_ratio_post = "Y-unique ratio from the filtered auxiliary metrics pass.",
      y_unique_reads_post = "Y-unique read count from the filtered auxiliary metrics pass.",
      total_nuclear_reads_post = "Total nuclear read denominator used for y_unique_ratio_post.",
      y_unique_regions_file = "Path to the Y-specific regions file used for Y-unique metrics.",
      metrics_pre_mapq = "MAPQ threshold used for the nonfiltered auxiliary metrics pass.",
      metrics_pre_require_flags = "Required SAM flag mask used for the nonfiltered auxiliary metrics pass.",
      metrics_pre_exclude_flags = "Excluded SAM flag mask used for the nonfiltered auxiliary metrics pass.",
      metrics_post_mapq = "MAPQ threshold used for the filtered auxiliary metrics pass.",
      metrics_post_require_flags = "Required SAM flag mask used for the filtered auxiliary metrics pass.",
      metrics_post_exclude_flags = "Excluded SAM flag mask used for the filtered auxiliary metrics pass.",
      seqff_source = "Provenance for SeqFF values, for example computed or file:<name>.",
      y_unique_source = "Provenance for Y-unique values, for example computed or file:<name>.",
      read_counts_source = "Provenance for alignment and binned read-count summary metrics.",
      nipter_bed_status = "NIPTeR BED write status for this sample.",
      nipter_bed_written = "Whether the NIPTeR BED artifact was written successfully.",
      nipter_bed_error = "Error text recorded when NIPTeR BED writing failed.",
      "Undocumented sample QC column."
    )
  }, character(1L))

  include_in_seqff_summary <- column_name %in% c(
    "sample_name", "bam",
    "seqff_pre", "enet_pre", "wrsc_pre", "fetal_fraction_pre",
    "seqff_post", "enet_post", "wrsc_post", "fetal_fraction_post",
    "metrics_pre_mapq", "metrics_pre_require_flags", "metrics_pre_exclude_flags",
    "metrics_post_mapq", "metrics_post_require_flags", "metrics_post_exclude_flags",
    "seqff_source"
  )

  include_in_y_unique_summary <- column_name %in% c(
    "sample_name", "bam",
    "y_unique_ratio_pre", "y_unique_reads_pre", "total_nuclear_reads_pre",
    "y_unique_ratio_post", "y_unique_reads_post", "total_nuclear_reads_post",
    "y_unique_regions_file",
    "metrics_pre_mapq", "metrics_pre_require_flags", "metrics_pre_exclude_flags",
    "metrics_post_mapq", "metrics_post_require_flags", "metrics_post_exclude_flags",
    "y_unique_source"
  )

  data.frame(
    column_name = column_name,
    layer = layer,
    description = description,
    required_on_success = required_on_success,
    include_in_sample_qc = TRUE,
    include_in_seqff_summary = include_in_seqff_summary,
    include_in_y_unique_summary = include_in_y_unique_summary
  )
})

.seqff_summary_schema <- function() {
  .sample_qc_column_dictionary$column_name[
    .sample_qc_column_dictionary$include_in_seqff_summary
  ]
}

.y_unique_summary_schema <- function() {
  .sample_qc_column_dictionary$column_name[
    .sample_qc_column_dictionary$include_in_y_unique_summary
  ]
}

.sample_qc_schema <- function() {
  .sample_qc_column_dictionary$column_name[
    .sample_qc_column_dictionary$include_in_sample_qc
  ]
}

.sample_qc_required_on_success <- function() {
  .sample_qc_column_dictionary$column_name[
    .sample_qc_column_dictionary$required_on_success
  ]
}

.sample_qc_column_dictionary_frame <- function() {
  .sample_qc_column_dictionary
}

.validate_sample_qc_row <- function(row, require_on_success = FALSE) {
  stopifnot(is.data.frame(row), nrow(row) == 1L)

  schema <- .sample_qc_schema()
  missing_cols <- setdiff(schema, names(row))
  if (length(missing_cols)) {
    stop(
      "Sample QC row is missing schema columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(require_on_success)) {
    required_cols <- .sample_qc_required_on_success()
    missing_vals <- required_cols[!vapply(required_cols, function(col) {
      .sample_qc_value_present(row[[col]])
    }, logical(1L))]
    if (length(missing_vals)) {
      stop(
        "Successful sample QC row is missing required values: ",
        paste(missing_vals, collapse = ", "),
        call. = FALSE
      )
    }
  }

  row[, schema, drop = FALSE]
}

.read_metrics_table <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))
  stopifnot(file.exists(path))

  first <- readLines(path, n = 1L, warn = FALSE)
  delim <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    ","
  } else if (length(first) && gregexpr(",", first, fixed = TRUE)[[1L]][1L] > 0L &&
             (!grepl("\t", first, fixed = TRUE) ||
                lengths(regmatches(first, gregexpr(",", first, fixed = TRUE))) >=
                  lengths(regmatches(first, gregexpr("\t", first, fixed = TRUE))))) {
    ","
  } else {
    "\t"
  }

  utils::read.delim(
    path,
    sep = delim,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.match_sample_col <- function(df,
                              candidates = c(
                                "sample_name", "sample", "Sample_name", "Sample",
                                "sample_id", "SampleID"
                              )) {
  hit <- candidates[candidates %in% names(df)][1L]
  if (is.na(hit) || !nzchar(hit)) NULL else hit
}

.match_metric_row <- function(path, sample_name) {
  df <- .read_metrics_table(path)
  sample_col <- .match_sample_col(df)
  if (is.null(sample_col)) {
    if (nrow(df) == 1L) {
      return(df[1L, , drop = FALSE])
    }
    stop(
      "Cannot infer sample column in metrics table: ", path,
      call. = FALSE
    )
  }

  hit <- which(as.character(df[[sample_col]]) == sample_name)
  if (!length(hit)) {
    stop(
      "Sample '", sample_name, "' not found in metrics table: ", path,
      call. = FALSE
    )
  }
  if (length(hit) > 1L) {
    stop(
      "Sample '", sample_name, "' appears multiple times in metrics table: ", path,
      call. = FALSE
    )
  }
  df[hit, , drop = FALSE]
}

.record_to_list <- function(record) {
  if (is.null(record)) {
    return(NULL)
  }
  if (is.list(record) && !length(record)) {
    return(list())
  }
  if (is.data.frame(record)) {
    if (nrow(record) != 1L) {
      stop("Metric records supplied as data frames must have exactly one row.", call. = FALSE)
    }
    out <- as.list(record[1L, , drop = FALSE])
    names(out) <- names(record)
    return(out)
  }
  if (is.atomic(record) && !is.null(names(record))) {
    return(as.list(record))
  }
  if (is.list(record) && !is.null(names(record))) {
    return(record)
  }
  stop("Metric record must be a named vector, named list, one-row data frame, or NULL.",
       call. = FALSE)
}

.metric_scalar <- function(record, candidates, type = c("numeric", "character", "logical")) {
  type <- match.arg(type)
  rec <- .record_to_list(record)
  if (is.null(rec) || !length(rec)) {
    return(NULL)
  }
  hit <- candidates[candidates %in% names(rec)][1L]
  if (is.na(hit) || !nzchar(hit)) {
    return(NULL)
  }
  value <- rec[[hit]]
  if (!length(value)) {
    return(NULL)
  }
  value <- value[[1L]]
  if (type == "numeric") {
    out <- suppressWarnings(as.numeric(value))
    if (!length(out)) return(NULL)
    out
  } else if (type == "logical") {
    if (is.logical(value)) {
      out <- value[[1L]]
    } else {
      chr <- trimws(tolower(as.character(value)))
      if (isTRUE(is.na(chr)) || !nzchar(chr) || identical(chr, "na")) {
        return(NULL)
      }
      if (chr %in% c("true", "t", "1")) {
        out <- TRUE
      } else if (chr %in% c("false", "f", "0")) {
        out <- FALSE
      } else {
        return(NULL)
      }
    }
    out
  } else {
    if (isTRUE(is.na(value))) {
      NA_character_
    } else {
      as.character(value)
    }
  }
}

.read_y_unique_regions_file <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))
  stopifnot(file.exists(path))

  df <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df) && all(c("Chromosome", "Start", "End") %in% names(df))) {
    if (!"GeneName" %in% names(df)) {
      df$GeneName <- NA_character_
    }
    return(df[, c("Chromosome", "Start", "End", "GeneName"), drop = FALSE])
  }

  df <- utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(df) < 3L) {
    stop("Y-unique regions file must have at least 3 columns: ", path, call. = FALSE)
  }
  df <- df[, seq_len(min(4L, ncol(df))), drop = FALSE]
  names(df)[1:3] <- c("Chromosome", "Start", "End")
  if (ncol(df) < 4L) {
    df$GeneName <- NA_character_
  } else {
    names(df)[4L] <- "GeneName"
  }
  df
}

.format_y_unique_genes <- function(regions) {
  if (is.null(regions) || !nrow(regions) || !"GeneName" %in% names(regions)) {
    return(NULL)
  }
  vals <- unique(as.character(regions$GeneName))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) NULL else paste(vals, collapse = ";")
}

.format_y_unique_intervals <- function(regions) {
  if (is.null(regions) || !nrow(regions)) {
    return(NULL)
  }
  gene <- if ("GeneName" %in% names(regions)) as.character(regions$GeneName) else rep(NA_character_, nrow(regions))
  paste(
    sprintf(
      "%s:%s-%s%s",
      as.character(regions$Chromosome),
      as.integer(regions$Start),
      as.integer(regions$End),
      ifelse(!is.na(gene) & nzchar(gene), paste0("(", gene, ")"), "")
    ),
    collapse = ";"
  )
}

.seqff_metric_record <- function(record = NULL,
                                 pre = NULL,
                                 post = NULL,
                                 source = NULL) {
  rec <- .record_to_list(record)
  pre_rec <- .record_to_list(pre)
  post_rec <- .record_to_list(post)

  out <- list(
    seqff_pre = .metric_scalar(rec, c("SeqFF_pre", "seqff_pre"), "numeric"),
    enet_pre = .metric_scalar(rec, c("Enet_pre", "enet_pre"), "numeric"),
    wrsc_pre = .metric_scalar(rec, c("WRSC_pre", "wrsc_pre"), "numeric"),
    seqff_post = .metric_scalar(rec, c("SeqFF_post", "seqff_post"), "numeric"),
    enet_post = .metric_scalar(rec, c("Enet_post", "enet_post"), "numeric"),
    wrsc_post = .metric_scalar(rec, c("WRSC_post", "wrsc_post"), "numeric"),
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("seqff_source", "metrics_seqff_source"), "character")
    )
  )

  if (!is.null(pre_rec)) {
    out$seqff_pre <- .metric_scalar(pre_rec, c("SeqFF", "seqff", "SeqFF_pre", "seqff_pre"), "numeric")
    out$enet_pre <- .metric_scalar(pre_rec, c("Enet", "enet", "Enet_pre", "enet_pre"), "numeric")
    out$wrsc_pre <- .metric_scalar(pre_rec, c("WRSC", "wrsc", "WRSC_pre", "wrsc_pre"), "numeric")
  }
  if (!is.null(post_rec)) {
    out$seqff_post <- .metric_scalar(post_rec, c("SeqFF", "seqff", "SeqFF_post", "seqff_post"), "numeric")
    out$enet_post <- .metric_scalar(post_rec, c("Enet", "enet", "Enet_post", "enet_post"), "numeric")
    out$wrsc_post <- .metric_scalar(post_rec, c("WRSC", "wrsc", "WRSC_post", "wrsc_post"), "numeric")
  }

  out
}

.y_unique_metric_record <- function(record = NULL,
                                    pre = NULL,
                                    post = NULL,
                                    regions = NULL,
                                    regions_file = NULL,
                                    source = NULL) {
  rec <- .record_to_list(record)
  pre_rec <- .record_to_list(pre)
  post_rec <- .record_to_list(post)

  out <- list(
    ratio_pre = .metric_scalar(rec, c("y_unique_ratio_pre"), "numeric"),
    reads_pre = .metric_scalar(rec, c("y_unique_reads_pre"), "numeric"),
    total_pre = .metric_scalar(rec, c("total_nuclear_reads_pre"), "numeric"),
    ratio_post = .metric_scalar(rec, c("y_unique_ratio_post"), "numeric"),
    reads_post = .metric_scalar(rec, c("y_unique_reads_post"), "numeric"),
    total_post = .metric_scalar(rec, c("total_nuclear_reads_post"), "numeric"),
    regions_file = .null_coalesce(
      regions_file,
      .metric_scalar(rec, c("y_unique_regions_file", "YUniqueRegionsFile"), "character")
    ),
    genes = .metric_scalar(rec, c("y_unique_genes", "YUniqueGenes"), "character"),
    intervals = .metric_scalar(rec, c("y_unique_regions", "YUniqueRegions"), "character"),
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("y_unique_source", "metrics_y_unique_source"), "character")
    )
  )

  if (!is.null(pre_rec)) {
    out$ratio_pre <- .metric_scalar(pre_rec, c("ratio", "y_unique_ratio", "YUniqueRatio", "Ratio"), "numeric")
    out$reads_pre <- .metric_scalar(pre_rec, c("y_unique_reads", "YUniqueReads"), "numeric")
    out$total_pre <- .metric_scalar(pre_rec, c("total_nuclear_reads", "TotalNuclearReads"), "numeric")
  }
  if (!is.null(post_rec)) {
    out$ratio_post <- .metric_scalar(post_rec, c("ratio", "y_unique_ratio", "YUniqueRatio", "Ratio"), "numeric")
    out$reads_post <- .metric_scalar(post_rec, c("y_unique_reads", "YUniqueReads"), "numeric")
    out$total_post <- .metric_scalar(post_rec, c("total_nuclear_reads", "TotalNuclearReads"), "numeric")
  }

  regions_df <- NULL
  if (!is.null(post_rec) && is.data.frame(post_rec$regions)) {
    regions_df <- post_rec$regions
  } else if (!is.null(pre_rec) && is.data.frame(pre_rec$regions)) {
    regions_df <- pre_rec$regions
  } else if (is.data.frame(regions)) {
    regions_df <- regions
  } else if (!is.null(out$regions_file) && nzchar(out$regions_file) && file.exists(out$regions_file)) {
    regions_df <- .read_y_unique_regions_file(out$regions_file)
  }

  if (!is.null(regions_df)) {
    out$genes <- .null_coalesce(out$genes, .format_y_unique_genes(regions_df))
    out$intervals <- .null_coalesce(out$intervals, .format_y_unique_intervals(regions_df))
    if (is.null(out$regions_file) && !is.null(regions_file) && file.exists(regions_file)) {
      out$regions_file <- normalizePath(regions_file, winslash = "/", mustWork = TRUE)
    }
  } else if (!is.null(regions_file) && file.exists(regions_file)) {
    out$regions_file <- normalizePath(regions_file, winslash = "/", mustWork = TRUE)
  }

  out
}

.read_counts_metric_record <- function(record = NULL,
                                       rows = NULL,
                                       source = NULL) {
  if (!is.null(rows)) {
    stopifnot(is.data.frame(rows))
    record <- .bam_read_counts_metadata(rows)
  }
  rec <- .record_to_list(record)

  out <- list(
    input_total = .metric_scalar(rec, c("read_counts_input_total"), "numeric"),
    mapped_total = .metric_scalar(rec, c("read_counts_mapped_total"), "numeric"),
    unmapped_total = .metric_scalar(rec, c("read_counts_unmapped_total"), "numeric"),
    count_pre = .metric_scalar(rec, c("read_counts_binned_pre_sum"), "numeric"),
    count_post = .metric_scalar(rec, c("read_counts_binned_post_sum"), "numeric"),
    count_fwd = .metric_scalar(rec, c("read_counts_binned_fwd_sum"), "numeric"),
    count_rev = .metric_scalar(rec, c("read_counts_binned_rev_sum"), "numeric"),
    n_nonzero_bins_post = .metric_scalar(
      rec,
      c("read_counts_binned_nonzero_bins_post"),
      "numeric"
    ),
    retained_fraction = .metric_scalar(rec, c("read_counts_binned_retained_fraction"), "numeric"),
    gc_pre = .metric_scalar(rec, c("gc_read_perc_pre"), "numeric"),
    gc_post = .metric_scalar(rec, c("gc_read_perc_post"), "numeric"),
    mean_mapq_post = .metric_scalar(rec, c("mean_mapq_post"), "numeric"),
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("read_counts_source"), "character")
    )
  )

  if (is.null(out$source) && !is.null(rows)) {
    out$source <- "bam_bin_counts(gc,mq)"
  }

  out
}

.samtools_stats_metric_record <- function(record = NULL,
                                          pre = NULL,
                                          post = NULL,
                                          source = NULL) {
  rec <- .record_to_list(record)
  pre_rec <- .record_to_list(pre)
  post_rec <- .record_to_list(post)

  parse_prefixed_stage <- function(stage) {
    out <- list()
    out[[paste0(stage, "_bookkeeping_mode")]] <- .metric_scalar(
      rec,
      c(
        paste0("samtools_stats_", stage, "_bookkeeping_mode"),
        paste0(stage, "_bookkeeping_mode")
      ),
      "character"
    )
    for (field in .samtools_stats_summary_fields) {
      out[[paste0(stage, "_", field)]] <- .metric_scalar(
        rec,
        c(
          paste0("samtools_stats_", stage, "_", field),
          paste0(stage, "_", field)
        ),
        "numeric"
      )
    }
    out
  }

  parse_summary_stage <- function(stage, stage_rec) {
    out <- list()
    out[[paste0(stage, "_bookkeeping_mode")]] <- if (isTRUE(.metric_scalar(
      stage_rec,
      c("report_filtered_stream_bookkeeping"),
      "logical"
    ))) {
      "filtered_stream"
    } else {
      "input_file"
    }
    for (field in .samtools_stats_summary_fields) {
      out[[paste0(stage, "_", field)]] <- .metric_scalar(
        stage_rec,
        c(field),
        "numeric"
      )
    }
    out
  }

  out <- list(
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("samtools_stats_source", "source"), "character")
    )
  )

  if (!is.null(pre_rec)) {
    out <- utils::modifyList(out, parse_summary_stage("pre", pre_rec))
  } else {
    out <- utils::modifyList(out, parse_prefixed_stage("pre"))
  }

  if (!is.null(post_rec)) {
    out <- utils::modifyList(out, parse_summary_stage("post", post_rec))
  } else {
    out <- utils::modifyList(out, parse_prefixed_stage("post"))
  }

  out
}

.read_mosdepth_summary_total <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))
  stopifnot(file.exists(path))

  df <- utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
  )
  required_cols <- c("chrom", "length", "bases", "mean", "min", "max")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols)) {
    stop(
      "Mosdepth summary is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\n  file: ", path,
      call. = FALSE
    )
  }

  total_idx <- match(c("total_region", "total"), as.character(df$chrom), nomatch = 0L)
  total_idx <- total_idx[total_idx > 0L][1L]
  if (!length(total_idx) || is.na(total_idx) || total_idx < 1L) {
    stop("Mosdepth summary does not contain a total row: ", path, call. = FALSE)
  }

  row <- df[total_idx, , drop = FALSE]
  length_bases <- suppressWarnings(as.numeric(row$length[[1L]]))
  bases <- suppressWarnings(as.numeric(row$bases[[1L]]))
  mean_depth <- suppressWarnings(as.numeric(row$mean[[1L]]))
  min_depth <- suppressWarnings(as.numeric(row$min[[1L]]))
  max_depth <- suppressWarnings(as.numeric(row$max[[1L]]))

  if (!all(is.finite(c(length_bases, bases, mean_depth, min_depth, max_depth)))) {
    stop("Mosdepth total row contains non-finite numeric values: ", path, call. = FALSE)
  }
  if (length_bases <= 0) {
    stop("Mosdepth total row has non-positive reference length: ", path, call. = FALSE)
  }

  list(
    summary_txt = normalizePath(path, winslash = "/", mustWork = TRUE),
    length_bases = length_bases,
    bases = bases,
    covered_fraction = bases / length_bases,
    mean_depth = mean_depth,
    min_depth = min_depth,
    max_depth = max_depth
  )
}

.coverage_metric_record <- function(record = NULL,
                                    pre_summary = NULL,
                                    post_summary = NULL,
                                    binsize = NULL,
                                    pre_mapq = NULL,
                                    pre_exclude_flags = NULL,
                                    post_mapq = NULL,
                                    post_exclude_flags = NULL,
                                    source = NULL) {
  rec <- .record_to_list(record)

  out <- list(
    binsize = .metric_scalar(rec, c("coverage_binsize", "binsize"), "numeric"),
    pre_mapq = .metric_scalar(rec, c("coverage_pre_mapq", "pre_mapq"), "numeric"),
    pre_exclude_flags = .metric_scalar(rec, c("coverage_pre_exclude_flags", "pre_exclude_flags"), "numeric"),
    pre_summary_txt = .metric_scalar(rec, c("coverage_pre_summary_txt", "pre_summary_txt"), "character"),
    pre_length_bases = .metric_scalar(rec, c("coverage_pre_length_bases", "pre_length_bases"), "numeric"),
    pre_bases = .metric_scalar(rec, c("coverage_pre_bases", "pre_bases"), "numeric"),
    pre_covered_fraction = .metric_scalar(rec, c("coverage_pre_covered_fraction", "pre_covered_fraction"), "numeric"),
    pre_mean_depth = .metric_scalar(rec, c("coverage_pre_mean_depth", "pre_mean_depth"), "numeric"),
    pre_min_depth = .metric_scalar(rec, c("coverage_pre_min_depth", "pre_min_depth"), "numeric"),
    pre_max_depth = .metric_scalar(rec, c("coverage_pre_max_depth", "pre_max_depth"), "numeric"),
    post_mapq = .metric_scalar(rec, c("coverage_post_mapq", "post_mapq"), "numeric"),
    post_exclude_flags = .metric_scalar(rec, c("coverage_post_exclude_flags", "post_exclude_flags"), "numeric"),
    post_summary_txt = .metric_scalar(rec, c("coverage_post_summary_txt", "post_summary_txt"), "character"),
    post_length_bases = .metric_scalar(rec, c("coverage_post_length_bases", "post_length_bases"), "numeric"),
    post_bases = .metric_scalar(rec, c("coverage_post_bases", "post_bases"), "numeric"),
    post_covered_fraction = .metric_scalar(rec, c("coverage_post_covered_fraction", "post_covered_fraction"), "numeric"),
    post_mean_depth = .metric_scalar(rec, c("coverage_post_mean_depth", "post_mean_depth"), "numeric"),
    post_min_depth = .metric_scalar(rec, c("coverage_post_min_depth", "post_min_depth"), "numeric"),
    post_max_depth = .metric_scalar(rec, c("coverage_post_max_depth", "post_max_depth"), "numeric"),
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("coverage_source", "source"), "character")
    )
  )

  if (!is.null(pre_summary)) {
    pre <- .read_mosdepth_summary_total(pre_summary)
    out$pre_summary_txt <- pre$summary_txt
    out$pre_length_bases <- pre$length_bases
    out$pre_bases <- pre$bases
    out$pre_covered_fraction <- pre$covered_fraction
    out$pre_mean_depth <- pre$mean_depth
    out$pre_min_depth <- pre$min_depth
    out$pre_max_depth <- pre$max_depth
  }
  if (!is.null(post_summary)) {
    post <- .read_mosdepth_summary_total(post_summary)
    out$post_summary_txt <- post$summary_txt
    out$post_length_bases <- post$length_bases
    out$post_bases <- post$bases
    out$post_covered_fraction <- post$covered_fraction
    out$post_mean_depth <- post$mean_depth
    out$post_min_depth <- post$min_depth
    out$post_max_depth <- post$max_depth
  }

  if (!is.null(binsize)) {
    out$binsize <- as.numeric(binsize)[[1L]]
  }
  if (!is.null(pre_mapq)) {
    out$pre_mapq <- as.numeric(pre_mapq)[[1L]]
  }
  if (!is.null(pre_exclude_flags)) {
    out$pre_exclude_flags <- as.numeric(pre_exclude_flags)[[1L]]
  }
  if (!is.null(post_mapq)) {
    out$post_mapq <- as.numeric(post_mapq)[[1L]]
  }
  if (!is.null(post_exclude_flags)) {
    out$post_exclude_flags <- as.numeric(post_exclude_flags)[[1L]]
  }
  if (is.null(out$source) && (!is.null(pre_summary) || !is.null(post_summary))) {
    out$source <- "native_mosdepth_dense"
  }

  out
}

.nipter_preprocess_qc_metric_record <- function(record = NULL) {
  rec <- .record_to_list(record)
  if (is.null(rec) || !length(rec)) {
    return(list())
  }

  list(
    gc_loess_valid_bins = .metric_scalar(rec, c("gc_loess_valid_bins"), "numeric"),
    gc_loess_total_bins = .metric_scalar(rec, c("gc_loess_total_bins"), "numeric"),
    gc_loess_invalid_bins = .metric_scalar(rec, c("gc_loess_invalid_bins"), "numeric"),
    gc_loess_valid_bin_fraction_pct = .metric_scalar(rec, c("gc_loess_valid_bin_fraction_pct"), "numeric"),
    gc_loess_invalid_out_of_chromosome_range = .metric_scalar(rec, c("gc_loess_invalid_out_of_chromosome_range"), "numeric"),
    gc_loess_invalid_gc_missing_or_nonfinite = .metric_scalar(rec, c("gc_loess_invalid_gc_missing_or_nonfinite"), "numeric"),
    gc_loess_invalid_gc_nonpositive = .metric_scalar(rec, c("gc_loess_invalid_gc_nonpositive"), "numeric"),
    gc_loess_invalid_raw_nonfinite = .metric_scalar(rec, c("gc_loess_invalid_raw_nonfinite"), "numeric"),
    gc_loess_invalid_raw_nonpositive = .metric_scalar(rec, c("gc_loess_invalid_raw_nonpositive"), "numeric"),
    gc_loess_invalid_corrected_nonfinite = .metric_scalar(rec, c("gc_loess_invalid_corrected_nonfinite"), "numeric"),
    gc_curve_has_valid_bins = .metric_scalar(rec, c("gc_curve_has_valid_bins"), "logical"),
    gc_curve_has_loess_support = .metric_scalar(rec, c("gc_curve_has_loess_support"), "logical"),
    nipter_autosomal_bin_cv_pre_gc_correction = .metric_scalar(rec, c("nipter_autosomal_bin_cv_pre_gc_correction"), "numeric"),
    nipter_autosomal_bin_cv_post_gc_correction = .metric_scalar(rec, c("nipter_autosomal_bin_cv_post_gc_correction"), "numeric"),
    nipter_gc_correction_bin_scale_mean = .metric_scalar(rec, c("nipter_gc_correction_bin_scale_mean"), "numeric"),
    nipter_gc_correction_bin_scale_sd = .metric_scalar(rec, c("nipter_gc_correction_bin_scale_sd"), "numeric"),
    nipter_gc_correction_bin_scale_cv = .metric_scalar(rec, c("nipter_gc_correction_bin_scale_cv"), "numeric"),
    nipter_gc_correlation_pre_gc_correction = .metric_scalar(rec, c("nipter_gc_correlation_pre_gc_correction"), "numeric"),
    nipter_gc_correlation_post_gc_correction = .metric_scalar(rec, c("nipter_gc_correlation_post_gc_correction"), "numeric"),
    gc_curve_plot = .metric_scalar(rec, c("nipter_gc_curve_png"), "character"),
    gc_curve_data_bgz = .metric_scalar(rec, c("nipter_gc_curve_data_bgz"), "character"),
    corrected_bin_ratio_genome_plot = .metric_scalar(rec, c("nipter_corrected_bin_ratio_genome_plot_png"), "character")
  )
}

.gc_correction_record <- function(record = NULL,
                                  applied = NULL,
                                  method = NULL,
                                  loess_span = NULL,
                                  include_sex = NULL,
                                  binsize = NULL,
                                  table_bgz = NULL,
                                  fasta = NULL) {
  rec <- .record_to_list(record)

  parsed_applied <- if (is.null(applied)) {
    .metric_scalar(
      rec,
      c("applied", "nipter_gc_correction_applied", "gc_correction_applied"),
      "logical"
    )
  } else {
    as.logical(applied)[[1L]]
  }
  parsed_method <- if (is.null(method)) {
    .metric_scalar(
      rec,
      c("method", "nipter_gc_correction_method", "gc_correction_method"),
      "character"
    )
  } else {
    as.character(method)[[1L]]
  }
  parsed_loess_span <- if (is.null(loess_span)) {
    .metric_scalar(
      rec,
      c("loess_span", "nipter_gc_correction_loess_span", "gc_correction_loess_span"),
      "numeric"
    )
  } else {
    as.numeric(loess_span)[[1L]]
  }
  parsed_include_sex <- if (is.null(include_sex)) {
    .metric_scalar(
      rec,
      c("include_sex", "nipter_gc_correction_include_sex", "gc_correction_include_sex"),
      "logical"
    )
  } else {
    as.logical(include_sex)[[1L]]
  }
  parsed_binsize <- if (is.null(binsize)) {
    .metric_scalar(
      rec,
      c("binsize", "nipter_gc_correction_binsize", "gc_correction_binsize"),
      "numeric"
    )
  } else {
    as.numeric(binsize)[[1L]]
  }
  parsed_table_bgz <- if (is.null(table_bgz)) {
    .metric_scalar(
      rec,
      c("table_bgz", "nipter_gc_correction_table_bgz", "gc_correction_table_bgz"),
      "character"
    )
  } else {
    as.character(table_bgz)[[1L]]
  }
  parsed_fasta <- if (is.null(fasta)) {
    .metric_scalar(
      rec,
      c("fasta", "nipter_gc_correction_fasta", "gc_correction_fasta"),
      "character"
    )
  } else {
    as.character(fasta)[[1L]]
  }

  parsed_applied <- if (is.null(parsed_applied) || is.na(parsed_applied)) FALSE else isTRUE(parsed_applied)
  parsed_method <- if (is.null(parsed_method) || !nzchar(parsed_method) || is.na(parsed_method)) {
    if (parsed_applied) "loess" else "none"
  } else {
    parsed_method
  }
  parsed_include_sex <- if (is.null(parsed_include_sex) || is.na(parsed_include_sex)) FALSE else isTRUE(parsed_include_sex)

  list(
    applied = parsed_applied,
    method = parsed_method,
    loess_span = if (is.null(parsed_loess_span) || is.na(parsed_loess_span)) NA_real_ else parsed_loess_span,
    include_sex = parsed_include_sex,
    binsize = if (is.null(parsed_binsize) || is.na(parsed_binsize)) NA_real_ else parsed_binsize,
    table_bgz = if (is.null(parsed_table_bgz) || is.na(parsed_table_bgz) || !nzchar(parsed_table_bgz)) NA_character_ else parsed_table_bgz,
    fasta = if (is.null(parsed_fasta) || is.na(parsed_fasta) || !nzchar(parsed_fasta)) NA_character_ else parsed_fasta
  )
}

.metric_filter_record <- function(record = NULL, stage = c("pre", "post")) {
  stage <- match.arg(stage)
  rec <- .record_to_list(record)
  if (is.null(rec) || !length(rec)) {
    return(list())
  }

  list(
    mapq = .metric_scalar(rec, c(paste0("metrics_", stage, "_mapq")), "numeric"),
    require_flags = .metric_scalar(rec, c(paste0("metrics_", stage, "_require_flags")), "numeric"),
    exclude_flags = .metric_scalar(rec, c(paste0("metrics_", stage, "_exclude_flags")), "numeric")
  )
}

.sample_qc_row <- function(sample_name,
                           bam = NULL,
                           seqff = NULL,
                           y_unique = NULL,
                           read_counts = NULL,
                           bam_stats = NULL,
                           coverage = NULL,
                           nipter_qc = NULL,
                           gc_correction = NULL,
                           filters_pre = NULL,
                           filters_post = NULL) {
  stopifnot(is.character(sample_name), length(sample_name) == 1L, nzchar(sample_name))
  seqff <- .null_coalesce(seqff, list())
  y_unique <- .null_coalesce(y_unique, list())
  read_counts <- .null_coalesce(read_counts, list())
  bam_stats <- .samtools_stats_metric_record(record = bam_stats)
  coverage <- .coverage_metric_record(record = coverage)
  nipter_qc <- .null_coalesce(nipter_qc, list())
  gc_correction <- .gc_correction_record(record = gc_correction)
  filters_pre <- .null_coalesce(.record_to_list(filters_pre), list())
  filters_post <- .null_coalesce(.record_to_list(filters_post), list())

  row <- data.frame(
    sample_name = sample_name,
    bam = if (is.null(bam)) NA_character_ else bam,
    sample_processing_status = NA_character_,
    sample_failure_stage = NA_character_,
    sample_failure_step = NA_character_,
    sample_failure_message = NA_character_,
    read_counts_input_total = unname(.null_coalesce(read_counts$input_total, NA_real_)),
    read_counts_mapped_total = unname(.null_coalesce(read_counts$mapped_total, NA_real_)),
    read_counts_unmapped_total = unname(.null_coalesce(read_counts$unmapped_total, NA_real_)),
    read_counts_binned_pre_sum = unname(.null_coalesce(read_counts$count_pre, NA_real_)),
    read_counts_binned_post_sum = unname(.null_coalesce(read_counts$count_post, NA_real_)),
    read_counts_binned_fwd_sum = unname(.null_coalesce(read_counts$count_fwd, NA_real_)),
    read_counts_binned_rev_sum = unname(.null_coalesce(read_counts$count_rev, NA_real_)),
    read_counts_binned_nonzero_bins_post = unname(.null_coalesce(read_counts$n_nonzero_bins_post, NA_real_)),
    read_counts_binned_retained_fraction = unname(.null_coalesce(read_counts$retained_fraction, NA_real_)),
    gc_read_perc_pre = unname(.null_coalesce(read_counts$gc_pre, NA_real_)),
    gc_read_perc_post = unname(.null_coalesce(read_counts$gc_post, NA_real_)),
    mean_mapq_post = unname(.null_coalesce(read_counts$mean_mapq_post, NA_real_)),
    coverage_binsize = unname(.null_coalesce(coverage$binsize, NA_real_)),
    coverage_pre_mapq = unname(.null_coalesce(coverage$pre_mapq, NA_real_)),
    coverage_pre_exclude_flags = unname(.null_coalesce(coverage$pre_exclude_flags, NA_real_)),
    coverage_pre_summary_txt = .null_coalesce(coverage$pre_summary_txt, NA_character_),
    coverage_pre_length_bases = unname(.null_coalesce(coverage$pre_length_bases, NA_real_)),
    coverage_pre_bases = unname(.null_coalesce(coverage$pre_bases, NA_real_)),
    coverage_pre_covered_fraction = unname(.null_coalesce(coverage$pre_covered_fraction, NA_real_)),
    coverage_pre_mean_depth = unname(.null_coalesce(coverage$pre_mean_depth, NA_real_)),
    coverage_pre_min_depth = unname(.null_coalesce(coverage$pre_min_depth, NA_real_)),
    coverage_pre_max_depth = unname(.null_coalesce(coverage$pre_max_depth, NA_real_)),
    coverage_post_mapq = unname(.null_coalesce(coverage$post_mapq, NA_real_)),
    coverage_post_exclude_flags = unname(.null_coalesce(coverage$post_exclude_flags, NA_real_)),
    coverage_post_summary_txt = .null_coalesce(coverage$post_summary_txt, NA_character_),
    coverage_post_length_bases = unname(.null_coalesce(coverage$post_length_bases, NA_real_)),
    coverage_post_bases = unname(.null_coalesce(coverage$post_bases, NA_real_)),
    coverage_post_covered_fraction = unname(.null_coalesce(coverage$post_covered_fraction, NA_real_)),
    coverage_post_mean_depth = unname(.null_coalesce(coverage$post_mean_depth, NA_real_)),
    coverage_post_min_depth = unname(.null_coalesce(coverage$post_min_depth, NA_real_)),
    coverage_post_max_depth = unname(.null_coalesce(coverage$post_max_depth, NA_real_)),
    coverage_source = .null_coalesce(coverage$source, NA_character_),
    nipter_gc_correction_applied = gc_correction$applied,
    nipter_gc_correction_method = gc_correction$method,
    nipter_gc_correction_loess_span = unname(.null_coalesce(gc_correction$loess_span, NA_real_)),
    nipter_gc_correction_include_sex = gc_correction$include_sex,
    nipter_gc_correction_binsize = unname(.null_coalesce(gc_correction$binsize, NA_real_)),
    nipter_gc_correction_table_bgz = .null_coalesce(gc_correction$table_bgz, NA_character_),
    nipter_gc_correction_fasta = .null_coalesce(gc_correction$fasta, NA_character_),
    gc_loess_valid_bins = unname(.null_coalesce(nipter_qc$gc_loess_valid_bins, NA_real_)),
    gc_loess_total_bins = unname(.null_coalesce(nipter_qc$gc_loess_total_bins, NA_real_)),
    gc_loess_invalid_bins = unname(.null_coalesce(nipter_qc$gc_loess_invalid_bins, NA_real_)),
    gc_loess_valid_bin_fraction_pct = unname(.null_coalesce(nipter_qc$gc_loess_valid_bin_fraction_pct, NA_real_)),
    gc_loess_invalid_out_of_chromosome_range = unname(.null_coalesce(nipter_qc$gc_loess_invalid_out_of_chromosome_range, NA_real_)),
    gc_loess_invalid_gc_missing_or_nonfinite = unname(.null_coalesce(nipter_qc$gc_loess_invalid_gc_missing_or_nonfinite, NA_real_)),
    gc_loess_invalid_gc_nonpositive = unname(.null_coalesce(nipter_qc$gc_loess_invalid_gc_nonpositive, NA_real_)),
    gc_loess_invalid_raw_nonfinite = unname(.null_coalesce(nipter_qc$gc_loess_invalid_raw_nonfinite, NA_real_)),
    gc_loess_invalid_raw_nonpositive = unname(.null_coalesce(nipter_qc$gc_loess_invalid_raw_nonpositive, NA_real_)),
    gc_loess_invalid_corrected_nonfinite = unname(.null_coalesce(nipter_qc$gc_loess_invalid_corrected_nonfinite, NA_real_)),
    gc_curve_has_valid_bins = .null_coalesce(nipter_qc$gc_curve_has_valid_bins, NA),
    gc_curve_has_loess_support = .null_coalesce(nipter_qc$gc_curve_has_loess_support, NA),
    nipter_autosomal_bin_cv_pre_gc_correction = unname(.null_coalesce(nipter_qc$nipter_autosomal_bin_cv_pre_gc_correction, NA_real_)),
    nipter_autosomal_bin_cv_post_gc_correction = unname(.null_coalesce(nipter_qc$nipter_autosomal_bin_cv_post_gc_correction, NA_real_)),
    nipter_gc_correction_bin_scale_mean = unname(.null_coalesce(nipter_qc$nipter_gc_correction_bin_scale_mean, NA_real_)),
    nipter_gc_correction_bin_scale_sd = unname(.null_coalesce(nipter_qc$nipter_gc_correction_bin_scale_sd, NA_real_)),
    nipter_gc_correction_bin_scale_cv = unname(.null_coalesce(nipter_qc$nipter_gc_correction_bin_scale_cv, NA_real_)),
    nipter_gc_correlation_pre_gc_correction = unname(.null_coalesce(nipter_qc$nipter_gc_correlation_pre_gc_correction, NA_real_)),
    nipter_gc_correlation_post_gc_correction = unname(.null_coalesce(nipter_qc$nipter_gc_correlation_post_gc_correction, NA_real_)),
    nipter_gc_curve_png = .null_coalesce(nipter_qc$gc_curve_plot, NA_character_),
    nipter_gc_curve_data_bgz = .null_coalesce(nipter_qc$gc_curve_data_bgz, NA_character_),
    nipter_corrected_bin_ratio_genome_plot_png = .null_coalesce(nipter_qc$corrected_bin_ratio_genome_plot, NA_character_),
    seqff_pre = unname(.null_coalesce(seqff$seqff_pre, NA_real_)),
    enet_pre = unname(.null_coalesce(seqff$enet_pre, NA_real_)),
    wrsc_pre = unname(.null_coalesce(seqff$wrsc_pre, NA_real_)),
    fetal_fraction_pre = unname(.null_coalesce(seqff$seqff_pre, NA_real_)),
    seqff_post = unname(.null_coalesce(seqff$seqff_post, NA_real_)),
    enet_post = unname(.null_coalesce(seqff$enet_post, NA_real_)),
    wrsc_post = unname(.null_coalesce(seqff$wrsc_post, NA_real_)),
    fetal_fraction_post = unname(.null_coalesce(seqff$seqff_post, NA_real_)),
    y_unique_ratio_pre = unname(.null_coalesce(y_unique$ratio_pre, NA_real_)),
    y_unique_reads_pre = unname(.null_coalesce(y_unique$reads_pre, NA_real_)),
    total_nuclear_reads_pre = unname(.null_coalesce(y_unique$total_pre, NA_real_)),
    y_unique_ratio_post = unname(.null_coalesce(y_unique$ratio_post, NA_real_)),
    y_unique_reads_post = unname(.null_coalesce(y_unique$reads_post, NA_real_)),
    total_nuclear_reads_post = unname(.null_coalesce(y_unique$total_post, NA_real_)),
    y_unique_regions_file = .null_coalesce(y_unique$regions_file, NA_character_),
    metrics_pre_mapq = unname(.null_coalesce(.metric_scalar(filters_pre, c("mapq"), "numeric"), NA_real_)),
    metrics_pre_require_flags = unname(.null_coalesce(.metric_scalar(filters_pre, c("require_flags"), "numeric"), NA_real_)),
    metrics_pre_exclude_flags = unname(.null_coalesce(.metric_scalar(filters_pre, c("exclude_flags"), "numeric"), NA_real_)),
    metrics_post_mapq = unname(.null_coalesce(.metric_scalar(filters_post, c("mapq"), "numeric"), NA_real_)),
    metrics_post_require_flags = unname(.null_coalesce(.metric_scalar(filters_post, c("require_flags"), "numeric"), NA_real_)),
    metrics_post_exclude_flags = unname(.null_coalesce(.metric_scalar(filters_post, c("exclude_flags"), "numeric"), NA_real_)),
    seqff_source = .null_coalesce(seqff$source, NA_character_),
    y_unique_source = .null_coalesce(y_unique$source, NA_character_),
    read_counts_source = .null_coalesce(read_counts$source, NA_character_),
    nipter_bed_status = NA_character_,
    nipter_bed_written = NA,
    nipter_bed_error = NA_character_
  )

  samtools_stats_row <- as.data.frame(
    c(
      list(
        samtools_stats_source = .null_coalesce(bam_stats$source, NA_character_),
        samtools_stats_pre_bookkeeping_mode = .null_coalesce(
          bam_stats$pre_bookkeeping_mode,
          NA_character_
        )
      ),
      stats::setNames(
        lapply(.samtools_stats_summary_fields, function(field) {
          unname(.null_coalesce(bam_stats[[paste0("pre_", field)]], NA_real_))
        }),
        paste0("samtools_stats_pre_", .samtools_stats_summary_fields)
      ),
      list(
        samtools_stats_post_bookkeeping_mode = .null_coalesce(
          bam_stats$post_bookkeeping_mode,
          NA_character_
        )
      ),
      stats::setNames(
        lapply(.samtools_stats_summary_fields, function(field) {
          unname(.null_coalesce(bam_stats[[paste0("post_", field)]], NA_real_))
        }),
        paste0("samtools_stats_post_", .samtools_stats_summary_fields)
      )
    ),
    check.names = FALSE
  )

  row <- data.frame(row, samtools_stats_row, check.names = FALSE)

  .validate_sample_qc_row(row)
}

.sample_qc_row_with_processing_status <- function(row,
                                                  status = NULL,
                                                  failure_stage = NULL,
                                                  failure_step = NULL,
                                                  failure_message = NULL) {
  stopifnot(is.data.frame(row), nrow(row) == 1L)

  schema <- .sample_qc_schema()
  missing_cols <- setdiff(schema, names(row))
  if (length(missing_cols)) {
    for (nm in missing_cols) {
      row[[nm]] <- NA
    }
  }
  row <- .validate_sample_qc_row(row)

  if (!is.null(status)) {
    status <- match.arg(as.character(status)[[1L]], c("ok", "failed"))
    row$sample_processing_status <- status
    if (identical(status, "ok")) {
      row$sample_failure_stage <- NA_character_
      row$sample_failure_step <- NA_character_
      row$sample_failure_message <- NA_character_
    }
  }
  if (!is.null(failure_stage)) {
    stage <- as.character(failure_stage)[[1L]]
    row$sample_failure_stage <- if (is.na(stage) || !nzchar(stage)) NA_character_ else stage
  }
  if (!is.null(failure_step)) {
    step <- as.character(failure_step)[[1L]]
    row$sample_failure_step <- if (is.na(step) || !nzchar(step)) NA_character_ else step
  }
  if (!is.null(failure_message)) {
    msg <- as.character(failure_message)[[1L]]
    row$sample_failure_message <- if (is.na(msg) || !nzchar(msg)) NA_character_ else msg
  }

  if (identical(row$sample_processing_status[[1L]], "failed") &&
      (is.na(row$sample_failure_stage[[1L]]) || !nzchar(row$sample_failure_stage[[1L]]))) {
    stop("Failed sample QC rows must record sample_failure_stage.", call. = FALSE)
  }
  if (identical(row$sample_processing_status[[1L]], "failed") &&
      (is.na(row$sample_failure_message[[1L]]) || !nzchar(row$sample_failure_message[[1L]]))) {
    stop("Failed sample QC rows must record sample_failure_message.", call. = FALSE)
  }

  .validate_sample_qc_row(
    row,
    require_on_success = identical(row$sample_processing_status[[1L]], "ok") &&
      isTRUE(row$nipter_bed_written[[1L]])
  )
}

.sample_qc_row_with_nipter_status <- function(row,
                                              written = NULL,
                                              error = NULL,
                                              status = NULL) {
  stopifnot(is.data.frame(row), nrow(row) == 1L)

  schema <- .sample_qc_schema()
  missing_cols <- setdiff(schema, names(row))
  if (length(missing_cols)) {
    for (nm in missing_cols) {
      row[[nm]] <- NA
    }
  }
  row <- .validate_sample_qc_row(row)

  if (!is.null(written)) {
    row$nipter_bed_written <- as.logical(written)[[1L]]
  }
  if (!is.null(error)) {
    err <- as.character(error)[[1L]]
    row$nipter_bed_error <- if (is.na(err) || !nzchar(err)) NA_character_ else err
  }

  if (!is.null(status)) {
    status <- match.arg(as.character(status)[[1L]], c("ok", "failed", "not_run"))
    row$nipter_bed_status <- status
    if (identical(status, "ok")) {
      row$nipter_bed_written <- TRUE
      row$nipter_bed_error <- NA_character_
    } else {
      row$nipter_bed_written <- FALSE
      if (identical(status, "not_run")) {
        row$nipter_bed_error <- NA_character_
      }
    }
  } else if (isTRUE(row$nipter_bed_written[[1L]])) {
    row$nipter_bed_status <- "ok"
    row$nipter_bed_error <- NA_character_
  } else if (identical(row$nipter_bed_written[[1L]], FALSE) ||
             (!is.na(row$nipter_bed_error[[1L]]) &&
                nzchar(row$nipter_bed_error[[1L]]))) {
    row$nipter_bed_status <- "failed"
  }

  .validate_sample_qc_row(
    row,
    require_on_success = identical(row$sample_processing_status[[1L]], "ok") &&
      isTRUE(row$nipter_bed_written[[1L]])
  )
}

.sample_qc_tabix_metadata <- function(seqff = NULL,
                                      y_unique = NULL,
                                      read_counts = NULL,
                                      bam_stats = NULL,
                                      coverage = NULL,
                                      nipter_qc = NULL,
                                      gc_correction = NULL,
                                      filters_pre = NULL,
                                      filters_post = NULL) {
  seqff <- .null_coalesce(seqff, list())
  y_unique <- .null_coalesce(y_unique, list())
  read_counts <- .null_coalesce(read_counts, list())
  bam_stats <- .samtools_stats_metric_record(record = bam_stats)
  coverage <- .coverage_metric_record(record = coverage)
  nipter_qc <- .null_coalesce(nipter_qc, list())
  gc_correction <- if (is.null(gc_correction)) NULL else .gc_correction_record(record = gc_correction)
  filters_pre <- .null_coalesce(.record_to_list(filters_pre), list())
  filters_post <- .null_coalesce(.record_to_list(filters_post), list())

  meta <- list(
    metrics_pre_mapq = .metric_scalar(filters_pre, c("mapq"), "numeric"),
    metrics_pre_require_flags = .metric_scalar(filters_pre, c("require_flags"), "numeric"),
    metrics_pre_exclude_flags = .metric_scalar(filters_pre, c("exclude_flags"), "numeric"),
    metrics_post_mapq = .metric_scalar(filters_post, c("mapq"), "numeric"),
    metrics_post_require_flags = .metric_scalar(filters_post, c("require_flags"), "numeric"),
    metrics_post_exclude_flags = .metric_scalar(filters_post, c("exclude_flags"), "numeric"),
    seqff_pre = seqff$seqff_pre,
    enet_pre = seqff$enet_pre,
    wrsc_pre = seqff$wrsc_pre,
    fetal_fraction_pre = seqff$seqff_pre,
    seqff_post = seqff$seqff_post,
    enet_post = seqff$enet_post,
    wrsc_post = seqff$wrsc_post,
    fetal_fraction_post = seqff$seqff_post,
    seqff_source = seqff$source,
    y_unique_ratio_pre = y_unique$ratio_pre,
    y_unique_reads_pre = y_unique$reads_pre,
    total_nuclear_reads_pre = y_unique$total_pre,
    y_unique_ratio_post = y_unique$ratio_post,
    y_unique_reads_post = y_unique$reads_post,
    total_nuclear_reads_post = y_unique$total_post,
    y_unique_source = y_unique$source,
    y_unique_regions_file = y_unique$regions_file,
    samtools_stats_source = bam_stats$source,
    samtools_stats_pre_bookkeeping_mode = bam_stats$pre_bookkeeping_mode,
    samtools_stats_post_bookkeeping_mode = bam_stats$post_bookkeeping_mode,
    gc_loess_valid_bins = nipter_qc$gc_loess_valid_bins,
    gc_loess_total_bins = nipter_qc$gc_loess_total_bins,
    gc_loess_invalid_bins = nipter_qc$gc_loess_invalid_bins,
    gc_loess_valid_bin_fraction_pct = nipter_qc$gc_loess_valid_bin_fraction_pct,
    gc_loess_invalid_out_of_chromosome_range = nipter_qc$gc_loess_invalid_out_of_chromosome_range,
    gc_loess_invalid_gc_missing_or_nonfinite = nipter_qc$gc_loess_invalid_gc_missing_or_nonfinite,
    gc_loess_invalid_gc_nonpositive = nipter_qc$gc_loess_invalid_gc_nonpositive,
    gc_loess_invalid_raw_nonfinite = nipter_qc$gc_loess_invalid_raw_nonfinite,
    gc_loess_invalid_raw_nonpositive = nipter_qc$gc_loess_invalid_raw_nonpositive,
    gc_loess_invalid_corrected_nonfinite = nipter_qc$gc_loess_invalid_corrected_nonfinite,
    gc_curve_has_valid_bins = nipter_qc$gc_curve_has_valid_bins,
    gc_curve_has_loess_support = nipter_qc$gc_curve_has_loess_support,
    nipter_autosomal_bin_cv_pre_gc_correction = nipter_qc$nipter_autosomal_bin_cv_pre_gc_correction,
    nipter_autosomal_bin_cv_post_gc_correction = nipter_qc$nipter_autosomal_bin_cv_post_gc_correction,
    nipter_gc_correction_bin_scale_mean = nipter_qc$nipter_gc_correction_bin_scale_mean,
    nipter_gc_correction_bin_scale_sd = nipter_qc$nipter_gc_correction_bin_scale_sd,
    nipter_gc_correction_bin_scale_cv = nipter_qc$nipter_gc_correction_bin_scale_cv,
    nipter_gc_correlation_pre_gc_correction = nipter_qc$nipter_gc_correlation_pre_gc_correction,
    nipter_gc_correlation_post_gc_correction = nipter_qc$nipter_gc_correlation_post_gc_correction,
    nipter_corrected_bin_ratio_genome_plot_png = nipter_qc$corrected_bin_ratio_genome_plot,
    read_counts_input_total = read_counts$input_total,
    read_counts_mapped_total = read_counts$mapped_total,
    read_counts_unmapped_total = read_counts$unmapped_total,
    gc_read_perc_pre = read_counts$gc_pre,
    gc_read_perc_post = read_counts$gc_post,
    mean_mapq_post = read_counts$mean_mapq_post,
    read_counts_binned_pre_sum = read_counts$count_pre,
    read_counts_binned_post_sum = read_counts$count_post,
    read_counts_binned_retained_fraction = read_counts$retained_fraction,
    read_counts_binned_fwd_sum = read_counts$count_fwd,
    read_counts_binned_rev_sum = read_counts$count_rev,
    read_counts_binned_nonzero_bins_post = read_counts$n_nonzero_bins_post,
    read_counts_source = read_counts$source,
    coverage_binsize = coverage$binsize,
    coverage_pre_mapq = coverage$pre_mapq,
    coverage_pre_exclude_flags = coverage$pre_exclude_flags,
    coverage_pre_summary_txt = coverage$pre_summary_txt,
    coverage_pre_length_bases = coverage$pre_length_bases,
    coverage_pre_bases = coverage$pre_bases,
    coverage_pre_covered_fraction = coverage$pre_covered_fraction,
    coverage_pre_mean_depth = coverage$pre_mean_depth,
    coverage_pre_min_depth = coverage$pre_min_depth,
    coverage_pre_max_depth = coverage$pre_max_depth,
    coverage_post_mapq = coverage$post_mapq,
    coverage_post_exclude_flags = coverage$post_exclude_flags,
    coverage_post_summary_txt = coverage$post_summary_txt,
    coverage_post_length_bases = coverage$post_length_bases,
    coverage_post_bases = coverage$post_bases,
    coverage_post_covered_fraction = coverage$post_covered_fraction,
    coverage_post_mean_depth = coverage$post_mean_depth,
    coverage_post_min_depth = coverage$post_min_depth,
    coverage_post_max_depth = coverage$post_max_depth,
    coverage_source = coverage$source
  )

  meta <- c(
    meta,
    stats::setNames(
      lapply(.samtools_stats_summary_fields, function(field) {
        bam_stats[[paste0("pre_", field)]]
      }),
      paste0("samtools_stats_pre_", .samtools_stats_summary_fields)
    ),
    stats::setNames(
      lapply(.samtools_stats_summary_fields, function(field) {
        bam_stats[[paste0("post_", field)]]
      }),
      paste0("samtools_stats_post_", .samtools_stats_summary_fields)
    )
  )

  if (!is.null(gc_correction)) {
    meta <- c(meta, list(
      nipter_gc_correction_applied = gc_correction$applied,
      nipter_gc_correction_method = gc_correction$method,
      nipter_gc_correction_loess_span = gc_correction$loess_span,
      nipter_gc_correction_include_sex = gc_correction$include_sex,
      nipter_gc_correction_binsize = gc_correction$binsize,
      nipter_gc_correction_table_bgz = gc_correction$table_bgz,
      nipter_gc_correction_fasta = gc_correction$fasta
    ))
  }

  meta <- meta[!vapply(meta, is.null, logical(1))]
  .normalize_tabix_metadata(meta)
}
