#' Compute samtools-stats-style BAM summary numbers
#'
#' Runs a native htslib-backed summary pass over a BAM or CRAM file and returns
#' a subset of the summary-number (`SN`) fields produced by `samtools stats`.
#' The core mismatch, duplicate-flag, and CIGAR-denominator logic is ported
#' from upstream `samtools` `stats.c`.
#'
#' @param path Path to the input BAM or CRAM file.
#' @param reference Optional reference FASTA path for CRAM input when required.
#' @param min_mapq Integer; additional minimum MAPQ filter applied before the
#'   samtools-style summary logic. This is a package-level extension for
#'   preprocessing workflows; native `samtools stats` does not expose a direct
#'   MAPQ threshold option.
#' @param require_flags Integer; required SAM flag mask, equivalent in spirit
#'   to samtools `-f`.
#' @param exclude_flags Integer; excluded SAM flag mask, equivalent in spirit
#'   to samtools `-F`.
#' @param decompression_threads Integer; number of htslib decompression worker
#'   threads to use. `0` disables worker threads.
#' @param report_filtered_stream_bookkeeping Logical. When `FALSE` (default),
#'   `raw_total_sequences` stays anchored to the original BAM/CRAM input and
#'   `filtered_sequences` reports the number excluded by `min_mapq` /
#'   flag-based filtering. When `TRUE`, these two fields instead follow the
#'   semantics of `samtools view ... | samtools stats -`, where the retained
#'   stream is treated as the input and `filtered_sequences` is therefore `0`.
#'
#' @details
#' `error_rate` follows upstream `samtools stats` semantics exactly:
#' `mismatches_from_nm / bases_mapped_cigar`, where mismatches are accumulated
#' from the `NM` auxiliary tag and `bases_mapped_cigar` counts `M`, `I`, `=`,
#' and `X` CIGAR operations.
#'
#' `reads_duplicated` and the derived duplicate fractions are duplicate-flag
#' summaries. They are **not** Picard `PERCENT_DUPLICATION`.
#'
#' For exact CLI comparison when using `min_mapq`, compare against a pipeline of
#' the form `samtools view -u -q ... | samtools stats -`.
#'
#' @return A one-row `data.frame`.
#' @export
bam_samtools_stats_summary <- function(path,
                                       reference = NULL,
                                       min_mapq = 0L,
                                       require_flags = 0L,
                                       exclude_flags = 0L,
                                       decompression_threads = 0L,
                                       report_filtered_stream_bookkeeping = FALSE) {
  stopifnot(is.character(path), length(path) == 1L, nzchar(path), file.exists(path))

  if (is.null(reference)) {
    reference <- ""
  }
  stopifnot(is.character(reference), length(reference) == 1L)
  if (nzchar(reference)) {
    stopifnot(file.exists(reference))
  }

  out <- samtools_stats_summary_cpp(
    path = path,
    reference = reference,
    min_mapq = as.integer(min_mapq),
    require_flags = as.integer(require_flags),
    exclude_flags = as.integer(exclude_flags),
    decompression_threads = as.integer(decompression_threads),
    report_filtered_stream_bookkeeping = isTRUE(report_filtered_stream_bookkeeping)
  )

  data.frame(out, check.names = FALSE)
}

.samtools_stats_sn_schema <- local({
  data.frame(
    field = c(
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
      "error_rate"
    ),
    label = c(
      "raw total sequences",
      "filtered sequences",
      "sequences",
      "1st fragments",
      "last fragments",
      "reads mapped",
      "reads mapped and paired",
      "reads unmapped",
      "reads properly paired",
      "reads paired",
      "reads duplicated",
      "reads MQ0",
      "reads QC failed",
      "non-primary alignments",
      "supplementary alignments",
      "pairs on different chromosomes",
      "total length",
      "total first fragment length",
      "total last fragment length",
      "bases mapped",
      "bases mapped (cigar)",
      "bases duplicated",
      "mismatches",
      "error rate"
    ),
    comment = c(
      "# excluding supplementary and secondary reads",
      "",
      "",
      "",
      "",
      "",
      "# paired-end technology bit set + both mates mapped",
      "",
      "# proper-pair bit set",
      "# paired-end technology bit set",
      "# PCR or optical duplicate bit set",
      "# mapped and MQ=0",
      "",
      "",
      "",
      "",
      "# ignores clipping",
      "# ignores clipping",
      "# ignores clipping",
      "# ignores clipping",
      "# more accurate",
      "",
      "# from NM fields",
      "# mismatches / bases mapped (cigar)"
    ),
    check.names = FALSE
  )
})

.format_samtools_stats_sn_value <- function(field, value) {
  stopifnot(length(value) == 1L)
  if (!is.finite(value)) {
    stop("Cannot write non-finite value for samtools stats field: ", field, call. = FALSE)
  }
  if (identical(field, "error_rate")) {
    return(sprintf("%e", as.numeric(value)))
  }
  sprintf("%.0f", as.numeric(value))
}

.samtools_stats_sn_lines <- function(x) {
  if (!is.data.frame(x) || nrow(x) != 1L) {
    stop("x must be a one-row data.frame returned by bam_samtools_stats_summary().", call. = FALSE)
  }

  schema <- .samtools_stats_sn_schema
  missing_fields <- setdiff(schema$field, names(x))
  if (length(missing_fields)) {
    stop(
      "Cannot write samtools-stats summary; missing fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  lines <- c(
    "# This file contains the RWisecondorX samtools-stats-compatible SN summary section.",
    "# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part."
  )

  body <- vapply(seq_len(nrow(schema)), function(i) {
    field <- schema$field[[i]]
    label <- schema$label[[i]]
    comment <- schema$comment[[i]]
    value_txt <- .format_samtools_stats_sn_value(field, x[[field]][[1L]])
    line <- paste("SN", paste0(label, ":"), value_txt, sep = "\t")
    if (nzchar(comment)) {
      line <- paste(line, comment, sep = "\t")
    }
    line
  }, character(1L))

  c(lines, body)
}

#' Write a samtools-stats-compatible `SN` summary section
#'
#' Writes the implemented summary-number subset from
#' [bam_samtools_stats_summary()] in a `samtools stats`-compatible `SN` text
#' format.
#'
#' @param x One-row `data.frame` returned by [bam_samtools_stats_summary()].
#' @param file Output path for the text file.
#'
#' @details
#' This writer emits the `SN` section only. It does not attempt to reproduce the
#' full `samtools stats` output surface such as insert-size histograms, GC-depth,
#' or cycle-level tables.
#'
#' @return `file`, invisibly.
#' @export
write_bam_samtools_stats_summary <- function(x, file) {
  stopifnot(is.character(file), length(file) == 1L, nzchar(file))
  base::writeLines(.samtools_stats_sn_lines(x), con = file, sep = "\n", useBytes = TRUE)
  invisible(file)
}
