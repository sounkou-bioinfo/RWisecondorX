library(RWisecondorX)

bam <- system.file("extdata", "fixture_mixed.bam", package = "Rduckhts")
has_samtools <- nzchar(Sys.which("samtools"))

expect_true(file.exists(bam), info = "Rduckhts fixture_mixed.bam is available")

if (!file.exists(bam) || !has_samtools) {
  exit_file("samtools CLI and RWisecondorX fixture BAM are required for bam-stats conformance checks")
}

.parse_samtools_sn <- function(lines) {
  sn <- lines[grepl("^SN\t", lines)]
  if (!length(sn)) {
    stop("No SN rows found in samtools stats output.")
  }
  fields <- strsplit(sn, "\t", fixed = TRUE)
  labels <- vapply(fields, function(x) sub(":$", "", x[[2L]]), character(1L))
  values <- vapply(fields, function(x) suppressWarnings(as.numeric(x[[3L]])), numeric(1L))
  stats::setNames(values, labels)
}

.native <- bam_samtools_stats_summary(bam)
.cli <- .parse_samtools_sn(system2("samtools", c("stats", bam), stdout = TRUE, stderr = TRUE))

expect_equal(.native$raw_total_sequences, .cli[["raw total sequences"]], info = "raw total sequences matches samtools stats")
expect_equal(.native$filtered_sequences, .cli[["filtered sequences"]], info = "filtered sequences matches samtools stats")
expect_equal(.native$sequences, .cli[["sequences"]], info = "sequences matches samtools stats")
expect_equal(.native$first_fragments, .cli[["1st fragments"]], info = "1st fragments matches samtools stats")
expect_equal(.native$last_fragments, .cli[["last fragments"]], info = "last fragments matches samtools stats")
expect_equal(.native$other_fragments, .native$sequences - .native$first_fragments - .native$last_fragments, info = "other fragments are algebraically consistent")
expect_equal(.native$reads_mapped, .cli[["reads mapped"]], info = "reads mapped matches samtools stats")
expect_equal(.native$reads_mapped_and_paired, .cli[["reads mapped and paired"]], info = "reads mapped and paired matches samtools stats")
expect_equal(.native$reads_unmapped, .cli[["reads unmapped"]], info = "reads unmapped matches samtools stats")
expect_equal(.native$reads_properly_paired, .cli[["reads properly paired"]], info = "reads properly paired matches samtools stats")
expect_equal(.native$reads_paired, .cli[["reads paired"]], info = "reads paired matches samtools stats")
expect_equal(.native$reads_duplicated, .cli[["reads duplicated"]], info = "reads duplicated matches samtools stats")
expect_equal(.native$reads_mq0, .cli[["reads MQ0"]], info = "reads MQ0 matches samtools stats")
expect_equal(.native$reads_qc_failed, .cli[["reads QC failed"]], info = "reads QC failed matches samtools stats")
expect_equal(.native$non_primary_alignments, .cli[["non-primary alignments"]], info = "non-primary alignments matches samtools stats")
expect_equal(.native$supplementary_alignments, .cli[["supplementary alignments"]], info = "supplementary alignments matches samtools stats")
expect_equal(.native$pairs_on_different_chromosomes, .cli[["pairs on different chromosomes"]], info = "pairs on different chromosomes matches samtools stats")
expect_equal(.native$total_length, .cli[["total length"]], info = "total length matches samtools stats")
expect_equal(.native$total_first_fragment_length, .cli[["total first fragment length"]], info = "total first fragment length matches samtools stats")
expect_equal(.native$total_last_fragment_length, .cli[["total last fragment length"]], info = "total last fragment length matches samtools stats")
expect_equal(.native$bases_mapped, .cli[["bases mapped"]], info = "bases mapped matches samtools stats")
expect_equal(.native$bases_mapped_cigar, .cli[["bases mapped (cigar)"]], info = "bases mapped (cigar) matches samtools stats")
expect_equal(.native$bases_duplicated, .cli[["bases duplicated"]], info = "bases duplicated matches samtools stats")
expect_equal(.native$mismatches_from_nm, .cli[["mismatches"]], info = "mismatches matches samtools stats")
expect_equal(.native$error_rate, .cli[["error rate"]], tolerance = 1e-12, info = "error rate matches samtools stats")
expect_equal(.native$duplicated_read_fraction, .native$reads_duplicated / .native$sequences, tolerance = 1e-15, info = "duplicated read fraction is algebraically consistent")
expect_equal(.native$duplicated_base_fraction, .native$bases_duplicated / .native$total_length, tolerance = 1e-15, info = "duplicated base fraction is algebraically consistent")

sn_path <- tempfile(fileext = ".samtools_stats.txt")
write_bam_samtools_stats_summary(.native, sn_path)
.written <- .parse_samtools_sn(readLines(sn_path, warn = FALSE))
expect_equal(.written[["raw total sequences"]], .cli[["raw total sequences"]], info = "writer preserves raw total sequences")
expect_equal(.written[["1st fragments"]], .cli[["1st fragments"]], info = "writer preserves first fragments")
expect_equal(.written[["last fragments"]], .cli[["last fragments"]], info = "writer preserves last fragments")
expect_equal(.written[["pairs on different chromosomes"]], .cli[["pairs on different chromosomes"]], info = "writer preserves cross-chromosome pair count")
expect_equal(.written[["error rate"]], .cli[["error rate"]], tolerance = 1e-7, info = "writer preserves error-rate text at samtools-compatible precision")

.native_filtered <- bam_samtools_stats_summary(bam, min_mapq = 1L, exclude_flags = 1024L)
.native_filtered_stream <- bam_samtools_stats_summary(
  bam,
  min_mapq = 1L,
  exclude_flags = 1024L,
  report_filtered_stream_bookkeeping = TRUE
)
cmd_filtered <- sprintf("samtools view -u -q 1 -F 1024 %s | samtools stats -", shQuote(bam))
.cli_filtered <- .parse_samtools_sn(system2("bash", c("-lc", shQuote(cmd_filtered)), stdout = TRUE, stderr = TRUE))

expect_equal(.native_filtered$sequences, .cli_filtered[["sequences"]], info = "filtered sequences retained by the pipeline match")
expect_equal(.native_filtered$other_fragments, .native_filtered$sequences - .native_filtered$first_fragments - .native_filtered$last_fragments, info = "filtered other fragments are algebraically consistent")
expect_equal(.native_filtered_stream$raw_total_sequences, .cli_filtered[["raw total sequences"]], info = "filtered-stream bookkeeping matches samtools raw total sequences")
expect_equal(.native_filtered_stream$filtered_sequences, .cli_filtered[["filtered sequences"]], info = "filtered-stream bookkeeping matches samtools filtered sequences")
expect_equal(.native_filtered$first_fragments, .cli_filtered[["1st fragments"]], info = "filtered first fragments match pipeline samtools stats")
expect_equal(.native_filtered$last_fragments, .cli_filtered[["last fragments"]], info = "filtered last fragments match pipeline samtools stats")
expect_equal(.native_filtered$reads_mapped, .cli_filtered[["reads mapped"]], info = "filtered reads mapped matches pipeline samtools stats")
expect_equal(.native_filtered$reads_duplicated, .cli_filtered[["reads duplicated"]], info = "filtered reads duplicated matches pipeline samtools stats")
expect_equal(.native_filtered$reads_mq0, .cli_filtered[["reads MQ0"]], info = "filtered reads MQ0 match pipeline samtools stats")
expect_equal(.native_filtered$reads_qc_failed, .cli_filtered[["reads QC failed"]], info = "filtered QC-failed reads match pipeline samtools stats")
expect_equal(.native_filtered$non_primary_alignments, .cli_filtered[["non-primary alignments"]], info = "filtered non-primary alignments match pipeline samtools stats")
expect_equal(.native_filtered$supplementary_alignments, .cli_filtered[["supplementary alignments"]], info = "filtered supplementary alignments match pipeline samtools stats")
expect_equal(.native_filtered$pairs_on_different_chromosomes, .cli_filtered[["pairs on different chromosomes"]], info = "filtered cross-chromosome pairs match pipeline samtools stats")
expect_equal(.native_filtered$total_length, .cli_filtered[["total length"]], info = "filtered total length matches pipeline samtools stats")
expect_equal(.native_filtered$total_first_fragment_length, .cli_filtered[["total first fragment length"]], info = "filtered first-fragment length matches pipeline samtools stats")
expect_equal(.native_filtered$total_last_fragment_length, .cli_filtered[["total last fragment length"]], info = "filtered last-fragment length matches pipeline samtools stats")
expect_equal(.native_filtered$bases_mapped, .cli_filtered[["bases mapped"]], info = "filtered bases mapped match pipeline samtools stats")
expect_equal(.native_filtered$bases_mapped_cigar, .cli_filtered[["bases mapped (cigar)"]], info = "filtered bases mapped (cigar) matches pipeline samtools stats")
expect_equal(.native_filtered$bases_duplicated, .cli_filtered[["bases duplicated"]], info = "filtered bases duplicated match pipeline samtools stats")
expect_equal(.native_filtered$mismatches_from_nm, .cli_filtered[["mismatches"]], info = "filtered mismatches matches pipeline samtools stats")
expect_equal(.native_filtered$error_rate, .cli_filtered[["error rate"]], tolerance = 1e-12, info = "filtered error rate matches pipeline samtools stats")
expect_equal(.native_filtered$duplicated_read_fraction, .native_filtered$reads_duplicated / .native_filtered$sequences, tolerance = 1e-15, info = "filtered duplicated read fraction is algebraically consistent")
expect_equal(.native_filtered$duplicated_base_fraction, .native_filtered$bases_duplicated / .native_filtered$total_length, tolerance = 1e-15, info = "filtered duplicated base fraction is algebraically consistent")
