library(tinytest)
library(RWisecondorX)

# ---------------------------------------------------------------------------
# CombinedStrands import modes from pre-fetched BED rows
# ---------------------------------------------------------------------------

combined_rows <- data.frame(
  chrom = c("1", "1", "X", "X"),
  start_pos = c(0L, 50000L, 0L, 50000L),
  end_pos = c(50000L, 100000L, 50000L, 100000L),
  count = c(10L, 20L, 5L, 7L),
  corrected_count = c(11, 22, 6, 8),
  stringsAsFactors = FALSE
)

combined_sample <- RWisecondorX:::.rows_to_nipter_combined(
  combined_rows,
  "combined",
  binsize = 50000L,
  autosomal_source = "corrected",
  sex_source = "raw"
)

expect_true(S7::S7_inherits(combined_sample, CombinedStrandsSample),
            info = "combined import mode returns CombinedStrandsSample")
expect_identical(combined_sample$correction_status_autosomal, "GC Corrected",
                 info = "combined import mode marks autosomes as GC-corrected")
expect_identical(combined_sample$correction_status_sex, "Uncorrected",
                 info = "combined import mode keeps sex raw")
expect_equal(unname(combined_sample$autosomal_chromosome_reads[[1L]]["1", 1:2]),
             c(11, 22),
             info = "combined import mode uses corrected autosomal values")
expect_identical(as.integer(combined_sample$sex_chromosome_reads[[1L]]["X", 1:2]),
                 c(5L, 7L),
                 info = "combined import mode uses raw sex values")


# ---------------------------------------------------------------------------
# SeparatedStrands import modes from pre-fetched BED rows
# ---------------------------------------------------------------------------

separated_rows <- data.frame(
  chrom = c("1", "1", "X", "X"),
  start_pos = c(0L, 50000L, 0L, 50000L),
  end_pos = c(50000L, 100000L, 50000L, 100000L),
  count = c(10L, 20L, 5L, 7L),
  count_fwd = c(6L, 11L, 3L, 4L),
  count_rev = c(4L, 9L, 2L, 3L),
  corrected_count = c(11, 22, 6, 8),
  corrected_fwd = c(7, 12, 4, 5),
  corrected_rev = c(4.5, 10, 2, 3.5),
  stringsAsFactors = FALSE
)

separated_sample <- RWisecondorX:::.rows_to_nipter_sep(
  separated_rows,
  "separated",
  binsize = 50000L,
  autosomal_source = "corrected",
  sex_source = "raw"
)

expect_true(S7::S7_inherits(separated_sample, SeparatedStrandsSample),
            info = "separated import mode returns SeparatedStrandsSample")
expect_identical(separated_sample$correction_status_autosomal, "GC Corrected",
                 info = "separated import mode marks autosomes as GC-corrected")
expect_identical(separated_sample$correction_status_sex, "Uncorrected",
                 info = "separated import mode keeps sex raw")
expect_equal(unname(separated_sample$autosomal_chromosome_reads[[1L]]["1F", 1:2]),
             c(7, 12),
             info = "separated import mode uses corrected forward autosomes")
expect_equal(unname(separated_sample$autosomal_chromosome_reads[[2L]]["1R", 1:2]),
             c(4.5, 10),
             info = "separated import mode uses corrected reverse autosomes")
expect_identical(as.integer(separated_sample$sex_chromosome_reads[[1L]]["XF", 1:2]),
                 c(3L, 4L),
                 info = "separated import mode uses raw forward sex values")
expect_identical(as.integer(separated_sample$sex_chromosome_reads[[2L]]["XR", 1:2]),
                 c(2L, 3L),
                 info = "separated import mode uses raw reverse sex values")


# ---------------------------------------------------------------------------
# Corrected import should error cleanly when corrected columns are absent
# ---------------------------------------------------------------------------

raw_only_rows <- data.frame(
  chrom = c("1", "X"),
  start_pos = c(0L, 0L),
  end_pos = c(50000L, 50000L),
  count = c(10L, 5L),
  corrected_count = c(NA_real_, NA_real_),
  stringsAsFactors = FALSE
)

expect_error(
  RWisecondorX:::.rows_to_nipter_combined(
    raw_only_rows,
    "raw_only",
    binsize = 50000L,
    autosomal_source = "corrected"
  ),
  pattern = "do not contain corrected autosomal values",
  info = "combined corrected import errors when corrected columns are absent"
)
