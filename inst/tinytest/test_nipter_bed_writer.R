library(tinytest)
library(RWisecondorX)

toy_binsize <- 5000L
toy_cols <- as.character(seq_len(4L))
toy_auto <- matrix(0L, nrow = 22L, ncol = 4L,
                   dimnames = list(as.character(1:22), toy_cols))
toy_sex <- matrix(0L, nrow = 2L, ncol = 4L,
                  dimnames = list(c("X", "Y"), toy_cols))
toy_auto["1", ] <- c(7L, 0L, 0L, 0L)
toy_auto["21", ] <- c(0L, 1L, 0L, 0L)

toy_no_lengths <- CombinedStrandsSample(
  sample_name = "toy_no_lengths",
  binsize = as.integer(toy_binsize),
  auto_matrix = toy_auto,
  sex_matrix_ = toy_sex
)
expect_error(
  nipter_sample_to_bed(toy_no_lengths, tempfile(fileext = ".bed.gz"), binsize = toy_binsize),
  pattern = "requires explicit chromosome lengths",
  info = "NIPTeR BED export refuses ambiguous padded matrices without chromosome lengths"
)

toy_lengths <- setNames(rep.int(as.integer(toy_binsize), 24L), c(as.character(1:22), "X", "Y"))
toy_lengths["1"] <- 12500L
toy_lengths["21"] <- 17500L
toy_with_lengths <- CombinedStrandsSample(
  sample_name = "toy_with_lengths",
  binsize = as.integer(toy_binsize),
  chrom_lengths = toy_lengths,
  auto_matrix = toy_auto,
  sex_matrix_ = toy_sex
)

toy_bed <- tempfile(fileext = ".bed.gz")
nipter_sample_to_bed(toy_with_lengths, toy_bed, binsize = toy_binsize)
toy_rt <- bed_to_nipter_sample(toy_bed)
expect_identical(
  unname(toy_rt$chrom_lengths[c("1", "21")]),
  unname(toy_lengths[c("1", "21")]),
  info = "NIPTeR BED round-trip preserves clipped chromosome lengths"
)

toy_lines <- readLines(gzfile(toy_bed))
toy_fields <- strsplit(toy_lines, "\t", fixed = TRUE)
toy_rows <- data.frame(
  chrom = vapply(toy_fields, `[[`, character(1L), 1L),
  start_pos = as.integer(vapply(toy_fields, `[[`, character(1L), 2L)),
  end_pos = as.integer(vapply(toy_fields, `[[`, character(1L), 3L)),
  count = as.integer(vapply(toy_fields, `[[`, character(1L), 4L)),
  stringsAsFactors = FALSE
)
toy_chr1 <- toy_rows[toy_rows[["chrom"]] == "1", ]
toy_chr21 <- toy_rows[toy_rows[["chrom"]] == "21", ]
expect_identical(nrow(toy_chr1), 3L,
                 info = "NIPTeR BED export emits only the required bins for chromosome 1")
expect_identical(max(toy_chr1[["end_pos"]]), 12500L,
                 info = "NIPTeR BED export clips chromosome 1 final interval to true length")
expect_identical(nrow(toy_chr21), 4L,
                 info = "NIPTeR BED export emits the required number of bins for chromosome 21")
expect_identical(max(toy_chr21[["end_pos"]]), 17500L,
                 info = "NIPTeR BED export clips chromosome 21 final interval to true length")
