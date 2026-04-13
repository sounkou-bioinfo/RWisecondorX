library(tinytest)

.fixture_path <- function(name) {
  candidates <- c(
    system.file("extdata", name, package = "RWisecondorX"),
    file.path("inst", "extdata", name),
    file.path("..", "..", "inst", "extdata", name)
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0L) return(NULL)
  candidates[[1L]]
}

.source_candidate <- function(path) {
  candidates <- c(path, file.path("..", "..", path))
  candidates[file.exists(candidates)][1L]
}

src_convert <- .source_candidate("R/convert.R")

if (!is.na(src_convert)) {
  source(src_convert)
} else if (requireNamespace("RWisecondorX", quietly = TRUE)) {
  bam_convert <- getExportedValue("RWisecondorX", "bam_convert")
} else {
  stop("Unable to locate R/convert.R for tests.", call. = FALSE)
}

if (!requireNamespace("Rduckhts", quietly = TRUE)) {
  exit_file("Rduckhts not available")
}

paired_bam <- .fixture_path("fixture_paired.bam")
single_bam <- .fixture_path("fixture_single.bam")
mixed_bam <- .fixture_path("fixture_mixed.bam")
mixed_cram <- .fixture_path("fixture_mixed.cram")
fixture_ref <- .fixture_path("fixture_ref.fa")

if (any(vapply(list(paired_bam, single_bam, mixed_bam, mixed_cram, fixture_ref), is.null, logical(1)))) {
  exit_file("Synthetic BAM/CRAM fixtures not available; run `make fixtures`")
}

paired_none <- bam_convert(paired_bam, binsize = 5000L, rmdup = "none")
paired_streaming <- bam_convert(paired_bam, binsize = 5000L, rmdup = "streaming")
paired_flag <- bam_convert(paired_bam, binsize = 5000L, rmdup = "flag")

expect_identical(sum(paired_none[["11"]]), 4L, info = "paired fixture keeps all reads with rmdup=none")
expect_identical(sum(paired_streaming[["11"]]), 2L, info = "paired fixture drops duplicate pair with rmdup=streaming")
expect_identical(sum(paired_flag[["11"]]), 4L, info = "paired fixture unchanged with rmdup=flag")

single_none <- bam_convert(single_bam, binsize = 5000L, rmdup = "none")
single_streaming <- bam_convert(single_bam, binsize = 5000L, rmdup = "streaming")
single_flag <- bam_convert(single_bam, binsize = 5000L, rmdup = "flag")

expect_identical(sum(single_none[["11"]]), 2L, info = "single-end fixture keeps both reads with rmdup=none")
expect_identical(sum(single_streaming[["11"]]), 1L, info = "single-end fixture drops consecutive duplicate with rmdup=streaming")
expect_identical(sum(single_flag[["11"]]), 2L, info = "single-end fixture unchanged with rmdup=flag")

mixed_none <- bam_convert(mixed_bam, binsize = 5000L, rmdup = "none")
mixed_streaming <- bam_convert(mixed_bam, binsize = 5000L, rmdup = "streaming")
mixed_flag <- bam_convert(mixed_bam, binsize = 5000L, rmdup = "flag")

expect_identical(sum(mixed_none[["11"]]), 10L, info = "mixed fixture total with rmdup=none")
expect_identical(sum(mixed_streaming[["11"]]), 5L, info = "mixed fixture total with rmdup=streaming")
expect_identical(sum(mixed_flag[["11"]]), 8L, info = "mixed fixture total with rmdup=flag")

expect_identical(as.integer(mixed_none[["11"]][1:3]), c(4L, 2L, 2L),
                 info = "mixed fixture bins for rmdup=none")
expect_identical(as.integer(mixed_streaming[["11"]][1:3]), c(2L, 2L, 1L),
                 info = "mixed fixture bins for rmdup=streaming")
expect_identical(as.integer(mixed_flag[["11"]][1:3]), c(4L, 0L, 2L),
                 info = "mixed fixture bins for rmdup=flag")

mixed_cram_streaming <- bam_convert(
  mixed_cram,
  reference = fixture_ref,
  binsize = 5000L,
  rmdup = "streaming"
)

expect_identical(mixed_cram_streaming[["11"]], mixed_streaming[["11"]],
                 info = "CRAM fixture matches BAM fixture when reference is supplied")
