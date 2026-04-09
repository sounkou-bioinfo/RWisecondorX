# Integration / conformance tests for bam_convert()
#
# These tests compare bam_convert() output against the official WisecondorX
# Python implementation invoked via condathis (bioconda). All tests skip
# cleanly when condathis, reticulate, or a suitable BAM file are unavailable.
#
# Conformance reference:
#   ../../duckhts/scripts/wisecondorx_convert_conformance.py achieves exact
#   bin-for-bin agreement (0 mismatches across 25,115 non-zero bins).

# ---------------------------------------------------------------------------
# Helper: find a usable test BAM
# ---------------------------------------------------------------------------

.find_test_bam <- function() {
  candidates <- c(
    system.file("extdata", "hg00106_chr11_fixture.bam", package = "RWisecondorX"),
    Sys.getenv("WISECONDORX_TEST_BAM", unset = NA_character_)
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0L) return(NULL)
  candidates[[1L]]
}

# ---------------------------------------------------------------------------
# Helper: run official wisecondorx convert via condathis
# ---------------------------------------------------------------------------

.wisecondorx_ref_bins <- function(bam, binsize = 5000L) {
  if (!requireNamespace("condathis", quietly = TRUE)) return(NULL)

  npz_out <- tempfile(fileext = ".npz")
  on.exit(unlink(npz_out), add = TRUE)

  status <- tryCatch(
    condathis::run(
      "wisecondorx",
      args    = c("convert", bam, npz_out, "--binsize", as.character(binsize)),
      channel = "bioconda"
    ),
    error = function(e) NULL
  )

  if (is.null(status) || !file.exists(npz_out)) return(NULL)

  # Load npz with reticulate/numpy and extract bin counts
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  np <- tryCatch(reticulate::import("numpy", convert = FALSE), error = function(e) NULL)
  if (is.null(np)) return(NULL)

  data <- tryCatch(
    np$load(npz_out, allow_pickle = TRUE),
    error = function(e) NULL
  )
  if (is.null(data)) return(NULL)

  # npz keys are chromosome numbers as strings
  keys <- reticulate::py_to_r(data$files)
  bins <- lapply(stats::setNames(keys, keys), function(k) {
    reticulate::py_to_r(data[k])
  })
  data$close()
  bins
}

# ---------------------------------------------------------------------------
# Basic unit test: bam_convert() returns expected structure
# ---------------------------------------------------------------------------

test_bam <- .find_test_bam()

if (is.null(test_bam)) {
  exit_file("No test BAM available; set WISECONDORX_TEST_BAM to run convert tests")
}

if (!requireNamespace("Rduckhts", quietly = TRUE)) {
  exit_file("Rduckhts not available")
}

bins <- bam_convert(test_bam, binsize = 5000L, rmdup = "streaming")

expect_true(is.list(bins), info = "bam_convert returns a list")
expect_identical(names(bins), as.character(1:24), info = "result has keys 1-24")

non_null <- Filter(Negate(is.null), bins)

# If the BAM has no human chromosomes (e.g. the bundled range.bam which uses
# CHROMOSOME_I/II/III/IV), there is nothing more to test for WisecondorX.
if (length(non_null) == 0L) {
  exit_file("Test BAM has no chr1-22/X/Y reads; set WISECONDORX_TEST_BAM to a human BAM")
}

# All non-null entries should be non-negative integer vectors
for (chr in names(non_null)) {
  v <- non_null[[chr]]
  expect_true(is.integer(v), info = paste("chr", chr, "is integer"))
  expect_true(all(v >= 0L), info = paste("chr", chr, "is non-negative"))
}

# rmdup="none" should return >= as many reads as rmdup="streaming"
bins_none    <- bam_convert(test_bam, binsize = 5000L, rmdup = "none")
total_streaming <- sum(unlist(Filter(Negate(is.null), bins)))
total_none      <- sum(unlist(Filter(Negate(is.null), bins_none)))
expect_true(total_none >= total_streaming,
            info = "rmdup=none keeps >= reads than rmdup=streaming")

# ---------------------------------------------------------------------------
# Conformance test: compare against official wisecondorx via condathis
# ---------------------------------------------------------------------------

ref_bins <- .wisecondorx_ref_bins(test_bam, binsize = 5000L)

if (is.null(ref_bins)) {
  exit_file("condathis/wisecondorx not available; skipping conformance test")
}

for (chr in intersect(names(bins), names(ref_bins))) {
  our_v   <- bins[[chr]]
  ref_v   <- ref_bins[[chr]]
  if (is.null(our_v) && is.null(ref_v)) next
  if (is.null(our_v)) our_v <- integer(length(ref_v))
  if (is.null(ref_v)) ref_v <- integer(length(our_v))

  # Lengths may differ by ±1 due to edge bins; compare the shared prefix
  n <- min(length(our_v), length(ref_v))
  expect_identical(our_v[seq_len(n)], ref_v[seq_len(n)],
                   info = paste("chr", chr, "bin counts match reference"))
}
