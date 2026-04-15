# Tests for bam_convert_npz()
#
# These tests run when reticulate + numpy are available; they skip cleanly
# otherwise.  Full pipeline (convert → newref → predict) conformance is covered
# in test_integration.R which additionally requires condathis.

library(tinytest)
library(RWisecondorX)

# ---------------------------------------------------------------------------
# Find a usable test BAM (same logic as test_integration.R)
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

.fixture_path <- function(name) {
  candidates <- c(
    system.file("extdata", name, package = "RWisecondorX"),
    file.path("inst", "extdata", name),
    file.path("..", "..", "inst", "extdata", name)
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0L) return(NULL)
  candidates[[1L]]
}

test_bam <- .find_test_bam()
if (is.null(test_bam)) {
  exit_file("No test BAM available")
}
Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")
if (!requireNamespace("reticulate", quietly = TRUE)) {
  exit_file("reticulate not available; skipping NPZ tests")
}
np <- tryCatch(reticulate::import("numpy", convert = FALSE), error = function(e) NULL)
if (is.null(np)) {
  exit_file("numpy not found in active Python environment; skipping NPZ tests")
}

# ---------------------------------------------------------------------------
# bam_convert_npz writes a readable NPZ file
# ---------------------------------------------------------------------------

# Pre-check: skip if BAM has no human chromosomes (e.g. bundled range.bam)
probe <- bam_convert(test_bam, binsize = 5000L, rmdup = "none")
if (length(Filter(Negate(is.null), probe)) == 0L) {
  exit_file("Test BAM has no chr1-22/X/Y reads; set WISECONDORX_TEST_BAM to a human BAM")
}

npz_out <- tempfile(fileext = ".npz")
on.exit(unlink(npz_out), add = TRUE)

bam_convert_npz(test_bam, npz_out, binsize = 5000L, rmdup = "streaming", np = np)

expect_true(file.exists(npz_out), info = "bam_convert_npz creates the output file")
expect_true(file.info(npz_out)$size > 0L, info = "output NPZ is not empty")

# Load back with numpy and verify structure
data <- np$load(npz_out, allow_pickle = TRUE)
npz_keys <- reticulate::py_to_r(data$files)

# Upstream format has three top-level keys: "sample", "binsize", "quality"
expect_true("sample" %in% npz_keys, info = "NPZ contains 'sample' key")
expect_true("binsize" %in% npz_keys, info = "NPZ contains 'binsize' key")
expect_true("quality" %in% npz_keys, info = "NPZ contains 'quality' key")

# binsize should round-trip correctly
bs <- reticulate::py_to_r(data["binsize"]$item())
expect_equal(as.integer(bs), 5000L, info = "binsize round-trips as 5000")

# sample should be a dict-like object extractable with .item()
sample <- reticulate::py_to_r(data["sample"]$item())
expect_true(is.list(sample), info = "sample.item() returns a list/dict")

keys <- names(sample)
expect_equal(length(keys), 24L, info = "sample dict has all 24 chromosome keys")
expect_true(setequal(keys, as.character(1:24)),
            info = "All sample keys are chromosome numbers 1-24")

# Each non-None entry should be int32 and non-negative; None entries are allowed
# (chromosomes without data in the BAM).
for (k in keys) {
  arr <- sample[[k]]
  if (is.null(arr)) next
  expect_true(is.integer(arr) || is.numeric(arr),
              info = paste("chr", k, "array is numeric"))
  expect_true(all(arr >= 0L),
              info = paste("chr", k, "values are non-negative"))
}
data$close()

# ---------------------------------------------------------------------------
# Round-trip: bam_convert() values match what we wrote to NPZ
# ---------------------------------------------------------------------------

bins <- bam_convert(test_bam, binsize = 5000L, rmdup = "streaming")
data2 <- np$load(npz_out, allow_pickle = TRUE)
sample2 <- reticulate::py_to_r(data2["sample"]$item())

for (k in keys) {
  r_vals  <- bins[[k]]
  np_vals <- sample2[[k]]
  if (is.null(r_vals) && is.null(np_vals)) next
  if (is.null(r_vals)) r_vals <- integer(length(np_vals))
  if (is.null(np_vals)) np_vals <- integer(length(r_vals))
  n <- min(length(r_vals), length(np_vals))
  expect_identical(as.integer(r_vals[seq_len(n)]), as.integer(np_vals[seq_len(n)]),
                   info = paste("chr", k, "round-trip identical"))
}
data2$close()

# ---------------------------------------------------------------------------
# CRAM path: bam_convert_npz() accepts a reference FASTA
# ---------------------------------------------------------------------------

test_cram <- .fixture_path("fixture_mixed.cram")
test_ref <- .fixture_path("fixture_ref.fa")

if (!is.null(test_cram) && !is.null(test_ref)) {
  npz_cram <- tempfile(fileext = ".npz")
  on.exit(unlink(npz_cram), add = TRUE)

  bam_convert_npz(
    bam = test_cram,
    reference = test_ref,
    npz = npz_cram,
    binsize = 5000L,
    rmdup = "streaming",
    np = np
  )

  expect_true(file.exists(npz_cram), info = "bam_convert_npz supports CRAM inputs with reference")
}

# ---------------------------------------------------------------------------
# bam_convert_npz forwards native filter arguments to bam_convert()
# ---------------------------------------------------------------------------

filter_bam <- .fixture_path("fixture_mixed.bam")
if (!is.null(filter_bam)) {
  npz_filtered <- tempfile(fileext = ".npz")
  on.exit(unlink(npz_filtered), add = TRUE)

  bam_convert_npz(
    bam = filter_bam,
    npz = npz_filtered,
    binsize = 50L,
    mapq = 0L,
    require_flags = 2L,
    exclude_flags = 1024L,
    rmdup = "none",
    np = np
  )

  filtered_bins <- bam_convert(
    filter_bam,
    binsize = 50L,
    mapq = 0L,
    require_flags = 2L,
    exclude_flags = 1024L,
    rmdup = "none"
  )

  filtered_npz <- np$load(npz_filtered, allow_pickle = TRUE)
  filtered_sample <- reticulate::py_to_r(filtered_npz["sample"]$item())

  for (k in names(filtered_sample)) {
    r_vals <- filtered_bins[[k]]
    np_vals <- filtered_sample[[k]]
    if (is.null(r_vals) && is.null(np_vals)) next
    if (is.null(r_vals)) r_vals <- integer(length(np_vals))
    if (is.null(np_vals)) np_vals <- integer(length(r_vals))
    n <- min(length(r_vals), length(np_vals))
    expect_identical(as.integer(r_vals[seq_len(n)]), as.integer(np_vals[seq_len(n)]),
                     info = paste("filtered NPZ round-trip identical for chr", k))
  }
  filtered_npz$close()
}
