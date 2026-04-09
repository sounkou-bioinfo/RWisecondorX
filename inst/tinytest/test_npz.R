# Tests for bam_convert_npz()
#
# These tests run when reticulate + numpy are available; they skip cleanly
# otherwise.  Full pipeline (convert → newref → predict) conformance is covered
# in test_integration.R which additionally requires condathis.

library(tinytest)

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
if (!requireNamespace("Rduckhts", quietly = TRUE)) {
  exit_file("Rduckhts not available")
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
keys <- reticulate::py_to_r(data$files)

expect_true(length(keys) > 0L, info = "NPZ contains at least one array")
expect_true(all(keys %in% as.character(1:24)),
            info = "All NPZ keys are valid chromosome numbers (1-24)")

# Each array should be int32 and non-negative
for (k in keys) {
  arr <- reticulate::py_to_r(data[k])
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

for (k in keys) {
  r_vals  <- bins[[k]]
  np_vals <- reticulate::py_to_r(data2[k])
  if (is.null(r_vals)) r_vals <- integer(length(np_vals))
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
