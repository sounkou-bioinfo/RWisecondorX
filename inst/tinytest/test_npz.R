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

.with_duckhts_con <- function(fun) {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  Rduckhts::rduckhts_load(con)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fun(con)
}

.bam_chr_lengths_test <- function(bam) {
  # tinytest evaluates top-level expressions one-by-one, so registering
  # on.exit() cleanup for a test-owned connection at file scope disconnects it
  # before later assertions run. Keep the full connection lifecycle inside this
  # helper and exercise the package's DuckDB-backed header path directly.
  .with_duckhts_con(function(con) RWisecondorX:::.bam_chr_lengths(con, bam))
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
# Exact NPZ shape + round-trip: native bins padded to upstream WisecondorX
# layout match what we wrote to NPZ
# ---------------------------------------------------------------------------

bins <- bam_convert(test_bam, binsize = 5000L, rmdup = "streaming")
chr_lengths <- .bam_chr_lengths_test(test_bam)
data2 <- np$load(npz_out, allow_pickle = TRUE)
sample2 <- reticulate::py_to_r(data2["sample"]$item())

for (k in keys) {
  chr_length <- chr_lengths[[k]]
  np_vals <- sample2[[k]]

  if (is.null(chr_length)) {
    expect_true(is.null(np_vals),
                info = paste("chr", k, "is absent from the BAM header and stays NULL"))
    next
  }

  expected_n_bins <- as.integer(as.double(chr_length) / 5000) + 1L
  native_vals <- bins[[k]]
  if (is.null(native_vals)) {
    native_vals <- integer(expected_n_bins)
  } else {
    native_vals <- as.integer(native_vals)
    expect_true(length(native_vals) <= expected_n_bins,
                info = paste("chr", k, "native dense vector does not exceed NPZ layout"))
    if (length(native_vals) < expected_n_bins) {
      native_vals <- c(native_vals, integer(expected_n_bins - length(native_vals)))
    }
  }

  expect_identical(length(np_vals), expected_n_bins,
                   info = paste("chr", k, "NPZ length matches WisecondorX header quirk"))
  expect_identical(as.integer(np_vals), native_vals,
                   info = paste("chr", k, "NPZ payload matches padded native bins"))
}
data2$close()
unlink(npz_out)

# ---------------------------------------------------------------------------
# CRAM path: bam_convert_npz() accepts a reference FASTA
# ---------------------------------------------------------------------------

test_cram <- .fixture_path("fixture_mixed.cram")
test_ref <- .fixture_path("fixture_ref.fa")

if (!is.null(test_cram) && !is.null(test_ref)) {
  npz_cram <- tempfile(fileext = ".npz")

  bam_convert_npz(
    bam = test_cram,
    reference = test_ref,
    npz = npz_cram,
    binsize = 5000L,
    rmdup = "streaming",
    np = np
  )

  expect_true(file.exists(npz_cram), info = "bam_convert_npz supports CRAM inputs with reference")
  unlink(npz_cram)
}

# ---------------------------------------------------------------------------
# bam_convert_npz forwards native filter arguments to bam_convert()
# ---------------------------------------------------------------------------

filter_bam <- .fixture_path("fixture_mixed.bam")
if (!is.null(filter_bam)) {
  npz_filtered <- tempfile(fileext = ".npz")

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
  unlink(npz_filtered)
}
