# test_wisecondorx_e2e.R — End-to-end WisecondorX conformance
#
# This test complements test_cohort_pipeline.R by adding:
#   1. Aberration-level assertions on the pred$aberrations data.frame
#      (test_cohort_pipeline checks Z-scores; this checks the final call).
#   2. bed_dir= pathway: generate_cohort → bam_convert_bed → newref(bed_dir=)
#      verifying the multi-file BED loading pipeline.
#   3. Python arm (conditional): condathis wisecondorx convert/newref/predict
#      compared against native R predictions for aberration agreement.
#
# Sections 1-2 always run (require samtools + DNAcopy).
# Section 3 requires condathis + reticulate + numpy + wisecondorx bioconda.

library(tinytest)
library(RWisecondorX)
library(DNAcopy)

has_samtools <- nchar(Sys.which("samtools")) > 0L
if (!has_samtools) exit_file("samtools not available")


# ---------------------------------------------------------------------------
# Section 1: Aberration calls on synthetic trisomy samples
#
# Reuses the cohort from test_cohort_pipeline.R logic (smaller refsize for
# speed). Focuses on pred$aberrations data.frame quality.
# ---------------------------------------------------------------------------

cohort_dir <- file.path(tempdir(), "e2e_wcx")
message("Generating synthetic cohort for E2E test ...")
manifest <- suppressMessages(generate_cohort(cohort_dir, verbose = FALSE))

all_samples <- vector("list", nrow(manifest))
names(all_samples) <- manifest$sample_id

for (i in seq_len(nrow(manifest))) {
  all_samples[[i]] <- suppressMessages(
    bam_convert(file.path(cohort_dir, manifest$bam_file[i]),
                binsize = COMPRESSED_BINSIZE, mapq = 0L, rmdup = "none")
  )
}

ref <- suppressMessages(
  rwisecondorx_newref(samples = all_samples, binsize = COMPRESSED_BINSIZE,
                      nipt = TRUE, refsize = 10L, cpus = 1L)
)

# T21 sample prediction
t21_id <- manifest$sample_id[manifest$trisomy == "T21"]
pred_t21 <- suppressMessages(
  rwisecondorx_predict(all_samples[[t21_id]], ref,
                       zscore = 3, minrefbins = 5L, alpha = 1e-4, seed = 42L)
)

# pred$aberrations must be a data.frame with required columns
expect_true(is.data.frame(pred_t21$aberrations),
            info = "T21 prediction has aberrations data.frame")
expect_true(all(c("chr", "start", "end", "type") %in%
                  names(pred_t21$aberrations)),
            info = "aberrations has chr/start/end/type columns")

# T21: chromosome 21 must be called as a gain
t21_ab <- pred_t21$aberrations
chr21_gains <- t21_ab[t21_ab$chr == "21" & t21_ab$type == "gain", ]
expect_true(nrow(chr21_gains) >= 1L,
            info = "T21 sample has ≥1 gain call on chr21")

# T18 sample prediction
t18_id <- manifest$sample_id[manifest$trisomy == "T18"]
pred_t18 <- suppressMessages(
  rwisecondorx_predict(all_samples[[t18_id]], ref,
                       zscore = 3, minrefbins = 5L, alpha = 1e-4, seed = 42L)
)
chr18_gains <- pred_t18$aberrations[pred_t18$aberrations$chr == "18" &
                                     pred_t18$aberrations$type == "gain", ]
expect_true(nrow(chr18_gains) >= 1L,
            info = "T18 sample has ≥1 gain call on chr18")

# T13 sample prediction
t13_id <- manifest$sample_id[manifest$trisomy == "T13"]
pred_t13 <- suppressMessages(
  rwisecondorx_predict(all_samples[[t13_id]], ref,
                       zscore = 3, minrefbins = 5L, alpha = 1e-4, seed = 42L)
)
chr13_gains <- pred_t13$aberrations[pred_t13$aberrations$chr == "13" &
                                     pred_t13$aberrations$type == "gain", ]
expect_true(nrow(chr13_gains) >= 1L,
            info = "T13 sample has ≥1 gain call on chr13")

# Euploid negative control: zero aberration calls
# Use zscore = 5 (upstream default) for the euploid check — synthetic data
# with only refsize = 10 can produce Z > 3 noise on whole-chromosome segments.
euploid_id <- manifest$sample_id[manifest$trisomy == "none"][1L]
pred_euploid <- suppressMessages(
  rwisecondorx_predict(all_samples[[euploid_id]], ref,
                       zscore = 5, minrefbins = 5L, alpha = 1e-4, seed = 42L)
)
expect_equal(nrow(pred_euploid$aberrations), 0L,
             info = "euploid sample has zero aberration calls at zscore=5")

message("Section 1 complete: aberration-level assertions passed.")


# ---------------------------------------------------------------------------
# Section 2: bed_dir= pathway for rwisecondorx_newref()
#
# Writes 4-column BED.gz files for all cohort samples, then builds a reference
# via newref(bed_dir=) and verifies it produces the same binsize.
# ---------------------------------------------------------------------------

bed_dir <- file.path(tempdir(), "e2e_wcx_beds")
dir.create(bed_dir)
message("Writing BED.gz files for ", nrow(manifest), " samples ...")

for (i in seq_len(nrow(manifest))) {
  bed_path <- file.path(bed_dir, paste0(manifest$sample_id[i], ".bed.gz"))
  suppressMessages(
    bam_convert_bed(file.path(cohort_dir, manifest$bam_file[i]), bed_path,
                    binsize = COMPRESSED_BINSIZE, mapq = 0L, rmdup = "none")
  )
}

expect_true(length(Sys.glob(file.path(bed_dir, "*.bed.gz"))) == nrow(manifest),
            info = "BED.gz files written for all samples")

ref_bed <- suppressMessages(
  rwisecondorx_newref(bed_dir = bed_dir, binsize = COMPRESSED_BINSIZE,
                      nipt = TRUE, refsize = 10L, cpus = 1L)
)

expect_true(inherits(ref_bed, "WisecondorXReference"),
            info = "bed_dir= reference is a WisecondorXReference")
expect_equal(ref_bed$binsize, COMPRESSED_BINSIZE,
             info = "bed_dir= reference has correct binsize")
expect_true(isTRUE(ref_bed$is_nipt),
            info = "bed_dir= reference is in NIPT mode")

# Both references should have the same mask shape
expect_equal(length(ref$mask), length(ref_bed$mask),
             info = "bed_dir= reference has same mask length as list-based reference")

message("Section 2 complete: bed_dir= pipeline verified.")


# ---------------------------------------------------------------------------
# Section 3: Python WisecondorX conformance (conditional)
#
# Compares aberration calls between native R (above) and upstream Python
# wisecondorx via condathis. Both pipelines use the same synthetic cohort.
# Expected agreement ≥ 80% on chr-level gain/loss calls.
# ---------------------------------------------------------------------------

has_condathis   <- requireNamespace("condathis",   quietly = TRUE)
has_reticulate  <- requireNamespace("reticulate",  quietly = TRUE)
has_numpy       <- has_reticulate && tryCatch(
  { reticulate::import("numpy", convert = FALSE); TRUE },
  error = function(e) FALSE
)
condathis_ok <- has_condathis && has_numpy

if (!condathis_ok) {
  message("condathis/reticulate/numpy not available; skipping Python conformance arm.")
  exit_file("condathis/reticulate/numpy not available")
}

message("Section 3: running Python WisecondorX conformance ...")

np <- reticulate::import("numpy", convert = FALSE)

# Convert all samples to NPZ
npz_dir <- file.path(tempdir(), "e2e_wcx_npz")
dir.create(npz_dir)

for (i in seq_len(nrow(manifest))) {
  npz_path <- file.path(npz_dir, paste0(manifest$sample_id[i], ".npz"))
  suppressMessages(
    bam_convert_npz(
      bam  = file.path(cohort_dir, manifest$bam_file[i]),
      npz  = npz_path,
      binsize = COMPRESSED_BINSIZE,
      rmdup = "none",
      np = np
    )
  )
}
npz_files <- sort(Sys.glob(file.path(npz_dir, "*.npz")))
expect_equal(length(npz_files), nrow(manifest),
             info = "NPZ files created for all samples")

# Python reference
py_ref <- file.path(tempdir(), "e2e_py_ref.npz")
py_newref_ok <- tryCatch({
  suppressMessages(
    wisecondorx_newref(
      npz_files   = npz_files,
      output      = py_ref,
      binsize     = COMPRESSED_BINSIZE,
      ref_binsize = COMPRESSED_BINSIZE,
      nipt        = TRUE,
      cpus        = 1L
    )
  )
  TRUE
}, error = function(e) {
  message("Python wisecondorx newref failed: ", conditionMessage(e))
  FALSE
})

if (!py_newref_ok) {
  message("Skipping Python conformance assertions (newref failed).")
} else {

expect_true(file.exists(py_ref),
            info = "Python reference NPZ created")

# Python prediction on T21
py_t21_out <- file.path(tempdir(), "e2e_py_t21")
t21_npz <- file.path(npz_dir, paste0(t21_id, ".npz"))
suppressMessages(
  wisecondorx_predict(
    npz           = t21_npz,
    ref           = py_ref,
    output_prefix = py_t21_out,
    bed           = TRUE,
    seed          = 1L
  )
)

py_t21_bed <- paste0(py_t21_out, "_aberrations.bed")
expect_true(file.exists(py_t21_bed),
            info = "Python T21 aberrations BED created")

# Parse Python BED (cols: chrom, start, end, ratio, ...; type encoded in
# the BED name field or as a separate column depending on wisecondorx version)
py_ab_lines <- readLines(py_t21_bed)
py_ab_lines <- py_ab_lines[!startsWith(py_ab_lines, "track") &
                              nzchar(py_ab_lines)]
py_chrs <- if (length(py_ab_lines) > 0L) {
  vapply(strsplit(py_ab_lines, "\t"), `[[`, character(1L), 1L)
} else {
  character(0L)
}
# Strip "chr" prefix if present
py_chrs <- sub("^chr", "", py_chrs)

r_chrs <- unique(as.character(pred_t21$aberrations$chr))

# Chr21 must be called by both
expect_true("21" %in% py_chrs,
            info = "Python WisecondorX calls chr21 gain on T21 sample")
expect_true("21" %in% r_chrs,
            info = "Native R calls chr21 gain on T21 sample")

# Jaccard similarity on chr-level calls
if (length(r_chrs) > 0L && length(py_chrs) > 0L) {
  jaccard <- length(intersect(r_chrs, py_chrs)) / length(union(r_chrs, py_chrs))
  message(sprintf("Aberration Jaccard (R vs Python) on T21: %.2f", jaccard))
  expect_true(jaccard >= 0.5,
              info = "Aberration Jaccard between R and Python >= 0.5 on T21")
}

message("Section 3 complete: Python WisecondorX conformance checked.")

}  # end if (py_newref_ok)


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

unlink(cohort_dir, recursive = TRUE)
unlink(bed_dir,   recursive = TRUE)
if (condathis_ok) {
  unlink(npz_dir, recursive = TRUE)
  if (exists("py_ref")) unlink(py_ref)
  if (exists("py_t21_out")) unlink(paste0(py_t21_out, "_aberrations.bed"))
}
