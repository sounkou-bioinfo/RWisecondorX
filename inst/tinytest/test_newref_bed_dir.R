# test_newref_bed_dir.R — Smoke test for rwisecondorx_newref(bed_dir=)
#
# Creates 10 copies of the fixture BAM's BED output, then verifies that
# rwisecondorx_newref() accepts bed_dir= and returns a WisecondorXReference.
# Uses nipt=TRUE so no male gonosomal reference is attempted (fixture is chr11
# only; fewer than 5 samples per gender would fail the gonosomal sub-ref).
#
# Skips cleanly if the fixture BAM is not present (make fixtures not yet run).

library(tinytest)
library(RWisecondorX)

fixture_bam <- system.file("extdata", "hg00106_chr11_fixture.bam",
                           package = "RWisecondorX")
if (!nzchar(fixture_bam)) {
  exit_file("hg00106_chr11_fixture.bam not available; skipping bed_dir tests")
}

# Write 10 BED.gz copies (minimum sample count for rwisecondorx_newref is 10)
bed_dir <- tempfile("newref_beds_")
dir.create(bed_dir)

bed_paths <- file.path(bed_dir, paste0("sample_", 1:10, ".bed.gz"))
for (p in bed_paths) {
  bam_convert_bed(fixture_bam, p, binsize = 5000L, rmdup = "streaming")
}

expect_true(all(file.exists(bed_paths)),
            info = "all 10 BED.gz files created")

# rwisecondorx_newref with bed_dir= should succeed and return a reference.
# yfrac is explicit because the fixture is chr11-only (no Y reads → all
# Y-fractions are zero → Mclust hangs on constant input).
ref <- rwisecondorx_newref(bed_dir = bed_dir, binsize = 5000L, nipt = TRUE,
                           refsize = 5L, yfrac = 0.001)

expect_true(inherits(ref, "WisecondorXReference"),
            info = "bed_dir= returns WisecondorXReference")
expect_equal(ref$binsize, 5000L,
             info = "reference binsize matches input")
expect_true(ref$is_nipt,
            info = "NIPT mode recorded in reference")
expect_false(ref$has_male,
             info = "no male gonosomal ref in NIPT mode")

# Supplying both samples= and bed_dir= must error
expect_error(
  rwisecondorx_newref(samples = list(), bed_dir = bed_dir),
  info = "samples + bed_dir together raises error"
)

# bed_dir with no matching files must error
empty_dir <- tempfile("empty_")
dir.create(empty_dir)
expect_error(
  rwisecondorx_newref(bed_dir = empty_dir, binsize = 5000L),
  info = "empty bed_dir raises error"
)
