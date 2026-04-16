library(tinytest)
library(RWisecondorX)

donor_bam <- system.file("extdata", "nipter_conformance_fixture.bam",
                         package = "RWisecondorX")
if (!nzchar(donor_bam)) {
  exit_file("Bundled NIPTeR conformance fixture is missing; run `make fixtures`.")
}

tmp_dir <- tempfile("bam_sim_", tmpdir = getwd())
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(tmp_dir)) {
  stop("Failed to create temporary BAM simulation directory.", call. = FALSE)
}

sim_t21 <- file.path(tmp_dir, "sim_t21.bam")
sim_t21_mosaic <- file.path(tmp_dir, "sim_t21_mosaic.bam")

res_t21 <- simulate_trisomy_bam(
  bam = donor_bam,
  out_bam = sim_t21,
  trisomy_chr = "chr21",
  fetal_fraction = 0.20,
  seed = 42L
)

res_t21_mosaic <- simulate_trisomy_bam(
  bam = donor_bam,
  out_bam = sim_t21_mosaic,
  trisomy_chr = "21",
  fetal_fraction = 0.20,
  mosaic = TRUE,
  seed = 42L
)

expect_true(file.exists(sim_t21), info = "simulate_trisomy_bam writes output BAM")
expect_true(file.exists(paste0(sim_t21, ".bai")),
            info = "simulate_trisomy_bam writes BAM index by default")
expect_true(file.exists(sim_t21_mosaic), info = "mosaic simulation writes output BAM")
expect_true(file.exists(paste0(sim_t21_mosaic, ".bai")),
            info = "mosaic simulation writes BAM index by default")

expect_true(res_t21$output_records < res_t21$input_records,
            info = "non-target thinning reduces total record count")
expect_true(res_t21_mosaic$output_records > res_t21$output_records,
            info = "mosaic simulation keeps more records than non-mosaic")
expect_equal(res_t21$target_chr, "21",
             info = "target chromosome is normalized in the result metadata")

orig_bins <- bam_convert(donor_bam, binsize = 50000L, rmdup = "none")
t21_bins <- bam_convert(sim_t21, binsize = 50000L, rmdup = "none")
t21_mosaic_bins <- bam_convert(sim_t21_mosaic, binsize = 50000L, rmdup = "none")

sum_chr <- function(x, chr) sum(x[[chr]])
sum_mapped <- function(x) sum(vapply(x[as.character(1:24)], sum, integer(1)))
sum_auto <- function(x) sum(vapply(x[as.character(1:22)], sum, integer(1)))

orig_chr21 <- sum_chr(orig_bins, "21")
t21_chr21 <- sum_chr(t21_bins, "21")
t21_mosaic_chr21 <- sum_chr(t21_mosaic_bins, "21")

expect_identical(t21_chr21, orig_chr21,
                 info = "target chromosome read count is preserved in non-mosaic simulation")
expect_identical(t21_mosaic_chr21, orig_chr21,
                 info = "target chromosome read count is preserved in mosaic simulation")

expect_true(sum_mapped(t21_bins) < sum_mapped(orig_bins),
            info = "non-mosaic simulation reduces mapped read count")
expect_true(sum_mapped(t21_mosaic_bins) > sum_mapped(t21_bins),
            info = "mosaic simulation retains more mapped reads than non-mosaic")

orig_frac21 <- orig_chr21 / sum_auto(orig_bins)
t21_frac21 <- t21_chr21 / sum_auto(t21_bins)
t21_mosaic_frac21 <- t21_mosaic_chr21 / sum_auto(t21_mosaic_bins)

expect_true(t21_frac21 > orig_frac21,
            info = "non-mosaic simulation increases chr21 autosomal fraction")
expect_true(t21_mosaic_frac21 > orig_frac21,
            info = "mosaic simulation increases chr21 autosomal fraction")
expect_true(t21_mosaic_frac21 < t21_frac21,
            info = "mosaic simulation gives weaker chr21 enrichment than non-mosaic")

unlink(tmp_dir, recursive = TRUE)
