library(tinytest)
library(RWisecondorX)
library(mclust)
library(DNAcopy)

.test_threads <- function(default = 4L) {
  threads <- suppressWarnings(
    as.integer(Sys.getenv("THREADS", unset = as.character(default)))
  )
  if (is.na(threads) || threads < 1L) default else threads
}

# ---------------------------------------------------------------------------
# End-to-end pipeline test using synthetic BAM cohort
#
# Generates 50 BAMs via generate_cohort(), bins them with bam_convert(),
# builds a reference with rwisecondorx_newref(), and runs predictions on
# the trisomy samples with rwisecondorx_predict().
#
# Trisomy simulation uses the Nguyen et al. 2023 fetal fraction model:
# reads are removed from non-target chromosomes with probability
# p = 0.5*f/(1+0.5*f) to simulate relative enrichment of the trisomy chr.
#
# We use fetal_fraction = 0.30 (30%), which is intentionally high but still
# biologically plausible (late third-trimester samples). Combined with
# reads_per_bin = 10 (~300k total reads), this gives sufficient signal-to-noise
# for the WisecondorX normalization pipeline to detect trisomy. At clinical
# depth (~10M reads), f = 0.05-0.10 is detectable.
# ---------------------------------------------------------------------------

test_ref_cpus <- .test_threads()
test_fetal_fraction <- 0.30
test_reads_per_bin <- 10
test_gc_bias_strength <- 1.0

# ---- Generate cohort ------------------------------------------------------

cohort_dir <- file.path(tempdir(), "cohort_pipeline_test")
on.exit(unlink(cohort_dir, recursive = TRUE), add = TRUE)
message("Generating synthetic BAM cohort into ", cohort_dir)
manifest <- suppressMessages(
  generate_cohort(cohort_dir, verbose = TRUE,
                  fetal_fraction = test_fetal_fraction,
                  reads_per_bin = test_reads_per_bin,
                  gc_bias_strength = test_gc_bias_strength)
)

expect_equal(nrow(manifest), 50L,
             info = "cohort has 50 samples")
expect_true(all(file.exists(file.path(cohort_dir, manifest$bam_file))),
            info = "all BAM files exist")
expect_true(all(file.exists(paste0(file.path(cohort_dir, manifest$bam_file), ".bai"))),
            info = "all BAM indexes exist")

# Verify manifest composition
expect_equal(sum(manifest$sex == "F"), 38L,
             info = "38 female samples (35 euploid + 3 trisomy)")
expect_equal(sum(manifest$sex == "M"), 12L,
             info = "12 male samples")
expect_equal(sum(manifest$trisomy != "none"), 3L,
             info = "3 trisomy samples")

# Verify fetal fraction is recorded for trisomy samples
trisomy_ff <- manifest$fetal_fraction[manifest$trisomy != "none"]
expect_true(all(!is.na(trisomy_ff)),
            info = "trisomy samples have non-NA fetal fraction")
expect_true(all(abs(trisomy_ff - test_fetal_fraction) < 1e-10),
            info = "trisomy fetal fractions match requested value")

# Verify trisomy samples have FEWER reads on average than euploid samples
# (Nguyen et al. model removes reads from non-target chromosomes)
euploid_f_reads <- manifest$n_reads[manifest$sex == "F" &
                                     manifest$trisomy == "none"]
trisomy_reads <- manifest$n_reads[manifest$trisomy != "none"]
expect_true(mean(trisomy_reads) < mean(euploid_f_reads),
            info = paste0("trisomy samples have fewer mean total reads than ",
                          "euploid females (Nguyen et al. removal model)"))


# ---- Bin all samples with bam_convert() -----------------------------------

message("Binning ", nrow(manifest), " BAMs with bam_convert() ...")
all_samples <- vector("list", nrow(manifest))
names(all_samples) <- manifest$sample_id

for (i in seq_len(nrow(manifest))) {
  bam_path <- file.path(cohort_dir, manifest$bam_file[i])
  all_samples[[i]] <- bam_convert(bam_path,
                                  binsize = COMPRESSED_BINSIZE,
                                  mapq = 0L,
                                  rmdup = "none")
  if (i %% 10L == 0L) message("  binned ", i, "/", nrow(manifest))
}

# Verify bin counts
n_chrs <- length(all_samples[[1]])
expect_true(n_chrs >= 23L,
            info = "binned samples have at least 23 chromosomes")

# Total reads per sample should be > 0
reads_per_sample <- vapply(all_samples, function(s) sum(vapply(s, sum, 0L)), 0L)
expect_true(all(reads_per_sample > 0L),
            info = "all samples have non-zero read counts")

# Female samples should have zero (or near-zero) Y-chr reads
female_euploid_idx <- which(manifest$sex == "F" & manifest$trisomy == "none")
y_counts_female <- vapply(all_samples[female_euploid_idx],
                          function(s) sum(s[["Y"]] %||% s[["24"]] %||% 0L),
                          0L)
expect_true(all(y_counts_female == 0L),
            info = "euploid female samples have zero Y-chr reads")


# ---- Build reference with rwisecondorx_newref() ---------------------------

message("Building reference with rwisecondorx_newref() ...")
ref <- rwisecondorx_newref(
    samples = all_samples,
    binsize = COMPRESSED_BINSIZE,
    nipt    = TRUE,
    refsize = 10L,
    cpus    = test_ref_cpus
  )

expect_true(is.list(ref),
            info = "newref returns a list-like reference object")
expect_true(S7::S7_inherits(ref, WisecondorXReference),
            info = "newref returns the typed WisecondorXReference S7 class")
expect_true(isTRUE(ref$is_nipt),
            info = "reference is in NIPT mode")
expect_equal(ref$binsize, COMPRESSED_BINSIZE,
             info = "reference binsize matches compressed binsize")
expect_true(!is.null(ref$mask),
            info = "reference has a bin mask")
expect_true(!is.null(ref$pca_components),
            info = "reference has PCA components")
expect_true(!is.null(ref$indexes),
            info = "reference has KNN indexes")
expect_true(!is.null(ref$null_ratios),
            info = "reference has null ratios")

message("Reference built successfully. Masked bins: ", sum(ref$mask))


# ---- Predict on trisomy samples -------------------------------------------

trisomy_samples <- manifest[manifest$trisomy != "none", ]
message("Running predictions on ", nrow(trisomy_samples), " trisomy samples ...")

for (i in seq_len(nrow(trisomy_samples))) {
  row <- trisomy_samples[i, ]
  sample_data <- all_samples[[row$sample_id]]
  trisomy_chr <- sub("T", "", row$trisomy)  # "T21" -> "21"

  pred <- suppressMessages(
    rwisecondorx_predict(
      sample    = sample_data,
      reference = ref,
      zscore    = 3,
      minrefbins = 5L,
      alpha     = 1e-4,
      seed      = 42L,
      parallel  = FALSE
    )
  )

  expect_true(is.list(pred),
              info = paste0(row$sample_id, ": predict returns a list"))
  expect_true(!is.null(pred$results_z),
              info = paste0(row$sample_id, ": prediction has Z-scores"))
  expect_true(!is.null(pred$results_r),
              info = paste0(row$sample_id, ": prediction has ratios"))

  # The trisomy chromosome should be detected via segment-level Z-scores.
  # WisecondorX's within-sample normalization absorbs whole-chromosome shifts
  # at the bin level, but CBS segmentation + null ratio comparison correctly
  # identifies the trisomy chromosome as a whole-chromosome aberration.
  #
  # We check two things:
  # 1. The per-chromosome statistics Z-score (from compute_statistics) > 3
  # 2. The trisomy chr Z-score exceeds the control chr Z-score

  chr_z <- pred$results_z[[trisomy_chr]]
  expect_true(!is.null(chr_z) && length(chr_z) > 0L,
              info = paste0(row$sample_id, ": chr", trisomy_chr,
                            " Z-score vector is non-NULL and non-empty"))

  # Per-chromosome statistics Z-score: computed from whole-chromosome segments
  # vs null ratio distributions. This is the primary trisomy detection metric.
  tri_chr_int <- as.integer(trisomy_chr)
  tri_stat <- pred$statistics[pred$statistics$chr == trisomy_chr, ]
  expect_true(nrow(tri_stat) == 1L,
              info = paste0(row$sample_id, ": statistics row for chr",
                            trisomy_chr, " exists"))

  stat_z <- tri_stat$zscore[1L]
  message(sprintf("  %s: chr%s stat_zscore = %.2f",
                  row$sample_id, trisomy_chr, stat_z))
  expect_true(stat_z > 3,
              info = paste0(row$sample_id, ": chr", trisomy_chr,
                            " statistics Z-score > 3 (trisomy detected)"))

  # Control chromosome: the Nguyen et al. removal model reduces reads on ALL
  # non-trisomy chromosomes, so their segment Z-scores may also be elevated.
  # The key distinction is the ratio DIRECTION: trisomy chr should have a
  # positive ratio (relative gain) while the control chr should have near-zero
  # or negative ratio (relative loss from read removal).
  control_chr <- if (trisomy_chr != "1") "1" else "2"
  ctrl_stat <- pred$statistics[pred$statistics$chr == control_chr, ]
  expect_true(nrow(ctrl_stat) == 1L,
              info = paste0(row$sample_id, ": statistics row for control chr",
                            control_chr, " exists"))

  tri_ratio <- tri_stat$ratio_mean[1L]
  ctrl_ratio <- ctrl_stat$ratio_mean[1L]
  expect_true(tri_ratio > ctrl_ratio,
              info = paste0(row$sample_id, ": trisomy chr", trisomy_chr,
                            " ratio (", round(tri_ratio, 4),
                            ") > control chr", control_chr,
                            " ratio (", round(ctrl_ratio, 4), ")"))

  # Check output writing works
  tmpdir <- tempfile(pattern = "wcx_out_")
  dir.create(tmpdir)
  outprefix <- file.path(tmpdir, row$sample_id)

  expect_true(exists("write_wisecondorx_output", mode = "function"),
              info = "write_wisecondorx_output is exported")
  result_path <- write_wisecondorx_output(pred, outprefix)
  expect_true(file.exists(paste0(outprefix, "_bins.bed")),
              info = paste0(row$sample_id, ": bins BED written"))
  expect_true(file.exists(paste0(outprefix, "_segments.bed")),
              info = paste0(row$sample_id, ": segments BED written"))

  unlink(tmpdir, recursive = TRUE)
}


# ---- Predict on a euploid sample (negative control) -----------------------

euploid_id <- manifest$sample_id[manifest$trisomy == "none"][1L]
message("Running euploid negative control: ", euploid_id)
pred_euploid <- suppressMessages(
  rwisecondorx_predict(
    sample    = all_samples[[euploid_id]],
    reference = ref,
    zscore    = 3,
    minrefbins = 5L,
    alpha     = 1e-4,
    seed      = 42L,
    parallel  = FALSE
  )
)

# Euploid sample should have no chromosome with extreme statistics Z-scores
euploid_stats <- pred_euploid$statistics
auto_stats <- euploid_stats[euploid_stats$chr %in% as.character(1:22), ]
extreme_chrs <- auto_stats$chr[abs(auto_stats$zscore) > 3]
expect_equal(length(extreme_chrs), 0L,
             info = paste0("euploid sample has no autosomes with |stat Z| > 3",
                           if (length(extreme_chrs) > 0)
                             paste0(" (offenders: ", paste(extreme_chrs, collapse = ","), ")")))

message("Cohort pipeline test complete.")
