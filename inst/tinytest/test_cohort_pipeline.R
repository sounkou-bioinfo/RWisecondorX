library(tinytest)
library(RWisecondorX)
library(mclust)
library(DNAcopy)

.test_threads <- function(default = 4L) {
  threads <- suppressWarnings(as.integer(Sys.getenv("THREADS", unset = as.character(default))))
  if (is.na(threads) || threads < 1L) {
    default
  } else {
    max(default, threads)
  }
}

# ---------------------------------------------------------------------------
# End-to-end pipeline test using synthetic BAM cohort
#
# Generates 50 BAMs via generate_cohort(), bins them with bam_convert(),
# builds a reference with rwisecondorx_newref(), and runs predictions on
# the trisomy samples with rwisecondorx_predict().
# ---------------------------------------------------------------------------

test_ref_cpus <- .test_threads()

# ---- Generate cohort ------------------------------------------------------

cohort_dir <- file.path(tempdir(), "cohort_pipeline_test")
message("Generating synthetic BAM cohort into ", cohort_dir)
manifest <- suppressMessages(
  generate_cohort(cohort_dir, verbose = TRUE)
)

expect_equal(nrow(manifest), 50L,
             info = "cohort has 50 samples")
expect_true(all(file.exists(file.path(cohort_dir, manifest$bam_file))),
            info = "all BAM files exist")

# Verify manifest composition
expect_equal(sum(manifest$sex == "F"), 38L,
             info = "38 female samples (35 euploid + 3 trisomy)")
expect_equal(sum(manifest$sex == "M"), 12L,
             info = "12 male samples")
expect_equal(sum(manifest$trisomy != "none"), 3L,
             info = "3 trisomy samples")


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

expect_true(inherits(ref, "WisecondorXReference") || is.list(ref),
            info = "newref returns a reference object")
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

  # The trisomy chromosome should have elevated Z-scores
  # (positive because it has 1.5x the expected reads)
  chr_z <- pred$results_z[[trisomy_chr]]
  if (!is.null(chr_z) && length(chr_z) > 0L) {
    mean_z <- mean(chr_z[chr_z != 0], na.rm = TRUE)
    message(sprintf("  %s: chr%s mean Z = %.2f",
                    row$sample_id, trisomy_chr, mean_z))
    expect_true(mean_z > 0,
                info = paste0(row$sample_id, ": chr", trisomy_chr,
                              " mean Z-score is positive (trisomy signal)"))
  }

  # Check a non-trisomy autosome has Z near zero
  control_chr <- if (trisomy_chr != "1") "1" else "2"
  chr_z_ctrl <- pred$results_z[[control_chr]]
  if (!is.null(chr_z_ctrl) && length(chr_z_ctrl) > 0L) {
    mean_z_ctrl <- mean(chr_z_ctrl[chr_z_ctrl != 0], na.rm = TRUE)
    expect_true(abs(mean_z_ctrl) < abs(mean_z),
                info = paste0(row$sample_id, ": control chr", control_chr,
                              " Z is smaller than trisomy chr"))
  }

  # Check output writing works
  tmpdir <- tempfile(pattern = "wcx_out_")
  dir.create(tmpdir)
  outprefix <- file.path(tmpdir, row$sample_id)

  if (exists("write_wisecondorx_output", mode = "function")) {
    result_path <- write_wisecondorx_output(pred, outprefix)
    expect_true(file.exists(paste0(outprefix, "_bins.bed")),
                info = paste0(row$sample_id, ": bins BED written"))
    expect_true(file.exists(paste0(outprefix, "_segments.bed")),
                info = paste0(row$sample_id, ": segments BED written"))
  }

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

# Euploid sample should have no chromosome with extreme Z-scores
chr_means <- vapply(pred_euploid$results_z[as.character(1:22)], function(z) {
  z <- z[z != 0]
  if (length(z) == 0) return(0)
  mean(z, na.rm = TRUE)
}, 0.0)

# No autosome should exceed |Z| > 3 on average for a euploid sample
extreme_chrs <- names(chr_means)[abs(chr_means) > 3]
expect_equal(length(extreme_chrs), 0L,
             info = paste0("euploid sample has no autosomes with |mean Z| > 3",
                           if (length(extreme_chrs) > 0)
                             paste0(" (offenders: ", paste(extreme_chrs, collapse = ","), ")")))

message("Cohort pipeline test complete.")


# ---- Cleanup ---------------------------------------------------------------

unlink(cohort_dir, recursive = TRUE)
