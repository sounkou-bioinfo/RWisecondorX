#!/usr/bin/env Rscript
# bench_500.R — End-to-end 500-sample pipeline benchmark with per-phase timing
#
# Usage: Rscript inst/scripts/bench_500.R [out_dir]

library(RWisecondorX)
library(mclust)
library(DNAcopy)

N <- 500L
out_dir <- if (length(commandArgs(TRUE)) > 0) {
  commandArgs(TRUE)[1]
} else {
  file.path(tempdir(), "bench_500")
}

cat(sprintf("=== RWisecondorX 500-sample benchmark ===\n"))
cat(sprintf("Output directory: %s\n", out_dir))
cat(sprintf("Samples: %d\n\n", N))

timings <- list()

# ---- Phase 1: Cohort generation -------------------------------------------
cat("Phase 1: Generating cohort...\n")
t0 <- proc.time()
manifest <- generate_cohort(out_dir, n_samples = N, verbose = FALSE)
t1 <- proc.time()
timings$cohort_gen <- (t1 - t0)["elapsed"]
cat(sprintf("  %d BAMs generated in %.1f s\n", nrow(manifest), timings$cohort_gen))
cat(sprintf("  Females: %d  Males: %d  Trisomy: %d\n",
            sum(manifest$sex == "F"), sum(manifest$sex == "M"),
            sum(manifest$trisomy != "none")))

# ---- Phase 2: Binning -----------------------------------------------------
cat("\nPhase 2: Binning all samples with bam_convert()...\n")
t0 <- proc.time()
all_samples <- vector("list", nrow(manifest))
names(all_samples) <- manifest$sample_id

for (i in seq_len(nrow(manifest))) {
  bam_path <- file.path(out_dir, manifest$bam_file[i])
  all_samples[[i]] <- bam_convert(bam_path,
                                  binsize = COMPRESSED_BINSIZE,
                                  mapq = 0L,
                                  rmdup = "none")
  if (i %% 100L == 0L) cat(sprintf("  binned %d/%d\n", i, nrow(manifest)))
}
t1 <- proc.time()
timings$binning <- (t1 - t0)["elapsed"]
cat(sprintf("  %d samples binned in %.1f s (%.2f s/sample)\n",
            N, timings$binning, timings$binning / N))

# ---- Phase 3: Reference building ------------------------------------------
cat("\nPhase 3: Building reference with rwisecondorx_newref()...\n")
t0 <- proc.time()
ref <- rwisecondorx_newref(
  samples = all_samples,
  binsize = COMPRESSED_BINSIZE,
  nipt    = TRUE,
  refsize = 10L,
  cpus    = 4L
)
t1 <- proc.time()
timings$newref <- (t1 - t0)["elapsed"]
cat(sprintf("  Reference built in %.1f s\n", timings$newref))
cat(sprintf("  Masked bins: %d / %d\n", sum(ref$mask), length(ref$mask)))

# ---- Phase 4: Prediction on trisomy samples --------------------------------
cat("\nPhase 4: Predicting trisomy samples...\n")
trisomy_rows <- manifest[manifest$trisomy != "none", ]
t0 <- proc.time()
for (i in seq_len(nrow(trisomy_rows))) {
  row <- trisomy_rows[i, ]
  sample_data <- all_samples[[row$sample_id]]
  trisomy_chr <- sub("T", "", row$trisomy)

  pred <- suppressMessages(
    rwisecondorx_predict(
      sample    = sample_data,
      reference = ref,
      zscore    = 3,
      alpha     = 1e-4,
      seed      = 42L
    )
  )

  chr_z <- pred$results_z[[trisomy_chr]]
  mean_z <- if (!is.null(chr_z) && length(chr_z) > 0L) {
    mean(chr_z[chr_z != 0], na.rm = TRUE)
  } else NA_real_
  detected <- !is.na(mean_z) && mean_z > 3
  cat(sprintf("  %s: chr%s mean Z = %.2f %s\n",
              row$sample_id, trisomy_chr, mean_z,
              if (detected) "[DETECTED]" else "[missed]"))
}
t1 <- proc.time()
timings$predict_trisomy <- (t1 - t0)["elapsed"]
cat(sprintf("  3 trisomy predictions in %.1f s\n", timings$predict_trisomy))

# ---- Phase 5: Prediction on a euploid negative control ---------------------
cat("\nPhase 5: Euploid negative control...\n")
euploid_id <- manifest$sample_id[manifest$trisomy == "none"][1L]
t0 <- proc.time()
pred_euploid <- suppressMessages(
  rwisecondorx_predict(
    sample    = all_samples[[euploid_id]],
    reference = ref,
    zscore    = 3,
    alpha     = 1e-4,
    seed      = 42L
  )
)
t1 <- proc.time()
timings$predict_euploid <- (t1 - t0)["elapsed"]

chr_means <- vapply(pred_euploid$results_z[as.character(1:22)], function(z) {
  z <- z[z != 0]
  if (length(z) == 0) return(0)
  mean(z, na.rm = TRUE)
}, 0.0)
extreme <- names(chr_means)[abs(chr_means) > 3]
cat(sprintf("  %s: %s\n", euploid_id,
            if (length(extreme) == 0L) "CLEAN (no |Z|>3 autosomes)"
            else paste0("WARN: extreme chrs: ", paste(extreme, collapse = ","))))
cat(sprintf("  Euploid prediction in %.1f s\n", timings$predict_euploid))

# ---- Summary ---------------------------------------------------------------
cat("\n=== TIMING SUMMARY ===\n")
total <- 0
for (nm in names(timings)) {
  cat(sprintf("  %-20s %7.1f s\n", nm, timings[[nm]]))
  total <- total + timings[[nm]]
}
cat(sprintf("  %-20s %7.1f s\n", "TOTAL", total))
cat(sprintf("\nDisk usage: %.1f MB\n",
            sum(file.size(file.path(out_dir, list.files(out_dir)))) / 1024 / 1024))

# ---- Cleanup ---------------------------------------------------------------
cat("\nCleaning up...\n")
unlink(out_dir, recursive = TRUE)
cat("Done.\n")
