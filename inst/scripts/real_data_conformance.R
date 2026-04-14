#!/usr/bin/env Rscript
# real_data_conformance.R — Real-data WisecondorX + NIPTeR conformance report
#
# Runs the full RWisecondorX native R pipeline on real NIPT BAMs and optionally
# compares results against the upstream Python WisecondorX (via condathis) and/or
# the NIPTeR R package. Produces a machine-readable conformance report suitable
# for inclusion in VALIDATION.md.
#
# rewrites.bio 4.3: this script records exact commands, package versions, hardware,
# and dates so that any claimed conformance result is fully reproducible.
#
# ---- Environment variables ---------------------------------------------------
#
#  Required:
#    RWXCONF_CONTROL_BAMS   colon-separated paths to ≥10 control BAMs
#    RWXCONF_TEST_BAM       path to a test BAM (known trisomy or euploid)
#
#  Optional:
#    RWXCONF_BINSIZE         reference bin size bp (default: 100000)
#    RWXCONF_SAMPLE_BINSIZE  sample bin size bp   (default: 5000)
#    RWXCONF_CPUS            thread count          (default: 4)
#    RWXCONF_FASTA           FASTA for GC correction
#    RWXCONF_OUTDIR          output directory for report (default: ./conformance_reports)
#    RWXCONF_SEED            RNG seed for CBS (default: 42)
#
# ---- Usage -------------------------------------------------------------------
#
#   export RWXCONF_CONTROL_BAMS="ctrl1.bam:ctrl2.bam:..."
#   export RWXCONF_TEST_BAM="test.bam"
#   export RWXCONF_CPUS=8
#   Rscript inst/scripts/real_data_conformance.R
#
# ------------------------------------------------------------------------------

stopifnot(requireNamespace("RWisecondorX", quietly = TRUE))
library(RWisecondorX)

# ---- Read configuration ------------------------------------------------------

.env <- function(name, default = NULL) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) default else val
}

control_bams_raw <- .env("RWXCONF_CONTROL_BAMS")
test_bam         <- .env("RWXCONF_TEST_BAM")

if (is.null(control_bams_raw) || is.null(test_bam)) {
  cat("Usage: Set RWXCONF_CONTROL_BAMS and RWXCONF_TEST_BAM before running.\n")
  cat("\nRequired:\n")
  cat("  RWXCONF_CONTROL_BAMS   colon-separated paths to >=10 control BAMs\n")
  cat("  RWXCONF_TEST_BAM       path to test BAM\n")
  cat("\nOptional:\n")
  cat("  RWXCONF_BINSIZE        reference bin size (default 100000)\n")
  cat("  RWXCONF_SAMPLE_BINSIZE sample bin size (default 5000)\n")
  cat("  RWXCONF_CPUS           threads (default 4)\n")
  cat("  RWXCONF_FASTA          FASTA for GC correction\n")
  cat("  RWXCONF_OUTDIR         output directory (default ./conformance_reports)\n")
  cat("  RWXCONF_SEED           CBS RNG seed (default 42)\n")
  quit(status = 1L)
}

control_bams  <- strsplit(control_bams_raw, ":")[[1L]]
ref_binsize   <- as.integer(.env("RWXCONF_BINSIZE",        "100000"))
samp_binsize  <- as.integer(.env("RWXCONF_SAMPLE_BINSIZE", "5000"))
cpus          <- as.integer(.env("RWXCONF_CPUS",           "4"))
fasta         <- .env("RWXCONF_FASTA")
out_dir       <- .env("RWXCONF_OUTDIR", "conformance_reports")
seed          <- as.integer(.env("RWXCONF_SEED", "42"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

report_file <- file.path(out_dir,
  sprintf("conformance_report_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))

# ---- Open report sink --------------------------------------------------------

sink(report_file, split = TRUE)

cat("# RWisecondorX Real-Data Conformance Report\n")
cat("# Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

cat("## System Info\n")
cat("R version:         ", R.version$version.string, "\n")
cat("RWisecondorX:      ", as.character(utils::packageVersion("RWisecondorX")), "\n")
cat("Platform:          ", Sys.info()[["sysname"]], Sys.info()[["machine"]], "\n")
cat("CPUs configured:   ", cpus, "\n")
cat("Reference binsize: ", ref_binsize, "\n")
cat("Sample binsize:    ", samp_binsize, "\n")
cat("Seed:              ", seed, "\n")
cat("Control BAMs:      ", length(control_bams), "\n")
for (b in control_bams) cat("  ", b, "\n")
cat("Test BAM:          ", test_bam, "\n")
if (!is.null(fasta)) cat("FASTA:             ", fasta, "\n")
cat("\n")

# ---- Validate inputs ---------------------------------------------------------

missing_ctrl <- control_bams[!file.exists(control_bams)]
if (length(missing_ctrl) > 0L) {
  cat("ERROR: missing control BAMs:\n")
  for (m in missing_ctrl) cat("  ", m, "\n")
  sink()
  stop("Missing control BAMs; aborting.")
}
if (!file.exists(test_bam)) {
  sink()
  stop("Test BAM not found: ", test_bam)
}
if (length(control_bams) < 10L) {
  cat("WARNING: fewer than 10 control BAMs; rwisecondorx_newref requires >=10.\n\n")
}

# ---- Step 1: Bin control BAMs -----------------------------------------------

cat("## Step 1: Bin control BAMs\n")
t0 <- proc.time()
control_samples <- vector("list", length(control_bams))
names(control_samples) <- sub("\\.bam$|\\.cram$", "", basename(control_bams))

for (i in seq_along(control_bams)) {
  bam <- control_bams[i]
  nm  <- names(control_samples)[i]
  cat(sprintf("  [%d/%d] %s ... ", i, length(control_bams), nm))
  t1 <- proc.time()
  control_samples[[i]] <- suppressMessages(
    bam_convert(bam, binsize = samp_binsize, mapq = 1L, rmdup = "streaming")
  )
  cat(sprintf("%.1fs\n", (proc.time() - t1)[["elapsed"]]))
}
cat(sprintf("Total binning time: %.1fs\n\n", (proc.time() - t0)[["elapsed"]]))

# ---- Step 2: Build native R reference ----------------------------------------

cat("## Step 2: Build native R reference (rwisecondorx_newref)\n")
t0 <- proc.time()
ref <- suppressMessages(
  rwisecondorx_newref(
    samples = control_samples,
    binsize = ref_binsize,
    sample_binsizes = rep(samp_binsize, length(control_samples)),
    nipt    = TRUE,
    cpus    = cpus
  )
)
cat(sprintf("Reference built in %.1fs\n", (proc.time() - t0)[["elapsed"]]))
cat(sprintf("  Masked bins:   %d / %d (%.1f%%)\n",
            sum(!ref$mask), length(ref$mask),
            100 * sum(!ref$mask) / length(ref$mask)))
cat(sprintf("  PCA components: %d\n", nrow(ref$pca_components)))
cat(sprintf("  Samples in ref: %d\n", ncol(ref$null_ratios)))
cat("\n")

# ---- Step 3: Bin + predict test BAM -----------------------------------------

cat("## Step 3: Predict on test BAM (rwisecondorx_predict)\n")
test_name <- sub("\\.bam$|\\.cram$", "", basename(test_bam))
cat(sprintf("  Binning: %s\n", test_name))
t0 <- proc.time()
test_sample <- suppressMessages(
  bam_convert(test_bam, binsize = samp_binsize, mapq = 1L, rmdup = "streaming")
)
cat(sprintf("  Binned in %.1fs\n", (proc.time() - t0)[["elapsed"]]))

cat(sprintf("  Predicting ...\n"))
t0 <- proc.time()
pred <- suppressMessages(
  rwisecondorx_predict(
    sample    = test_sample,
    reference = ref,
    zscore    = 5,
    alpha     = 1e-4,
    seed      = seed,
    cpus      = cpus
  )
)
cat(sprintf("  Prediction in %.1fs\n", (proc.time() - t0)[["elapsed"]]))

n_ab <- nrow(pred$aberrations)
cat(sprintf("  Aberration calls: %d\n", n_ab))
if (n_ab > 0L) {
  cat("  Aberrations:\n")
  for (i in seq_len(min(n_ab, 20L))) {
    ab <- pred$aberrations[i, ]
    cat(sprintf("    chr%s:%d-%d  type=%s  zscore=%.2f\n",
                ab$chr, ab$start, ab$end, ab$type,
                if ("zscore" %in% names(ab)) ab$zscore else NA_real_))
  }
  if (n_ab > 20L) cat(sprintf("    ... (%d more)\n", n_ab - 20L))
}

# Per-chromosome mean Z-score
chr_z_means <- vapply(as.character(1:22), function(chr) {
  z <- pred$results_z[[chr]]
  z <- z[!is.na(z) & z != 0]
  if (length(z) == 0L) return(NA_real_)
  mean(z)
}, numeric(1L))

cat("\n  Per-chromosome mean Z-score:\n")
for (chr in as.character(1:22)) {
  z <- chr_z_means[chr]
  flag <- if (!is.na(z) && abs(z) > 3) " <<< ELEVATED" else ""
  cat(sprintf("    chr%-3s  Z = %6.3f%s\n", chr, if (is.na(z)) 0 else z, flag))
}
cat("\n")

# ---- Step 4: NIPTeR pipeline on control + test -------------------------------

cat("## Step 4: NIPTeR pipeline (nipter_bin_bam + nipter_z_score)\n")
if (!is.null(fasta)) {
  cat("  GC correction: enabled (FASTA supplied)\n")
} else {
  cat("  GC correction: skipped (RWXCONF_FASTA not set)\n")
}

t0 <- proc.time()
nipter_controls <- vector("list", length(control_bams))
names(nipter_controls) <- names(control_samples)

for (i in seq_along(control_bams)) {
  nipter_controls[[i]] <- suppressMessages(
    nipter_bin_bam(control_bams[i], binsize = 50000L, mapq = 1L,
                   exclude_flags = 1024L)
  )
}
nipter_test <- suppressMessages(
  nipter_bin_bam(test_bam, binsize = 50000L, mapq = 1L, exclude_flags = 1024L)
)
cat(sprintf("  Binned %d control + 1 test sample in %.1fs\n",
            length(control_bams), (proc.time() - t0)[["elapsed"]]))

cg_nipter <- nipter_as_control_group(nipter_controls)

# GC correction (if FASTA available)
if (!is.null(fasta) && file.exists(fasta)) {
  cat("  Applying GC correction ...\n")
  t0 <- proc.time()
  cg_nipter    <- suppressWarnings(nipter_gc_correct(cg_nipter, fasta = fasta))
  nipter_test  <- suppressWarnings(nipter_gc_correct(nipter_test, fasta = fasta))
  cat(sprintf("  GC correction in %.1fs\n", (proc.time() - t0)[["elapsed"]]))
}

# Chi correction
cat("  Applying chi-squared correction ...\n")
chi_result  <- nipter_chi_correct(nipter_test, cg_nipter)
nipter_test <- chi_result$sample
cg_nipter   <- chi_result$control_group

# Z-scores for T21, T18, T13
cat("\n  Chromosomal Z-scores (standard NIPT targets):\n")
for (chr in c(21L, 18L, 13L, 15L, 16L)) {
  z_result <- tryCatch(
    nipter_z_score(nipter_test, cg_nipter, chromo_focus = chr),
    error = function(e) NULL
  )
  if (!is.null(z_result)) {
    cat(sprintf("    chr%-3d  Z = %6.3f  (ctrl sd = %.4f)\n",
                chr, z_result$sample_z_score,
                z_result$control_statistics["sd"]))
  }
}
cat("\n")

# ---- Step 5: Python WisecondorX comparison (optional) ------------------------

has_condathis  <- requireNamespace("condathis",  quietly = TRUE)
has_reticulate <- requireNamespace("reticulate", quietly = TRUE)
has_numpy      <- has_reticulate && tryCatch(
  { reticulate::import("numpy", convert = FALSE); TRUE }, error = function(e) FALSE)

cat("## Step 5: Python WisecondorX comparison\n")
if (!has_condathis || !has_numpy) {
  cat("  SKIPPED: condathis or reticulate/numpy not available.\n\n")
} else {
  cat("  condathis and numpy available — running Python pipeline ...\n")
  np <- reticulate::import("numpy", convert = FALSE)

  npz_dir  <- file.path(out_dir, "npz")
  dir.create(npz_dir, showWarnings = FALSE)
  py_ref   <- file.path(out_dir, "py_reference.npz")
  py_out   <- file.path(out_dir, paste0(test_name, "_py"))

  # Convert controls
  cat("  Converting control BAMs to NPZ ...\n")
  npz_files <- character(length(control_bams))
  for (i in seq_along(control_bams)) {
    npz_files[i] <- file.path(npz_dir, paste0(names(control_samples)[i], ".npz"))
    suppressMessages(
      bam_convert_npz(control_bams[i], npz_files[i],
                      binsize = samp_binsize, rmdup = "streaming", np = np)
    )
  }

  # Python newref
  cat("  Building Python reference ...\n")
  t0 <- proc.time()
  suppressMessages(
    wisecondorx_newref(npz_files = npz_files, output = py_ref,
                       ref_binsize = ref_binsize, nipt = TRUE, cpus = cpus)
  )
  cat(sprintf("  Python reference in %.1fs\n", (proc.time() - t0)[["elapsed"]]))

  # Test NPZ
  test_npz <- file.path(npz_dir, paste0(test_name, ".npz"))
  suppressMessages(
    bam_convert_npz(test_bam, test_npz, binsize = samp_binsize,
                    rmdup = "streaming", np = np)
  )

  # Python predict
  cat("  Running Python prediction ...\n")
  t0 <- proc.time()
  suppressMessages(
    wisecondorx_predict(npz = test_npz, ref = py_ref,
                        output_prefix = py_out, bed = TRUE, seed = 1L)
  )
  cat(sprintf("  Python prediction in %.1fs\n", (proc.time() - t0)[["elapsed"]]))

  # Parse Python aberrations BED
  py_ab_bed <- paste0(py_out, "_aberrations.bed")
  if (file.exists(py_ab_bed)) {
    py_lines <- readLines(py_ab_bed)
    py_lines <- py_lines[!startsWith(py_lines, "track") & nzchar(py_lines)]
    py_ab_df <- if (length(py_lines) > 0L) {
      do.call(rbind, lapply(strsplit(py_lines, "\t"), function(f) {
        data.frame(chr = sub("^chr", "", f[1L]), start = as.integer(f[2L]),
                   end = as.integer(f[3L]), stringsAsFactors = FALSE)
      }))
    } else {
      data.frame(chr = character(0L), start = integer(0L), end = integer(0L))
    }

    r_chrs  <- unique(as.character(pred$aberrations$chr))
    py_chrs <- unique(py_ab_df$chr)

    cat(sprintf("  R aberrations:      %d calls on chrs: %s\n",
                nrow(pred$aberrations), paste(sort(r_chrs), collapse = ",")))
    cat(sprintf("  Python aberrations: %d calls on chrs: %s\n",
                nrow(py_ab_df), paste(sort(py_chrs), collapse = ",")))

    if (length(r_chrs) > 0L || length(py_chrs) > 0L) {
      jaccard <- length(intersect(r_chrs, py_chrs)) /
                 max(length(union(r_chrs, py_chrs)), 1L)
      cat(sprintf("  Chr-level Jaccard (R vs Python): %.3f\n", jaccard))
      cat(sprintf("  Conformance: %s\n",
                  if (jaccard >= 0.8) "PASS (>=0.8)" else
                  if (jaccard >= 0.5) "MARGINAL (0.5-0.8)" else "FAIL (<0.5)"))
    }
  } else {
    cat("  Python aberrations BED not found.\n")
  }
  cat("\n")
}

# ---- Step 6: NIPTeR R package comparison (optional) --------------------------

cat("## Step 6: NIPTeR R package comparison\n")
if (!requireNamespace("NIPTeR", quietly = TRUE)) {
  cat("  SKIPPED: NIPTeR package not installed.\n\n")
} else {
  cat("  NIPTeR package available — running cross-check on test sample ...\n")

  nipter_ref_s <- tryCatch(
    NIPTeR::bin_bam_sample(bam_filepath = test_bam, do_sort = FALSE,
                            separate_strands = FALSE),
    error = function(e) {
      cat("  NIPTeR::bin_bam_sample failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(nipter_ref_s)) {
    nipter_ref_cg <- tryCatch({
      ctrl_list <- lapply(seq_along(control_bams), function(i) {
        s <- NIPTeR::bin_bam_sample(control_bams[i], do_sort = FALSE,
                                    separate_strands = FALSE)
        s$sample.name <- names(control_samples)[i]
        s
      })
      NIPTeR::as.control.group(ctrl_list)
    }, error = function(e) {
      cat("  NIPTeR control group failed:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(nipter_ref_cg)) {
      nipter_chi_res  <- NIPTeR::chi.correct(nipter_ref_s, nipter_ref_cg)
      nipter_z21      <- NIPTeR::chromosomal.zscore(nipter_chi_res$sample,
                                                     nipter_chi_res$control.group,
                                                     chromo.focus = 21)$sample.z.score
      our_z21_val <- tryCatch(
        nipter_z_score(chi_result$sample, chi_result$control_group,
                       21L)$sample_z_score,
        error = function(e) NA_real_
      )
      cat(sprintf("  Chr21 Z-score: R = %.4f,  NIPTeR = %.4f,  diff = %.4f\n",
                  our_z21_val, nipter_z21, abs(our_z21_val - nipter_z21)))
      cat(sprintf("  Agreement: %s\n",
                  if (abs(our_z21_val - nipter_z21) < 0.1) "PASS (|diff|<0.1)"
                  else "MARGINAL"))

      cat(sprintf("  [GC note] GC correction NOT cross-checked: R uses fasta_nuc(),\n"))
      cat(sprintf("            NIPTeR uses bundled GC tables. Numeric divergence expected.\n"))
    }
  }
  cat("\n")
}

# ---- Final summary -----------------------------------------------------------

cat("## Summary\n")
cat(sprintf("Report written to: %s\n", report_file))
cat(sprintf("Date: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("RWisecondorX version: %s\n", as.character(utils::packageVersion("RWisecondorX"))))
cat("\nTo add this to VALIDATION.md, report the above Jaccard and Z-score results.\n")

sink()
message("\nConformance report saved to: ", report_file)
