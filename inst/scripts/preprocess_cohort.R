#!/usr/bin/env Rscript
#
# preprocess_cohort.R
#
# Manifest-driven preprocessing for real BAM/CRAM cohorts. This script:
#   1. stages a BAM manifest under --out-root/manifests
#   2. precomputes a NIPTeR GC table once
#   3. converts every BAM to WisecondorX BED.gz
#   4. optionally writes WisecondorX NPZ files for later upstream conformance
#   5. converts every BAM to NIPTeR BED.gz with separated strands and GC-corrected columns
#
# It intentionally does not build references or score samples.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with install.packages('optparse').",
       call. = FALSE)
}

library(optparse)

.read_file_list <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) >= 1L && identical(lines[[1L]], "bam")) {
    lines <- lines[-1L]
  }
  missing <- lines[!file.exists(lines)]
  if (length(missing) > 0L) {
    stop("Manifest contains missing BAM/CRAM paths:\n  ",
         paste(head(missing, 10L), collapse = "\n  "),
         if (length(missing) > 10L) {
           paste0("\n  ... and ", length(missing) - 10L, " more")
         },
         call. = FALSE)
  }
  lines
}

.write_manifest <- function(paths, out_path) {
  utils::write.table(
    data.frame(bam = paths, stringsAsFactors = FALSE),
    file = out_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

.ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("Failed to create directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.run_one <- function(expr, label) {
  cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), label))
  force(expr)
}

option_list <- list(
  make_option("--bam-list", type = "character", default = NULL,
              help = "Text file with one BAM/CRAM path per line [required]"),
  make_option("--out-root", type = "character", default = NULL,
              help = "Working root for manifests/BEDs/NPZs/logs [required]"),
  make_option("--fasta", type = "character", default = NULL,
              help = "Indexed reference FASTA for GC precompute [required]"),
  make_option("--threads", type = "integer", default = 20L,
              help = "Reserved thread budget for downstream work [default: %default]"),
  make_option("--wcx-binsize", type = "integer", default = 100000L,
              help = "WisecondorX convert binsize [default: %default]"),
  make_option("--nipter-binsize", type = "integer", default = 50000L,
              help = "NIPTeR binsize [default: %default]"),
  make_option("--nipter-mapq", type = "integer", default = 40L,
              help = "NIPTeR MAPQ filter [default: %default]"),
  make_option("--nipter-exclude-flags", type = "integer", default = 1024L,
              help = "NIPTeR exclude flags [default: %default]"),
  make_option("--nipter-separate-strands", action = "store_true", default = TRUE,
              help = "Write NIPTeR 9-column separated-strand BEDs [default: %default]"),
  make_option("--wcx-write-npz", action = "store_true", default = FALSE,
              help = "Also write WisecondorX NPZ files for later upstream conformance [default: %default]"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Overwrite existing outputs [default: %default]")
)

parser <- OptionParser(
  usage = "%prog --bam-list cohort.txt --out-root /mnt/data/niptseq/testdata --fasta ref.fa [options]",
  option_list = option_list
)
opts <- parse_args(parser)

if (is.null(opts$`bam-list`) || !file.exists(opts$`bam-list`)) {
  stop("--bam-list is required and must exist.", call. = FALSE)
}
if (is.null(opts$`out-root`) || !nzchar(opts$`out-root`)) {
  stop("--out-root is required.", call. = FALSE)
}
if (is.null(opts$fasta) || !file.exists(opts$fasta)) {
  stop("--fasta is required and must exist.", call. = FALSE)
}

bams <- .read_file_list(opts$`bam-list`)

out_root <- .ensure_dir(opts$`out-root`)
dirs <- list(
  manifests = .ensure_dir(file.path(out_root, "manifests")),
  gc = .ensure_dir(file.path(out_root, "gc")),
  wcx_beds = .ensure_dir(file.path(out_root, "wcx_beds")),
  wcx_npz = .ensure_dir(file.path(out_root, "wcx_npz")),
  nipter_beds = .ensure_dir(file.path(out_root, "nipter_beds")),
  logs = .ensure_dir(file.path(out_root, "logs"))
)

staged_manifest <- file.path(dirs$manifests, basename(opts$`bam-list`))
if (file.exists(staged_manifest) && !isTRUE(opts$overwrite)) {
  stop("Staged manifest already exists: ", staged_manifest,
       ". Use --overwrite to replace it.", call. = FALSE)
}
.write_manifest(bams, staged_manifest)

gc_table <- file.path(
  dirs$gc,
  sprintf("gc_%s.tsv.bgz", as.integer(opts$`nipter-binsize`))
)

library(RWisecondorX)

.run_one({
  nipter_gc_precompute(
    fasta = opts$fasta,
    out = gc_table,
    binsize = as.integer(opts$`nipter-binsize`)
  )
}, sprintf("Precompute GC table → %s", gc_table))

.run_one({
  for (i in seq_along(bams)) {
    bam <- bams[[i]]
    stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
    out_bed <- file.path(dirs$wcx_beds, paste0(stem, ".bed.gz"))
    if (file.exists(out_bed) && !isTRUE(opts$overwrite)) {
      next
    }
    cat(sprintf("  [%d/%d] WisecondorX BED %s\n", i, length(bams), stem))
    bam_convert_bed(
      bam = bam,
      bed = out_bed,
      binsize = as.integer(opts$`wcx-binsize`),
      mapq = 1L,
      rmdup = "streaming"
    )
  }
}, sprintf("Convert %d BAMs to WisecondorX BED.gz", length(bams)))

if (isTRUE(opts$`wcx-write-npz`)) {
  .run_one({
    for (i in seq_along(bams)) {
      bam <- bams[[i]]
      stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
      out_npz <- file.path(dirs$wcx_npz, paste0(stem, ".npz"))
      if (file.exists(out_npz) && !isTRUE(opts$overwrite)) {
        next
      }
      cat(sprintf("  [%d/%d] WisecondorX NPZ %s\n", i, length(bams), stem))
      bam_convert_npz(
        bam = bam,
        npz = out_npz,
        binsize = as.integer(opts$`wcx-binsize`),
        mapq = 1L,
        rmdup = "streaming"
      )
    }
  }, sprintf("Convert %d BAMs to WisecondorX NPZ", length(bams)))
}

.run_one({
  for (i in seq_along(bams)) {
    bam <- bams[[i]]
    stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
    out_bed <- file.path(dirs$nipter_beds, paste0(stem, ".bed.gz"))
    if (file.exists(out_bed) && !isTRUE(opts$overwrite)) {
      next
    }
    cat(sprintf("  [%d/%d] NIPTeR BED %s\n", i, length(bams), stem))
    raw_sample <- nipter_bin_bam(
      bam = bam,
      binsize = as.integer(opts$`nipter-binsize`),
      mapq = as.integer(opts$`nipter-mapq`),
      exclude_flags = as.integer(opts$`nipter-exclude-flags`),
      rmdup = "none",
      separate_strands = isTRUE(opts$`nipter-separate-strands`)
    )
    corrected_sample <- nipter_gc_correct(
      raw_sample,
      gc_table = gc_table
    )
    nipter_sample_to_bed(
      sample = raw_sample,
      corrected = corrected_sample,
      bed = out_bed,
      binsize = as.integer(opts$`nipter-binsize`)
    )
  }
}, sprintf("Convert %d BAMs to NIPTeR BED.gz", length(bams)))

cat("\nPreprocessing layout ready under:\n")
cat("  out_root:    ", out_root, "\n", sep = "")
cat("  manifest:    ", staged_manifest, "\n", sep = "")
cat("  gc_table:    ", gc_table, "\n", sep = "")
cat("  wcx_beds:    ", dirs$wcx_beds, "\n", sep = "")
cat("  wcx_npz:     ", dirs$wcx_npz, "\n", sep = "")
cat("  nipter_beds: ", dirs$nipter_beds, "\n", sep = "")
