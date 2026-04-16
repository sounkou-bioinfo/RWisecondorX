#!/usr/bin/env Rscript
#
# preprocess_cohort.R
#
# Manifest-driven preprocessing for real BAM/CRAM cohorts. This script:
#   1. stages a BAM manifest under --out-root/manifests
#   2. precomputes a NIPTeR GC table once
#   3. optionally computes SeqFF fetal-fraction estimates from BAMs
#   4. converts every BAM to native RWisecondorX BED.gz
#   5. optionally writes upstream WisecondorX NPZ files via the Python CLI
#   6. converts every BAM to NIPTeR BED.gz with separated strands and GC-corrected columns
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

.log_sample <- function(i, n, label, stem) {
  cat(sprintf("[%s] [%d/%d] %s %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              i, n, label, stem))
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
  make_option("--seqff", action = "store_true", default = TRUE,
              help = "Compute SeqFF fetal-fraction estimates from BAMs [default: %default]"),
  make_option("--seqff-samtools-bin", type = "character", default = "samtools",
              help = "Samtools executable used by SeqFF [default: %default]"),
  make_option("--wcx-write-npz", action = "store_true", default = FALSE,
              help = "Also write upstream WisecondorX NPZ files via the Python CLI [default: %default]"),
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
  seqff = .ensure_dir(file.path(out_root, "seqff")),
  rwcx_beds = .ensure_dir(file.path(out_root, "rwcx_beds")),
  wisecondorx_npz = .ensure_dir(file.path(out_root, "wisecondorx_npz")),
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

if (isTRUE(opts$seqff)) {
  .run_one({
    seqff_rows <- vector("list", length(bams))
    for (i in seq_along(bams)) {
      bam <- bams[[i]]
      stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
      out_tsv <- file.path(dirs$seqff, paste0(stem, ".seqff.tsv"))
      if (file.exists(out_tsv) && !isTRUE(opts$overwrite)) {
        seqff_rows[[i]] <- utils::read.delim(
          out_tsv,
          sep = "\t",
          header = TRUE,
          stringsAsFactors = FALSE
        )
        next
      }
      .log_sample(i, length(bams), "SeqFF", stem)
      ff <- seqff_predict(
        input = bam,
        input_type = "bam",
        samtools_bin = opts$`seqff-samtools-bin`,
        samtools_exclude_flags = as.integer(opts$`nipter-exclude-flags`),
        samtools_min_mapq = as.integer(opts$`nipter-mapq`)
      )
      row <- data.frame(
        sample = stem,
        bam = bam,
        SeqFF = unname(ff[["SeqFF"]]),
        Enet = unname(ff[["Enet"]]),
        WRSC = unname(ff[["WRSC"]]),
        stringsAsFactors = FALSE
      )
      utils::write.table(
        row,
        file = out_tsv,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      )
      seqff_rows[[i]] <- row
    }
    seqff_summary <- do.call(rbind, seqff_rows)
    utils::write.table(
      seqff_summary,
      file = file.path(dirs$seqff, "seqff_summary.tsv"),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }, sprintf("Compute SeqFF fetal-fraction estimates for %d BAMs", length(bams)))
}

.run_one({
  for (i in seq_along(bams)) {
    bam <- bams[[i]]
    stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
    out_bed <- file.path(dirs$rwcx_beds, paste0(stem, ".bed.gz"))
    if (file.exists(out_bed) && !isTRUE(opts$overwrite)) {
      next
    }
    .log_sample(i, length(bams), "RWisecondorX BED", stem)
    bam_convert_bed(
      bam = bam,
      bed = out_bed,
      binsize = as.integer(opts$`wcx-binsize`),
      mapq = 1L,
      rmdup = "streaming"
    )
  }
}, sprintf("Convert %d BAMs to native RWisecondorX BED.gz", length(bams)))

if (isTRUE(opts$`wcx-write-npz`)) {
  .run_one({
    for (i in seq_along(bams)) {
      bam <- bams[[i]]
      stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
      out_npz <- file.path(dirs$wisecondorx_npz, paste0(stem, ".npz"))
      if (file.exists(out_npz) && !isTRUE(opts$overwrite)) {
        next
      }
      .log_sample(i, length(bams), "WisecondorX CLI NPZ", stem)
      wisecondorx_convert(
        bam = bam,
        npz = out_npz,
        binsize = as.integer(opts$`wcx-binsize`),
        normdup = FALSE
      )
    }
  }, sprintf("Convert %d BAMs to upstream WisecondorX NPZ via Python CLI", length(bams)))
}

.run_one({
  for (i in seq_along(bams)) {
    bam <- bams[[i]]
    stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
    out_bed <- file.path(dirs$nipter_beds, paste0(stem, ".bed.gz"))
    if (file.exists(out_bed) && !isTRUE(opts$overwrite)) {
      next
    }
    .log_sample(i, length(bams), "NIPTeR BED", stem)
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
cat("  seqff:       ", dirs$seqff, "\n", sep = "")
cat("  rwcx_beds:   ", dirs$rwcx_beds, "\n", sep = "")
cat("  wisecondorx_npz: ", dirs$wisecondorx_npz, "\n", sep = "")
cat("  nipter_beds: ", dirs$nipter_beds, "\n", sep = "")
