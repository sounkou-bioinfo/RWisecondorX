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
#   7. writes native mosdepth-compatible 50 kb coverage outputs for raw-style
#      and filtered-style NIPT coverage
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

.compact_metadata <- function(x) {
  keep <- !vapply(x, is.null, logical(1L))
  x[keep]
}

.ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("Failed to create directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.normalize_legacy_wisecondorx_npz <- function(path) {
  legacy_path <- paste0(path, ".npz")
  if (!file.exists(path) && file.exists(legacy_path)) {
    ok <- file.rename(legacy_path, path)
    if (!isTRUE(ok) || !file.exists(path)) {
      stop(
        "Failed to rename legacy WisecondorX NPZ output:\n  ",
        legacy_path, "\n  -> ", path,
        call. = FALSE
      )
    }
  }
  invisible(path)
}

.run_one <- function(expr, label) {
  cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), label))
  force(expr)
  invisible(NULL)
}

.log_sample <- function(i, n, label, stem) {
  cat(sprintf("[%s] [%d/%d] %s %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              i, n, label, stem))
}

.run_stage_samples <- function(bams, jobs, worker) {
  jobs <- as.integer(jobs)
  stopifnot(length(jobs) == 1L, jobs >= 1L)
  if (length(bams) == 0L) {
    return(list())
  }
  res <- if (.Platform$OS.type == "unix" && jobs > 1L) {
    parallel::mclapply(
      X = seq_along(bams),
      FUN = function(i) worker(i, bams[[i]]),
      mc.cores = jobs,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_along(bams), function(i) worker(i, bams[[i]]))
  }

  bad <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(bad)) {
    msgs <- vapply(res[bad], as.character, character(1))
    stop(
      "Parallel stage failed:\n",
      paste(sprintf("[[%d]] %s", which(bad), msgs), collapse = "\n"),
      call. = FALSE
    )
  }

  res
}

.with_duckhts_con <- function(expr, envir = parent.frame()) {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  Rduckhts::rduckhts_load(con)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  eval(
    substitute(expr),
    envir = list2env(list(con = con), parent = envir)
  )
}

.per_worker_threads <- function(total_threads, jobs) {
  total_threads <- as.integer(total_threads)
  jobs <- as.integer(jobs)
  max(1L, total_threads %/% max(1L, jobs))
}

.mosdepth_thread_args <- function(total_threads, jobs) {
  worker_threads <- .per_worker_threads(total_threads, jobs)
  if (worker_threads <= 1L) {
    return(list(threads = 1L, processing_threads = 1L))
  }
  list(
    threads = worker_threads - 1L,
    processing_threads = 1L
  )
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
  make_option("--jobs", type = "integer", default = 4L,
              help = "Number of samples to process in parallel per stage [default: %default]"),
  make_option("--wcx-binsize", type = "integer", default = 100000L,
              help = "WisecondorX convert binsize [default: %default]"),
  make_option("--nipter-binsize", type = "integer", default = 50000L,
              help = "NIPTeR binsize [default: %default]"),
  make_option("--nipter-mapq", type = "integer", default = 40L,
              help = "NIPTeR MAPQ filter [default: %default]"),
  make_option("--nipter-exclude-flags", type = "integer", default = 1024L,
              help = "NIPTeR exclude flags [default: %default]"),
  make_option("--nipter-gc-include-sex", action = "store_true", default = FALSE,
              help = "Also GC-correct NIPTeR X/Y bins before BED export [default: %default]"),
  make_option("--nipter-separate-strands", action = "store_true", default = TRUE,
              help = "Write NIPTeR 9-column separated-strand BEDs [default: %default]"),
  make_option("--seqff", action = "store_true", default = TRUE,
              help = "Compute SeqFF fetal-fraction estimates from BAMs [default: %default]"),
  make_option("--wcx-write-npz", action = "store_true", default = FALSE,
              help = "Also write upstream WisecondorX NPZ files via the Python CLI [default: %default]"),
  make_option("--tabix-metadata", action = "store_true", default = FALSE,
              help = "Write optional ##RWX_<key>=<value> provenance lines into BED/tabix outputs [default: %default]"),
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
if (is.na(opts$jobs) || as.integer(opts$jobs) < 1L) {
  stop("--jobs must be >= 1.", call. = FALSE)
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
  mosdepth_raw_50k = .ensure_dir(file.path(out_root, "mosdepth_raw_50k")),
  mosdepth_filtered_50k = .ensure_dir(file.path(out_root, "mosdepth_filtered_50k")),
  logs = .ensure_dir(file.path(out_root, "logs"))
)

staged_manifest <- file.path(dirs$manifests, basename(opts$`bam-list`))
if (!file.exists(staged_manifest) || isTRUE(opts$overwrite)) {
  .write_manifest(bams, staged_manifest)
}

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
    seqff_rows <- .run_stage_samples(bams, opts$jobs, function(i, bam) {
      stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
      out_tsv <- file.path(dirs$seqff, paste0(stem, ".seqff.tsv"))
      if (file.exists(out_tsv) && !isTRUE(opts$overwrite)) {
        return(utils::read.delim(
          out_tsv,
          sep = "\t",
          header = TRUE,
          stringsAsFactors = FALSE
        ))
      }
      .log_sample(i, length(bams), "SeqFF", stem)
      ff <- seqff_predict(
        input = bam,
        input_type = "bam",
        mapq = as.integer(opts$`nipter-mapq`),
        exclude_flags = as.integer(opts$`nipter-exclude-flags`),
        reference = opts$fasta
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
      row
    })
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
  .run_stage_samples(bams, opts$jobs, function(i, bam) {
    stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
    out_bed <- file.path(dirs$rwcx_beds, paste0(stem, ".bed.gz"))
    if (file.exists(out_bed) && !isTRUE(opts$overwrite)) {
      return(invisible(NULL))
    }
    .log_sample(i, length(bams), "RWisecondorX BED", stem)
    bam_convert_bed(
      bam = bam,
      bed = out_bed,
      binsize = as.integer(opts$`wcx-binsize`),
      mapq = 1L,
      rmdup = "streaming",
      reference = opts$fasta,
      metadata = if (isTRUE(opts$`tabix-metadata`)) .compact_metadata(list(
        format = "rwisecondorx_bed",
        schema = "count_v1",
        binsize = as.integer(opts$`wcx-binsize`),
        mapq = 1L,
        require_flags = 0L,
        exclude_flags = 0L,
        rmdup = "streaming",
        reference = normalizePath(opts$fasta, winslash = "/", mustWork = TRUE)
      )) else NULL
    )
    invisible(NULL)
  })
}, sprintf("Convert %d BAMs to native RWisecondorX BED.gz", length(bams)))

if (isTRUE(opts$`wcx-write-npz`)) {
  .run_one({
    RWisecondorX:::.ensure_wisecondorx_env("wisecondorx")
    .run_stage_samples(bams, opts$jobs, function(i, bam) {
      stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
      out_npz <- file.path(dirs$wisecondorx_npz, paste0(stem, ".npz"))
      .normalize_legacy_wisecondorx_npz(out_npz)
      if (file.exists(out_npz) && !isTRUE(opts$overwrite)) {
        return(invisible(NULL))
      }
      .log_sample(i, length(bams), "WisecondorX CLI NPZ", stem)
      wisecondorx_convert(
        bam = bam,
        npz = out_npz,
        binsize = as.integer(opts$`wcx-binsize`),
        reference = opts$fasta,
        normdup = FALSE
      )
      invisible(NULL)
    })
  }, sprintf("Convert %d BAMs to upstream WisecondorX NPZ via Python CLI", length(bams)))
}

.run_one({
  .run_stage_samples(bams, opts$jobs, function(i, bam) {
    stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)
    out_bed <- file.path(dirs$nipter_beds, paste0(stem, ".bed.gz"))
    if (file.exists(out_bed) && !isTRUE(opts$overwrite)) {
      return(invisible(NULL))
    }
    .log_sample(i, length(bams), "NIPTeR BED", stem)
    raw_sample <- nipter_bin_bam(
      bam = bam,
      binsize = as.integer(opts$`nipter-binsize`),
      mapq = as.integer(opts$`nipter-mapq`),
      exclude_flags = as.integer(opts$`nipter-exclude-flags`),
      rmdup = "none",
      separate_strands = isTRUE(opts$`nipter-separate-strands`),
      reference = opts$fasta
    )
    corrected_sample <- nipter_gc_correct(
      raw_sample,
      gc_table = gc_table,
      include_sex = isTRUE(opts$`nipter-gc-include-sex`)
    )
    nipter_sample_to_bed(
      sample = raw_sample,
      corrected = corrected_sample,
      bed = out_bed,
      binsize = as.integer(opts$`nipter-binsize`),
      metadata = if (isTRUE(opts$`tabix-metadata`)) .compact_metadata(list(
        format = "nipter_bed",
        schema = if (isTRUE(opts$`nipter-separate-strands`)) "separated_v1" else "combined_v1",
        binsize = as.integer(opts$`nipter-binsize`),
        mapq = as.integer(opts$`nipter-mapq`),
        require_flags = 0L,
        exclude_flags = as.integer(opts$`nipter-exclude-flags`),
        rmdup = "none",
        corrected_columns = TRUE,
        gc_method = "loess",
        gc_include_sex = isTRUE(opts$`nipter-gc-include-sex`),
        reference = normalizePath(opts$fasta, winslash = "/", mustWork = TRUE),
        gc_table = normalizePath(gc_table, winslash = "/", mustWork = TRUE)
      )) else NULL
    )
    invisible(NULL)
  })
}, sprintf("Convert %d BAMs to NIPTeR BED.gz", length(bams)))

.run_one({
  mosdepth_threads <- .mosdepth_thread_args(opts$threads, opts$jobs)
  .run_stage_samples(bams, opts$jobs, function(i, bam) {
    stem <- sub("\\.(bam|cram)$", "", basename(bam), ignore.case = TRUE)

    raw_prefix <- file.path(dirs$mosdepth_raw_50k, stem)
    raw_summary <- paste0(raw_prefix, ".mosdepth.summary.txt")
    if (!file.exists(raw_summary) || isTRUE(opts$overwrite)) {
      .log_sample(i, length(bams), "Mosdepth raw 50k", stem)
      .with_duckhts_con({
        Rduckhts::rduckhts_mosdepth(
          con = con,
          prefix = raw_prefix,
          path = bam,
          by = "50000",
          fasta = opts$fasta,
          no_per_base = TRUE,
          threads = mosdepth_threads$threads,
          processing_threads = mosdepth_threads$processing_threads,
          flag = 1796L,
          include_flag = 0L,
          fast_mode = TRUE,
          mapq = 0L,
          overwrite = isTRUE(opts$overwrite)
        )
      })
    }

    filtered_prefix <- file.path(dirs$mosdepth_filtered_50k, stem)
    filtered_summary <- paste0(filtered_prefix, ".mosdepth.summary.txt")
    if (!file.exists(filtered_summary) || isTRUE(opts$overwrite)) {
      .log_sample(i, length(bams), "Mosdepth filtered 50k", stem)
      .with_duckhts_con({
        Rduckhts::rduckhts_mosdepth(
          con = con,
          prefix = filtered_prefix,
          path = bam,
          by = "50000",
          fasta = opts$fasta,
          no_per_base = TRUE,
          threads = mosdepth_threads$threads,
          processing_threads = mosdepth_threads$processing_threads,
          flag = bitwOr(1796L, 1024L),
          include_flag = 0L,
          fast_mode = TRUE,
          mapq = as.integer(opts$`nipter-mapq`),
          overwrite = isTRUE(opts$overwrite)
        )
      })
    }

    invisible(NULL)
  })
}, sprintf("Write native mosdepth-compatible 50 kb outputs for %d BAMs", length(bams)))

cat("\nPreprocessing layout ready under:\n")
cat("  out_root:    ", out_root, "\n", sep = "")
cat("  manifest:    ", staged_manifest, "\n", sep = "")
cat("  jobs:        ", as.integer(opts$jobs), "\n", sep = "")
cat("  gc_table:    ", gc_table, "\n", sep = "")
cat("  seqff:       ", dirs$seqff, "\n", sep = "")
cat("  rwcx_beds:   ", dirs$rwcx_beds, "\n", sep = "")
cat("  wisecondorx_npz: ", dirs$wisecondorx_npz, "\n", sep = "")
cat("  nipter_beds: ", dirs$nipter_beds, "\n", sep = "")
cat("  mosdepth_raw_50k: ", dirs$mosdepth_raw_50k, "\n", sep = "")
cat("  mosdepth_filtered_50k: ", dirs$mosdepth_filtered_50k, "\n", sep = "")
