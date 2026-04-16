#!/usr/bin/env Rscript
# inst/scripts/build_reference.R
#
# CLI script for building a WisecondorX or NIPTeR reference from BAM/CRAM
# files or pre-binned BED.gz files.
#
# Input can be specified as:
#   --bam-dir / --bed-dir   : all matching files in a directory
#   --bam-list / --bed-list : a text file with one path per line
#
# Usage:
#   Rscript build_reference.R --mode wisecondorx --bam-dir /path/to/bams --out ref.rds
#   Rscript build_reference.R --mode nipter --bam-list samples.txt --out control.rds
#   Rscript build_reference.R --mode nipter --bam-list samples.txt --fasta hg38.fa --out control.rds
#   Rscript build_reference.R --mode wisecondorx --bed-dir /path/to/beds --out ref.rds
#   Rscript build_reference.R --mode nipter --bed-list beds.txt --gc-table hg38_gc.tsv.bgz --out control.rds

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with: install.packages('optparse')",
       call. = FALSE)
}

library(optparse)

option_list <- list(
  make_option("--mode", type = "character", default = "wisecondorx",
              help = "Reference type: 'wisecondorx' or 'nipter' [default: %default]"),

  # --- Input sources (exactly one BAM source + optional BED source) ---
  make_option("--bam-dir", type = "character", default = NULL,
              help = "Directory of BAM/CRAM files to bin"),
  make_option("--bam-list", type = "character", default = NULL,
              help = "Text file with one BAM/CRAM path per line"),
  make_option("--bed-dir", type = "character", default = NULL,
              help = "Directory of pre-binned BED.gz files (skip binning)"),
  make_option("--bed-list", type = "character", default = NULL,
              help = "Text file with one BED.gz path per line (skip binning)"),
  make_option("--out", type = "character", default = NULL,
              help = "Output path for the reference RDS file [required]"),

  # --- File discovery ---
  make_option("--bam-pattern", type = "character", default = "*.bam",
              help = "Glob pattern for BAM files in --bam-dir [default: %default]"),
  make_option("--bed-pattern", type = "character", default = "*.bed.gz",
              help = "Glob pattern for BED files in --bed-dir [default: %default]"),

  # --- Binning parameters ---
  make_option("--binsize", type = "integer", default = NULL,
              help = paste0("Bin size in bp. Defaults: 5000 (wisecondorx), ",
                            "50000 (nipter) [default: mode-dependent]")),
  make_option("--ref-binsize", type = "integer", default = NULL,
              help = paste0("Reference binsize for wisecondorx newref (rescale target). ",
                            "Ignored in nipter mode [default: same as --binsize]")),
  make_option("--mapq", type = "integer", default = NULL,
              help = "Minimum MAPQ for binning [default: 1 wisecondorx, 0 nipter]"),
  make_option("--rmdup", type = "character", default = NULL,
              help = paste0("Dedup mode: 'streaming', 'flag', 'none'. ",
                            "Defaults: 'streaming' (wisecondorx), 'none' (nipter)")),
  make_option("--exclude-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask to exclude (samtools -F) [default: %default]"),
  make_option("--require-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask to require (samtools -f) [default: %default]"),
  make_option("--reference", type = "character", default = NULL,
              help = "FASTA reference for CRAM decoding"),

  # --- Binning output location ---
  make_option("--bed-out-dir", type = "character", default = NULL,
              help = paste0("Directory for binned BED.gz output (used with --bam-dir/--bam-list). ",
                            "Defaults to <bam-dir>/beds or a temp directory")),

  # --- NIPTeR GC correction ---
  make_option("--gc-table", type = "character", default = NULL,
              help = paste0("Precomputed GC table (TSV.bgz from nipter_gc_precompute). ",
                            "NIPTeR mode: GC-correct each sample before writing BED. ",
                            "Requires binning from BAMs (ignored with --bed-dir/--bed-list)")),
  make_option("--fasta", type = "character", default = NULL,
              help = paste0("FASTA reference for on-the-fly GC correction (alternative to --gc-table). ",
                            "NIPTeR mode only. Slower than --gc-table for large cohorts")),

  # --- WisecondorX newref parameters ---
  make_option("--refsize", type = "integer", default = 300L,
              help = "KNN neighbors for wisecondorx newref [default: %default]"),
  make_option("--nipt", action = "store_true", default = FALSE,
              help = "Enable NIPT mode in wisecondorx newref [default: %default]"),
  make_option("--cpus", type = "integer", default = 4L,
              help = "Thread count for KNN/CBS [default: %default]"),

  # --- NIPTeR parameters ---
  make_option("--separate-strands", action = "store_true", default = FALSE,
              help = "NIPTeR mode: produce SeparatedStrands samples [default: %default]"),
  make_option("--description", type = "character", default = "General control group",
              help = "Description for the NIPTeR control group [default: %default]")
)

parser <- OptionParser(
  usage = "usage: %prog [options]",
  description = paste(
    "Build a WisecondorX reference or NIPTeR control group from BAM/CRAM or BED.gz files.",
    "",
    "Input files can be specified via --bam-dir/--bed-dir (directory + glob pattern)",
    "or --bam-list/--bed-list (text file with one path per line).",
    "",
    "Examples:",
    "  # WisecondorX: bin BAMs at 5kb, build reference at 100kb",
    "  %prog --mode wisecondorx --bam-dir bams/ --binsize 5000 --ref-binsize 100000 --out ref.rds",
    "",
    "  # WisecondorX: use a file list of pre-binned BEDs",
    "  %prog --mode wisecondorx --bed-list beds.txt --out ref.rds",
    "",
    "  # NIPTeR: bin BAMs with GC correction",
    "  %prog --mode nipter --bam-dir bams/ --gc-table hg38_gc_50k.tsv.bgz --out control.rds",
    "",
    "  # NIPTeR: bin BAMs with on-the-fly GC correction from FASTA",
    "  %prog --mode nipter --bam-dir bams/ --fasta hg38.fa --out control.rds",
    "",
    "  # NIPTeR: pre-filtered BAMs from a file list",
    "  %prog --mode nipter --bam-list samples.txt --mapq 40 --exclude-flags 1024 --out control.rds",
    sep = "\n"
  ),
  option_list = option_list
)

opts <- parse_args(parser)

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------

mode <- tolower(opts$mode)
if (!mode %in% c("wisecondorx", "nipter")) {
  stop("--mode must be 'wisecondorx' or 'nipter', got: ", opts$mode, call. = FALSE)
}

if (is.null(opts$out)) {
  stop("--out is required (path for the output RDS file)", call. = FALSE)
}

# Resolve BAM input sources
has_bam_dir  <- !is.null(opts$`bam-dir`)
has_bam_list <- !is.null(opts$`bam-list`)
has_bed_dir  <- !is.null(opts$`bed-dir`)
has_bed_list <- !is.null(opts$`bed-list`)

has_bam_input <- has_bam_dir || has_bam_list
has_bed_input <- has_bed_dir || has_bed_list

if (!has_bam_input && !has_bed_input) {
  stop("At least one input source is required: --bam-dir, --bam-list, --bed-dir, or --bed-list",
       call. = FALSE)
}

if (has_bam_input && has_bed_input) {
  stop("Provide BAM input (--bam-dir/--bam-list) OR BED input (--bed-dir/--bed-list), not both",
       call. = FALSE)
}

if (has_bam_dir && has_bam_list) {
  stop("--bam-dir and --bam-list are mutually exclusive", call. = FALSE)
}

if (has_bed_dir && has_bed_list) {
  stop("--bed-dir and --bed-list are mutually exclusive", call. = FALSE)
}

if (!is.null(opts$`gc-table`) && !is.null(opts$fasta)) {
  stop("--gc-table and --fasta are mutually exclusive; provide one or the other", call. = FALSE)
}

if (!is.null(opts$fasta) && mode != "nipter") {
  stop("--fasta GC correction is only valid in nipter mode", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Resolve file lists
# ---------------------------------------------------------------------------

.read_file_list <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  missing <- lines[!file.exists(lines)]
  if (length(missing) > 0L) {
    stop("Files not found (listed in ", path, "):\n  ",
         paste(head(missing, 10), collapse = "\n  "),
         if (length(missing) > 10) paste0("\n  ... and ", length(missing) - 10, " more"),
         call. = FALSE)
  }
  lines
}

if (has_bam_dir) {
  stopifnot(dir.exists(opts$`bam-dir`))
  bam_files <- sort(Sys.glob(file.path(opts$`bam-dir`, opts$`bam-pattern`)))
  if (length(bam_files) == 0L) {
    stop("No files matching '", opts$`bam-pattern`, "' in ", opts$`bam-dir`, call. = FALSE)
  }
} else if (has_bam_list) {
  bam_files <- .read_file_list(opts$`bam-list`)
  if (length(bam_files) == 0L) {
    stop("No files listed in ", opts$`bam-list`, call. = FALSE)
  }
} else {
  bam_files <- NULL
}

if (has_bed_dir) {
  stopifnot(dir.exists(opts$`bed-dir`))
  bed_files <- sort(Sys.glob(file.path(opts$`bed-dir`, opts$`bed-pattern`)))
  if (length(bed_files) == 0L) {
    stop("No files matching '", opts$`bed-pattern`, "' in ", opts$`bed-dir`, call. = FALSE)
  }
} else if (has_bed_list) {
  bed_files <- .read_file_list(opts$`bed-list`)
  if (length(bed_files) == 0L) {
    stop("No files listed in ", opts$`bed-list`, call. = FALSE)
  }
} else {
  bed_files <- NULL
}

# ---------------------------------------------------------------------------
# Set mode-dependent defaults
# ---------------------------------------------------------------------------

if (mode == "wisecondorx") {
  binsize <- if (is.null(opts$binsize)) 5000L else opts$binsize
  mapq    <- if (is.null(opts$mapq)) 1L else opts$mapq
  rmdup   <- if (is.null(opts$rmdup)) "streaming" else opts$rmdup
} else {
  binsize <- if (is.null(opts$binsize)) 50000L else opts$binsize
  mapq    <- if (is.null(opts$mapq)) 0L else opts$mapq
  rmdup   <- if (is.null(opts$rmdup)) "none" else opts$rmdup
}

ref_binsize <- opts$`ref-binsize`
gc_table    <- opts$`gc-table`
gc_fasta    <- opts$fasta

library(RWisecondorX)

# ---------------------------------------------------------------------------
# Step 1: Bin BAMs → BED.gz (if BAM input provided)
# ---------------------------------------------------------------------------

if (!is.null(bam_files)) {
  # Determine output directory for BED files
  bed_out_dir <- opts$`bed-out-dir`
  if (is.null(bed_out_dir)) {
    if (has_bam_dir) {
      bed_out_dir <- file.path(opts$`bam-dir`, "beds")
    } else {
      bed_out_dir <- tempfile("rwx_beds_")
    }
  }
  if (!dir.exists(bed_out_dir)) dir.create(bed_out_dir, recursive = TRUE)

  cat(sprintf("Binning %d BAM files into %s (binsize=%d, mapq=%d, rmdup=%s) ...\n",
              length(bam_files), bed_out_dir, binsize, mapq, rmdup))

  bed_files <- character(length(bam_files))

  # Create a single DuckDB connection for the entire binning loop (NIPTeR GC

  # correction path). Avoids leaking one connection per BAM when on.exit()
  # handlers only fire at script exit.
  gc_con <- NULL
  if (mode == "nipter" && (!is.null(gc_table) || !is.null(gc_fasta))) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    gc_con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(gc_con)
    on.exit(DBI::dbDisconnect(gc_con, shutdown = TRUE), add = TRUE)
  }

  for (i in seq_along(bam_files)) {
    bam <- bam_files[i]
    bed_name <- paste0(tools::file_path_sans_ext(basename(bam)), ".bed.gz")
    bed_path <- file.path(bed_out_dir, bed_name)
    bed_files[i] <- bed_path

    if (file.exists(bed_path)) {
      cat(sprintf("  [%d/%d] %s (exists, skipping)\n", i, length(bam_files), basename(bam)))
      next
    }

    cat(sprintf("  [%d/%d] %s\n", i, length(bam_files), basename(bam)))

    if (mode == "wisecondorx") {
      bam_convert_bed(
        bam           = bam,
        bed           = bed_path,
        binsize       = binsize,
        mapq          = mapq,
        require_flags = opts$`require-flags`,
        exclude_flags = opts$`exclude-flags`,
        rmdup         = rmdup,
        reference     = opts$reference
      )
    } else {
      if (!is.null(gc_table) || !is.null(gc_fasta)) {
        raw_sample <- nipter_bin_bam(
          bam              = bam,
          binsize          = binsize,
          mapq             = mapq,
          require_flags    = opts$`require-flags`,
          exclude_flags    = opts$`exclude-flags`,
          rmdup            = rmdup,
          separate_strands = opts$`separate-strands`,
          con              = gc_con,
          reference        = opts$reference
        )
        corrected_sample <- if (!is.null(gc_table)) {
          nipter_gc_correct(raw_sample, gc_table = gc_table, con = gc_con)
        } else {
          nipter_gc_correct(raw_sample, fasta = gc_fasta, binsize = binsize, con = gc_con)
        }

        nipter_sample_to_bed(
          sample    = raw_sample,
          bed       = bed_path,
          binsize   = binsize,
          corrected = corrected_sample,
          con       = gc_con
        )
      } else {
        nipter_bin_bam_bed(
          bam              = bam,
          bed              = bed_path,
          binsize          = binsize,
          mapq             = mapq,
          require_flags    = opts$`require-flags`,
          exclude_flags    = opts$`exclude-flags`,
          rmdup            = rmdup,
          separate_strands = opts$`separate-strands`,
          reference        = opts$reference
        )
      }
    }
  }

  cat("Binning complete.\n")
}

# ---------------------------------------------------------------------------
# Step 2: Build reference from BED files
# ---------------------------------------------------------------------------

if (mode == "wisecondorx") {
  cat(sprintf("Building WisecondorX reference from %d BED files ...\n", length(bed_files)))

  # Load samples from BED files
  samples <- lapply(bed_files, function(f) bed_to_sample(f))
  cat(sprintf("  Loaded %d samples.\n", length(samples)))

  ref <- rwisecondorx_newref(
    samples     = samples,
    binsize     = if (is.null(ref_binsize)) binsize else ref_binsize,
    nipt        = opts$nipt,
    refsize     = opts$refsize,
    cpus        = opts$cpus
  )

  saveRDS(ref, opts$out)
  cat(sprintf("Reference saved to %s\n", opts$out))
  cat(sprintf("  Masked autosomal bins: %d / %d\n", sum(ref$mask), length(ref$mask)))
  cat(sprintf("  Gender model: %s\n",
              if (is.null(ref$gender_model)) "none" else "trained"))

} else {
  cat(sprintf("Building NIPTeR control group from %d BED files ...\n", length(bed_files)))

  # Load each BED into a NIPTeRSample
  samples <- lapply(bed_files, function(f) bed_to_nipter_sample(f))
  ctrl <- nipter_as_control_group(samples, description = opts$description)

  # If BEDs already have corrected counts (from GC correction during binning),
  # they were loaded by bed_to_nipter_sample(). If the user wants to GC-correct
  # from pre-existing BEDs, they should re-bin with --gc-table.
  saveRDS(ctrl, opts$out)
  cat(sprintf("Control group saved to %s\n", opts$out))
  cat(sprintf("  Samples: %d\n", length(ctrl$samples)))
  cat(sprintf("  Description: %s\n", ctrl$description))
}
