#!/usr/bin/env Rscript
# inst/scripts/build_reference.R
#
# CLI script for building a WisecondorX or NIPTeR reference from BAM/CRAM
# files or pre-binned BED.gz files.
#
# Usage:
#   Rscript build_reference.R --mode wisecondorx --bam-dir /path/to/bams --out ref.rds
#   Rscript build_reference.R --mode nipter --bam-dir /path/to/bams --out control.rds
#   Rscript build_reference.R --mode wisecondorx --bed-dir /path/to/beds --out ref.rds
#
# The --bam-dir option bins BAMs into BED.gz files (written to --bed-dir or a
# temp directory), then builds the reference.  The --bed-dir option skips
# binning and loads pre-existing BED.gz files directly.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with: install.packages('optparse')",
       call. = FALSE)
}

library(optparse)

option_list <- list(
  make_option("--mode", type = "character", default = "wisecondorx",
              help = "Reference type: 'wisecondorx' or 'nipter' [default: %default]"),
  make_option("--bam-dir", type = "character", default = NULL,
              help = "Directory of BAM/CRAM files to bin (mutually exclusive with --bed-dir)"),
  make_option("--bed-dir", type = "character", default = NULL,
              help = "Directory of pre-binned BED.gz files (skip binning step)"),
  make_option("--out", type = "character", default = NULL,
              help = "Output path for the reference RDS file [required]"),
  make_option("--bam-pattern", type = "character", default = "*.bam",
              help = "Glob pattern for BAM files in --bam-dir [default: %default]"),
  make_option("--bed-pattern", type = "character", default = "*.bed.gz",
              help = "Glob pattern for BED files in --bed-dir [default: %default]"),
  make_option("--binsize", type = "integer", default = NULL,
              help = paste0("Bin size in bp. Defaults: 5000 (wisecondorx binning), ",
                            "50000 (nipter binning). For newref this is the reference ",
                            "binsize used after rescaling [default: mode-dependent]")),
  make_option("--ref-binsize", type = "integer", default = NULL,
              help = paste0("Reference binsize for wisecondorx newref (rescale target). ",
                            "If set, samples are rescaled from --binsize to --ref-binsize. ",
                            "Ignored in nipter mode [default: same as --binsize]")),
  make_option("--mapq", type = "integer", default = NULL,
              help = "Minimum MAPQ for binning [default: 1 wisecondorx, 0 nipter]"),
  make_option("--rmdup", type = "character", default = NULL,
              help = paste0("Dedup mode for binning: 'streaming', 'flag', 'none'. ",
                            "Defaults: 'streaming' (wisecondorx), 'none' (nipter)")),
  make_option("--exclude-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask for reads to exclude (samtools -F) [default: %default]"),
  make_option("--require-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask for reads to require (samtools -f) [default: %default]"),
  make_option("--reference", type = "character", default = NULL,
              help = "FASTA reference for CRAM inputs"),
  make_option("--refsize", type = "integer", default = 300L,
              help = "KNN neighbors for wisecondorx newref [default: %default]"),
  make_option("--nipt", action = "store_true", default = FALSE,
              help = "Enable NIPT mode in wisecondorx newref [default: %default]"),
  make_option("--cpus", type = "integer", default = 4L,
              help = "Thread count for KNN/CBS [default: %default]"),
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
    "Examples:",
    "  # WisecondorX: bin BAMs at 5kb, build reference at 100kb",
    "  %prog --mode wisecondorx --bam-dir bams/ --binsize 5000 --ref-binsize 100000 --out ref.rds",
    "",
    "  # WisecondorX: use pre-binned BEDs",
    "  %prog --mode wisecondorx --bed-dir beds/ --out ref.rds",
    "",
    "  # NIPTeR: bin BAMs at 50kb, build control group",
    "  %prog --mode nipter --bam-dir bams/ --out control.rds",
    "",
    "  # NIPTeR: pre-filtered BAMs (MAPQ >= 40, skip duplicates)",
    "  %prog --mode nipter --bam-dir bams/ --mapq 40 --exclude-flags 1024 --out control.rds",
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

if (is.null(opts$`bam-dir`) && is.null(opts$`bed-dir`)) {
  stop("At least one of --bam-dir or --bed-dir is required", call. = FALSE)
}

if (!is.null(opts$`bam-dir`) && !is.null(opts$`bed-dir`)) {
  stop("--bam-dir and --bed-dir are mutually exclusive", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Set mode-dependent defaults
# ---------------------------------------------------------------------------

if (mode == "wisecondorx") {
  binsize <- opts$binsize %||% 5000L
  mapq    <- opts$mapq    %||% 1L
  rmdup   <- opts$rmdup   %||% "streaming"
} else {
  binsize <- opts$binsize %||% 50000L
  mapq    <- opts$mapq    %||% 0L
  rmdup   <- opts$rmdup   %||% "none"
}

ref_binsize <- opts$`ref-binsize`

library(RWisecondorX)

# ---------------------------------------------------------------------------
# Step 1: Bin BAMs → BED.gz (if --bam-dir provided)
# ---------------------------------------------------------------------------

bed_dir <- opts$`bed-dir`

if (!is.null(opts$`bam-dir`)) {
  bam_dir <- opts$`bam-dir`
  stopifnot(dir.exists(bam_dir))

  bam_files <- Sys.glob(file.path(bam_dir, opts$`bam-pattern`))
  if (length(bam_files) == 0L) {
    stop("No files matching '", opts$`bam-pattern`, "' found in ", bam_dir, call. = FALSE)
  }

  # Write BEDs alongside BAMs if no --bed-dir, or into --bed-dir
  if (is.null(bed_dir)) {
    bed_dir <- file.path(bam_dir, "beds")
  }
  if (!dir.exists(bed_dir)) dir.create(bed_dir, recursive = TRUE)

  cat(sprintf("Binning %d BAM files into %s (binsize=%d, mapq=%d, rmdup=%s) ...\n",
              length(bam_files), bed_dir, binsize, mapq, rmdup))

  for (i in seq_along(bam_files)) {
    bam <- bam_files[i]
    bed_name <- paste0(tools::file_path_sans_ext(basename(bam)), ".bed.gz")
    bed_path <- file.path(bed_dir, bed_name)

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

  cat("Binning complete.\n")
}

stopifnot(dir.exists(bed_dir))

# ---------------------------------------------------------------------------
# Step 2: Build reference
# ---------------------------------------------------------------------------

if (mode == "wisecondorx") {
  cat(sprintf("Building WisecondorX reference from %s ...\n", bed_dir))

  ref <- rwisecondorx_newref(
    bed_dir     = bed_dir,
    bed_pattern = opts$`bed-pattern`,
    binsize     = ref_binsize %||% binsize,
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
  cat(sprintf("Building NIPTeR control group from %s ...\n", bed_dir))

  ctrl <- nipter_control_group_from_beds(
    bed_dir     = bed_dir,
    pattern     = opts$`bed-pattern`,
    description = opts$description
  )

  saveRDS(ctrl, opts$out)
  cat(sprintf("Control group saved to %s\n", opts$out))
  cat(sprintf("  Samples: %d\n", length(ctrl$samples)))
  cat(sprintf("  Description: %s\n", ctrl$description))
}
