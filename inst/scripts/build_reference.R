#!/usr/bin/env Rscript
# inst/scripts/build_reference.R
#
# Build a reference/control object from preprocessed cohort artifacts only.
#
# Supported modes:
#   - rwisecondorx : native RWisecondorX reference from 4-column BED.gz files
#   - wisecondorx  : upstream WisecondorX reference NPZ via condathis wrapper
#   - nipter       : NIPTeR control group from 5/9-column BED.gz files
#
# This script intentionally does not accept BAM/CRAM input. Preprocessing
# belongs in convert_sample.R or preprocess_cohort.R.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with: install.packages('optparse')",
       call. = FALSE)
}

library(optparse)

.read_file_list <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  missing <- lines[!file.exists(lines)]
  if (length(missing) > 0L) {
    stop("Files not found (listed in ", path, "):\n  ",
         paste(head(missing, 10L), collapse = "\n  "),
         if (length(missing) > 10L) {
           paste0("\n  ... and ", length(missing) - 10L, " more")
         },
         call. = FALSE)
  }
  lines
}

option_list <- list(
  make_option("--mode", type = "character", default = "rwisecondorx",
              help = "Reference type: 'rwisecondorx', 'wisecondorx', or 'nipter' [default: %default]"),
  make_option("--bed-dir", type = "character", default = NULL,
              help = "Directory of pre-binned BED.gz files"),
  make_option("--bed-list", type = "character", default = NULL,
              help = "Text file with one BED.gz path per line"),
  make_option("--npz-dir", type = "character", default = NULL,
              help = "Directory of WisecondorX NPZ files"),
  make_option("--npz-list", type = "character", default = NULL,
              help = "Text file with one WisecondorX NPZ path per line"),
  make_option("--out", type = "character", default = NULL,
              help = "Output path [required]"),
  make_option("--bed-pattern", type = "character", default = "*.bed.gz",
              help = "Glob pattern for BED files in --bed-dir [default: %default]"),
  make_option("--npz-pattern", type = "character", default = "*.npz",
              help = "Glob pattern for NPZ files in --npz-dir [default: %default]"),
  make_option("--binsize", type = "integer", default = NULL,
              help = "Input bin size. Defaults: 100000 (rwisecondorx/wisecondorx), 50000 (nipter)"),
  make_option("--ref-binsize", type = "integer", default = NULL,
              help = "Reference output binsize for wisecondorx/rwisecondorx [default: same as --binsize for rwisecondorx, 100000 for wisecondorx]"),
  make_option("--refsize", type = "integer", default = 300L,
              help = "KNN neighbors for wisecondorx/rwisecondorx reference building [default: %default]"),
  make_option("--nipt", action = "store_true", default = FALSE,
              help = "Enable NIPT mode in wisecondorx/rwisecondorx reference building [default: %default]"),
  make_option("--cpus", type = "integer", default = 4L,
              help = "Thread count for reference building [default: %default]"),
  make_option("--qc-json", type = "character", default = NULL,
              help = "Optional QC JSON output path for reference modes. Defaults to <out stem>_qc.json for rwisecondorx; upstream wisecondorx writes this automatically."),
  make_option("--description", type = "character", default = "General control group",
              help = "Description for the NIPTeR control group [default: %default]")
)

parser <- OptionParser(
  usage = paste(
    "%prog --mode rwisecondorx --bed-dir wcx_beds/ --out rwx_ref.rds",
    "%prog --mode wisecondorx --npz-dir wcx_npz/ --out wcx_ref.npz",
    "%prog --mode nipter --bed-dir nipter_beds/ --out nipter_cg.rds",
    sep = "\n"
  ),
  description = paste(
    "Build a reference/control object from preprocessed cohort artifacts only.",
    "",
    "Modes:",
    "  rwisecondorx : native RWisecondorX reference from BED.gz",
    "  wisecondorx  : upstream WisecondorX reference from NPZ",
    "  nipter       : NIPTeR control group from BED.gz",
    "",
    "Examples:",
    "  %prog --mode rwisecondorx --bed-dir wcx_beds/ --binsize 100000 --out rwx_ref.rds",
    "  %prog --mode wisecondorx --npz-dir wcx_npz/ --binsize 100000 --ref-binsize 100000 --out wcx_ref.npz",
    "  %prog --mode nipter --bed-dir nipter_beds/ --out nipter_cg.rds",
    sep = "\n"
  ),
  option_list = option_list
)

opts <- parse_args(parser)
mode <- tolower(opts$mode)
if (!mode %in% c("rwisecondorx", "wisecondorx", "nipter")) {
  stop("--mode must be 'rwisecondorx', 'wisecondorx', or 'nipter'", call. = FALSE)
}
if (is.null(opts$out) || !nzchar(opts$out)) {
  stop("--out is required", call. = FALSE)
}

has_bed_dir <- !is.null(opts$`bed-dir`)
has_bed_list <- !is.null(opts$`bed-list`)
has_npz_dir <- !is.null(opts$`npz-dir`)
has_npz_list <- !is.null(opts$`npz-list`)

if (mode %in% c("rwisecondorx", "nipter")) {
  if ((has_bed_dir || has_bed_list) == FALSE) {
    stop(mode, " mode requires --bed-dir or --bed-list", call. = FALSE)
  }
  if (has_npz_dir || has_npz_list) {
    stop(mode, " mode does not accept NPZ input", call. = FALSE)
  }
  if (has_bed_dir && has_bed_list) {
    stop("--bed-dir and --bed-list are mutually exclusive", call. = FALSE)
  }
} else {
  if ((has_npz_dir || has_npz_list) == FALSE) {
    stop("wisecondorx mode requires --npz-dir or --npz-list", call. = FALSE)
  }
  if (has_bed_dir || has_bed_list) {
    stop("wisecondorx mode does not accept BED input", call. = FALSE)
  }
  if (has_npz_dir && has_npz_list) {
    stop("--npz-dir and --npz-list are mutually exclusive", call. = FALSE)
  }
}

if (mode == "rwisecondorx") {
  binsize <- if (is.null(opts$binsize)) 100000L else as.integer(opts$binsize)
  ref_binsize <- if (is.null(opts$`ref-binsize`)) binsize else as.integer(opts$`ref-binsize`)
} else if (mode == "wisecondorx") {
  binsize <- if (is.null(opts$binsize)) 100000L else as.integer(opts$binsize)
  ref_binsize <- if (is.null(opts$`ref-binsize`)) 100000L else as.integer(opts$`ref-binsize`)
} else {
  binsize <- if (is.null(opts$binsize)) 50000L else as.integer(opts$binsize)
  ref_binsize <- NA_integer_
}

if (mode %in% c("rwisecondorx", "wisecondorx") && ref_binsize %% binsize != 0L) {
  stop("--ref-binsize must be a multiple of --binsize", call. = FALSE)
}

.default_qc_json <- function(path) {
  if (grepl("\\.npz$", path, ignore.case = TRUE)) {
    sub("\\.npz$", "_qc.json", path, ignore.case = TRUE)
  } else {
    sub("\\.rds$", "_qc.json", path, ignore.case = TRUE)
  }
}

qc_json <- opts$`qc-json`
if (is.null(qc_json) && mode == "rwisecondorx") {
  qc_json <- .default_qc_json(opts$out)
}

if (has_bed_dir) {
  stopifnot(dir.exists(opts$`bed-dir`))
  input_files <- sort(Sys.glob(file.path(opts$`bed-dir`, opts$`bed-pattern`)))
  if (length(input_files) == 0L) {
    stop("No files matching '", opts$`bed-pattern`, "' in ", opts$`bed-dir`, call. = FALSE)
  }
} else if (has_bed_list) {
  input_files <- .read_file_list(opts$`bed-list`)
} else if (has_npz_dir) {
  stopifnot(dir.exists(opts$`npz-dir`))
  input_files <- sort(Sys.glob(file.path(opts$`npz-dir`, opts$`npz-pattern`)))
  if (length(input_files) == 0L) {
    stop("No files matching '", opts$`npz-pattern`, "' in ", opts$`npz-dir`, call. = FALSE)
  }
} else {
  input_files <- .read_file_list(opts$`npz-list`)
}

library(RWisecondorX)

if (mode == "rwisecondorx") {
  cat(sprintf("Building native rwisecondorx reference from %d BED files ...\n", length(input_files)))
  ref <- rwisecondorx_newref(
    bed_dir = if (has_bed_dir) opts$`bed-dir` else NULL,
    samples = if (has_bed_list) lapply(input_files, bed_to_sample) else NULL,
    binsize = ref_binsize,
    nipt = opts$nipt,
    refsize = as.integer(opts$refsize),
    cpus = as.integer(opts$cpus)
  )
  saveRDS(ref, opts$out)
  qc <- rwisecondorx_ref_qc(ref, output_json = qc_json)
  cat(sprintf("Native reference saved to %s\n", opts$out))
  if (!is.null(qc_json)) {
    cat(sprintf("Native reference QC written to %s\n", qc_json))
  }
  cat(sprintf("Native reference QC verdict: %s\n", qc$overall_verdict))
} else if (mode == "wisecondorx") {
  cat(sprintf("Building upstream wisecondorx reference from %d NPZ files ...\n", length(input_files)))
  wisecondorx_newref(
    npz_files = input_files,
    output = opts$out,
    binsize = binsize,
    ref_binsize = ref_binsize,
    nipt = opts$nipt,
    refsize = as.integer(opts$refsize),
    cpus = as.integer(opts$cpus)
  )
  cat(sprintf("Upstream reference saved to %s\n", opts$out))
  cat(sprintf("Upstream reference QC written to %s\n", .default_qc_json(opts$out)))
} else {
  cat(sprintf("Building NIPTeR control group from %d BED files ...\n", length(input_files)))
  ctrl <- if (has_bed_dir) {
    nipter_control_group_from_beds(
      bed_dir = opts$`bed-dir`,
      pattern = opts$`bed-pattern`,
      binsize = binsize,
      description = opts$description
    )
  } else {
    nipter_as_control_group(
      lapply(input_files, bed_to_nipter_sample),
      description = opts$description
    )
  }
  saveRDS(ctrl, opts$out)
  cat(sprintf("Control group saved to %s\n", opts$out))
}
