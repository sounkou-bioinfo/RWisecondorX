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

.read_name_list <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  unique(lines[nzchar(lines) & !startsWith(lines, "#")])
}

.parse_name_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  unique(parts[nzchar(parts)])
}

.sanitize_inline_name_arg <- function(values) {
  values <- unique(values[nzchar(values)])
  if (length(values) == 1L && file.exists(values)) {
    return(character())
  }
  values
}

.read_delim_table <- function(path) {
  stopifnot(file.exists(path))
  first <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (length(first) && grepl("\t", first, fixed = TRUE)) "\t" else ","
  utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.match_col <- function(df, candidates, label) {
  nms <- names(df)
  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) {
    stop(
      "Could not find a ", label, " column. Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  nms[hit[[1L]]]
}

.read_sample_sex_map <- function(path) {
  df <- .read_delim_table(path)
  sample_col <- .match_col(
    df,
    c("sample_name", "sample", "Sample_name", "Sample", "sample_id", "SampleID"),
    "sample-name"
  )
  sex_col <- .match_col(
    df,
    c("sample_sex", "SampleSex", "sex", "Sex", "consensus_gender", "ConsensusGender"),
    "sample-sex"
  )

  vals <- trimws(tolower(as.character(df[[sex_col]])))
  vals[vals %in% c("f", "female", "xx")] <- "female"
  vals[vals %in% c("m", "male", "xy")] <- "male"
  vals[vals %in% c("u", "unknown", "na", "n/a")] <- "unknown"
  vals[vals %in% c("ambiguous", "undetermined")] <- "ambiguous"
  bad <- !vals %in% c("female", "male", "ambiguous", "unknown")
  if (any(bad)) {
    stop(
      "Invalid sample sex values in ", path, ": ",
      paste(unique(vals[bad]), collapse = ", "),
      call. = FALSE
    )
  }

  out <- stats::setNames(vals, as.character(df[[sample_col]]))
  dup <- duplicated(names(out))
  if (any(dup)) {
    stop(
      "Duplicate sample names in sample-sex table: ",
      paste(unique(names(out)[dup]), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

.read_y_unique_map <- function(path) {
  df <- .read_delim_table(path)
  sample_col <- .match_col(
    df,
    c("sample_name", "sample", "Sample_name", "Sample", "sample_id", "SampleID"),
    "sample-name"
  )
  value_col <- .match_col(
    df,
    c("y_unique_ratio_post"),
    "Y-unique ratio"
  )
  vals <- suppressWarnings(as.numeric(df[[value_col]]))
  if (any(!is.finite(vals))) {
    stop("Non-numeric Y-unique ratio values found in ", path, call. = FALSE)
  }
  out <- stats::setNames(vals, as.character(df[[sample_col]]))
  dup <- duplicated(names(out))
  if (any(dup)) {
    stop(
      "Duplicate sample names in Y-unique table: ",
      paste(unique(names(out)[dup]), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

.sample_stem <- function(path) {
  sub("\\.bed(\\.gz)?$|\\.tsv(\\.bgz)?$|\\.npz$", "",
      basename(path),
      ignore.case = TRUE)
}

.apply_sample_filters <- function(input_files,
                                  include_samples = character(),
                                  exclude_samples = character()) {
  sample_names <- vapply(input_files, .sample_stem, character(1L))
  dup <- duplicated(sample_names)
  if (any(dup)) {
    stop(
      "Duplicate sample stems in input files: ",
      paste(unique(sample_names[dup]), collapse = ", "),
      call. = FALSE
    )
  }

  all_names <- unique(sample_names)
  include_samples <- unique(include_samples[nzchar(include_samples)])
  exclude_samples <- unique(exclude_samples[nzchar(exclude_samples)])

  if (length(include_samples)) {
    missing_include <- setdiff(include_samples, all_names)
    if (length(missing_include)) {
      stop(
        "Requested include-sample names not found: ",
        paste(missing_include, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (length(exclude_samples)) {
    missing_exclude <- setdiff(exclude_samples, all_names)
    if (length(missing_exclude)) {
      stop(
        "Requested exclude-sample names not found: ",
        paste(missing_exclude, collapse = ", "),
        call. = FALSE
      )
    }
  }

  keep <- rep(TRUE, length(input_files))
  reason <- rep(NA_character_, length(input_files))
  if (length(include_samples)) {
    keep <- keep & sample_names %in% include_samples
    reason[!sample_names %in% include_samples] <- "not_in_include_list"
  }
  if (length(exclude_samples)) {
    reason[sample_names %in% exclude_samples] <- "user_excluded"
    keep[sample_names %in% exclude_samples] <- FALSE
  }

  list(
    files = input_files[keep],
    membership = data.frame(
      sample_name = sample_names,
      input_path = normalizePath(input_files, winslash = "/", mustWork = FALSE),
      kept_initially = keep,
      initial_exclusion_reason = reason,
      stringsAsFactors = FALSE
    )
  )
}

.default_model_rds <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    sub("\\.rds$", "_model.rds", path, ignore.case = TRUE)
  } else {
    paste0(path, "_model.rds")
  }
}

.default_qc_prefix <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    sub("\\.rds$", "_qc", path, ignore.case = TRUE)
  } else {
    paste0(path, "_qc")
  }
}

.ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

.write_tsv <- function(df, path) {
  .ensure_parent_dir(path)
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

.write_lines <- function(lines, path) {
  .ensure_parent_dir(path)
  writeLines(lines, path, useBytes = TRUE)
}

.qc_settings_df <- function(settings) {
  keys <- names(settings)
  values <- vapply(settings, function(x) {
    if (is.null(x)) {
      return(NA_character_)
    }
    if (length(x) == 1L) {
      return(as.character(x))
    }
    key_names <- names(x)
    if (is.null(key_names)) {
      key_names <- as.character(seq_along(x))
    }
    paste(sprintf("%s=%s", key_names, as.character(x)),
          collapse = ";")
  }, character(1L))
  data.frame(
    key = keys,
    value = values,
    stringsAsFactors = FALSE
  )
}

.write_nipter_qc_bundle <- function(qc, prefix) {
  .ensure_parent_dir(paste0(prefix, ".rds"))
  saveRDS(qc, paste0(prefix, ".rds"))
  .write_tsv(qc$chromosome_summary, paste0(prefix, ".chromosome_summary.tsv"))
  .write_tsv(qc$sample_summary, paste0(prefix, ".sample_summary.tsv"))
  if (!is.null(qc$sex_summary)) {
    .write_tsv(qc$sex_summary, paste0(prefix, ".sex_summary.tsv"))
  }
  if (!is.null(qc$sex_model_summary)) {
    .write_tsv(qc$sex_model_summary, paste0(prefix, ".sex_model_summary.tsv"))
  }
  if (!is.null(qc$bin_summary)) {
    .write_tsv(qc$bin_summary, paste0(prefix, ".bin_summary.tsv"))
  }
  .write_tsv(.qc_settings_df(qc$settings), paste0(prefix, ".settings.tsv"))
}

.reference_sex_map <- function(reference_model) {
  frame <- as.data.frame(reference_model$reference_frame, stringsAsFactors = FALSE)
  if (!all(c("Sample_name", "ConsensusGender") %in% names(frame))) {
    return(NULL)
  }
  stats::setNames(as.character(frame$ConsensusGender), frame$Sample_name)
}

.reference_sex_outliers <- function(reference_model) {
  frame <- as.data.frame(reference_model$reference_frame, stringsAsFactors = FALSE)
  if (!all(c("Sample_name", "IsRefSexOutlier") %in% names(frame))) {
    return(character())
  }
  as.character(frame$Sample_name[isTRUE(frame$IsRefSexOutlier) | frame$IsRefSexOutlier %in% TRUE])
}

.normalize_named_by_samples <- function(values, sample_names) {
  if (is.null(values)) {
    return(NULL)
  }
  stopifnot(all(sample_names %in% names(values)))
  values[sample_names]
}

.append_reason <- function(current, new_reason) {
  ifelse(
    is.na(current) | !nzchar(current),
    new_reason,
    paste(current, new_reason, sep = ";")
  )
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
  make_option("--model-rds", type = "character", default = NULL,
              help = "Optional NIPTeR reference-model RDS output path in nipter mode. Defaults to <out stem>_model.rds."),
  make_option("--control-qc-prefix", type = "character", default = NULL,
              help = "Optional output prefix for NIPTeR control-group QC files in nipter mode. Defaults to <out stem>_qc."),
  make_option("--include-bin-qc", action = "store_true", default = FALSE,
              help = "Include per-bin chi/CV QC output for NIPTeR control groups [default: %default]"),
  make_option("--write-qc-plots", action = "store_true", default = FALSE,
              help = "Write PNG plots for the NIPTeR QC bundle and sex-model reference spaces [default: %default]"),
  make_option("--qc-plot-dpi", type = "integer", default = 150L,
              help = "PNG DPI for --write-qc-plots in nipter mode [default: %default]"),
  make_option("--description", type = "character", default = "General control group",
              help = "Description for the NIPTeR control group [default: %default]"),
  make_option("--include-samples", type = "character", default = NULL,
              help = "Comma-separated sample names to retain before reference building"),
  make_option("--include-samples-file", type = "character", default = NULL,
              help = "Text file with one sample name per line to retain before reference building"),
  make_option("--exclude-samples", type = "character", default = NULL,
              help = "Comma-separated sample names to exclude before reference building"),
  make_option("--exclude-samples-file", type = "character", default = NULL,
              help = "Text file with one sample name per line to exclude before reference building"),
  make_option("--retained-samples-out", type = "character", default = NULL,
              help = "Optional text output path for the final retained sample names. Defaults to <control qc prefix>.retained_samples.txt in nipter mode."),
  make_option("--excluded-samples-out", type = "character", default = NULL,
              help = "Optional TSV output path for excluded sample names and reasons. Defaults to <control qc prefix>.excluded_samples.tsv in nipter mode."),
  make_option("--sample-sex-tsv", type = "character", default = NULL,
              help = "Optional TSV/CSV with sample sex annotations for nipter mode. Must contain sample and sex columns."),
  make_option("--sample-sex-source", type = "character", default = NULL,
              help = "Optional provenance string for sample sex annotations in nipter mode."),
  make_option("--y-unique-tsv", type = "character", default = NULL,
              help = "Optional TSV/CSV with per-sample Y-unique ratios for nipter mode. Must contain sample_name/sample and y_unique_ratio_post columns."),
  make_option("--sample-qc-tsv", type = "character", default = NULL,
              help = "Optional TSV/CSV with per-sample QC metrics for nipter mode. Used for hard filtering by read depth and/or GC."),
  make_option("--sample-qc-sample-col", type = "character", default = NULL,
              help = "Optional sample-name column in --sample-qc-tsv. When omitted, common names such as sample_name/Sample are inferred."),
  make_option("--sample-qc-total-unique-reads-col", type = "character", default = NULL,
              help = "Optional total-unique-reads column in --sample-qc-tsv. When omitted, common names such as TotalUniqueReads are inferred."),
  make_option("--sample-qc-gc-col", type = "character", default = NULL,
              help = "Optional GC column in --sample-qc-tsv. When omitted, common names such as GCPCTAfterFiltering are inferred."),
  make_option("--min-total-unique-reads", type = "double", default = NULL,
              help = "nipter mode: drop controls below this TotalUniqueReads threshold."),
  make_option("--max-total-unique-reads", type = "double", default = NULL,
              help = "nipter mode: drop controls above this TotalUniqueReads threshold."),
  make_option("--gc-mad-cutoff", type = "double", default = NULL,
              help = "nipter mode: drop GC outlier controls beyond this many MADs from the cohort median."),
  make_option("--nipter-autosomal-source", type = "character", default = "auto",
              help = "nipter mode: import autosomal counts from 'auto', 'raw', or 'corrected' BED columns [default: %default]"),
  make_option("--nipter-sex-counts", type = "character", default = "match",
              help = "nipter mode: import sex counts from 'match', 'raw', or 'corrected' BED columns [default: %default]"),
  make_option("--prune-outliers", action = "store_true", default = FALSE,
              help = "nipter mode: iteratively chi-correct, diagnose, and drop autosomal outlier controls before model building [default: %default]"),
  make_option("--prune-chi-cutoff", type = "double", default = 3.5,
              help = "nipter mode: chi cutoff for iterative outlier pruning [default: %default]"),
  make_option("--prune-z-cutoff", type = "double", default = 3,
              help = "nipter mode: absolute z-score cutoff for iterative outlier pruning [default: %default]"),
  make_option("--prune-outlier-rule", type = "character", default = "any_aberrant_score",
              help = "nipter mode: outlier pruning rule: 'any_aberrant_score' or 'bidirectional_or_multichromosome' [default: %default]"),
  make_option("--prune-max-aberrant-chromosomes", type = "integer", default = 2L,
              help = "nipter mode: max distinct aberrant chromosomes allowed before dropping a control [default: %default]"),
  make_option("--prune-min-controls", type = "integer", default = 10L,
              help = "nipter mode: minimum retained controls required during pruning [default: %default]"),
  make_option("--prune-max-iterations", type = "integer", default = 100L,
              help = "nipter mode: maximum pruning iterations [default: %default]"),
  make_option("--prune-collapse-strands", action = "store_true", default = FALSE,
              help = "nipter mode: collapse separated strands during outlier diagnosis instead of using strand-resolved diagnostics [default: %default]"),
  make_option("--drop-ref-sex-outliers", action = "store_true", default = FALSE,
              help = "nipter mode: after a provisional reference build, drop ConsensusGender sex outliers and rebuild the final model on the same retained set [default: %default]")
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
autosomal_source <- match.arg(tolower(opts$`nipter-autosomal-source`), c("auto", "raw", "corrected"))
sex_counts <- match.arg(tolower(opts$`nipter-sex-counts`), c("match", "raw", "corrected"))

include_samples <- unique(c(
  .sanitize_inline_name_arg(.parse_name_arg(opts$`include-samples`)),
  if (!is.null(opts$`include-samples-file`)) .read_name_list(opts$`include-samples-file`) else character()
))
exclude_samples <- unique(c(
  .sanitize_inline_name_arg(.parse_name_arg(opts$`exclude-samples`)),
  if (!is.null(opts$`exclude-samples-file`)) .read_name_list(opts$`exclude-samples-file`) else character()
))

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

filter_info <- .apply_sample_filters(
  input_files,
  include_samples = include_samples,
  exclude_samples = exclude_samples
)
input_files <- filter_info$files
membership <- filter_info$membership
if (length(input_files) < 2L) {
  stop("Fewer than 2 input files remain after sample filtering.", call. = FALSE)
}

library(RWisecondorX)

if (mode == "rwisecondorx") {
  cat(sprintf("Building native rwisecondorx reference from %d BED files ...\n", length(input_files)))
  ref <- rwisecondorx_newref(
    bed_dir = if (has_bed_dir && !length(include_samples) && !length(exclude_samples)) opts$`bed-dir` else NULL,
    samples = if (has_bed_dir && !length(include_samples) && !length(exclude_samples)) NULL else lapply(input_files, bed_to_sample),
    binsize = ref_binsize,
    nipt = opts$nipt,
    refsize = as.integer(opts$refsize),
    cpus = as.integer(opts$cpus)
  )
  .ensure_parent_dir(opts$out)
  saveRDS(ref, opts$out)
  qc <- rwisecondorx_ref_qc(ref, output_json = qc_json)
  cat(sprintf("Native reference saved to %s\n", opts$out))
  if (!is.null(qc_json)) {
    cat(sprintf("Native reference QC written to %s\n", qc_json))
  }
  cat(sprintf("Native reference QC verdict: %s\n", qc$overall_verdict))
} else if (mode == "wisecondorx") {
  cat(sprintf("Building upstream wisecondorx reference from %d NPZ files ...\n", length(input_files)))
  .ensure_parent_dir(opts$out)
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
  sample_sex <- if (!is.null(opts$`sample-sex-tsv`)) {
    .read_sample_sex_map(opts$`sample-sex-tsv`)
  } else {
    NULL
  }
  sex_source <- if (!is.null(sample_sex)) {
    if (!is.null(opts$`sample-sex-source`) && nzchar(opts$`sample-sex-source`)) {
      opts$`sample-sex-source`
    } else {
      paste0("file:", basename(opts$`sample-sex-tsv`))
    }
  } else {
    NULL
  }
  y_unique_ratios <- if (!is.null(opts$`y-unique-tsv`)) {
    .read_y_unique_map(opts$`y-unique-tsv`)
  } else {
    NULL
  }
  sample_qc <- if (!is.null(opts$`sample-qc-tsv`)) {
    .read_delim_table(opts$`sample-qc-tsv`)
  } else {
    NULL
  }

  cat(sprintf("Building NIPTeR control group from %d BED files ...\n", length(input_files)))
  ctrl <- if (has_bed_dir && !length(include_samples) && !length(exclude_samples)) {
    nipter_control_group_from_beds(
      bed_dir = opts$`bed-dir`,
      pattern = opts$`bed-pattern`,
      binsize = binsize,
      autosomal_source = autosomal_source,
      sex_counts = sex_counts,
      description = opts$description,
      sample_sex = sample_sex,
      sex_source = sex_source
    )
  } else {
    nipter_as_control_group(
      lapply(
        input_files,
        function(path) {
          bed_to_nipter_sample(
            path,
            binsize = binsize,
            autosomal_source = autosomal_source,
            sex_source = sex_counts
          )
        }
      ),
      description = opts$description,
      sample_sex = sample_sex,
      sex_source = sex_source
    )
  }

  qc_filter_result <- NULL
  if (!is.null(sample_qc) ||
      !is.null(opts$`min-total-unique-reads`) ||
      !is.null(opts$`max-total-unique-reads`) ||
      !is.null(opts$`gc-mad-cutoff`)) {
    qc_filter_result <- nipter_filter_control_group_qc(
      ctrl,
      sample_qc = sample_qc,
      sample_col = opts$`sample-qc-sample-col`,
      total_unique_reads_col = opts$`sample-qc-total-unique-reads-col`,
      gc_col = opts$`sample-qc-gc-col`,
      min_total_unique_reads = opts$`min-total-unique-reads`,
      max_total_unique_reads = opts$`max-total-unique-reads`,
      gc_mad_cutoff = opts$`gc-mad-cutoff`
    )
    ctrl <- qc_filter_result$control_group
  }

  prune_result <- NULL
  if (isTRUE(opts$`prune-outliers`)) {
    prune_result <- nipter_prune_control_group_outliers(
      ctrl,
      chi_cutoff = opts$`prune-chi-cutoff`,
      collapse_strands = isTRUE(opts$`prune-collapse-strands`),
      z_cutoff = opts$`prune-z-cutoff`,
      outlier_rule = opts$`prune-outlier-rule`,
      max_aberrant_chromosomes = as.integer(opts$`prune-max-aberrant-chromosomes`),
      min_controls = as.integer(opts$`prune-min-controls`),
      max_iterations = as.integer(opts$`prune-max-iterations`),
      verbose = TRUE
    )
    ctrl <- prune_result$control_group
    if (!isTRUE(prune_result$converged)) {
      stop(
        "NIPTeR control pruning did not converge (reason: ",
        prune_result$stop_reason,
        "). Refuse to build a fallback reference.",
        call. = FALSE
      )
    }
  }

  current_names <- control_names(ctrl)
  current_y_unique <- .normalize_named_by_samples(y_unique_ratios, current_names)

  .ensure_parent_dir(opts$out)
  saveRDS(ctrl, opts$out)
  model_rds <- if (is.null(opts$`model-rds`)) {
    .default_model_rds(opts$out)
  } else {
    opts$`model-rds`
  }
  ref_model <- nipter_build_reference(
    ctrl,
    y_unique_ratios = current_y_unique,
    sample_qc = sample_qc,
    sample_qc_sample_col = opts$`sample-qc-sample-col`,
    sample_qc_total_unique_reads_col = opts$`sample-qc-total-unique-reads-col`,
    sample_qc_gc_col = opts$`sample-qc-gc-col`,
    min_total_unique_reads = opts$`min-total-unique-reads`,
    max_total_unique_reads = opts$`max-total-unique-reads`,
    gc_mad_cutoff = opts$`gc-mad-cutoff`,
    build_params = list(
      control_group_rds = normalizePath(opts$out, winslash = "/", mustWork = FALSE),
      source_type = if (has_bed_dir) "bed_dir" else "bed_list",
      autosomal_source = autosomal_source,
      sex_counts = sex_counts,
      prune_outliers = isTRUE(opts$`prune-outliers`),
      prune_outlier_rule = opts$`prune-outlier-rule`
    )
  )
  ctrl <- ref_model$control_group

  sex_outliers <- character()
  if (isTRUE(opts$`drop-ref-sex-outliers`)) {
    sex_outliers <- .reference_sex_outliers(ref_model)
    if (length(sex_outliers)) {
      if ((n_controls(ctrl) - length(sex_outliers)) < as.integer(opts$`prune-min-controls`)) {
        stop(
          "Dropping reference sex outliers would leave fewer than --prune-min-controls samples.",
          call. = FALSE
        )
      }
      ctrl <- nipter_drop_control_group_samples(ctrl, sex_outliers)
      current_names <- control_names(ctrl)
      current_y_unique <- .normalize_named_by_samples(y_unique_ratios, current_names)
      .ensure_parent_dir(opts$out)
      saveRDS(ctrl, opts$out)
      ref_model <- nipter_build_reference(
        ctrl,
        y_unique_ratios = current_y_unique,
        sample_qc = sample_qc,
        sample_qc_sample_col = opts$`sample-qc-sample-col`,
        sample_qc_total_unique_reads_col = opts$`sample-qc-total-unique-reads-col`,
        sample_qc_gc_col = opts$`sample-qc-gc-col`,
        min_total_unique_reads = opts$`min-total-unique-reads`,
        max_total_unique_reads = opts$`max-total-unique-reads`,
        gc_mad_cutoff = opts$`gc-mad-cutoff`,
        build_params = list(
          control_group_rds = normalizePath(opts$out, winslash = "/", mustWork = FALSE),
          source_type = if (has_bed_dir) "bed_dir" else "bed_list",
          autosomal_source = autosomal_source,
          sex_counts = sex_counts,
          prune_outliers = isTRUE(opts$`prune-outliers`),
          prune_outlier_rule = opts$`prune-outlier-rule`,
          dropped_reference_sex_outliers = sex_outliers
        )
      )
      ctrl <- ref_model$control_group
    }
  }

  ref_model <- nipter_build_gaunosome_models(ref_model)
  .ensure_parent_dir(opts$out)
  saveRDS(ctrl, opts$out)
  .ensure_parent_dir(model_rds)
  saveRDS(ref_model, model_rds)
  qc_prefix <- if (is.null(opts$`control-qc-prefix`)) {
    .default_qc_prefix(opts$out)
  } else {
    opts$`control-qc-prefix`
  }
  reference_sex <- .reference_sex_map(ref_model)
  qc <- nipter_control_group_qc(
    ctrl,
    sample_sex = if (is.null(sample_sex)) reference_sex else sample_sex[control_names(ctrl)],
    reference_model = ref_model,
    z_cutoff = opts$`prune-z-cutoff`,
    collapse_strands = isTRUE(opts$`prune-collapse-strands`),
    max_aberrant_chromosomes = as.integer(opts$`prune-max-aberrant-chromosomes`),
    outlier_rule = opts$`prune-outlier-rule`,
    include_bins = isTRUE(opts$`include-bin-qc`)
  )
  qc$settings <- c(
    qc$settings,
    list(
      autosomal_source = autosomal_source,
      sex_counts = sex_counts,
      min_total_unique_reads = opts$`min-total-unique-reads`,
      max_total_unique_reads = opts$`max-total-unique-reads`,
      gc_mad_cutoff = opts$`gc-mad-cutoff`,
      n_qc_filtered_samples = if (is.null(qc_filter_result)) 0L else length(qc_filter_result$excluded_samples),
      prune_outliers = isTRUE(opts$`prune-outliers`),
      prune_converged = if (is.null(prune_result)) NA else isTRUE(prune_result$converged),
      prune_stop_reason = if (is.null(prune_result)) NA_character_ else prune_result$stop_reason,
      n_pruned_autosomal_outliers = if (is.null(prune_result)) 0L else length(prune_result$dropped_samples),
      n_dropped_reference_sex_outliers = length(sex_outliers)
    )
  )
  .write_nipter_qc_bundle(qc, qc_prefix)
  if (isTRUE(opts$`write-qc-plots`)) {
    write_nipter_reference_plots(
      qc = qc,
      reference = ref_model,
      outprefix = qc_prefix,
      dpi = as.integer(opts$`qc-plot-dpi`)
    )
  }

  retained_out <- if (is.null(opts$`retained-samples-out`)) {
    paste0(qc_prefix, ".retained_samples.txt")
  } else {
    opts$`retained-samples-out`
  }
  excluded_out <- if (is.null(opts$`excluded-samples-out`)) {
    paste0(qc_prefix, ".excluded_samples.tsv")
  } else {
    opts$`excluded-samples-out`
  }

  if (!is.null(qc_filter_result) && length(qc_filter_result$excluded_samples)) {
    membership$initial_exclusion_reason[membership$sample_name %in% qc_filter_result$excluded_samples] <-
      .append_reason(
        membership$initial_exclusion_reason[membership$sample_name %in% qc_filter_result$excluded_samples],
        qc_filter_result$exclusion_table$exclusion_reason[
          match(
            membership$sample_name[membership$sample_name %in% qc_filter_result$excluded_samples],
            qc_filter_result$exclusion_table$sample_name
          )
        ]
      )
  }
  if (!is.null(prune_result) && length(prune_result$dropped_samples)) {
    membership$initial_exclusion_reason[membership$sample_name %in% prune_result$dropped_samples] <-
      .append_reason(
        membership$initial_exclusion_reason[membership$sample_name %in% prune_result$dropped_samples],
        "autosomal_outlier"
      )
  }
  if (length(sex_outliers)) {
    membership$initial_exclusion_reason[membership$sample_name %in% sex_outliers] <-
      .append_reason(
        membership$initial_exclusion_reason[membership$sample_name %in% sex_outliers],
        "sex_outlier"
      )
  }
  membership$kept_final <- membership$sample_name %in% control_names(ctrl)

  .write_lines(control_names(ctrl), retained_out)
  excluded_df <- membership[!membership$kept_final, c("sample_name", "input_path", "initial_exclusion_reason"), drop = FALSE]
  if (!nrow(excluded_df)) {
    excluded_df <- data.frame(
      sample_name = character(),
      input_path = character(),
      initial_exclusion_reason = character(),
      stringsAsFactors = FALSE
    )
  }
  .write_tsv(excluded_df, excluded_out)
  .write_tsv(
    as.data.frame(ref_model$reference_frame, stringsAsFactors = FALSE),
    paste0(qc_prefix, ".reference_frame.tsv")
  )
  if (!is.null(prune_result)) {
    .write_tsv(prune_result$iteration_log, paste0(qc_prefix, ".prune_iteration_log.tsv"))
  }

  cat(sprintf("Control group saved to %s\n", opts$out))
  cat(sprintf("NIPTeR reference model saved to %s\n", model_rds))
  cat(sprintf("NIPTeR control-group QC files written with prefix %s\n", qc_prefix))
  cat(sprintf("Retained sample list written to %s\n", retained_out))
  cat(sprintf("Excluded sample table written to %s\n", excluded_out))
}
