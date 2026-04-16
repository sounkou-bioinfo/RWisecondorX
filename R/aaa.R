# Suppress R CMD check NOTEs for symbols used inside evalq() in the mclust
# namespace (nipter_sex.R and rwisecondorx_utils.R).
utils::globalVariables(c("Mclust", "emControl", "data"))


# %||% is base R since 4.4.0; define fallback for R 4.1-4.3
if (getRversion() < "4.4.0") {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
}


# ---- Shared chromosome name normalization ---------------------------------
#
# Two conventions coexist: the WisecondorX layer uses "23"/"24" as internal
# keys for X/Y, while the NIPTeR layer uses "X"/"Y" natively. Both need to
# strip the "chr" prefix from external sources (BAM headers, BED files,
# FASTA output). This shared helper handles both steps.

#' Normalize chromosome names: strip "chr" prefix and optionally map X/Y to
#' 23/24.
#'
#' Works on both scalar and vector inputs. Used throughout the package to
#' convert external chromosome names (from BAM headers, BED files, FASTA
#' nucleotide tables) to internal representation.
#'
#' @param x Character vector of chromosome names.
#' @param xy_to_numeric Logical; if `TRUE` (default), map "X"/"x" to "23" and
#'   "Y"/"y" to "24". Set to `FALSE` for the NIPTeR layer which uses "X"/"Y"
#'   natively.
#' @return Character vector with normalized names.
#' @noRd
.normalize_chr_name <- function(x, xy_to_numeric = TRUE) {
  x <- sub("^[Cc][Hh][Rr]", "", x)
  if (xy_to_numeric) {
    x[x == "X" | x == "x"] <- "23"
    x[x == "Y" | x == "y"] <- "24"
  }
  x
}

#' SRA Run Metadata Utilities
#'
#' These helpers standardize how `RWisecondorX` stages SRA run metadata for
#' conformance fixtures. They use the NCBI SRA backend `runinfo` endpoint and
#' default to writing project metadata into `inst/extdata/` in the source tree.
#'
#' The SRA Run Selector and backend `runinfo` endpoints can expose slightly
#' different columns. For reproducible package fixtures, `RWisecondorX`
#' standardizes on the backend `runinfo` CSV output.
#'
#' @param accession A study or project accession such as `PRJNA400134`.
#' @param dest Output CSV path. Defaults to `inst/extdata/<accession>_runinfo.csv`
#'   in the current source checkout.
#' @param overwrite Logical; overwrite an existing file.
#' @param quiet Logical; passed through to the downloader.
#' @param mode Retrieval mode. `"acc"` queries the backend with `acc=<accession>`;
#'   `"term"` queries with `term=<accession>`.
#' @param downloader Function used to retrieve the URL. Defaults to
#'   [utils::download.file()]. This is injectable so tests can avoid network
#'   calls.
#'
#' @return `sra_runinfo_url()` returns a length-1 character URL.
#'   `download_sra_runinfo()` returns the output path, invisibly.
#'   `read_sra_runinfo()` returns a `data.frame`.
#'
#' @examples
#' sra_runinfo_url("PRJNA400134")
#'
#' \dontrun{
#' download_sra_runinfo("PRJNA400134")
#' metadata <- read_sra_runinfo("PRJNA400134")
#' }
#'
#' @name sra_runinfo
NULL

#' @rdname sra_runinfo
#' @export
sra_runinfo_url <- function(accession, mode = c("acc", "term")) {
  mode <- match.arg(mode)
  stopifnot(is.character(accession), length(accession) == 1L, nzchar(accession))

  query <- utils::URLencode(accession, reserved = TRUE)
  parameter <- if (identical(mode, "term")) "term" else "acc"

  paste0(
    "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?",
    parameter,
    "=",
    query
  )
}

#' @rdname sra_runinfo
#' @export
download_sra_runinfo <- function(
  accession,
  dest = file.path("inst", "extdata", paste0(accession, "_runinfo.csv")),
  overwrite = FALSE,
  quiet = TRUE,
  mode = c("acc", "term"),
  downloader = utils::download.file
) {
  mode <- match.arg(mode)
  stopifnot(is.character(dest), length(dest) == 1L, nzchar(dest))

  if (file.exists(dest) && !isTRUE(overwrite)) {
    stop(
      "Refusing to overwrite existing file: ",
      dest,
      call. = FALSE
    )
  }

  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  url <- sra_runinfo_url(accession, mode = mode)
  status <- downloader(url = url, destfile = dest, quiet = quiet)

  if (!file.exists(dest) || (is.numeric(status) && !identical(status, 0L))) {
    stop("Failed to download SRA run metadata for ", accession, call. = FALSE)
  }

  invisible(dest)
}

#' @rdname sra_runinfo
#' @param path Path to the run metadata CSV. Defaults to the bundled file under
#'   `inst/extdata/<accession>_runinfo.csv`.
#' @export
read_sra_runinfo <- function(
  accession,
  path = system.file(
    "extdata",
    paste0(accession, "_runinfo.csv"),
    package = "RWisecondorX"
  )
) {
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "Run metadata file not found for ", accession,
      ". Supply `path` explicitly or stage the CSV under inst/extdata/.",
      call. = FALSE
    )
  }

  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}


.onLoad <- function(libname, pkgname) {
  S7::methods_register()

  ns <- asNamespace(pkgname)
  register_compat_method <- function(generic, cls, fun) {
    registerS3method(generic, paste0(pkgname, "::", cls), fun, envir = ns)
  }

  register_compat_method("$", "NIPTSample", .nipt_sample_dollar)
  register_compat_method("$", "CombinedStrandsSample", .nipt_sample_dollar)
  register_compat_method("$", "SeparatedStrandsSample", .nipt_sample_dollar)
  register_compat_method("[[", "NIPTSample", .nipt_sample_subset2)
  register_compat_method("[[", "CombinedStrandsSample", .nipt_sample_subset2)
  register_compat_method("[[", "SeparatedStrandsSample", .nipt_sample_subset2)

  register_compat_method("$", "NIPTControlGroup", .nipt_control_group_dollar)
  register_compat_method("$", "CombinedControlGroup", .nipt_control_group_dollar)
  register_compat_method("$", "SeparatedControlGroup", .nipt_control_group_dollar)
  register_compat_method("[[", "NIPTControlGroup", .nipt_control_group_subset2)
  register_compat_method("[[", "CombinedControlGroup", .nipt_control_group_subset2)
  register_compat_method("[[", "SeparatedControlGroup", .nipt_control_group_subset2)

}
