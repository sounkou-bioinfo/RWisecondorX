# Suppress R CMD check NOTEs for symbols used inside evalq() in the mclust
# namespace (nipter_sex.R and rwisecondorx_utils.R).
utils::globalVariables(c("Mclust", "emControl", "data"))

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
