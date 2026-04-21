# wisecondorx_npz.R — WisecondorX-compatible NPZ serialization
#
.wcx_upstream_n_bins <- function(chr_length, binsize) {
  as.integer(as.double(chr_length) / as.double(binsize)) + 1L
}

.compat_npz_writer_script <- function() {
  normalizePath(
    system.file("extdata", "write_compat_npz.py",
                package = "RWisecondorX", mustWork = TRUE),
    winslash = "/",
    mustWork = TRUE
  )
}

#' Convert BAM/CRAM to WisecondorX NPZ format
#'
#' Runs [bam_convert()] and serialises the resulting bin-count list to a
#' `.npz` file that is byte-compatible with the file written by
#' `wisecondorx convert`.  The NPZ must be created by numpy so that
#' `wisecondorx newref` can load it; this function therefore requires
#' `reticulate` and a Python environment with `numpy` installed; any numpy
#' version from 1.16 onward works.
#'
#' This function exists for Python WisecondorX CLI conformance and
#' interoperability. The native R pipeline uses [bam_convert_bed()] together
#' with [rwisecondorx_newref()] and [rwisecondorx_predict()].
#'
#' The resulting NPZ has three top-level keys:
#' \describe{
#'   \item{`sample`}{0-d object array wrapping a dict mapping `"1"`..`"24"`
#'     to int32 arrays (chromosome bin counts).}
#'   \item{`binsize`}{Scalar int (bin size in bp).}
#'   \item{`quality`}{0-d object array wrapping a dict of QC counters
#'     (populated with zeros since the native binning kernel does not track
#'     per-read filter stats).}
#' }
#'
#' When the writer's numpy is >= 2.0 the internal pickle references
#' `numpy._core`, which older numpy (< 2.0) cannot unpickle.  This function
#' patches the pickle bytestream to use `numpy.core` so the NPZ is readable
#' by both numpy 1.x and 2.x.
#'
#' Native [bam_convert()] returns dense vectors covering the chromosome span in
#' the BAM header. Upstream `wisecondorx convert` has one additional historical
#' quirk: it allocates `int(length / binsize + 1)` bins for every chromosome in
#' the header, so chromosomes whose lengths are exact multiples of `binsize`
#' carry one extra trailing all-zero bin. This function pads to that upstream
#' NPZ layout during serialisation so Python `wisecondorx newref/predict`
#' remains byte-compatible, while the native R list/BED paths keep the cleaner
#' header-span contract.
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param reference Optional FASTA reference path for CRAM inputs. Passed to
#'   [bam_convert()].
#' @param npz Path for the output `.npz` file (created or overwritten).
#' @param binsize Bin size in base pairs (default 5000).
#' @param mapq Minimum mapping quality passed to [bam_convert()].
#' @param require_flags Integer bitmask of required SAM flags passed to
#'   [bam_convert()].
#' @param exclude_flags Integer bitmask of excluded SAM flags passed to
#'   [bam_convert()].
#' @param rmdup Duplicate-removal strategy passed to [bam_convert()].
#' @param con Optional open DBI connection with duckhts already loaded.
#'   Passed through to [bam_convert()].
#' @param np Optional numpy module imported via `reticulate::import("numpy")`.
#'   If `NULL` (default) it is imported automatically.
#'
#' @return `npz` (invisibly).
#'
#' @seealso [bam_convert()], [bam_convert_bed()], [rwisecondorx_newref()],
#'   [rwisecondorx_predict()], [wisecondorx_newref()], [wisecondorx_predict()]
#'
#' @examples
#' \dontrun{
#' bam_convert_npz("sample.bam", "sample.npz", binsize = 5000, rmdup = "streaming")
#' }
#'
#' @export
bam_convert_npz <- function(bam,
                            npz,
                            binsize = 5000L,
                            mapq = 1L,
                            require_flags = 0L,
                            exclude_flags = 0L,
                            rmdup   = c("streaming", "none", "flag"),
                            con     = NULL,
                            np      = NULL,
                            reference = NULL) {
  # This function is needed ONLY for Python WisecondorX CLI conformance.
  # For the native R pipeline, use bam_convert_bed() + rwisecondorx_newref().

  rmdup <- match.arg(rmdup)
  stopifnot(is.character(npz), length(npz) == 1L, nzchar(npz))
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  stopifnot(is.numeric(mapq), length(mapq) == 1L, mapq >= 0L)
  stopifnot(is.numeric(require_flags), length(require_flags) == 1L, require_flags >= 0L)
  stopifnot(is.numeric(exclude_flags), length(exclude_flags) == 1L, exclude_flags >= 0L)
  binsize <- as.integer(binsize)
  mapq <- as.integer(mapq)
  require_flags <- as.integer(require_flags)
  exclude_flags <- as.integer(exclude_flags)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to write NPZ files. Install it with: install.packages('reticulate')",
         call. = FALSE)
  }

  if (is.null(np)) {
    np <- tryCatch(
      reticulate::import("numpy", convert = FALSE),
      error = function(e) {
        stop("numpy not found in the active Python environment. ",
             "Install it with: reticulate::py_install('numpy')",
             call. = FALSE)
      }
    )
  }

  # We need the DuckDB connection for both bam_convert() and the BAM header

  # query, so manage it here rather than letting bam_convert() create its own.
  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  bins <- bam_convert(
    bam = bam,
    reference = reference,
    binsize = binsize,
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    rmdup = rmdup,
    con = con
  )

  if (all(vapply(bins, is.null, logical(1L)))) {
    stop("bam_convert returned no chr1-22/X/Y bins for ", bam, call. = FALSE)
  }

  # Upstream convert_reads() allocates int(length / binsize + 1) bins for every
  # chromosome in the BAM header, even chromosomes with 0 reads. Mirror that
  # quirk here so Python wisecondorx newref/predict can consume our NPZs
  # without shape mismatches, while bam_convert() itself keeps the denser
  # native header-span contract.
  chr_lengths <- .bam_chr_lengths(con, bam)
  for (nm in names(bins)) {
    chr_length <- chr_lengths[[nm]]
    if (!is.null(chr_length)) {
      n_bins <- .wcx_upstream_n_bins(chr_length, binsize)
      if (is.null(bins[[nm]])) {
        bins[[nm]] <- integer(n_bins)
      } else if (length(bins[[nm]]) < n_bins) {
        bins[[nm]] <- c(as.integer(bins[[nm]]),
                        integer(n_bins - length(bins[[nm]])))
      } else if (length(bins[[nm]]) > n_bins) {
        stop("bam_convert produced more bins than the WisecondorX NPZ layout for chr ",
             nm, call. = FALSE)
      }
    } else if (is.null(bins[[nm]])) {
      next
    } else {
      stop("bam_convert returned bins for chromosome ", nm,
           " that is absent from the BAM header", call. = FALSE)
    }
  }

  for (nm in names(bins)) {
    if (!is.null(bins[[nm]])) {
      bins[[nm]] <- as.integer(bins[[nm]])
    }
  }

  # Chromosomes still NULL are absent from the BAM header entirely and should
  # remain Python None in the serialized dict.
  py_builtins <- reticulate::import_builtins(convert = FALSE)
  sample_dict <- py_builtins$dict()
  for (nm in names(bins)) {
    if (is.null(bins[[nm]])) {
      sample_dict[nm] <- py_builtins$None
    } else {
      sample_dict[nm] <- np$array(bins[[nm]], dtype = np$int32)
    }
  }

  # Quality info -- upstream stores read-level QC stats.  We populate with
  # zeros since the native binning kernel does not track per-read filter counters.
  # These fields are never read by newref or predict.
  quality_dict <- py_builtins$dict()
  for (qk in c("mapped", "unmapped", "no_coordinate", "filter_rmdup",
               "filter_mapq", "pre_retro", "post_retro", "pair_fail")) {
    quality_dict[qk] <- 0L
  }

  # Write the NPZ using a packaged Python helper that produces pickle
  # bytestreams compatible with both numpy 1.x and 2.x. numpy >= 2.0 changed
  # its internal module layout (numpy._core vs numpy.core); when
  # savez_compressed pickles a dict-valued 0-d object array, the pickle
  # references the writer's module path. A numpy 1.x reader cannot unpickle
  # numpy._core references, so the helper patches the pickle bytes to
  # numpy.core for cross-version readability.
  reticulate::py_run_file(.compat_npz_writer_script(), local = FALSE)

  writer <- reticulate::py_eval("_write_compat_npz", convert = FALSE)
  writer(npz, sample_dict, quality_dict, as.integer(binsize))

  invisible(npz)
}
