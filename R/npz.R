#' Convert BAM/CRAM to WisecondorX NPZ format
#'
#' Runs [bam_convert()] and serialises the resulting bin-count list to a
#' `.npz` file that is byte-compatible with the file written by
#' `wisecondorx convert`.  The NPZ must be created by numpy so that
#' `wisecondorx newref` can load it; this function therefore requires
#' `reticulate` and a Python environment with `numpy` installed (any version
#' ≥ 1.16 works — the format is stable).
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param npz Path for the output `.npz` file (created or overwritten).
#' @param binsize Bin size in base pairs (default 5000).
#' @param rmdup Duplicate-removal strategy passed to [bam_convert()].
#' @param con Optional open DBI connection with duckhts already loaded.
#'   Passed through to [bam_convert()].
#' @param np Optional numpy module imported via `reticulate::import("numpy")`.
#'   If `NULL` (default) it is imported automatically.
#'
#' @return `npz` (invisibly).
#'
#' @seealso [bam_convert()], [wisecondorx_newref()], [wisecondorx_predict()]
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
                            rmdup   = c("streaming", "none", "flag"),
                            con     = NULL,
                            np      = NULL) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.character(npz), length(npz) == 1L, nzchar(npz))

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

  bins <- bam_convert(bam, binsize = binsize, rmdup = rmdup, con = con)

  # Keep only chromosomes that have data; match wisecondorx key format (string)
  arrays <- Filter(Negate(is.null), bins)
  if (length(arrays) == 0L) {
    stop("bam_convert returned no data for ", bam, call. = FALSE)
  }

  # Convert each integer vector to a numpy int32 array (wisecondorx expects int32)
  np_arrays <- lapply(arrays, function(v) {
    np$array(as.integer(v), dtype = np$int32)
  })

  # savez_compressed matches wisecondorx's output: compressed npz, keyword args
  # are the array names (chromosome keys "1"-"24").
  do.call(np$savez_compressed, c(list(npz), np_arrays))

  invisible(npz)
}
