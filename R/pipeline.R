#' Build a WisecondorX reference panel
#'
#' Calls `wisecondorx newref` via [condathis::run()] on a set of NPZ files
#' produced by [bam_convert_npz()].  The NPZ files must all use the same
#' `binsize`; the reference is built at `ref_binsize` (must be a multiple of
#' `binsize`).
#'
#' @param npz_files Character vector of paths to `.npz` files (one per sample).
#' @param output Path for the output reference `.npz` file.
#' @param binsize Convert-step bin size in base pairs (default 5000).
#'   Must match the `binsize` used when calling [bam_convert_npz()].
#' @param ref_binsize Reference bin size in base pairs (default 50000).
#'   Must be a multiple of `binsize`.
#' @param cpus Number of CPUs to pass to `wisecondorx newref` (default 1).
#' @param env_name Name of the conda environment containing `wisecondorx`
#'   (default `"wisecondorx"`).  Created automatically by condathis on first
#'   use via the `bioconda` channel.
#' @param extra_args Character vector of additional arguments passed verbatim
#'   to `wisecondorx newref` (e.g. `c("--yfrac", "0.4")`).
#'
#' @return `output` (invisibly).
#'
#' @seealso [bam_convert_npz()], [wisecondorx_predict()]
#'
#' @examples
#' \dontrun{
#' # Build a reference from 30 NIPT controls
#' wisecondorx_newref(
#'   npz_files = list.files("controls/", "\\.npz$", full.names = TRUE),
#'   output    = "reference.npz",
#'   ref_binsize = 50000L
#' )
#' }
#'
#' @export
wisecondorx_newref <- function(npz_files,
                               output,
                               binsize     = 5000L,
                               ref_binsize = 50000L,
                               cpus        = 1L,
                               env_name    = "wisecondorx",
                               extra_args  = character(0)) {
  .check_condathis()
  stopifnot(is.character(npz_files), length(npz_files) >= 1L)
  stopifnot(all(file.exists(npz_files)))
  stopifnot(is.character(output), length(output) == 1L, nzchar(output))
  stopifnot(is.numeric(ref_binsize), ref_binsize %% binsize == 0)

  .ensure_wisecondorx_env(env_name)

  # condathis::run does not expand relative paths; use absolute paths
  npz_abs <- normalizePath(npz_files, mustWork = TRUE)
  out_abs <- normalizePath(output, mustWork = FALSE)

  args <- c("newref", npz_abs, out_abs,
            "--binsize", as.character(as.integer(ref_binsize)),
            "--cpus", as.character(as.integer(cpus)),
            extra_args)

  do.call(condathis::run, c(list("wisecondorx"), as.list(args),
                            list(env_name = env_name)))
  invisible(output)
}


#' Predict copy-number aberrations with WisecondorX
#'
#' Calls `wisecondorx predict` via [condathis::run()] on a single-sample NPZ
#' file against a reference panel built by [wisecondorx_newref()].
#'
#' @param npz Path to the sample `.npz` file (from [bam_convert_npz()]).
#' @param ref Path to the reference `.npz` file (from [wisecondorx_newref()]).
#' @param output_prefix Output prefix for all `wisecondorx predict` output
#'   files (e.g. `"results/sample1"`).
#' @param ref_binsize Reference bin size in base pairs (default 50000).
#'   Must match the `ref_binsize` used when building the reference.
#' @param zscore Z-score threshold for aberration calling (default 5).
#' @param bed Logical; also write BED output files (default `FALSE`).
#' @param gender_file Optional path to a gender model file produced by
#'   `wisecondorx gender` (passed as `--gender`).
#' @param env_name Name of the conda environment containing `wisecondorx`
#'   (default `"wisecondorx"`).
#' @param extra_args Character vector of additional arguments passed verbatim
#'   to `wisecondorx predict`.
#'
#' @return `output_prefix` (invisibly).
#'
#' @seealso [bam_convert_npz()], [wisecondorx_newref()]
#'
#' @examples
#' \dontrun{
#' wisecondorx_predict(
#'   npz    = "sample.npz",
#'   ref    = "reference.npz",
#'   output_prefix = "results/sample"
#' )
#' }
#'
#' @export
wisecondorx_predict <- function(npz,
                                ref,
                                output_prefix,
                                ref_binsize = 50000L,
                                zscore      = 5,
                                bed         = FALSE,
                                gender_file = NULL,
                                env_name    = "wisecondorx",
                                extra_args  = character(0)) {
  .check_condathis()
  stopifnot(is.character(npz), length(npz) == 1L, file.exists(npz))
  stopifnot(is.character(ref), length(ref) == 1L, file.exists(ref))
  stopifnot(is.character(output_prefix), length(output_prefix) == 1L, nzchar(output_prefix))

  .ensure_wisecondorx_env(env_name)

  npz_abs <- normalizePath(npz, mustWork = TRUE)
  ref_abs <- normalizePath(ref, mustWork = TRUE)
  out_abs <- normalizePath(output_prefix, mustWork = FALSE)

  args <- c("predict", npz_abs, ref_abs, out_abs,
            "--binsize", as.character(as.integer(ref_binsize)),
            "--zscore", as.character(zscore))

  if (isTRUE(bed)) args <- c(args, "--bed")
  if (!is.null(gender_file)) {
    gender_abs <- normalizePath(gender_file, mustWork = TRUE)
    args <- c(args, "--gender", gender_abs)
  }
  args <- c(args, extra_args)

  do.call(condathis::run, c(list("wisecondorx"), as.list(args),
                            list(env_name = env_name)))
  invisible(output_prefix)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.check_condathis <- function() {
  if (!requireNamespace("condathis", quietly = TRUE)) {
    stop(
      "condathis is required for wisecondorx_newref() and wisecondorx_predict(). ",
      "Install it with: install.packages('condathis')",
      call. = FALSE
    )
  }
}

.ensure_wisecondorx_env <- function(env_name) {
  # Create the environment if it does not exist yet.
  # condathis::create_env() is idempotent if the env already exists.
  tryCatch(
    condathis::create_env(
      packages  = "wisecondorx",
      channels  = c("bioconda", "conda-forge"),
      env_name  = env_name
    ),
    error = function(e) {
      stop(
        "Failed to create conda environment '", env_name, "' with wisecondorx.\n",
        "Make sure condathis can reach the bioconda channel.\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}
