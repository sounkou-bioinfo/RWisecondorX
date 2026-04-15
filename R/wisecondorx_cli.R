#' Convert BAM/CRAM to WisecondorX NPZ format (upstream CLI wrapper)
#'
#' Calls `wisecondorx convert` via [condathis::run()] to convert an aligned
#' BAM or CRAM file to a `.npz` file for downstream use with
#' [wisecondorx_newref()] or [wisecondorx_predict()].
#'
#' This wrapper exposes the upstream `wisecondorx convert` CLI flags:
#' `--reference`, `--binsize`, and `--normdup`.
#'
#' For a fully native R implementation (no Python dependency) that uses
#' Rduckhts instead of pysam, see [bam_convert()] and [bam_convert_npz()].
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param npz Path for the output `.npz` file (created or overwritten).
#' @param reference Optional path to a FASTA reference file. Required when
#'   `bam` is a CRAM file. Passed to upstream `--reference`.
#' @param binsize Bin size in base pairs (default 5000). The reference bin
#'   size should be a multiple of this value. Passed to upstream `--binsize`.
#' @param normdup Logical; when `TRUE`, passes `--normdup` to skip duplicate
#'   removal. Recommended for NIPT data where read depth is low.
#' @param env_name Name of the conda environment containing `wisecondorx`
#'   (default `"wisecondorx"`). Created automatically by condathis on first
#'   use via the `bioconda` channel.
#' @param extra_args Character vector of additional arguments passed verbatim
#'   after the mapped CLI flags. For forward compatibility with future upstream
#'   WisecondorX releases.
#'
#' @return `npz` (invisibly).
#'
#' @seealso [bam_convert()], [bam_convert_npz()], [wisecondorx_newref()],
#'   [wisecondorx_predict()]
#'
#' @examples
#' \dontrun{
#' wisecondorx_convert(
#'   bam = "sample.bam",
#'   npz = "sample.npz",
#'   binsize = 5000L,
#'   normdup = FALSE
#' )
#' }
#'
#' @export
wisecondorx_convert <- function(bam,
                                npz,
                                reference  = NULL,
                                binsize    = 5000L,
                                normdup    = FALSE,
                                env_name   = "wisecondorx",
                                extra_args = character(0)) {
  .check_condathis()
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  stopifnot(is.character(npz), length(npz) == 1L, nzchar(npz))
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  stopifnot(is.logical(normdup), length(normdup) == 1L)
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }

  .ensure_wisecondorx_env(env_name)

  args <- .wisecondorx_convert_args(
    bam        = bam,
    npz        = npz,
    reference  = reference,
    binsize    = binsize,
    normdup    = normdup,
    extra_args = extra_args
  )

  .run_wisecondorx_cli(args, env_name)
  invisible(npz)
}


#' Build a WisecondorX reference panel
#'
#' Calls `wisecondorx newref` via [condathis::run()] on a set of NPZ files
#' produced by [bam_convert_npz()]. The NPZ files must all use the same
#' `binsize`; the reference is built at `ref_binsize` (must be a multiple of
#' `binsize`).
#'
#' This wrapper exposes the current upstream `wisecondorx newref` CLI flags
#' documented in the WisecondorX README: `--nipt`, `--binsize`, `--refsize`,
#' `--yfrac`, `--plotyfrac`, and `--cpus`.
#'
#' @param npz_files Character vector of paths to `.npz` files (one per sample).
#' @param output Path for the output reference `.npz` file.
#' @param binsize Convert-step bin size in base pairs (default 5000).
#'   Must match the `binsize` used when calling [bam_convert_npz()].
#' @param ref_binsize Reference bin size in base pairs (default 100000).
#'   Passed to upstream `--binsize`. Must be a multiple of `binsize`.
#' @param nipt Logical; pass upstream `--nipt`.
#' @param refsize Number of reference locations per target bin. Passed to
#'   upstream `--refsize`.
#' @param yfrac Optional numeric Y-read fraction cutoff. Passed to upstream
#'   `--yfrac`.
#' @param plotyfrac Optional output path for the Y-fraction histogram and
#'   mixture-model plot. Passed to upstream `--plotyfrac`.
#' @param cpus Number of CPUs to pass to `wisecondorx newref` (default 1).
#' @param env_name Name of the conda environment containing `wisecondorx`
#'   (default `"wisecondorx"`). Created automatically by condathis on first
#'   use via the `bioconda` channel.
#' @param extra_args Character vector of additional arguments passed verbatim
#'   after the mapped CLI flags. Keep this for forward compatibility with future
#'   upstream WisecondorX releases.
#'
#' @return `output` (invisibly).
#'
#' @seealso [bam_convert_npz()], [wisecondorx_predict()]
#'
#' @examples
#' \dontrun{
#' wisecondorx_newref(
#'   npz_files = list.files("controls/", "\\.npz$", full.names = TRUE),
#'   output = "reference.npz",
#'   binsize = 5000L,
#'   ref_binsize = 100000L,
#'   nipt = TRUE,
#'   refsize = 300L,
#'   yfrac = 0.05,
#'   cpus = 1L
#' )
#' }
#'
#' @export
wisecondorx_newref <- function(npz_files,
                               output,
                               binsize     = 5000L,
                               ref_binsize = 100000L,
                               nipt        = FALSE,
                               refsize     = 300L,
                               yfrac       = NULL,
                               plotyfrac   = NULL,
                               cpus        = 1L,
                               env_name    = "wisecondorx",
                               extra_args  = character(0)) {
  .check_condathis()
  stopifnot(is.character(npz_files), length(npz_files) >= 1L)
  stopifnot(all(file.exists(npz_files)))
  stopifnot(is.character(output), length(output) == 1L, nzchar(output))
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  stopifnot(is.numeric(ref_binsize), length(ref_binsize) == 1L, ref_binsize >= 1L)
  stopifnot(ref_binsize %% binsize == 0)

  .ensure_wisecondorx_env(env_name)

  args <- .wisecondorx_newref_args(
    npz_files = npz_files,
    output = output,
    ref_binsize = ref_binsize,
    nipt = nipt,
    refsize = refsize,
    yfrac = yfrac,
    plotyfrac = plotyfrac,
    cpus = cpus,
    extra_args = extra_args
  )

  .run_wisecondorx_cli(args, env_name)
  invisible(output)
}


#' Predict copy-number aberrations with WisecondorX
#'
#' Calls `wisecondorx predict` via [condathis::run()] on a single-sample NPZ
#' file against a reference panel built by [wisecondorx_newref()].
#'
#' This wrapper exposes the current upstream `wisecondorx predict` CLI flags
#' documented in the WisecondorX README: `--minrefbins`, `--maskrepeats`,
#' `--zscore`, `--alpha`, `--beta`, `--blacklist`, `--gender`, `--bed`,
#' `--plot`, `--add-plot-title`, `--regions`, `--ylim`, `--cairo`, and
#' `--seed`.
#'
#' @param npz Path to the sample `.npz` file (from [bam_convert_npz()]).
#' @param ref Path to the reference `.npz` file (from [wisecondorx_newref()]).
#' @param output_prefix Output prefix for all `wisecondorx predict` output
#'   files (e.g. `"results/sample1"`).
#' @param minrefbins Minimum number of sensible reference bins per target bin.
#'   Passed to upstream `--minrefbins`.
#' @param maskrepeats Number of repeat-masking cycles. Passed to upstream
#'   `--maskrepeats`.
#' @param zscore Z-score threshold for aberration calling. Passed to upstream
#'   `--zscore`.
#' @param alpha P-value cutoff for circular binary segmentation breakpoints.
#'   Passed to upstream `--alpha`.
#' @param beta Optional ratio cutoff for aberration calling. Passed to upstream
#'   `--beta`. When set, upstream ignores `zscore`.
#' @param blacklist Optional path to a headerless BED blacklist file. Passed to
#'   upstream `--blacklist`.
#' @param gender Optional forced gender, `"F"` or `"M"`. Passed to upstream
#'   `--gender`.
#' @param bed Logical; also write BED output files. Passed to upstream `--bed`.
#' @param plot Logical; also write plot output files. Passed to upstream
#'   `--plot`.
#' @param regions Optional path to a headerless BED file with regions to mark on
#'   the plot. Passed to upstream `--regions`.
#' @param ylim Optional numeric vector of length 2 giving the y-axis interval,
#'   e.g. `c(-2, 2)`. Passed to upstream `--ylim` as `[a,b]`.
#' @param cairo Logical; use Cairo bitmap output. Passed to upstream `--cairo`.
#' @param seed Optional integer random seed for segmentation. Passed to
#'   upstream `--seed`.
#' @param add_plot_title Logical; when `TRUE`, adds the output basename as the
#'   plot title. Passed to upstream `--add-plot-title`. Only effective when
#'   `plot = TRUE`.
#' @param env_name Name of the conda environment containing `wisecondorx`
#'   (default `"wisecondorx"`).
#' @param extra_args Character vector of additional arguments passed verbatim
#'   after the mapped CLI flags. Keep this for forward compatibility with future
#'   upstream WisecondorX releases.
#'
#' @return `output_prefix` (invisibly).
#'
#' @seealso [bam_convert_npz()], [wisecondorx_newref()]
#'
#' @examples
#' \dontrun{
#' wisecondorx_predict(
#'   npz = "sample.npz",
#'   ref = "reference.npz",
#'   output_prefix = "results/sample",
#'   minrefbins = 150L,
#'   maskrepeats = 5L,
#'   zscore = 5,
#'   alpha = 1e-4,
#'   beta = NULL,
#'   blacklist = NULL,
#'   gender = "F",
#'   bed = TRUE,
#'   plot = TRUE,
#'   add_plot_title = TRUE,
#'   regions = NULL,
#'   ylim = c(-2, 2),
#'   cairo = FALSE,
#'   seed = 1L
#' )
#' }
#'
#' @export
wisecondorx_predict <- function(npz,
                                ref,
                                output_prefix,
                                minrefbins     = 150L,
                                maskrepeats    = 5L,
                                zscore         = 5,
                                alpha          = 1e-4,
                                beta           = NULL,
                                blacklist      = NULL,
                                gender         = NULL,
                                bed            = FALSE,
                                plot           = FALSE,
                                regions        = NULL,
                                ylim           = NULL,
                                cairo          = FALSE,
                                seed           = NULL,
                                add_plot_title = FALSE,
                                env_name       = "wisecondorx",
                                extra_args     = character(0)) {
  .check_condathis()
  stopifnot(is.character(npz), length(npz) == 1L, file.exists(npz))
  stopifnot(is.character(ref), length(ref) == 1L, file.exists(ref))
  stopifnot(is.character(output_prefix), length(output_prefix) == 1L, nzchar(output_prefix))

  .ensure_wisecondorx_env(env_name)

  args <- .wisecondorx_predict_args(
    npz = npz,
    ref = ref,
    output_prefix = output_prefix,
    minrefbins = minrefbins,
    maskrepeats = maskrepeats,
    zscore = zscore,
    alpha = alpha,
    beta = beta,
    blacklist = blacklist,
    gender = gender,
    bed = bed,
    plot = plot,
    regions = regions,
    ylim = ylim,
    cairo = cairo,
    seed = seed,
    add_plot_title = add_plot_title,
    extra_args = extra_args
  )

  .run_wisecondorx_cli(args, env_name)
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

.with_condathis_cache <- function(expr) {
  # condathis resets most conda/mamba env vars itself, but libmamba still uses
  # HOME/XDG cache roots for proc-lock state. Isolating both avoids stale
  # ~/.cache/mamba/proc/proc.lock failures from unrelated sessions.
  old_xdg_cache <- Sys.getenv("XDG_CACHE_HOME", unset = NA_character_)
  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_xdg_data <- Sys.getenv("XDG_DATA_HOME", unset = NA_character_)
  old_xdg_config <- Sys.getenv("XDG_CONFIG_HOME", unset = NA_character_)
  resolved_xdg_data <- if (!is.na(old_xdg_data) && nzchar(old_xdg_data)) {
    old_xdg_data
  } else {
    file.path(path.expand("~"), ".local", "share")
  }
  resolved_xdg_config <- if (!is.na(old_xdg_config) && nzchar(old_xdg_config)) {
    old_xdg_config
  } else {
    file.path(path.expand("~"), ".config")
  }
  cache_home <- file.path(tempdir(), "rwisecondorx-condathis-home")
  cache_dir <- file.path(cache_home, ".cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(
    HOME = cache_home,
    XDG_CACHE_HOME = cache_dir,
    XDG_DATA_HOME = resolved_xdg_data,
    XDG_CONFIG_HOME = resolved_xdg_config
  )
  on.exit(
    {
      if (is.na(old_home)) Sys.unsetenv("HOME") else Sys.setenv(HOME = old_home)
      if (is.na(old_xdg_cache)) Sys.unsetenv("XDG_CACHE_HOME") else Sys.setenv(XDG_CACHE_HOME = old_xdg_cache)
      if (is.na(old_xdg_data)) Sys.unsetenv("XDG_DATA_HOME") else Sys.setenv(XDG_DATA_HOME = old_xdg_data)
      if (is.na(old_xdg_config)) Sys.unsetenv("XDG_CONFIG_HOME") else Sys.setenv(XDG_CONFIG_HOME = old_xdg_config)
    },
    add = TRUE
  )
  force(expr)
}

.ensure_wisecondorx_env <- function(env_name) {
  # Create the environment if it does not exist yet.
  # condathis::create_env() is idempotent if the env already exists.
  tryCatch(
    .with_condathis_cache(
      condathis::create_env(
        packages  = "wisecondorx",
        channels  = c("bioconda", "conda-forge"),
        env_name  = env_name
      )
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

.run_wisecondorx_cli <- function(args, env_name) {
  .with_condathis_cache(
    do.call(
      condathis::run,
      c(list("wisecondorx"), as.list(args), list(env_name = env_name))
    )
  )
}

.wisecondorx_convert_args <- function(bam,
                                      npz,
                                      reference,
                                      binsize,
                                      normdup,
                                      extra_args = character(0)) {
  bam_abs <- normalizePath(bam, mustWork = TRUE)
  npz_abs <- normalizePath(npz, mustWork = FALSE)

  args <- c("convert", bam_abs, npz_abs)
  args <- .append_optional_path(args, "--reference", reference, must_work = TRUE)
  args <- .append_value(args, "--binsize", as.integer(binsize))
  args <- .append_flag(args, "--normdup", isTRUE(normdup))
  c(args, extra_args)
}

.wisecondorx_newref_args <- function(npz_files,
                                     output,
                                     ref_binsize,
                                     nipt,
                                     refsize,
                                     yfrac,
                                     plotyfrac,
                                     cpus,
                                     extra_args = character(0)) {
  npz_abs <- normalizePath(npz_files, mustWork = TRUE)
  out_abs <- normalizePath(output, mustWork = FALSE)

  stopifnot(is.logical(nipt), length(nipt) == 1L)
  stopifnot(is.numeric(refsize), length(refsize) == 1L, refsize >= 1L)
  stopifnot(is.numeric(cpus), length(cpus) == 1L, cpus >= 1L)

  args <- c("newref", npz_abs, out_abs)
  args <- .append_flag(args, "--nipt", isTRUE(nipt))
  args <- .append_value(args, "--binsize", as.integer(ref_binsize))
  args <- .append_value(args, "--refsize", as.integer(refsize))
  args <- .append_optional_scalar(args, "--yfrac", yfrac)
  args <- .append_optional_path(args, "--plotyfrac", plotyfrac, must_work = FALSE)
  args <- .append_value(args, "--cpus", as.integer(cpus))
  c(args, extra_args)
}

.wisecondorx_predict_args <- function(npz,
                                      ref,
                                      output_prefix,
                                      minrefbins,
                                      maskrepeats,
                                      zscore,
                                      alpha,
                                      beta,
                                      blacklist,
                                      gender,
                                      bed,
                                      plot,
                                      regions,
                                      ylim,
                                      cairo,
                                      seed,
                                      add_plot_title,
                                      extra_args = character(0)) {
  npz_abs <- normalizePath(npz, mustWork = TRUE)
  ref_abs <- normalizePath(ref, mustWork = TRUE)
  out_abs <- normalizePath(output_prefix, mustWork = FALSE)

  stopifnot(is.numeric(minrefbins), length(minrefbins) == 1L, minrefbins >= 1L)
  stopifnot(is.numeric(maskrepeats), length(maskrepeats) == 1L, maskrepeats >= 0L)
  stopifnot(is.numeric(zscore), length(zscore) == 1L)
  stopifnot(is.numeric(alpha), length(alpha) == 1L)
  stopifnot(is.logical(bed), length(bed) == 1L)
  stopifnot(is.logical(plot), length(plot) == 1L)
  stopifnot(is.logical(cairo), length(cairo) == 1L)
  stopifnot(is.logical(add_plot_title), length(add_plot_title) == 1L)

  if (!is.null(gender)) {
    stopifnot(is.character(gender), length(gender) == 1L, nzchar(gender))
    gender <- match.arg(toupper(gender), c("F", "M"))
  }

  if (!is.null(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim) == 2L)
  }

  args <- c("predict", npz_abs, ref_abs, out_abs)
  args <- .append_value(args, "--minrefbins", as.integer(minrefbins))
  args <- .append_value(args, "--maskrepeats", as.integer(maskrepeats))
  args <- .append_value(args, "--zscore", zscore)
  args <- .append_value(args, "--alpha", alpha)
  args <- .append_optional_scalar(args, "--beta", beta)
  args <- .append_optional_path(args, "--blacklist", blacklist, must_work = TRUE)
  args <- .append_optional_scalar(args, "--gender", gender)
  args <- .append_flag(args, "--bed", isTRUE(bed))
  args <- .append_flag(args, "--plot", isTRUE(plot))
  args <- .append_flag(args, "--add-plot-title", isTRUE(add_plot_title))
  args <- .append_optional_path(args, "--regions", regions, must_work = TRUE)
  args <- .append_optional_scalar(args, "--ylim", .format_ylim(ylim))
  args <- .append_flag(args, "--cairo", isTRUE(cairo))
  args <- .append_optional_scalar(args, "--seed", seed)
  c(args, extra_args)
}

.append_value <- function(args, flag, value) {
  c(args, flag, as.character(value))
}

.append_flag <- function(args, flag, include) {
  if (isTRUE(include)) c(args, flag) else args
}

.append_optional_scalar <- function(args, flag, value) {
  if (is.null(value)) return(args)
  c(args, flag, as.character(value))
}

.append_optional_path <- function(args, flag, path, must_work) {
  if (is.null(path)) return(args)
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))
  c(args, flag, normalizePath(path, mustWork = must_work))
}

.format_ylim <- function(ylim) {
  if (is.null(ylim)) return(NULL)
  sprintf("[%s,%s]", format(ylim[[1L]], trim = TRUE), format(ylim[[2L]], trim = TRUE))
}
