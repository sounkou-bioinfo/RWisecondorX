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
  if (!file.exists(npz)) {
    stop(
      "wisecondorx convert did not create the expected output file: ",
      npz,
      call. = FALSE
    )
  }
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
#' @param cpus Number of CPUs to pass to `wisecondorx newref` (default 4).
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
#'   cpus = 4L
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
                               cpus        = 4L,
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

  .validate_wisecondorx_predict_threshold(ref, minrefbins)
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

.wisecondorx_env_state <- new.env(parent = emptyenv())

.condathis_safe_wd <- function() {
  wd <- tryCatch(getwd(), error = function(e) NULL)
  wd_ok <- is.character(wd) &&
    length(wd) == 1L &&
    !is.na(wd) &&
    nzchar(wd) &&
    isTRUE(tryCatch(dir.exists(wd), error = function(e) FALSE))
  if (wd_ok) {
    return(wd)
  }

  fallback <- tempdir()
  dir.create(fallback, recursive = TRUE, showWarnings = FALSE)
  fallback
}

.with_condathis_cache <- function(expr) {
  # condathis resets most conda/mamba env vars itself, but libmamba still uses
  # HOME/XDG cache roots for proc-lock state. Isolating both avoids stale
  # ~/.cache/mamba/proc/proc.lock failures from unrelated sessions.
  old_xdg_cache <- Sys.getenv("XDG_CACHE_HOME", unset = NA_character_)
  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_xdg_data <- Sys.getenv("XDG_DATA_HOME", unset = NA_character_)
  old_xdg_config <- Sys.getenv("XDG_CONFIG_HOME", unset = NA_character_)
  old_wd <- tryCatch(getwd(), error = function(e) NA_character_)
  safe_wd <- .condathis_safe_wd()
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

  # Register cleanup BEFORE mutating state, so an interrupt between Sys.setenv
  # and on.exit does not permanently corrupt the process environment.
  on.exit(
    {
      if (is.na(old_home)) Sys.unsetenv("HOME") else Sys.setenv(HOME = old_home)
      if (is.na(old_xdg_cache)) Sys.unsetenv("XDG_CACHE_HOME") else Sys.setenv(XDG_CACHE_HOME = old_xdg_cache)
      if (is.na(old_xdg_data)) Sys.unsetenv("XDG_DATA_HOME") else Sys.setenv(XDG_DATA_HOME = old_xdg_data)
      if (is.na(old_xdg_config)) Sys.unsetenv("XDG_CONFIG_HOME") else Sys.setenv(XDG_CONFIG_HOME = old_xdg_config)
      if (!is.na(old_wd) && nzchar(old_wd) && dir.exists(old_wd)) {
        setwd(old_wd)
      }
    },
    add = TRUE
  )

  Sys.setenv(
    HOME = cache_home,
    XDG_CACHE_HOME = cache_dir,
    XDG_DATA_HOME = resolved_xdg_data,
    XDG_CONFIG_HOME = resolved_xdg_config
  )
  if (is.na(old_wd) || !identical(old_wd, safe_wd)) {
    setwd(safe_wd)
  }
  force(expr)
}

.bundled_wisecondorx_cbs_script <- function() {
  path <- system.file("extdata", "wisecondorx_CBS.R", package = "RWisecondorX")
  if (is.character(path) && length(path) == 1L && nzchar(path) && file.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  ns_path <- tryCatch(
    getNamespaceInfo(asNamespace("RWisecondorX"), "path"),
    error = function(e) NULL
  )
  if (is.character(ns_path) && length(ns_path) == 1L && nzchar(ns_path)) {
    cand <- file.path(ns_path, "inst", "extdata", "wisecondorx_CBS.R")
    if (file.exists(cand)) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Bundled WisecondorX CBS.R helper is missing from the RWisecondorX package.",
    call. = FALSE
  )
}

.wisecondorx_site_package_dir <- function(env_name) {
  env_dir <- condathis::get_env_dir(env_name)
  hits <- Sys.glob(file.path(
    env_dir,
    "lib",
    "python*",
    "site-packages",
    "wisecondorx",
    "__init__.py"
  ))
  if (!length(hits)) {
    stop(
      "Could not locate the installed wisecondorx Python package under conda env '",
      env_name, "'.",
      call. = FALSE
    )
  }
  normalizePath(dirname(hits[[1L]]), winslash = "/", mustWork = TRUE)
}

.prepare_wisecondorx_predict_shadow <- function(env_name) {
  key <- paste0("predict_shadow__", env_name)
  cached <- get0(key, envir = .wisecondorx_env_state, inherits = FALSE)
  if (is.character(cached) && length(cached) == 1L &&
      nzchar(cached) && dir.exists(file.path(cached, "wisecondorx"))) {
    return(cached)
  }

  pkg_dir <- .wisecondorx_site_package_dir(env_name)
  cbs_src <- .bundled_wisecondorx_cbs_script()
  shadow_root <- file.path(tempdir(), sprintf("rwisecondorx-wcx-shadow-%s", env_name))
  shadow_pkg_root <- file.path(shadow_root, "wisecondorx")

  if (dir.exists(shadow_root)) {
    unlink(shadow_root, recursive = TRUE, force = TRUE)
  }
  dir.create(shadow_root, recursive = TRUE, showWarnings = FALSE)

  ok <- file.copy(pkg_dir, shadow_root, recursive = TRUE)
  if (!isTRUE(ok) || !dir.exists(shadow_pkg_root)) {
    stop(
      "Failed to create temporary shadow copy of the wisecondorx Python package.",
      call. = FALSE
    )
  }

  include_dir <- file.path(shadow_pkg_root, "include")
  dir.create(include_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(cbs_src, file.path(include_dir, "CBS.R"), overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to stage bundled CBS.R into the temporary wisecondorx shadow package.",
         call. = FALSE)
  }

  assign(key, shadow_root, envir = .wisecondorx_env_state)
  shadow_root
}

.with_wisecondorx_predict_patch <- function(expr, env_name) {
  shadow_root <- .prepare_wisecondorx_predict_shadow(env_name)
  mpl_dir <- file.path(tempdir(), sprintf("rwisecondorx-mpl-%s", env_name))
  dir.create(mpl_dir, recursive = TRUE, showWarnings = FALSE)

  old_pythonpath <- Sys.getenv("PYTHONPATH", unset = NA_character_)
  old_mplconfig <- Sys.getenv("MPLCONFIGDIR", unset = NA_character_)
  on.exit(
    {
      if (is.na(old_pythonpath)) {
        Sys.unsetenv("PYTHONPATH")
      } else {
        Sys.setenv(PYTHONPATH = old_pythonpath)
      }
      if (is.na(old_mplconfig)) {
        Sys.unsetenv("MPLCONFIGDIR")
      } else {
        Sys.setenv(MPLCONFIGDIR = old_mplconfig)
      }
    },
    add = TRUE
  )

  pythonpath <- shadow_root
  if (!is.na(old_pythonpath) && nzchar(old_pythonpath)) {
    pythonpath <- paste(shadow_root, old_pythonpath, sep = .Platform$path.sep)
  }
  Sys.setenv(
    PYTHONPATH = pythonpath,
    MPLCONFIGDIR = mpl_dir
  )
  force(expr)
}

.ensure_wisecondorx_env <- function(env_name) {
  if (isTRUE(get0(env_name, envir = .wisecondorx_env_state, inherits = FALSE))) {
    return(invisible(NULL))
  }

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
  assign(env_name, TRUE, envir = .wisecondorx_env_state)
  invisible(NULL)
}

.run_wisecondorx_cli <- function(args, env_name) {
  runner <- function() {
    .with_condathis_cache(
      do.call(
        condathis::run,
        c(list("wisecondorx"), as.list(args), list(env_name = env_name))
      )
    )
  }

  if (length(args) >= 1L && identical(args[[1L]], "predict")) {
    .with_wisecondorx_predict_patch(runner(), env_name = env_name)
  } else {
    runner()
  }
}

.wisecondorx_convert_args <- function(bam,
                                      npz,
                                      reference,
                                      binsize,
                                      normdup,
                                      extra_args = character(0)) {
  bam_abs <- normalizePath(bam, mustWork = TRUE)
  npz_abs <- normalizePath(npz, mustWork = FALSE)
  npz_prefix <- sub("\\.npz$", "", npz_abs, ignore.case = TRUE)

  args <- c("convert", bam_abs, npz_prefix)
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

.validate_wisecondorx_predict_threshold <- function(ref, minrefbins) {
  stopifnot(is.numeric(minrefbins), length(minrefbins) == 1L, minrefbins >= 1L)

  refsize <- .wisecondorx_npz_refsize(ref)
  if (is.na(refsize) || minrefbins <= refsize) {
    return(invisible(NULL))
  }

  stop(
    sprintf(
      "Invalid WisecondorX predict settings: minrefbins=%d exceeds the reference refsize=%d. ",
      as.integer(minrefbins),
      as.integer(refsize)
    ),
    "No target bin can retain enough reference bins under this combination. ",
    "Lower 'minrefbins' or rebuild the reference with a larger 'refsize'.",
    call. = FALSE
  )
}

.wisecondorx_npz_refsize <- function(ref) {
  stopifnot(is.character(ref), length(ref) == 1L, file.exists(ref))

  entries <- utils::unzip(ref, list = TRUE)
  if (!nrow(entries)) {
    return(NA_integer_)
  }

  hit <- entries$Name[grepl("^indexes(\\.[FM])?\\.npy$", basename(entries$Name))]
  if (!length(hit)) {
    return(NA_integer_)
  }

  dims <- vapply(hit, function(name) {
    con <- unz(ref, name, open = "rb")
    on.exit(close(con), add = TRUE)
    shape <- .npy_shape(con)
    if (!length(shape)) {
      return(NA_integer_)
    }
    if (length(shape) == 1L) {
      return(as.integer(shape[[1L]]))
    }
    as.integer(shape[[2L]])
  }, integer(1L))

  dims <- dims[is.finite(dims) & dims > 0L]
  if (!length(dims)) {
    return(NA_integer_)
  }
  max(dims)
}

.npy_shape <- function(con) {
  magic <- readBin(con, what = "raw", n = 6L)
  if (length(magic) != 6L ||
      !identical(magic[[1L]], as.raw(0x93)) ||
      !identical(rawToChar(magic[2:6]), "NUMPY")) {
    return(integer())
  }

  ver <- readBin(con, what = "integer", size = 1L, n = 2L,
                 signed = FALSE, endian = "little")
  if (length(ver) != 2L) {
    return(integer())
  }
  major <- ver[[1L]]
  header_len_size <- if (major <= 1L) 2L else 4L
  header_len <- readBin(con, what = "integer", size = header_len_size, n = 1L,
                        signed = FALSE, endian = "little")
  if (!length(header_len) || !is.finite(header_len) || header_len <= 0L) {
    return(integer())
  }

  header <- rawToChar(readBin(con, what = "raw", n = header_len), multiple = FALSE)
  m <- regexec("shape[^\\(]*\\(([^\\)]*)\\)", header)
  reg <- regmatches(header, m)[[1L]]
  if (length(reg) < 2L) {
    return(integer())
  }

  parts <- strsplit(reg[[2L]], ",", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (!length(parts)) {
    return(integer())
  }

  suppressWarnings(as.integer(parts))
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
