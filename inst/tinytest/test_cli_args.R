library(tinytest)
library(RWisecondorX)

# ---------------------------------------------------------------------------
# CLI argument builder tests
# ---------------------------------------------------------------------------

# Expose internal helpers for unit testing
.wisecondorx_convert_args <- RWisecondorX:::.wisecondorx_convert_args
.wisecondorx_newref_args <- RWisecondorX:::.wisecondorx_newref_args
.wisecondorx_predict_args <- RWisecondorX:::.wisecondorx_predict_args
.format_ylim <- RWisecondorX:::.format_ylim
.condathis_safe_wd <- RWisecondorX:::.condathis_safe_wd

# ---------------------------------------------------------------------------
# convert CLI argument mapping
# ---------------------------------------------------------------------------

tmp_bam <- tempfile(tmpdir = tempdir(), fileext = ".bam")
tmp_ref <- tempfile(tmpdir = tempdir(), fileext = ".fa")
tmp_npz_out <- file.path(tempdir(), "sample.npz")

file.create(tmp_bam)
file.create(tmp_ref)

# Minimal: no reference, no normdup
convert_args_basic <- .wisecondorx_convert_args(
  bam        = tmp_bam,
  npz        = tmp_npz_out,
  reference  = NULL,
  binsize    = 5000L,
  normdup    = FALSE,
  extra_args = character(0)
)

expect_identical(
  convert_args_basic,
  c(
    "convert",
    normalizePath(tmp_bam, mustWork = TRUE),
    normalizePath(tmp_npz_out, mustWork = FALSE),
    "--binsize", "5000"
  )
)

# With reference and normdup
convert_args_full <- .wisecondorx_convert_args(
  bam        = tmp_bam,
  npz        = tmp_npz_out,
  reference  = tmp_ref,
  binsize    = 10000L,
  normdup    = TRUE,
  extra_args = c("--extra", "val")
)

expect_identical(
  convert_args_full,
  c(
    "convert",
    normalizePath(tmp_bam, mustWork = TRUE),
    normalizePath(tmp_npz_out, mustWork = FALSE),
    "--reference", normalizePath(tmp_ref, mustWork = TRUE),
    "--binsize", "10000",
    "--normdup",
    "--extra", "val"
  )
)

# ---------------------------------------------------------------------------
# newref CLI argument mapping
# ---------------------------------------------------------------------------

tmp_dir <- tempdir()
npz1 <- tempfile(tmpdir = tmp_dir, fileext = ".npz")
npz2 <- tempfile(tmpdir = tmp_dir, fileext = ".npz")
plotyfrac <- tempfile(tmpdir = tmp_dir, fileext = ".png")
output_ref <- file.path(tmp_dir, "reference.npz")

file.create(npz1)
file.create(npz2)

newref_args <- .wisecondorx_newref_args(
  npz_files = c(npz1, npz2),
  output = output_ref,
  ref_binsize = 50000L,
  nipt = TRUE,
  refsize = 300L,
  yfrac = 0.05,
  plotyfrac = plotyfrac,
  cpus = 2L,
  extra_args = c("--extra-flag", "value")
)

expect_identical(
  newref_args,
  c(
    "newref",
    normalizePath(npz1, mustWork = TRUE),
    normalizePath(npz2, mustWork = TRUE),
    normalizePath(output_ref, mustWork = FALSE),
    "--nipt",
    "--binsize", "50000",
    "--refsize", "300",
    "--yfrac", "0.05",
    "--plotyfrac", normalizePath(plotyfrac, mustWork = FALSE),
    "--cpus", "2",
    "--extra-flag", "value"
  )
)

expect_identical(
  formals(wisecondorx_newref)$cpus,
  4L,
  info = "wisecondorx_newref defaults to 4 CPUs"
)

.with_deleted_wd <- function(code) {
  old_wd <- getwd()
  doomed_wd <- tempfile(tmpdir = tempdir())
  dir.create(doomed_wd, recursive = TRUE, showWarnings = FALSE)
  setwd(doomed_wd)
  unlink(doomed_wd, recursive = TRUE)
  on.exit(setwd(old_wd), add = TRUE)
  force(code)
}

fallback_wd <- .with_deleted_wd(.condathis_safe_wd())
expect_true(dir.exists(fallback_wd),
            info = ".condathis_safe_wd falls back to an existing directory")

# ---------------------------------------------------------------------------
# predict CLI argument mapping
# ---------------------------------------------------------------------------

blacklist <- tempfile(tmpdir = tmp_dir, fileext = ".bed")
regions <- tempfile(tmpdir = tmp_dir, fileext = ".bed")
sample_npz <- tempfile(tmpdir = tmp_dir, fileext = ".npz")
reference_npz <- tempfile(tmpdir = tmp_dir, fileext = ".npz")
output_prefix <- file.path(tmp_dir, "results", "sample")

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
file.create(blacklist)
file.create(regions)
file.create(sample_npz)
file.create(reference_npz)

predict_args <- .wisecondorx_predict_args(
  npz = sample_npz,
  ref = reference_npz,
  output_prefix = output_prefix,
  minrefbins = 150L,
  maskrepeats = 5L,
  zscore = 5,
  alpha = 1e-4,
  beta = 0.2,
  blacklist = blacklist,
  gender = "m",
  bed = TRUE,
  plot = TRUE,
  regions = regions,
  ylim = c(-2, 2),
  cairo = TRUE,
  seed = 7L,
  add_plot_title = TRUE,
  extra_args = c("--future-flag", "future-value")
)

expect_identical(
  predict_args,
  c(
    "predict",
    normalizePath(sample_npz, mustWork = TRUE),
    normalizePath(reference_npz, mustWork = TRUE),
    normalizePath(output_prefix, mustWork = FALSE),
    "--minrefbins", "150",
    "--maskrepeats", "5",
    "--zscore", "5",
    "--alpha", "1e-04",
    "--beta", "0.2",
    "--blacklist", normalizePath(blacklist, mustWork = TRUE),
    "--gender", "M",
    "--bed",
    "--plot",
    "--add-plot-title",
    "--regions", normalizePath(regions, mustWork = TRUE),
    "--ylim", "[-2,2]",
    "--cairo",
    "--seed", "7",
    "--future-flag", "future-value"
  )
)

expect_identical(.format_ylim(c(-2, 2)), "[-2,2]")
expect_identical(.format_ylim(NULL), NULL)

# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------

# Invalid gender → match.arg error
expect_error(
  .wisecondorx_predict_args(
    npz           = sample_npz,
    ref           = reference_npz,
    output_prefix = output_prefix,
    minrefbins    = 150L,
    maskrepeats   = 5L,
    zscore        = 5,
    alpha         = 1e-4,
    beta          = NULL,
    blacklist     = NULL,
    gender        = "X",
    bed           = FALSE,
    plot          = FALSE,
    regions       = NULL,
    ylim          = NULL,
    cairo         = FALSE,
    seed          = NULL,
    add_plot_title = FALSE,
    extra_args    = character(0)
  ),
  info = "invalid gender 'X' errors via match.arg"
)

# ylim of wrong length → stopifnot error
expect_error(
  .wisecondorx_predict_args(
    npz           = sample_npz,
    ref           = reference_npz,
    output_prefix = output_prefix,
    minrefbins    = 150L,
    maskrepeats   = 5L,
    zscore        = 5,
    alpha         = 1e-4,
    beta          = NULL,
    blacklist     = NULL,
    gender        = NULL,
    bed           = FALSE,
    plot          = FALSE,
    regions       = NULL,
    ylim          = c(-2, 0, 2),
    cairo         = FALSE,
    seed          = NULL,
    add_plot_title = FALSE,
    extra_args    = character(0)
  ),
  info = "ylim of length 3 errors"
)

# cpus < 1 → stopifnot error
expect_error(
  .wisecondorx_newref_args(
    npz_files   = c(npz1, npz2),
    output      = output_ref,
    ref_binsize = 50000L,
    nipt        = FALSE,
    refsize     = 300L,
    yfrac       = NULL,
    plotyfrac   = NULL,
    cpus        = 0L,
    extra_args  = character(0)
  ),
  info = "cpus = 0 errors"
)

# NULL gender → --gender absent from args vector
predict_args_no_gender <- .wisecondorx_predict_args(
  npz           = sample_npz,
  ref           = reference_npz,
  output_prefix = output_prefix,
  minrefbins    = 150L,
  maskrepeats   = 5L,
  zscore        = 5,
  alpha         = 1e-4,
  beta          = NULL,
  blacklist     = NULL,
  gender        = NULL,
  bed           = FALSE,
  plot          = FALSE,
  regions       = NULL,
  ylim          = NULL,
  cairo         = FALSE,
  seed          = NULL,
  add_plot_title = FALSE,
  extra_args    = character(0)
)
expect_false("--gender" %in% predict_args_no_gender,
             info = "NULL gender omits --gender flag")
