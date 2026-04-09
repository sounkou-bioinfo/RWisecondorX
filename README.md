
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWisecondorX

<!-- badges: start -->

[![R-CMD-check](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`RWisecondorX` exposes WisecondorX-style copy number analysis workflows
from R on top of `Rduckhts` and DuckDB.

The package is intentionally split into two layers:

- `bam_convert()` runs the convert step in R/SQL and does not require
  Python at runtime.
- `bam_convert_npz()` optionally uses `reticulate` plus `numpy` to write
  WisecondorX-compatible `.npz` files.
- `wisecondorx_newref()` and `wisecondorx_predict()` optionally call the
  official `wisecondorx` CLI through `condathis`.

The convert implementation is designed for exact conformance with
upstream WisecondorX duplicate handling, including the streaming `larp`
/ `larp2` behaviour used during bin counting.

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("sounkou-bioinfo/RWisecondorX")
```

`RWisecondorX` imports `Rduckhts`, `DBI`, and `duckdb`. Optional
features use:

- `reticulate` for writing `.npz` files
- `condathis` for creating a reproducible conformance environment with
  the official bioconda `wisecondorx` package

## What The Package Covers

The current workflow surface is:

1.  Convert BAM/CRAM alignments to WisecondorX-style per-bin counts with
    `bam_convert()`.
2.  Optionally serialise those counts to `.npz` with
    `bam_convert_npz()`.
3.  Optionally build a reference panel with `wisecondorx_newref()`.
4.  Optionally run predictions with `wisecondorx_predict()`.

This keeps the core counting step in R while still allowing conformance
checks and downstream use of the official Python CLI when needed.

## Convert A BAM To Bin Counts

The package ships with a small chromosome 11 BAM fixture for examples
and tests.

``` r
library(RWisecondorX)

fixture_bam <- system.file(
  "extdata",
  "hg00106_chr11_fixture.bam",
  package = "RWisecondorX"
)

bins <- bam_convert(
  fixture_bam,
  binsize = 5000L,
  rmdup = "streaming"
)

names(Filter(Negate(is.null), bins))
#> [1] "11"

length(bins[["11"]])
#> [1] 1198

sum(bins[["11"]])
#> [1] 982
```

The `rmdup` argument supports the three intended duplicate-handling
modes:

- `"streaming"`: exact replication of upstream WisecondorX `larp` /
  `larp2`
- `"none"`: no duplicate removal, corresponding to
  `wisecondorx convert --normdup`
- `"flag"`: use the SAM duplicate flag (`0x400`) for pre-marked BAMs

## Optional NPZ And CLI Workflow

When you want to continue into the original WisecondorX reference and
prediction steps, the package provides an optional NPZ bridge via
`reticulate` and CLI wrappers via `condathis`.

``` r
library(RWisecondorX)

fixture_npz <- tempfile(fileext = ".npz")
Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")

if (!requireNamespace("reticulate", quietly = TRUE)) {
  "reticulate not available; skipping NPZ example"
} else {
  np <- tryCatch(
    reticulate::import("numpy", convert = FALSE),
    error = function(e) NULL
  )

  if (is.null(np)) {
    "numpy not available in the active Python environment; skipping NPZ example"
  } else {
    bam_convert_npz(
      bam = fixture_bam,
      npz = fixture_npz,
      binsize = 5000L,
      rmdup = "streaming",
      np = np
    )

    c(
      npz_created = file.exists(fixture_npz),
      npz_size = file.info(fixture_npz)$size
    )
  }
}
#> npz_created    npz_size 
#>           1         392
```

The CLI wrappers expose the upstream arguments documented in the
WisecondorX README. The chunk below is evaluated, but it only runs the
external CLI calls when `RWISECONDORX_EVAL_CLI_EXAMPLES=1` is set in the
environment.

``` r
library(RWisecondorX)

cli_enabled <- identical(Sys.getenv("RWISECONDORX_EVAL_CLI_EXAMPLES"), "1")

if (!cli_enabled) {
  "Set RWISECONDORX_EVAL_CLI_EXAMPLES=1 to run the condathis-backed CLI example"
} else if (!file.exists(fixture_npz)) {
  "NPZ example did not create an input file; skipping CLI example"
} else if (!requireNamespace("condathis", quietly = TRUE)) {
  "condathis not available; skipping CLI example"
} else {
  reference_npz <- tempfile(fileext = ".npz")
  plotyfrac_png <- tempfile(fileext = ".png")
  output_prefix <- tempfile(pattern = "wisecondorx-")

  newref_status <- tryCatch(
    {
      wisecondorx_newref(
        npz_files = c(fixture_npz, fixture_npz),
        output = reference_npz,
        binsize = 5000L,
        ref_binsize = 50000L,
        nipt = TRUE,
        refsize = 300L,
        yfrac = NULL,
        plotyfrac = plotyfrac_png,
        cpus = 1L
      )
      "ok"
    },
    error = function(e) conditionMessage(e)
  )

  predict_status <- if (identical(newref_status, "ok") && file.exists(reference_npz)) {
    tryCatch(
      {
        wisecondorx_predict(
          npz = fixture_npz,
          ref = reference_npz,
          output_prefix = output_prefix,
          ref_binsize = 50000L,
          minrefbins = 150L,
          maskrepeats = 5L,
          zscore = 5,
          alpha = 1e-4,
          beta = NULL,
          blacklist = NULL,
          gender = NULL,
          bed = TRUE,
          plot = FALSE,
          regions = NULL,
          ylim = NULL,
          cairo = FALSE,
          seed = 1L
        )
        "ok"
      },
      error = function(e) conditionMessage(e)
    )
  } else {
    "predict skipped because newref did not produce a reference"
  }

  c(
    newref_status = newref_status,
    reference_exists = file.exists(reference_npz),
    predict_status = predict_status
  )
}
#> [1] "Set RWISECONDORX_EVAL_CLI_EXAMPLES=1 to run the condathis-backed CLI example"
```

## SRA Metadata Helpers

The package also includes small helpers for staging SRA run metadata
used by conformance fixtures.

``` r
library(RWisecondorX)

sra_runinfo_url("PRJNA400134")
#> [1] "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA400134"
```

## Development

`README.Rmd` is the editable source for this document. Regenerate
`README.md` with the repository `Makefile`:

``` sh
make readme
make fixtures
make rd
```
