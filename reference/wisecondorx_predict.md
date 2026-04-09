# Predict copy-number aberrations with WisecondorX

Calls `wisecondorx predict` via
[`condathis::run()`](https://luciorq.github.io/condathis/reference/run.html)
on a single-sample NPZ file against a reference panel built by
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md).

## Usage

``` r
wisecondorx_predict(
  npz,
  ref,
  output_prefix,
  ref_binsize = 50000L,
  zscore = 5,
  bed = FALSE,
  gender_file = NULL,
  env_name = "wisecondorx",
  extra_args = character(0)
)
```

## Arguments

- npz:

  Path to the sample `.npz` file (from
  [`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md)).

- ref:

  Path to the reference `.npz` file (from
  [`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md)).

- output_prefix:

  Output prefix for all `wisecondorx predict` output files (e.g.
  `"results/sample1"`).

- ref_binsize:

  Reference bin size in base pairs (default 50000). Must match the
  `ref_binsize` used when building the reference.

- zscore:

  Z-score threshold for aberration calling (default 5).

- bed:

  Logical; also write BED output files (default `FALSE`).

- gender_file:

  Optional path to a gender model file produced by `wisecondorx gender`
  (passed as `--gender`).

- env_name:

  Name of the conda environment containing `wisecondorx` (default
  `"wisecondorx"`).

- extra_args:

  Character vector of additional arguments passed verbatim to
  `wisecondorx predict`.

## Value

`output_prefix` (invisibly).

## See also

[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wisecondorx_predict(
  npz    = "sample.npz",
  ref    = "reference.npz",
  output_prefix = "results/sample"
)
} # }
```
