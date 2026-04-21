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
  minrefbins = 150L,
  maskrepeats = 5L,
  zscore = 5,
  alpha = 1e-04,
  beta = NULL,
  blacklist = NULL,
  gender = NULL,
  bed = FALSE,
  plot = FALSE,
  regions = NULL,
  ylim = NULL,
  cairo = FALSE,
  seed = NULL,
  add_plot_title = FALSE,
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

- minrefbins:

  Minimum number of sensible reference bins per target bin. Passed to
  upstream `--minrefbins`.

- maskrepeats:

  Number of repeat-masking cycles. Passed to upstream `--maskrepeats`.

- zscore:

  Z-score threshold for aberration calling. Passed to upstream
  `--zscore`.

- alpha:

  P-value cutoff for circular binary segmentation breakpoints. Passed to
  upstream `--alpha`.

- beta:

  Optional ratio cutoff for aberration calling. Passed to upstream
  `--beta`. When set, upstream ignores `zscore`.

- blacklist:

  Optional path to a headerless BED blacklist file. Passed to upstream
  `--blacklist`.

- gender:

  Optional forced gender, `"F"` or `"M"`. Passed to upstream `--gender`.

- bed:

  Logical; also write BED output files. Passed to upstream `--bed`.

- plot:

  Logical; also write plot output files. Passed to upstream `--plot`.

- regions:

  Optional path to a headerless BED file with regions to mark on the
  plot. Passed to upstream `--regions`.

- ylim:

  Optional numeric vector of length 2 giving the y-axis interval, e.g.
  `c(-2, 2)`. Passed to upstream `--ylim` as `[a,b]`.

- cairo:

  Logical; use Cairo bitmap output. Passed to upstream `--cairo`.

- seed:

  Optional integer random seed for segmentation. Passed to upstream
  `--seed`.

- add_plot_title:

  Logical; when `TRUE`, adds the output basename as the plot title.
  Passed to upstream `--add-plot-title`. Only effective when
  `plot = TRUE`.

- env_name:

  Name of the conda environment containing `wisecondorx` (default
  `"wisecondorx"`).

- extra_args:

  Character vector of additional arguments passed verbatim after the
  mapped CLI flags. Keep this for forward compatibility with future
  upstream WisecondorX releases.

## Value

`output_prefix` (invisibly).

## Details

This wrapper exposes the current upstream `wisecondorx predict` CLI
flags documented in the WisecondorX README: `--minrefbins`,
`--maskrepeats`, `--zscore`, `--alpha`, `--beta`, `--blacklist`,
`--gender`, `--bed`, `--plot`, `--add-plot-title`, `--regions`,
`--ylim`, `--cairo`, and `--seed`.

Conda execution is run inside an isolated HOME/XDG sandbox so stale
libmamba lock files in the user's global cache do not leak into package
calls. Unlike
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md),
this wrapper does not expose the native `parallel` argument used for
`ParDNAcopy`; segmentation behavior is delegated entirely to the
upstream Python WisecondorX CLI.

## See also

[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wisecondorx_predict(
  npz = "sample.npz",
  ref = "reference.npz",
  output_prefix = "results/sample",
  minrefbins = 150L,
  maskrepeats = 5L,
  zscore = 5,
  alpha = 1e-4,
  beta = NULL,
  blacklist = NULL,
  gender = "F",
  bed = TRUE,
  plot = TRUE,
  add_plot_title = TRUE,
  regions = NULL,
  ylim = c(-2, 2),
  cairo = FALSE,
  seed = 1L
)
} # }
```
