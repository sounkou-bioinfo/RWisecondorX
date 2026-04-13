# Build a WisecondorX reference panel

Calls `wisecondorx newref` via
[`condathis::run()`](https://luciorq.github.io/condathis/reference/run.html)
on a set of NPZ files produced by
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md).
The NPZ files must all use the same `binsize`; the reference is built at
`ref_binsize` (must be a multiple of `binsize`).

## Usage

``` r
wisecondorx_newref(
  npz_files,
  output,
  binsize = 5000L,
  ref_binsize = 100000L,
  nipt = FALSE,
  refsize = 300L,
  yfrac = NULL,
  plotyfrac = NULL,
  cpus = 1L,
  env_name = "wisecondorx",
  extra_args = character(0)
)
```

## Arguments

- npz_files:

  Character vector of paths to `.npz` files (one per sample).

- output:

  Path for the output reference `.npz` file.

- binsize:

  Convert-step bin size in base pairs (default 5000). Must match the
  `binsize` used when calling
  [`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md).

- ref_binsize:

  Reference bin size in base pairs (default 100000). Passed to upstream
  `--binsize`. Must be a multiple of `binsize`.

- nipt:

  Logical; pass upstream `--nipt`.

- refsize:

  Number of reference locations per target bin. Passed to upstream
  `--refsize`.

- yfrac:

  Optional numeric Y-read fraction cutoff. Passed to upstream `--yfrac`.

- plotyfrac:

  Optional output path for the Y-fraction histogram and mixture-model
  plot. Passed to upstream `--plotyfrac`.

- cpus:

  Number of CPUs to pass to `wisecondorx newref` (default 1).

- env_name:

  Name of the conda environment containing `wisecondorx` (default
  `"wisecondorx"`). Created automatically by condathis on first use via
  the `bioconda` channel.

- extra_args:

  Character vector of additional arguments passed verbatim after the
  mapped CLI flags. Keep this for forward compatibility with future
  upstream WisecondorX releases.

## Value

`output` (invisibly).

## Details

This wrapper exposes the current upstream `wisecondorx newref` CLI flags
documented in the WisecondorX README: `--nipt`, `--binsize`,
`--refsize`, `--yfrac`, `--plotyfrac`, and `--cpus`.

## See also

[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wisecondorx_newref(
  npz_files = list.files("controls/", "\\.npz$", full.names = TRUE),
  output = "reference.npz",
  binsize = 5000L,
  ref_binsize = 100000L,
  nipt = TRUE,
  refsize = 300L,
  yfrac = 0.05,
  cpus = 1L
)
} # }
```
