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
  ref_binsize = 50000L,
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

  Reference bin size in base pairs (default 50000). Must be a multiple
  of `binsize`.

- cpus:

  Number of CPUs to pass to `wisecondorx newref` (default 1).

- env_name:

  Name of the conda environment containing `wisecondorx` (default
  `"wisecondorx"`). Created automatically by condathis on first use via
  the `bioconda` channel.

- extra_args:

  Character vector of additional arguments passed verbatim to
  `wisecondorx newref` (e.g. `c("--yfrac", "0.4")`).

## Value

`output` (invisibly).

## See also

[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Build a reference from 30 NIPT controls
wisecondorx_newref(
  npz_files = list.files("controls/", "\\.npz$", full.names = TRUE),
  output    = "reference.npz",
  ref_binsize = 50000L
)
} # }
```
