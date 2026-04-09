# Convert BAM/CRAM to WisecondorX NPZ format

Runs
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
and serialises the resulting bin-count list to a `.npz` file that is
byte-compatible with the file written by `wisecondorx convert`. The NPZ
must be created by numpy so that `wisecondorx newref` can load it; this
function therefore requires `reticulate` and a Python environment with
`numpy` installed (any version ≥ 1.16 works — the format is stable).

## Usage

``` r
bam_convert_npz(
  bam,
  npz,
  binsize = 5000L,
  rmdup = c("streaming", "none", "flag"),
  con = NULL,
  np = NULL
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- npz:

  Path for the output `.npz` file (created or overwritten).

- binsize:

  Bin size in base pairs (default 5000).

- rmdup:

  Duplicate-removal strategy passed to
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

- con:

  Optional open DBI connection with duckhts already loaded. Passed
  through to
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

- np:

  Optional numpy module imported via `reticulate::import("numpy")`. If
  `NULL` (default) it is imported automatically.

## Value

`npz` (invisibly).

## See also

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md),
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bam_convert_npz("sample.bam", "sample.npz", binsize = 5000, rmdup = "streaming")
} # }
```
