# Convert BAM/CRAM to WisecondorX NPZ format (upstream CLI wrapper)

Calls `wisecondorx convert` via
[`condathis::run()`](https://luciorq.github.io/condathis/reference/run.html)
to convert an aligned BAM or CRAM file to a `.npz` file for downstream
use with
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md)
or
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md).

## Usage

``` r
wisecondorx_convert(
  bam,
  npz,
  reference = NULL,
  binsize = 5000L,
  normdup = FALSE,
  env_name = "wisecondorx",
  extra_args = character(0)
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- npz:

  Path for the output `.npz` file (created or overwritten).

- reference:

  Optional path to a FASTA reference file. Required when `bam` is a CRAM
  file. Passed to upstream `--reference`.

- binsize:

  Bin size in base pairs (default 5000). The reference bin size should
  be a multiple of this value. Passed to upstream `--binsize`.

- normdup:

  Logical; when `TRUE`, passes `--normdup` to skip duplicate removal.
  Recommended for NIPT data where read depth is low.

- env_name:

  Name of the conda environment containing `wisecondorx` (default
  `"wisecondorx"`). Created automatically by condathis on first use via
  the `bioconda` channel.

- extra_args:

  Character vector of additional arguments passed verbatim after the
  mapped CLI flags. For forward compatibility with future upstream
  WisecondorX releases.

## Value

`npz` (invisibly).

## Details

This wrapper exposes the upstream `wisecondorx convert` CLI flags:
`--reference`, `--binsize`, and `--normdup`.

For a fully native R implementation (no Python dependency) that uses
Rduckhts instead of pysam, see
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
and
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md).

## See also

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md),
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wisecondorx_convert(
  bam = "sample.bam",
  npz = "sample.npz",
  binsize = 5000L,
  normdup = FALSE
)
} # }
```
