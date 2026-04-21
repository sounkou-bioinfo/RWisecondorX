# Convert BAM/CRAM to WisecondorX NPZ format

Runs
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
and serialises the resulting bin-count list to a `.npz` file that is
byte-compatible with the file written by `wisecondorx convert`. The NPZ
must be created by numpy so that `wisecondorx newref` can load it; this
function therefore requires `reticulate` and a Python environment with
`numpy` installed; any numpy version from 1.16 onward works.

## Usage

``` r
bam_convert_npz(
  bam,
  npz,
  binsize = 5000L,
  mapq = 1L,
  require_flags = 0L,
  exclude_flags = 0L,
  rmdup = c("streaming", "none", "flag"),
  con = NULL,
  np = NULL,
  reference = NULL
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- npz:

  Path for the output `.npz` file (created or overwritten).

- binsize:

  Bin size in base pairs (default 5000).

- mapq:

  Minimum mapping quality passed to
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

- require_flags:

  Integer bitmask of required SAM flags passed to
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

- exclude_flags:

  Integer bitmask of excluded SAM flags passed to
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

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

- reference:

  Optional FASTA reference path for CRAM inputs. Passed to
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

## Value

`npz` (invisibly).

## Details

This function exists for Python WisecondorX CLI conformance and
interoperability. The native R pipeline uses
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md)
together with
[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md)
and
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md).

The resulting NPZ has three top-level keys:

- `sample`:

  0-d object array wrapping a dict mapping `"1"`..`"24"` to int32 arrays
  (chromosome bin counts).

- `binsize`:

  Scalar int (bin size in bp).

- `quality`:

  0-d object array wrapping a dict of QC counters (populated with zeros
  since the native binning kernel does not track per-read filter stats).

When the writer's numpy is \>= 2.0 the internal pickle references
`numpy._core`, which older numpy (\< 2.0) cannot unpickle. This function
patches the pickle bytestream to use `numpy.core` so the NPZ is readable
by both numpy 1.x and 2.x.

Native
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
returns dense vectors covering the chromosome span in the BAM header.
Upstream `wisecondorx convert` has one additional historical quirk: it
allocates `int(length / binsize + 1)` bins for every chromosome in the
header, so chromosomes whose lengths are exact multiples of `binsize`
carry one extra trailing all-zero bin. This function pads to that
upstream NPZ layout during serialisation so Python
`wisecondorx newref/predict` remains byte-compatible, while the native R
list/BED paths keep the cleaner header-span contract.

## See also

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md),
[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md),
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bam_convert_npz("sample.bam", "sample.npz", binsize = 5000, rmdup = "streaming")
} # }
```
