# Read a WisecondorX-format BED file into a sample list

Reads a bgzipped BED file and returns a named list of integer vectors
suitable for
[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md),
[`scale_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/scale_sample.md),
or
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md).

## Usage

``` r
bed_to_sample(bed, binsize = NULL, con = NULL)
```

## Arguments

- bed:

  Path to a bgzipped (or plain) BED file with a `.tbi` index.

- binsize:

  Optional integer; bin size in base pairs. If `NULL` (default),
  inferred from the first row of the BED file.

- con:

  Optional open DBI connection with duckhts already loaded.

## Value

A named list with one integer vector per chromosome key (`"1"`–`"22"`,
`"23"` for X, `"24"` for Y). Each vector has length `max_bin + 1` for
that chromosome. Chromosomes absent from the BED file are `NULL`. This
is the same format returned by
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

## Details

Only the first 4 tab-delimited columns are used: `chrom`, `start`,
`end`, `count` (no column header row). Optional leading
`##RWX_<key>=<value>` metadata lines are ignored automatically. Any
trailing columns are ignored, so this reader stays strand-agnostic and
can also consume NIPTeR BEDs by taking their total-count column.
Coordinates are 0-based half-open intervals. The bin size is inferred
from the first row (`end - start`) unless explicitly provided.

## See also

[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md),
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md),
[`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Write bin counts to BED, then read them back
bam_convert_bed("sample.bam", "sample.bed.gz", binsize = 5000L)
bins <- bed_to_sample("sample.bed.gz")

# Use directly with the native WisecondorX pipeline
samples <- lapply(bed_files, bed_to_sample)
ref <- rwisecondorx_newref(samples, binsize = 100000L, nipt = TRUE)
} # }
```
