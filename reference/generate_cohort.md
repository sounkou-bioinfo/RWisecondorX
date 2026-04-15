# Generate a synthetic BAM cohort for testing

Creates coordinate-sorted, indexed BAM files in a specified directory
using compressed chromosome lengths: each 100kb GRCh37 bin is
represented as a 100bp region, yielding identical bin-count structure at
a fraction of the file size (~435KB per BAM).

## Usage

``` r
generate_cohort(out_dir, n_samples = 50L, verbose = TRUE)
```

## Arguments

- out_dir:

  Directory to write BAM files into (created if needed).

- n_samples:

  Total number of samples to generate (default 50L). Must be at least 6
  (3 trisomy + 1 male + 2 female minimum).

- verbose:

  Logical; emit progress messages via
  [`message()`](https://rdrr.io/r/base/message.html).

## Value

A data frame (manifest) with columns: `sample_id`, `sex`, `trisomy`,
`n_reads`, `bam_file`.

## Details

The cohort contains approximately 75% euploid females, 25% euploid
males, and 3 trisomy females (T21, T18, T13). Each sample uses
approximately 3 reads per bin (~91k reads total). Results are
deterministic (sample *i* uses `set.seed(42 + i)`).

The generated BAMs are NIPTeR-compatible (no unmapped reads, unique
positions per chromosome) and can be used with both the NIPTeR binning
layer
([`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md))
and the WisecondorX native pipeline
([`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md))
when binned with `binsize = COMPRESSED_BINSIZE` (100).

The manifest is also written to `manifest.tsv` in `out_dir`.

## Examples

``` r
if (FALSE) { # \dontrun{
manifest <- generate_cohort(tempdir())
head(manifest)
} # }
```
