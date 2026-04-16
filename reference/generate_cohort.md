# Generate a synthetic BAM cohort for testing

Creates coordinate-sorted, indexed BAM files in a specified directory
using compressed chromosome lengths: each 100kb GRCh37 bin is
represented as a 100bp region, yielding identical bin-count structure at
a fraction of the file size (~435KB per BAM).

## Usage

``` r
generate_cohort(
  out_dir,
  n_samples = 50L,
  fetal_fraction = 0.1,
  reads_per_bin = READS_PER_BIN_DEFAULT,
  gc_bias_strength = 0.3,
  mosaic = FALSE,
  verbose = TRUE
)
```

## Arguments

- out_dir:

  Directory to write BAM files into (created if needed).

- n_samples:

  Total number of samples to generate (default 50L). Must be at least 6
  (3 trisomy + 1 male + 2 female minimum).

- fetal_fraction:

  Fetal DNA fraction for trisomy samples (default 0.10). Higher values
  produce stronger trisomy signal. Typical NIPT range: 0.05-0.20.

- reads_per_bin:

  Average reads per bin for euploid autosomes (default 3). Higher values
  produce larger BAMs but better signal-to-noise for trisomy detection.
  At `reads_per_bin = 3`, BAMs are ~435KB; at 10, ~1.3MB; at 30, ~4.4MB.

- gc_bias_strength:

  Numeric scaling factor for simulated GC bias amplitude (default 0.3).
  Each sample gets a random quadratic GC bias curve; higher values
  produce more between-sample variance. Set to 0 for uniform coverage
  (no GC effect).

- mosaic:

  Logical; if `TRUE`, simulate 50% placental mosaicism (halved trisomy
  signal strength). Default `FALSE`.

- verbose:

  Logical; emit progress messages via
  [`message()`](https://rdrr.io/r/base/message.html).

## Value

A data frame (manifest) with columns: `sample_id`, `sex`, `trisomy`,
`fetal_fraction`, `n_reads`, `bam_file`.

## Details

The cohort contains approximately 75% euploid females, 25% euploid
males, and 3 trisomy females (T21, T18, T13). Each sample uses
approximately 3 reads per bin (~91k reads total). Results are
deterministic (sample *i* uses `set.seed(42 + i)`).

Trisomy simulation follows Nguyen et al. 2023
(doi:10.1101/2023.11.24.568620): instead of naively multiplying reads on
the trisomy chromosome, reads are stochastically removed from non-target
chromosomes so that the relative coverage of the trisomy chromosome
reflects an extra copy in the fetal genome. The removal probability is
`p = k*f / (1 + k*f)` where `k = 0.5` (non-mosaic) or `k = 0.25`
(mosaic).

The generated BAMs are NIPTeR-compatible (no unmapped reads, unique
positions per chromosome) and can be used with both the NIPTeR binning
layer
([`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md))
and the WisecondorX native pipeline
([`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md))
when binned with `binsize = COMPRESSED_BINSIZE` (100).

These BAMs are structurally valid and fast to generate, but the
underlying coverage model is still synthetic. Use real negative cohorts
for calibration, threshold selection, and any real conformance work.

The manifest is also written to `manifest.tsv` in `out_dir`.

## Examples

``` r
if (FALSE) { # \dontrun{
manifest <- generate_cohort(tempdir())
head(manifest)
} # }
```
