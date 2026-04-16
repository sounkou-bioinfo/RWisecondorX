# Simulate a trisomy BAM by thinning a donor alignment

Creates a BAM by copying all reads on the target chromosome unchanged
and randomly removing reads from all other mapped chromosomes according
to the Nguyen et al. fetal-fraction model. This preserves the donor
alignment's read names, flags, CIGAR strings, mapping qualities, and
coordinate order while enriching the target chromosome relative to the
rest of the genome.

## Usage

``` r
simulate_trisomy_bam(
  bam,
  out_bam,
  trisomy_chr,
  fetal_fraction = 0.1,
  mosaic = FALSE,
  seed = 1L,
  threads = 1L,
  reference = NULL,
  index = TRUE,
  overwrite = FALSE
)
```

## Arguments

- bam:

  Input BAM or CRAM path.

- out_bam:

  Output BAM path.

- trisomy_chr:

  Target chromosome to enrich, e.g. `"21"` or `"chr21"`.

- fetal_fraction:

  Numeric fetal fraction in `(0, 1]`. Default `0.10`.

- mosaic:

  Logical; if `TRUE`, use the 50% mosaic enrichment model (`k = 0.25`)
  instead of the non-mosaic model (`k = 0.5`).

- seed:

  Integer RNG seed for deterministic thinning.

- threads:

  Integer thread count for htslib I/O and BAM indexing.

- reference:

  Optional FASTA path required when the input is CRAM.

- index:

  Logical; if `TRUE` (default), create a `.bai` index for `out_bam`.

- overwrite:

  Logical; overwrite `out_bam` and its `.bai` if present.

## Value

A named list with input/output record counts and thinning metadata.

## Details

This is a BAM-level simulation primitive intended for realistic fixture
and benchmark generation from real negative donor alignments. It is more
faithful than the compressed synthetic cohort generator in
[`generate_cohort()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/generate_cohort.md),
but it is still a simulation and does not replace real-data validation.

## See also

[`generate_cohort()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/generate_cohort.md),
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)

## Examples

``` r
if (FALSE) { # \dontrun{
simulate_trisomy_bam(
  bam = "negative_sample.bam",
  out_bam = "simulated_t21.bam",
  trisomy_chr = "21",
  fetal_fraction = 0.12,
  seed = 42L
)
} # }
```
