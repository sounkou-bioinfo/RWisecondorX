# Simulate a trisomy cohort from donor BAMs/CRAMs

Generates a set of positive BAMs from one or more negative donor
alignments by applying
[`simulate_trisomy_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/simulate_trisomy_bam.md)
across a simulation grid. The result is a manifest that ties every
simulated BAM back to its donor, chromosome, fetal fraction, and mosaic
setting.

## Usage

``` r
simulate_trisomy_cohort(
  bams = NULL,
  bam_dir = NULL,
  out_dir,
  donor_ids = NULL,
  trisomy_chromosomes = c("21", "18", "13"),
  fetal_fraction = 0.1,
  mosaic = FALSE,
  seed = 1L,
  threads = 1L,
  reference = NULL,
  index = TRUE,
  include_donors = TRUE,
  pattern = "\\.(bam|cram)$",
  overwrite = FALSE
)
```

## Arguments

- bams:

  Optional character vector of donor BAM/CRAM paths.

- bam_dir:

  Optional directory containing donor BAM/CRAM files. Supply exactly one
  of `bams` or `bam_dir`.

- out_dir:

  Output directory for simulated positive BAMs and the manifest.

- donor_ids:

  Optional donor IDs, one per input BAM/CRAM. Defaults to filename stems
  with duplicates made unique.

- trisomy_chromosomes:

  Character vector of target chromosomes. Default `c("21", "18", "13")`.

- fetal_fraction:

  Numeric vector of fetal fractions. A separate simulation is produced
  for each distinct value.

- mosaic:

  Logical vector of mosaic settings. A separate simulation is produced
  for each distinct value. Default `FALSE`.

- seed:

  Integer base seed. Individual simulations use deterministic seeds
  derived from this base.

- threads:

  Integer thread count for htslib I/O and indexing.

- reference:

  Optional FASTA path applied to all CRAM inputs.

- index:

  Logical; create `.bai` indexes for simulated BAMs. Default `TRUE`.

- include_donors:

  Logical; include the donor negatives as rows in the manifest. Default
  `TRUE`.

- pattern:

  File-matching regex used only when `bam_dir` is supplied.

- overwrite:

  Logical; overwrite existing simulated BAMs and manifest.

## Value

A data frame manifest. The same table is written to
`<out_dir>/manifest.tsv`.

## Details

The donor BAMs themselves serve as the negative cohort; by default they
are recorded in the manifest but not copied into `out_dir`.

## See also

[`simulate_trisomy_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/simulate_trisomy_bam.md),
[`generate_cohort()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/generate_cohort.md)

## Examples

``` r
if (FALSE) { # \dontrun{
manifest <- simulate_trisomy_cohort(
  bam_dir = "negatives/",
  out_dir = "simulated_positives/",
  trisomy_chromosomes = c("21", "18"),
  fetal_fraction = c(0.08, 0.12)
)
} # }
```
