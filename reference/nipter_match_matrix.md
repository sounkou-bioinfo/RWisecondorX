# Compute the full pairwise SSD matrix for a control group

Returns the symmetric N×N matrix of sum-of-squared-differences between
all pairs of control samples' chromosomal fractions. The diagonal is
zero. Row means of this matrix are the per-sample "matching score" used
for QC in the production matching loop (see *Details*).

## Usage

``` r
nipter_match_matrix(
  control_group,
  exclude_chromosomes = c(13L, 18L, 21L),
  include_chromosomes = NULL,
  cpus = 1L
)
```

## Arguments

- control_group:

  A `NIPTeRControlGroup` object.

- exclude_chromosomes:

  Integer vector of chromosomes to exclude from the distance calculation
  (default `c(13, 18, 21)`).

- include_chromosomes:

  Integer vector of chromosomes to include. If `NULL` (default), uses
  all autosomal chromosomes minus `exclude_chromosomes`.

- cpus:

  Integer; OpenMP threads. Default `1L`.

## Value

A numeric N×N matrix with sample names as row and column names.

## Details

The production NIPT pipeline uses this matrix to identify outlier
controls before scoring: each sample's mean SSD against all others is
computed, and samples with mean SSD more than 3 SD above the group mean
are iteratively removed. This function replaces the `lapply` over
`match_control_group` calls in `CoverageProjectionSCA_Reports.R` with a
single vectorized Rcpp kernel.

## See also

[`nipter_match_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_control_group.md)
