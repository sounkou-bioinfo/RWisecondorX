# Select best-matching controls for a sample

Ranks control samples by similarity of chromosomal fractions to a test
sample. Uses sum-of-squared differences of collapsed chromosomal
fractions (chromosomes 1-12, 14-17, 19-20, 22 by default — excluding
trisomy chromosomes 13, 18, 21). The distance computation is accelerated
by an Rcpp+OpenMP kernel.

## Usage

``` r
nipter_match_control_group(
  sample,
  control_group,
  n,
  mode = c("subset", "report"),
  exclude_chromosomes = c(13L, 18L, 21L),
  include_chromosomes = NULL,
  cpus = 1L
)
```

## Arguments

- sample:

  A `NIPTeRSample` object to match against.

- control_group:

  A `NIPTeRControlGroup` object.

- n:

  Integer; number of best-matching controls to return.

- mode:

  `"subset"` (default) returns a new `NIPTeRControlGroup`; `"report"`
  returns a named numeric vector of sum-of-squares scores.

- exclude_chromosomes:

  Integer vector of chromosomes to exclude from the distance calculation
  (default `c(13, 18, 21)`).

- include_chromosomes:

  Integer vector of chromosomes to include. If `NULL` (default), uses
  all autosomal chromosomes minus `exclude_chromosomes`.

- cpus:

  Integer; number of OpenMP threads for the SSD computation. Default
  `1L`.

## Value

A `NIPTeRControlGroup` (if `mode = "subset"`) or a named numeric vector
of sum-of-squares distances (if `mode = "report"`).

## See also

[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md),
[`nipter_match_matrix()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_matrix.md)
