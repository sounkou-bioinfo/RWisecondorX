# Normalised Chromosome Value (NCV) score

Computes the NCV score for a test sample. The NCV method selects an
optimal set of denominator chromosomes that minimises the coefficient of
variation (CV) in the control group, then uses the resulting ratio to
Z-score the test sample.

## Usage

``` r
nipter_ncv_score(
  sample,
  control_group,
  chromo_focus,
  max_elements = 5L,
  exclude_chromosomes = c(13L, 18L, 21L),
  include_chromosomes = NULL
)
```

## Arguments

- sample:

  A `NIPTeRSample` object.

- control_group:

  A `NIPTeRControlGroup` object.

- chromo_focus:

  Integer; the target chromosome (1-22).

- max_elements:

  Maximum number of denominator chromosomes to try (default 5). The
  algorithm searches all combinations of 1 to `max_elements` from the
  candidate pool.

- exclude_chromosomes:

  Integer vector of chromosomes to exclude from the denominator pool
  (default `c(13, 18, 21)` — the trisomy chromosomes).

- include_chromosomes:

  Integer vector of chromosomes to force-include in the candidate pool.
  Default `NULL`.

## Value

A list of class `"NIPTeRNCV"` with elements:

- sample_score:

  The NCV score for the test sample (numeric).

- focus_chromosome:

  Character; the tested chromosome.

- denominators:

  Integer vector of selected denominator chromosomes.

- control_statistics:

  Named numeric vector with `mean`, `sd`, and `shapiro_p_value`.

- control_z_scores:

  Named numeric vector of control NCV scores.

- best_cv:

  The minimum CV achieved.

- correction_status:

  Character vector.

- sample_name:

  Character.

## Details

The candidate denominator pool is: chromosomes 1-12, 14-17, 19-20, 22
(the "control chromosomes") minus `exclude_chromosomes` minus
`chromo_focus`, plus any `include_chromosomes`.

For each combination size from 1 to `max_elements`, all
`choose(n_candidates, size)` combinations are evaluated. The ratio
`reads[focus] / sum(reads[denominators])` is computed per control
sample, and the combination minimising the coefficient of variation
(`sd/mean`) is selected.

## See also

[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md),
[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ncv21 <- nipter_ncv_score(sample, cg, chromo_focus = 21, max_elements = 5)
ncv21$sample_score
ncv21$denominators
} # }
```
