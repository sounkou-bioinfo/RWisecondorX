# Diagnose a NIPTeR control group

Computes per-chromosome Z-scores across all samples in a control group
and flags outliers (\|Z\| \> 3). The Shapiro-Wilk test is applied to
each chromosome's fraction distribution.

## Usage

``` r
nipter_diagnose_control_group(
  control_group,
  collapse_strands = TRUE,
  z_cutoff = 3
)
```

## Arguments

- control_group:

  A `NIPTeRControlGroup` object.

- collapse_strands:

  Logical; for `SeparatedStrands` control groups, collapse forward and
  reverse strand fractions into 22 autosomal fractions (`TRUE`, default)
  or keep the upstream-style 44-row strand-resolved diagnostics
  (`FALSE`). Ignored for `CombinedStrands`.

- z_cutoff:

  Absolute Z-score cutoff used to flag aberrant rows. Default `3`.

## Value

A list with three elements:

- z_scores:

  Chromosome-by-sample matrix of Z-scores (rows = chromosomes 1-22, or
  strand-resolved rows 1F-22F/1R-22R when `collapse_strands = FALSE` on
  a `SeparatedStrands` control group; columns = samples).

- aberrant_scores:

  A `data.frame` with columns `chromosome`, `sample_name`, `z_score` for
  all `|Z| > z_cutoff`, or `NULL` if none.

- statistics:

  A matrix with rows per chromosome and columns `mean`, `SD`,
  `shapiro_p_value`.

## See also

[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md)

## Examples

``` r
if (FALSE) { # \dontrun{
diag <- nipter_diagnose_control_group(cg)
diag$aberrant_scores
} # }
```
