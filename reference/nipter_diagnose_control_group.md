# Diagnose a NIPTeR control group

Computes per-chromosome Z-scores across all samples in a control group
and flags outliers (\|Z\| \> 3). The Shapiro-Wilk test is applied to
each chromosome's fraction distribution.

## Usage

``` r
nipter_diagnose_control_group(control_group)
```

## Arguments

- control_group:

  A `NIPTeRControlGroup` object.

## Value

A list with three elements:

- z_scores:

  Chromosome-by-sample matrix of Z-scores (rows = chromosomes 1-22,
  columns = samples).

- aberrant_scores:

  A `data.frame` with columns `chromosome`, `sample_name`, `z_score` for
  all \|Z\| \> 3, or `NULL` if none.

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
