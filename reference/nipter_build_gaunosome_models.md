# Build both gaunosome model families on a typed reference

Convenience wrapper that attaches both sex-chromosome NCV models and
sex-chromosome regression models to a
[`NIPTReferenceModel`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTReferenceModel.md).
This is the package-owned precomputation step for downstream X/Y
scoring.

## Usage

``` r
nipter_build_gaunosome_models(
  reference,
  candidate_chromosomes = c(1:12, 14:16, 20, 22),
  ncv_min_elements = 6L,
  ncv_max_elements = 9L,
  regression_n_models = 4L,
  regression_n_predictors = 4L,
  regression_extra_predictors = "GCPCTAfterFiltering",
  focus_chromosomes = c("X", "Y")
)
```

## Arguments

- reference:

  A `NIPTReferenceModel`.

- candidate_chromosomes:

  Integer vector of autosomes allowed in the denominator/predictor pool.
  Defaults to the production-style set `c(1:12, 14:16, 20, 22)`.

- ncv_min_elements:

  Minimum denominator-set size for the NCV search.

- ncv_max_elements:

  Maximum denominator-set size for the NCV search.

- regression_n_models:

  Number of regression models to build per sex/focus combination.

- regression_n_predictors:

  Maximum predictors per regression model.

- regression_extra_predictors:

  Optional character vector of additional numeric columns to use when
  present in `reference$reference_frame`.

- focus_chromosomes:

  Character vector; any subset of `c("X", "Y")`.

## Value

The input `NIPTReferenceModel`, enriched with both `$sex_ncv_models` and
`$sex_regression_models`.
