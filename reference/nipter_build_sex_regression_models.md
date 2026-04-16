# Build sex-chromosome regression models on a typed reference

Fits multiple sex-matched X/Y ratio models on the non-outlier subsets of
a
[`NIPTReferenceModel`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTReferenceModel.md).
The returned reference gains a typed `$sex_regression_models` field.

## Usage

``` r
nipter_build_sex_regression_models(
  reference,
  candidate_chromosomes = c(1:12, 14:16, 20, 22),
  n_models = 4L,
  n_predictors = 4L,
  extra_predictors = "GCPCTAfterFiltering",
  focus_chromosomes = c("X", "Y")
)
```

## Arguments

- reference:

  A `NIPTReferenceModel`.

- candidate_chromosomes:

  Integer vector of autosomes allowed as ratio predictors. Defaults to
  the production-style set `c(1:12, 14:16, 20, 22)`.

- n_models:

  Number of models to build per sex/focus combination.

- n_predictors:

  Maximum predictors per model.

- extra_predictors:

  Optional character vector of additional numeric columns to use when
  present in `reference$reference_frame`. The default
  `"GCPCTAfterFiltering"` keeps room for explicit QC metadata without
  hard-wiring that requirement into the core sample classes.

- focus_chromosomes:

  Character vector; any subset of `c("X", "Y")`.

## Value

The input `NIPTReferenceModel`, with a typed `$sex_regression_models`
field containing nested `female`/`male` and `X`/`Y` model lists.
