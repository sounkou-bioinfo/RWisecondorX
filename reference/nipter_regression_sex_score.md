# Score sex chromosomes with regression models from a typed reference

Applies the sex-matched X/Y regression models attached to a
`NIPTReferenceModel`. Each per-model score is the standardized ratio of
observed to predicted sex-chromosome fraction, mirroring the production
pipeline's regression-based z-score idea.

## Usage

``` r
nipter_regression_sex_score(
  sample,
  reference,
  focus_chromosome = c("X", "Y"),
  y_unique_ratio = NULL,
  sample_predictors = NULL
)
```

## Arguments

- sample:

  A `NIPTeRSample` or typed `NIPTSample`.

- reference:

  A `NIPTReferenceModel` with `$sex_regression_models` already built by
  [`nipter_build_sex_regression_models`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_sex_regression_models.md).

- focus_chromosome:

  Character scalar, `"X"` or `"Y"`.

- y_unique_ratio:

  Optional scalar passed through to
  [`nipter_predict_sex`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
  if the reference carries a `"y_unique"` sex model.

- sample_predictors:

  Optional named list of extra predictor values for the sample, used
  when the fitted regression models include extra columns such as
  `GCPCTAfterFiltering`.

## Value

A typed `NIPTSexRegressionScore`.
