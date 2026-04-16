# Score gaunosomes from a typed reference model

Computes the package-level X/Y scoring bundle against a prepared
[`NIPTReferenceModel`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTReferenceModel.md):
sex-matched z-scores, sex-matched NCV scores, and sex-matched regression
scores. The result is returned as one validated object with a compact
per-chromosome summary table.

## Usage

``` r
nipter_gaunosome_score(
  sample,
  reference,
  y_unique_ratio = NULL,
  min_controls = 2L,
  sample_predictors = NULL,
  focus_chromosomes = c("X", "Y")
)
```

## Arguments

- sample:

  A `NIPTeRSample` or typed `NIPTSample`.

- reference:

  A `NIPTReferenceModel` prepared with
  [`nipter_build_gaunosome_models`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_gaunosome_models.md).

- y_unique_ratio:

  Optional numeric scalar passed through to
  [`nipter_predict_sex`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
  when the reference includes a `"y_unique"` sex model.

- min_controls:

  Minimum number of non-outlier same-sex controls required for the
  z-score component. Default `2L`.

- sample_predictors:

  Optional named list of extra predictor values for the sample, used
  when the regression models include extra columns such as
  `GCPCTAfterFiltering`.

- focus_chromosomes:

  Character vector; any subset of `c("X", "Y")`.

## Value

A typed `NIPTGaunosomeScore` object containing the component sex, NCV,
and regression scores plus a compact summary data frame.
