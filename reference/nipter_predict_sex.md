# Predict fetal sex for a NIPTeR sample

Uses one or more `NIPTeRSexModel` objects to classify a sample as male
or female. When multiple models are supplied, a majority vote determines
the consensus call.

## Usage

``` r
nipter_predict_sex(sample, ..., y_unique_ratio = NULL)
```

## Arguments

- sample:

  A `NIPTeRSample` object.

- ...:

  One or more `NIPTeRSexModel` objects, a single list of such models, or
  one `NIPTReferenceModel` containing pre-built `sex_models`.

- y_unique_ratio:

  Optional numeric scalar; a pre-computed Y-unique ratio (from
  [`nipter_y_unique_ratio`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_y_unique_ratio.md))
  for the sample. This is only used when one of the models has
  `method = "y_unique"`. If a `"y_unique"` model is present and this
  argument is `NULL`, that model is skipped with a warning.

## Value

A list-like `NIPTeRSexPrediction` S7 object with elements:

- prediction:

  Character; `"male"` or `"female"` (consensus when multiple models).

- model_predictions:

  Named character vector of per-model predictions (names are the method
  names).

- fractions:

  Named numeric vector with `x_fraction` and `y_fraction` for the
  sample.

- sample_name:

  Character; the sample name.

## Details

For each model, the sample's sex chromosome fractions are computed and
classified using `predict(model$model, newdata = ...)`. The
classification is mapped to `"male"`/`"female"` using the model's
`male_cluster`.

For `"y_unique"` models, the `y_unique_ratio` argument is used as the
input feature instead of chromosome fractions derived from binned read
counts.

With multiple models, the consensus is the label that appears most
frequently (majority vote). In case of a tie, `"female"` is returned
(conservative default for NIPT, as false-male calls could mask sex
chromosome aneuploidies).

## See also

[`nipter_sex_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md),
[`nipter_sex_model_y_unique()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model_y_unique.md),
[`nipter_y_unique_ratio()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_y_unique_ratio.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single model
pred <- nipter_predict_sex(test_sample, sex_model)
pred$prediction

# Consensus from three models (majority vote)
yr <- nipter_y_unique_ratio("sample.bam")
pred <- nipter_predict_sex(test_sample, model_y, model_xy, model_yu,
                           y_unique_ratio = yr$ratio)
pred$prediction
} # }
```
