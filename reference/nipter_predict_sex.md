# Predict fetal sex for a NIPTeR sample

Uses one or more `NIPTeRSexModel` objects to classify a sample as male
or female. When multiple models are supplied, a majority vote determines
the consensus call.

## Usage

``` r
nipter_predict_sex(sample, ...)
```

## Arguments

- sample:

  A `NIPTeRSample` object.

- ...:

  One or more `NIPTeRSexModel` objects, or a single list of models.

## Value

A list of class `"NIPTeRSexPrediction"` with elements:

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

With multiple models, the consensus is the label that appears most
frequently (majority vote). In case of a tie, `"female"` is returned
(conservative default for NIPT, as false-male calls could mask sex
chromosome aneuploidies).

## See also

[`nipter_sex_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single model
pred <- nipter_predict_sex(test_sample, sex_model)
pred$prediction

# Consensus from two models (majority vote)
pred <- nipter_predict_sex(test_sample, model_y, model_xy)
pred$prediction
} # }
```
