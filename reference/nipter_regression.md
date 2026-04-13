# Regression-based Z-score for trisomy prediction

Implements NIPTeR's stepwise linear regression approach. Builds multiple
regression models that predict the focus chromosome's fraction from
other chromosomal fractions, then computes Z-scores based on the
prediction residuals.

## Usage

``` r
nipter_regression(
  sample,
  control_group,
  chromo_focus,
  n_models = 4L,
  n_predictors = 4L,
  exclude_chromosomes = c(13L, 18L, 21L),
  include_chromosomes = NULL,
  train_fraction = 0.6,
  overdispersion_rate = 1.15,
  force_practical_cv = FALSE,
  seed = NULL
)
```

## Arguments

- sample:

  A `NIPTeRSample` object.

- control_group:

  A `NIPTeRControlGroup` object.

- chromo_focus:

  Integer; the target chromosome (1-22).

- n_models:

  Number of regression models to build (default 4).

- n_predictors:

  Maximum number of predictor chromosomes per model (default 4).

- exclude_chromosomes:

  Integer vector of chromosomes to exclude from the predictor pool
  (default `c(13, 18, 21)`).

- include_chromosomes:

  Integer vector of chromosomes to force-include.

- train_fraction:

  Fraction of control samples used for training (default 0.6).

- overdispersion_rate:

  Overdispersion multiplier for theoretical CV (default 1.15).

- force_practical_cv:

  Logical; always use practical CV regardless of theoretical CV? Default
  `FALSE`.

- seed:

  Random seed for the train/test split. Default `NULL` (no seed set).

## Value

A list of class `"NIPTeRRegression"` with elements:

- models:

  A list of `n_models` sublists, each containing: `z_score` (numeric),
  `cv` (selected CV), `cv_type` (`"practical"` or `"theoretical"`),
  `predictors` (character vector of predictor chromosomes),
  `shapiro_p_value`, `control_z_scores` (named numeric vector).

- focus_chromosome:

  Character.

- correction_status:

  Character vector.

- sample_name:

  Character.

## Details

For `SeparatedStrands` samples the predictor pool is doubled (each
autosomal chromosome contributes both a forward and reverse fraction),
and complementary-strand exclusion within each model prevents selecting
both strands of the same chromosome as predictors in the same model.
Across models, only the exact predictor string is excluded; the
complementary strand remains available.

The algorithm:

1.  Split the control group 60/40 into train and test sets.

2.  Compute chromosomal fractions for all samples.

3.  For each of `n_models` models, greedily select `n_predictors`
    chromosomes by maximising adjusted R-squared (forward stepwise
    selection). Predictors already used in previous models are excluded.

4.  For each model, fit `lm(focus ~ predictors)` on the training set,
    predict on the test set, and compute:

    - Practical CV: `sd(obs/pred) / mean(obs/pred)` on the test set.

    - Theoretical CV: `overdispersion_rate / sqrt(total_reads_focus)`
      for the test sample.

    - The larger CV is used (unless `force_practical_cv = TRUE`).

    - Sample Z-score: `(sample_ratio - 1) / cv`.

## See also

[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md),
[`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
reg21 <- nipter_regression(sample, cg, chromo_focus = 21)
reg21$models[[1]]$z_score
} # }
```
