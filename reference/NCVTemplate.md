# Pre-computed NCV denominator template

Returned by `nipter_prepare_ncv()` and consumed by
[`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md).

## Usage

``` r
NCVTemplate(
  focus_chromosome = integer(0),
  denominators = character(0),
  ctrl_mean = numeric(0),
  ctrl_sd = numeric(0),
  ctrl_cv = numeric(0),
  shapiro_p = numeric(0),
  test_z_scores = numeric(0),
  test_sample_names = character(0),
  train_sample_names = character(0)
)
```

## Arguments

- focus_chromosome:

  Integer chromosome identifier of the numerator.

- denominators:

  Character vector of denominator chromosome labels.

- ctrl_mean:

  Mean control ratio for the chosen denominator set.

- ctrl_sd:

  Standard deviation of the control ratios.

- ctrl_cv:

  Coefficient of variation of the control ratios.

- shapiro_p:

  Shapiro-Wilk p-value for the control-ratio distribution.

- test_z_scores:

  Numeric vector of held-out control z-scores.

- test_sample_names:

  Names corresponding to `test_z_scores`.

- train_sample_names:

  Control-sample names used to fit the template.
