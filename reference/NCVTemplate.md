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
