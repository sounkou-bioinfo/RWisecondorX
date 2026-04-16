# Score sex chromosomes with NCV models from a typed reference

Computes sex-matched X or Y NCV scores from the `$sex_ncv_models`
attached to a `NIPTReferenceModel`.

## Usage

``` r
nipter_ncv_sex_score(
  sample,
  reference,
  focus_chromosome = c("X", "Y"),
  y_unique_ratio = NULL
)
```

## Arguments

- sample:

  A `NIPTeRSample` or typed `NIPTSample`.

- reference:

  A `NIPTReferenceModel` with `$sex_ncv_models` already built by
  [`nipter_build_sex_ncv_models`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_sex_ncv_models.md).

- focus_chromosome:

  Character scalar, `"X"` or `"Y"`.

- y_unique_ratio:

  Optional scalar passed through to
  [`nipter_predict_sex`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
  if the reference carries a `"y_unique"` sex model.

## Value

A typed `NIPTSexNCVScore`.
