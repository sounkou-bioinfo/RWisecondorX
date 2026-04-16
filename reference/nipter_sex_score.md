# Score sex chromosomes against a typed NIPT reference model

Computes sex-matched X/Y z-scores from the enriched
[`NIPTReferenceModel`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTReferenceModel.md)
built by
[`nipter_build_reference`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md).
The scoring mirrors the production pipeline's gaunosome z-score stage:
predicted fetal sex is used to select the relevant non-outlier reference
subset, and the sample's X/Y fractions are standardized against that
subset.

## Usage

``` r
nipter_sex_score(sample, reference, y_unique_ratio = NULL, min_controls = 2L)
```

## Arguments

- sample:

  A `NIPTeRSample` or typed `NIPTSample`.

- reference:

  A `NIPTReferenceModel` built by
  [`nipter_build_reference`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md).

- y_unique_ratio:

  Optional numeric scalar passed through to
  [`nipter_predict_sex`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
  when the reference model includes a `"y_unique"` sex model.

- min_controls:

  Minimum number of non-outlier reference samples required in a sex
  class to compute a score. Default `2L`.

## Value

A list-like `NIPTSexScore` S7 object with:

- sample_name:

  The sample identifier.

- predicted_sex:

  Consensus male/female call from the reference model's sex classifiers.

- sex_prediction:

  The underlying `NIPTeRSexPrediction` object.

- sample_metrics:

  Named numeric vector containing `FrChrReads_X`, `FrChrReads_Y`,
  `RR_X`, and `RR_Y`.

- z_scores:

  Named numeric vector containing the selected, XX-reference, and
  XY-reference X/Y z-scores.

- cv:

  Named numeric vector containing the corresponding coefficients of
  variation.

- reference_sizes:

  Named numeric vector with female, male, and selected same-sex
  reference counts after outlier exclusion.

- reference_sample_names:

  Character vector of the same-sex control samples actually used for the
  selected score.

## Details

The returned object also includes the alternate XX and XY
reference-group scores so downstream code can inspect both hypotheses
without rebuilding the reference statistics.

## See also

[`nipter_build_reference()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md),
[`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
