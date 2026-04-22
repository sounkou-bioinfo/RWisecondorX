# Score multiple samples against a typed gaunosome reference

Batch wrapper around
[`nipter_gaunosome_score`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gaunosome_score.md).
Scores each sample against one
[`NIPTReferenceModel`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTReferenceModel.md)
and returns a typed cohort-level report with both the per-sample objects
and one flattened summary table.

## Usage

``` r
nipter_gaunosome_report(
  samples,
  reference,
  y_unique_ratios = NULL,
  min_controls = 2L,
  sample_predictors = NULL,
  focus_chromosomes = c("X", "Y")
)
```

## Arguments

- samples:

  A non-empty list of `NIPTeRSample`/`NIPTSample` objects, or a
  `NIPTControlGroup` whose `$samples` should be scored.

- reference:

  A `NIPTReferenceModel` prepared with
  [`nipter_build_gaunosome_models`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_gaunosome_models.md).

- y_unique_ratios:

  Optional named numeric vector keyed by sample name, used when the
  reference includes a `"y_unique"` sex model.

- min_controls:

  Minimum number of non-outlier same-sex controls required for the
  z-score component. Default `2L`.

- sample_predictors:

  Optional named list keyed by sample name. Each value should be a named
  list of extra regression predictors for that sample, used when the
  fitted models include extra columns such as `gc_read_perc_post`.

- focus_chromosomes:

  Character vector; any subset of `c("X", "Y")`.

## Value

A typed `NIPTGaunosomeReport`.
