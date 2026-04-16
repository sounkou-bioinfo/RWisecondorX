# NIPTeR control group for SeparatedStrands samples

NIPTeR control group for SeparatedStrands samples

## Usage

``` r
SeparatedControlGroup(
  samples = list(),
  description = "General control group",
  sample_sex = NULL,
  sex_source = NULL,
  .cache = new.env(parent = emptyenv())
)
```

## Arguments

- samples:

  List of `NIPTeRSample` or `NIPTSample` objects.

- description:

  Human-readable label for the control group.

- sample_sex:

  Optional named character vector keyed by sample name with values in
  `female`, `male`, `ambiguous`, or `unknown`.

- sex_source:

  Optional scalar string describing the origin of `sample_sex`.

- .cache:

  Internal environment used for lazy cached summaries.
