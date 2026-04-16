# Build a typed NIPT reference model

Packages the control group, chromosome-level training frame, and
optional sex-prediction models into one validated reference object. The
embedded `reference_frame` is enriched with the derived sex-scoring
columns used by the production pipeline: `ConsensusGender`, `RR_X`,
`RR_Y`, `RR_X_SexClassMAD`, `RR_Y_SexClassMAD`, and `IsRefSexOutlier`.
When `y_unique_ratios` are supplied, the frame also carries
`YUniqueRatio`.

## Usage

``` r
nipter_build_reference(
  control_group,
  sample_sex = NULL,
  sex_source = NULL,
  sex_methods = c("y_fraction", "xy_fraction"),
  y_unique_ratios = NULL,
  build_params = list()
)
```

## Arguments

- control_group:

  A `NIPTControlGroup`.

- sample_sex:

  Optional character vector overriding the sex labels stored on
  `control_group`. Accepted values are `"female"`, `"male"`,
  `"ambiguous"`, and `"unknown"`.

- sex_source:

  Optional scalar string describing where `sample_sex` came from when
  supplied here, for example `"explicit"` or `"laboratory_lims"`.

- sex_methods:

  Character vector of fraction-based sex model methods to build. Allowed
  values are `"y_fraction"` and `"xy_fraction"`. Default builds both.

- y_unique_ratios:

  Optional named numeric vector of Y-unique ratios for building an
  additional `"y_unique"` sex model and storing `YUniqueRatio` in the
  reference frame.

- build_params:

  Optional named list of provenance parameters to attach to the returned
  reference model.

## Value

A list-like `NIPTReferenceModel` S7 object containing:

- control_group:

  The input control group, optionally re-annotated with explicit sample
  sex labels.

- reference_frame:

  A typed `NIPTReferenceFrame` with one row per control sample,
  including chromosome counts/fractions and the derived sex-aware
  reference columns used for X/Y scoring.

- sex_models:

  Named list of `NIPTeRSexModel` objects keyed by method.

- sample_sex_source:

  The provenance string for the control-group sex labels, if available.

- build_date:

  UTC timestamp for the build.

- build_params:

  Caller-supplied provenance metadata.
