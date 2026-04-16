# Build a chromosome-level reference frame from a NIPTeR control group

Produces the per-sample count and fraction table that downstream
sex-aware NIPT model building actually needs. This keeps
application-level gaunosome modelling data out of the raw `NIPTeRSample`
class while providing a stable training frame for future sex-chromosome
Z-score, NCV, and regression models.

## Usage

``` r
nipter_reference_frame(control_group, sample_sex = NULL)
```

## Arguments

- control_group:

  A `NIPTeRControlGroup` object.

- sample_sex:

  Optional character vector overriding the sex labels stored on
  `control_group`. Accepted values are `"female"`, `"male"`,
  `"ambiguous"`, and `"unknown"`.

## Value

A `data.frame` with one row per sample and columns: `Sample_name`,
optional `SampleSex`, `NChrReads_*`, and `FrChrReads_*` for chromosomes
`1:22`, `X`, and `Y`.
