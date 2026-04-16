# Build a typed QC report for a NIPTeR control group

Summarises the stability of a `NIPTControlGroup` at the chromosome,
sample, and optional autosomal-bin level. The chromosome summary exposes
the mean, standard deviation, coefficient of variation, and Shapiro-Wilk
normality statistic for each autosome. The sample summary adds SSD-based
matching scores so aberrant controls can be identified before scoring.
When sex labels are available, the report also includes sex-stratified
`RR_X`/`RR_Y` spread metrics relevant for gaunosome model readiness.
Optional bin-level output exposes the scaled-count CV and chi-correction
profile that underlies
[`nipter_chi_correct`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md).

## Usage

``` r
nipter_control_group_qc(
  control_group,
  sample_sex = NULL,
  chi_cutoff = 3.5,
  include_bins = FALSE
)
```

## Arguments

- control_group:

  A `NIPTControlGroup`.

- sample_sex:

  Optional character vector overriding the sex labels stored on
  `control_group`. Accepted values are `"female"`, `"male"`,
  `"ambiguous"`, and `"unknown"`.

- chi_cutoff:

  Numeric scalar used to flag overdispersed bins in the optional chi
  profile. Default `3.5`.

- include_bins:

  Logical; when `TRUE`, compute the autosomal bin-level scaled-count CV
  and chi profile. Default `FALSE`.

## Value

A typed `NIPTControlGroupQC` object.

## See also

[`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md),
[`nipter_match_matrix()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_matrix.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)
