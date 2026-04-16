# Build a NIPTeR control group from a list of binned samples

Collects multiple `NIPTeRSample` objects into a `NIPTeRControlGroup`,
the reference cohort used by
[`nipter_z_score`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md),
[`nipter_ncv_score`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md),
[`nipter_chi_correct`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md),
and
[`nipter_regression`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md).

## Usage

``` r
nipter_as_control_group(
  samples,
  description = "General control group",
  sample_sex = NULL,
  sex_source = NULL
)
```

## Arguments

- samples:

  A list of `NIPTeRSample` objects.

- description:

  Label for the group (default `"General control group"`).

- sample_sex:

  Optional character vector of known sex labels for the control samples.
  Accepted values are `"female"`, `"male"`, `"ambiguous"`, and
  `"unknown"`. When unnamed, values are matched in sample order; when
  named, names must match `sample_name`.

- sex_source:

  Optional string describing where `sample_sex` came from, e.g.
  `"explicit"`, `"consensus_gmm"`, or `"laboratory_lims"`.

## Value

An object of class `c("NIPTeRControlGroup", <strand_type>)`.

## Details

All samples must share the same strand type (`"CombinedStrands"` or
`"SeparatedStrands"`). Duplicate sample names are silently removed.

## See also

[`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md),
[`nipter_match_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_control_group.md),
[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)

## Examples

``` r
if (FALSE) { # \dontrun{
samples <- lapply(bam_files, nipter_bin_bam)
cg <- nipter_as_control_group(samples)
} # }
```
