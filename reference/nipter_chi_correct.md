# Chi-squared correction for overdispersed bins

Identifies bins in the control group with excess variance
(overdispersion) using a chi-squared test and downweights them. This
reduces the influence of noisy bins on downstream Z-scores and NCV
calculations.

## Usage

``` r
nipter_chi_correct(
  sample,
  control_group,
  chi_cutoff = 3.5,
  include_sex = FALSE
)
```

## Arguments

- sample:

  A `NIPTeRSample` object (the test sample).

- control_group:

  A `NIPTeRControlGroup` object.

- chi_cutoff:

  Normalised chi-squared threshold. Bins with
  `(chi - df) / sqrt(2*df) > chi_cutoff` are corrected. Default `3.5`.

- include_sex:

  Logical; correct sex chromosomes as well? Default `FALSE`. Sex
  chromosome bins use the same chi-squared weights derived from
  autosomes.

## Value

A list with two elements:

- sample:

  The corrected `NIPTeRSample`.

- control_group:

  The corrected `NIPTeRControlGroup`.

## Details

The correction is applied simultaneously to both the test sample and all
control group samples, maintaining consistency.

The algorithm follows NIPTeR's chi-squared correction:

1.  For each control sample, scale autosomal bin counts so that total
    reads match the overall mean across all control samples.

2.  Compute the expected count per bin (mean of scaled counts).

3.  Compute chi-squared per bin: \\\chi^2 = \sum_i (expected -
    scaled_i)^2 / expected\\.

4.  Normalise: \\z = (\chi^2 - df) / \sqrt{2 \cdot df}\\ where \\df =
    n\_{controls} - 1\\.

5.  For bins where \\z \> \text{chi\\cutoff}\\: divide reads by \\\chi^2
    / df\\.

## See also

[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md),
[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md),
[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- nipter_chi_correct(sample, control_group)
sample_corrected <- result$sample
cg_corrected     <- result$control_group
} # }
```
