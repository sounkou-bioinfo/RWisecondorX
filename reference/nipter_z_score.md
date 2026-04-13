# Chromosomal Z-score

Computes the Z-score for a given chromosome of a test sample against a
control group. The Z-score measures how many standard deviations the
sample's chromosomal fraction deviates from the control mean.

## Usage

``` r
nipter_z_score(sample, control_group, chromo_focus)
```

## Arguments

- sample:

  A `NIPTeRSample` object (the test sample).

- control_group:

  A `NIPTeRControlGroup` object.

- chromo_focus:

  Integer; the chromosome to test (1-22).

## Value

A list of class `"NIPTeRZScore"` with elements:

- sample_z_score:

  The test sample's Z-score (numeric scalar).

- focus_chromosome:

  Character; the tested chromosome.

- control_statistics:

  Named numeric vector with `mean`, `sd`, and `shapiro_p_value`.

- control_z_scores:

  Named numeric vector of Z-scores for each control sample.

- correction_status:

  Character vector of correction statuses.

- sample_name:

  Character; the test sample name.

## See also

[`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md),
[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md),
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)

## Examples

``` r
if (FALSE) { # \dontrun{
z21 <- nipter_z_score(sample, cg, chromo_focus = 21)
z21$sample_z_score
} # }
```
