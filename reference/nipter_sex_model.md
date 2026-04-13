# Build a sex prediction model from a NIPTeR control group

Fits a two-component Gaussian mixture model (GMM) on sex chromosome
fractions derived from a `NIPTeRControlGroup`. The model distinguishes
male and female samples based on Y-chromosome read fraction (univariate)
or X+Y chromosome fractions (bivariate).

## Usage

``` r
nipter_sex_model(control_group, method = c("y_fraction", "xy_fraction"))
```

## Arguments

- control_group:

  A `NIPTeRControlGroup` object. Ideally chi-squared corrected via
  [`nipter_chi_correct`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md).

- method:

  Character; the feature space for the GMM:

  `"y_fraction"`

  :   Univariate: Y-chromosome read fraction relative to total autosomal
      reads.

  `"xy_fraction"`

  :   Bivariate: X and Y chromosome fractions.

## Value

An object of class `"NIPTeRSexModel"` with elements:

- model:

  The
  [`mclust::Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html)
  fitted object.

- method:

  Character; the method used.

- male_cluster:

  Integer (1 or 2); which cluster is male (higher Y fraction).

- classifications:

  Named character vector of `"male"`/ `"female"` labels for each control
  sample.

- fractions:

  Matrix of the input features used for fitting (samples as rows).
  Columns are `"y_fraction"` or `c("x_fraction", "y_fraction")`.

## Details

This mirrors the approach used in clinical NIPT pipelines (see
*Details*), ported to R using
[`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html)
in place of Python's `sklearn.GaussianMixture`.

The algorithm:

1.  Compute sex chromosome read fractions for every sample in the
    control group. Y fraction = `sum(Y bins) / sum(autosomal bins)`; X
    fraction = `sum(X bins) / sum(autosomal bins)`.

2.  Fit a two-component Gaussian mixture with equal mixing proportions
    (`mclust::Mclust(data, G = 2, control = mclust::emControl(equalPro = TRUE))`).

3.  Identify the male cluster as the component with the higher median Y
    fraction.

This follows the clinical NIPT pipeline pattern of building sex models
from control cohorts. The user's pipeline builds three models (Y-unique
ratio from samtools, XY fractions, Y fraction) and takes a majority
vote;
[`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
implements the consensus when given multiple models.

## See also

[`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Build sex model from chi-corrected control group
sex_model <- nipter_sex_model(chi_cg, method = "y_fraction")

# Bivariate X+Y model
sex_model_xy <- nipter_sex_model(chi_cg, method = "xy_fraction")
} # }
```
