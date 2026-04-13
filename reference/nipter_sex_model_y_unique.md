# Build a sex prediction model from Y-unique ratios

Fits a two-component Gaussian mixture model (GMM) on Y-unique region
read ratios, which are computed from BAM files by
[`nipter_y_unique_ratio`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_y_unique_ratio.md).

## Usage

``` r
nipter_sex_model_y_unique(ratios)
```

## Arguments

- ratios:

  Named numeric vector of Y-unique ratios (one per sample). Typically
  obtained by calling
  [`nipter_y_unique_ratio`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_y_unique_ratio.md)
  on each BAM in the control cohort.

## Value

An object of class `"NIPTeRSexModel"` with elements:

- model:

  The
  [`mclust::Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html)
  fitted object.

- method:

  `"y_unique"`.

- male_cluster:

  Integer (1 or 2); which cluster is male (higher Y-unique ratio).

- classifications:

  Named character vector of `"male"`/ `"female"` labels for each input
  sample.

- fractions:

  The input `ratios` vector (named).

## Details

This is a companion to
[`nipter_sex_model`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md),
which operates on binned `NIPTeRSample` fractions. The Y-unique model
operates at the BAM level and does not require prior binning. The
resulting model object is compatible with
[`nipter_predict_sex`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
for majority-vote consensus.

The algorithm:

1.  Fit a two-component Gaussian mixture with equal mixing proportions
    (`mclust::Mclust(ratios, G = 2, control = mclust::emControl(equalPro = TRUE))`).

2.  Identify the male cluster as the component with the higher median
    Y-unique ratio.

## See also

[`nipter_y_unique_ratio()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_y_unique_ratio.md),
[`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md),
[`nipter_sex_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute Y-unique ratios for a cohort
bams <- list.files("bams/", pattern = "\\.bam$", full.names = TRUE)
ratios <- vapply(bams, function(b) nipter_y_unique_ratio(b)$ratio,
                 numeric(1L))
names(ratios) <- basename(bams)

# Build sex model
model_yu <- nipter_sex_model_y_unique(ratios)
model_yu$classifications
} # }
```
