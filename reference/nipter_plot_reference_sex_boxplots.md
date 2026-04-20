# Plot reference sex-model distributions

Draws the per-sex X/Y fraction distributions used by the reference sex
models. When the reference frame carries `YUniqueRatio`, it is included
as an additional facet.

## Usage

``` r
nipter_plot_reference_sex_boxplots(reference)
```

## Arguments

- reference:

  A `NIPTReferenceModel` from
  [`nipter_build_reference`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md).

## Value

A `ggplot` object.
