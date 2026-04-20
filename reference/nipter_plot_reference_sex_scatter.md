# Plot reference sex-model scatter spaces

Draws either the X/Y fraction space or the RR_X/RR_Y ratio space used by
the reference sex and gaunosome models.

## Usage

``` r
nipter_plot_reference_sex_scatter(reference, space = c("fraction", "ratio"))
```

## Arguments

- reference:

  A `NIPTReferenceModel` from
  [`nipter_build_reference`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md).

- space:

  Either `"fraction"` (default) or `"ratio"`.

## Value

A `ggplot` object.
