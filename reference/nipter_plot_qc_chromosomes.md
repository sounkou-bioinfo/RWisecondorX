# Plot chromosome-level NIPTeR control-group QC

Draws the chromosome-level coefficient of variation profile from a
[`NIPTControlGroupQC`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTControlGroupQC.md)
object.

## Usage

``` r
nipter_plot_qc_chromosomes(
  qc,
  theme = c("minimal", "light", "bw", "classic"),
  base_size = 11,
  method_palette = NULL,
  shapiro_shapes = NULL,
  cv_step = 0.05
)
```

## Arguments

- qc:

  A `NIPTControlGroupQC` from
  [`nipter_control_group_qc`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_control_group_qc.md).

- theme:

  Theme name. One of `"minimal"` (default), `"light"`, `"bw"`, or
  `"classic"`.

- base_size:

  Base font size passed to the selected ggplot2 theme.

- method_palette:

  Optional named color vector keyed by method labels.

- shapiro_shapes:

  Optional named shape vector keyed by Shapiro labels `"<0.05"`,
  `">=0.05"`, and `"NA"`.

- cv_step:

  Y-axis increment for CV facets. Default `0.05`.

## Value

A `ggplot` object.
