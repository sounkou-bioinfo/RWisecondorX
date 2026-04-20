# Plot autosomal bin-level NIPTeR control-group QC

Draws either the scaled-count CV or chi-normalised overdispersion
profile across autosomal bins from a
[`NIPTControlGroupQC`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTControlGroupQC.md)
object built with `include_bins = TRUE`.

## Usage

``` r
nipter_plot_qc_bins(qc, metric = c("cv_scaled", "chi_z"))
```

## Arguments

- qc:

  A `NIPTControlGroupQC` from
  [`nipter_control_group_qc`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_control_group_qc.md)
  with `$bin_summary`.

- metric:

  Either `"cv_scaled"` (default) or `"chi_z"`.

## Value

A `ggplot` object.
