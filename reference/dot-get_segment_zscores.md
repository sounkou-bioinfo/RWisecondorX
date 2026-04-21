# Compute between-sample Z-scores for CBS segments

For each segment, computes a Z-score by comparing the observed segment
ratio to a null distribution derived from the stored null ratios.
Mirrors `overall_tools.get_z_score()`.

## Usage

``` r
.get_segment_zscores(
  cbs_result,
  results_nr,
  results_r,
  results_w,
  zscore_cap = 1000
)
```

## Arguments

- cbs_result:

  Data frame with chr, start, end, ratio.

- results_nr:

  List of per-chromosome null-ratio matrices.

- results_r:

  List of per-chromosome log2-ratio vectors.

- results_w:

  List of per-chromosome weight vectors.

- zscore_cap:

  Numeric; absolute cap applied to finite segment z-scores.

## Value

Numeric vector of Z-scores, one per segment.
