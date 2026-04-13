# Split segments that span large NA gaps

Exact port of the upstream CBS.R segment-splitting logic.

## Usage

``` r
.split_na_segments(seg_df, full_df, na_threshold)
```

## Arguments

- seg_df:

  Data frame with chr, s, e, r.

- full_df:

  Data frame with chromosome, x, y, w.

- na_threshold:

  Consecutive NA stretch that triggers a split.

## Value

Modified segment data frame.
