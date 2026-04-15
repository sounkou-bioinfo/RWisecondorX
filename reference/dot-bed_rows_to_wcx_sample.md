# Convert BED rows (from rduckhts_tabix_multi) to WisecondorX sample format

Convert BED rows (from rduckhts_tabix_multi) to WisecondorX sample
format

## Usage

``` r
.bed_rows_to_wcx_sample(rows)
```

## Arguments

- rows:

  data.frame with columns chrom, start_pos, end_pos, count.

## Value

Named list of integer vectors keyed by chromosome ("1"–"24").
