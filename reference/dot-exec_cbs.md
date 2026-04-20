# Execute CBS on WisecondorX results

Runs circular binary segmentation on log2-ratios using weights, then
splits segments spanning large NA gaps and recalculates segmental
ratios. Exact port of upstream `CBS.R` logic.

## Usage

``` r
.exec_cbs(
  results_r,
  results_w,
  ref_gender,
  alpha,
  binsize,
  seed,
  parallel = TRUE,
  cpus = 1L
)
```

## Arguments

- results_r:

  List of per-chromosome log2-ratio vectors.

- results_w:

  List of per-chromosome weight vectors.

- ref_gender:

  `"F"` or `"M"` — determines whether chrY is included.

- alpha:

  CBS breakpoint p-value threshold.

- binsize:

  Reference bin size in bp.

- seed:

  Optional RNG seed.

- parallel:

  Logical; use ParDNAcopy when `TRUE` (default). This requires the
  `ParDNAcopy` package to be installed. Set `parallel = FALSE` to use
  [`DNAcopy::segment()`](https://rdrr.io/pkg/DNAcopy/man/segment.html)
  explicitly.

- cpus:

  Integer; number of threads passed to `parSegment()`. Only used when
  `parallel = TRUE` and ParDNAcopy is available.

## Value

Data frame with columns `chr` (integer, 1-based), `start` (integer,
0-based bin index), `end` (integer, exclusive bin index), `ratio`
(numeric, weighted mean segment ratio).
