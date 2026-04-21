# Predict copy number aberrations using WisecondorX

Native R implementation of the WisecondorX `predict` pipeline. Takes a
binned sample and a reference (from
[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md))
and detects copy number aberrations using within-sample normalization,
PCA correction, and circular binary segmentation (CBS).

## Usage

``` r
rwisecondorx_predict(
  sample,
  reference,
  sample_binsize = NULL,
  outprefix = NULL,
  minrefbins = 150L,
  maskrepeats = 5L,
  alpha = 1e-04,
  zscore = 5,
  optimal_cutoff_sd_multiplier = 3,
  within_sample_mask_iterations = 3L,
  within_sample_mask_quantile = 0.99,
  cbs_split_min_gap_bp = 2000000L,
  segment_zscore_cap = 1000,
  beta = NULL,
  blacklist = NULL,
  gender = NULL,
  seed = NULL,
  parallel = TRUE,
  cpus = 4L
)
```

## Arguments

- sample:

  Named list of integer vectors keyed by chromosome (`"1"`–`"24"`), as
  returned by
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

- reference:

  A `WisecondorXReference` object from
  [`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
  or a list with equivalent structure.

- sample_binsize:

  Integer; bin size of the input sample. Used for rescaling to the
  reference bin size.

- outprefix:

  Character; path prefix for output files. BED and statistics files will
  be written as `<outprefix>_bins.bed`, `<outprefix>_segments.bed`,
  `<outprefix>_aberrations.bed`, `<outprefix>_statistics.txt`. If `NULL`
  (default), no files are written and results are returned as a list.

- minrefbins:

  Integer; minimum reference bins per target. Bins with fewer references
  are zeroed out. Default `150L`.

- maskrepeats:

  Integer; number of iterative distance-masking cycles. Default `5L`.

- alpha:

  Numeric; CBS breakpoint p-value threshold. Default `1e-4`.

- zscore:

  Numeric; Z-score cutoff for aberration calling. Default `5`.

- optimal_cutoff_sd_multiplier:

  Numeric; number of population standard deviations added to the mean
  distance when iteratively deriving the optimal within-reference
  cutoff. Default `3`.

- within_sample_mask_iterations:

  Integer; number of iterative within-sample masking passes. Default
  `3L`.

- within_sample_mask_quantile:

  Numeric scalar in `(0, 1)`; bins with
  `|z| >= qnorm(within_sample_mask_quantile)` are excluded from the next
  within-sample normalization pass. Default `0.99`.

- cbs_split_min_gap_bp:

  Integer; NA stretches spanning more than this many base pairs trigger
  CBS segment splitting. Default `2000000L`.

- segment_zscore_cap:

  Numeric; absolute cap applied to segment z-scores after they are
  derived from the null-ratio distributions. Default `1000`.

- beta:

  Optional numeric; if given, ratio-based cutoff is used instead of
  Z-score. Should approximate purity (0, 1\]. Default `NULL`.

- blacklist:

  Optional character; path to a headerless BED file of regions to mask.

- gender:

  Optional character; force gender (`"F"` or `"M"`).

- seed:

  Optional integer; RNG seed for CBS reproducibility.

- parallel:

  Logical; use
  [`ParDNAcopy::parSegment()`](https://rdrr.io/pkg/ParDNAcopy/man/parSegment.html)
  for CBS when `TRUE` (default). Requires the `ParDNAcopy` package. Set
  `parallel = FALSE` to use
  [`DNAcopy::segment()`](https://rdrr.io/pkg/DNAcopy/man/segment.html)
  explicitly.

- cpus:

  Integer; number of threads for parallel CBS (`parSegment`) and any
  other OpenMP-accelerated steps. Default `4L`.

## Value

A list with class `"WisecondorXPrediction"` containing:

- results_r:

  List of per-chromosome log2-ratio vectors.

- results_z:

  List of per-chromosome Z-score vectors.

- results_w:

  List of per-chromosome weight vectors.

- results_c:

  Data frame of CBS segments (`chr`, `start`, `end`, `zscore`, `ratio`).

- aberrations:

  Data frame of called aberrations.

- statistics:

  Data frame of per-chromosome statistics.

- gender:

  Predicted (or forced) gender.

- algorithm_params:

  Named list of native prediction parameters used for this result.

- n_reads:

  Total read count.

- binsize:

  Reference bin size.

## Details

CBS is performed using the DNAcopy Bioconductor package (or optionally
ParDNAcopy for parallel segmentation). Both must be available in
`Suggests`. Unlike
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md),
this function's `parallel` argument controls native R CBS execution
through
[`ParDNAcopy::parSegment()`](https://rdrr.io/pkg/ParDNAcopy/man/parSegment.html);
the upstream Python CLI wrapper delegates segmentation to upstream
WisecondorX and does not expose this native switch.

## See also

[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`scale_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/scale_sample.md)
