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
  for CBS when available. Default `TRUE`. Falls back to
  [`DNAcopy::segment()`](https://rdrr.io/pkg/DNAcopy/man/segment.html)
  with a message if ParDNAcopy is not installed.

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

- n_reads:

  Total read count.

- binsize:

  Reference bin size.

## Details

CBS is performed using the DNAcopy Bioconductor package (or optionally
ParDNAcopy for parallel segmentation). Both must be available in
`Suggests`.

## See also

[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`scale_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/scale_sample.md)
