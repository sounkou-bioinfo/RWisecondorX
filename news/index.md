# Changelog

## RWisecondorX (development version)

### Sex prediction via Gaussian mixture models

- New
  [`nipter_sex_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md)
  fits a 2-component GMM on sex chromosome fractions from a
  `NIPTeRControlGroup` using
  [`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html).
  Supports `"y_fraction"` (univariate on Y-chromosome fraction) and
  `"xy_fraction"` (bivariate on X + Y fractions) methods. The male
  cluster is identified as the component with higher median Y fraction.

- New
  [`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
  classifies a `NIPTeRSample` as male or female given one or more
  `NIPTeRSexModel` objects. Multiple models use majority vote consensus
  (tie defaults to “female” — conservative for NIPT).

- `mclust` added to `Suggests` in DESCRIPTION.

### `nipter_bin_bam_bed()` SeparatedStrands output

- [`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md)
  gains a `separate_strands` parameter. When `TRUE`, outputs a 7-column
  BED (`chrom`, `start`, `end`, `count`, `count_fwd`, `count_rev`,
  `corrected_count`) where `count = count_fwd + count_rev`. When `FALSE`
  (default), the 5-column BED format is unchanged.

### SeparatedStrands support

- [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
  gains a `separate_strands` parameter. When `TRUE`, returns a list with
  `fwd` (forward strand) and `rev` (reverse strand) data frames.

- [`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
  gains `separate_strands = TRUE` support, producing `NIPTeRSample`
  objects with class `c("NIPTeRSample", "SeparatedStrands")`. Autosomal
  chromosome reads are stored as a list of two matrices (forward and
  reverse) with rownames `"1F".."22F"` and `"1R".."22R"`.

- All NIPTeR statistical functions now support SeparatedStrands samples:

  - [`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md):
    LOESS/bin-weight fitted on summed strand counts, corrections applied
    independently to each strand matrix.
  - [`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md):
    chi-squared computed on summed strands, correction applied
    per-strand.
  - [`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md)
    and
    [`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md):
    use collapsed (F+R) fractions.
  - [`nipter_regression()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md):
    doubles the predictor pool (44 candidates: `"1F".."22F"`,
    `"1R".."22R"`) with complementary exclusion within each model
    (selecting `"5F"` excludes both `"5F"` and `"5R"`).

### WisecondorX CLI wrapper fixes

- [`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md):
  removed hallucinated `ref_binsize` / `--binsize` parameter that does
  not exist in the upstream `predict` subcommand.

- [`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md):
  fixed `ref_binsize` default from `50000L` to `100000L` to match
  upstream default.

- [`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md):
  added `add_plot_title` parameter mapping to upstream
  `--add-plot-title`.
