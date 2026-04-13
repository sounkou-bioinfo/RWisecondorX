# Changelog

## RWisecondorX (development version)

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
