# RWisecondorX (development version)

## BED.gz reader functions — close the round-trip

* New `bed_to_sample()` reads a 4-column BED.gz file (written by
  `bam_convert_bed()`) back into the named-list-of-integer-vectors format
  expected by `rwisecondorx_newref()` and `rwisecondorx_predict()`. This
  closes the WisecondorX round-trip: bin once with `bam_convert_bed()`, store
  the BED.gz, and reload for analysis without re-reading the BAM.

* New `bed_to_nipter_sample()` reads a 5-column (CombinedStrands) or 7-column
  (SeparatedStrands) BED.gz file (written by `nipter_bin_bam_bed()`) into a
  `NIPTeRSample` object compatible with all NIPTeR statistical functions. Column
  count is auto-detected: if the `strand` field in `read_bed()` contains
  non-null integers it is a 7-column SeparatedStrands file, otherwise 5-column
  CombinedStrands. Handles literal `"NA"` strings in the `corrected_count` field
  via `TRY_CAST`. Sample name is inferred from the filename or set explicitly.

* Both functions accept an optional DuckDB connection for reuse across multiple
  files, creating one internally (with `allow_unsigned_extensions = "true"`)
  when none is supplied.

* New `inst/tinytest/test_bed_reader.R` — 34 assertions covering WisecondorX
  4-column round-trip, NIPTeR CombinedStrands 5-column round-trip,
  SeparatedStrands 7-column round-trip (all four matrices), sample name
  inference, and integration with `scale_sample()`.

* Total test count: 379 assertions, all passing.

## Native WisecondorX implementation

* New `rwisecondorx_newref()` — pure-R/Rcpp implementation of the WisecondorX
  `newref` pipeline. Takes a list of binned samples (from `bam_convert()`) and
  builds a PCA-based reference: gender model training (2-component GMM on
  Y-fractions), global bin masking, per-partition normalize/PCA/filter/KNN
  reference building, and null ratio computation. Supports NIPT mode, custom
  Y-fraction cutoff, and multi-threaded KNN via OpenMP.

* New `rwisecondorx_predict()` — pure-R implementation of the WisecondorX
  `predict` pipeline. Coverage normalization, PCA projection, iterative
  within-sample normalization with aberration masking, gonosomal normalization,
  result inflation, log-transformation, optional blacklist masking, CBS
  segmentation (via DNAcopy or ParDNAcopy), segment Z-scoring, and aberration
  calling. Supports both Z-score and beta/ratio calling modes.

* New `rwisecondorx_write_bins_bed()`, `rwisecondorx_write_segments_bed()`,
  `rwisecondorx_write_aberrations_bed()`, and `rwisecondorx_write_statistics()`
  for writing prediction results to BED and statistics files.

* New `scale_sample()` for rescaling binned samples between different bin sizes
  (e.g. 5kb → 100kb).

* `R/rwisecondorx_utils.R` — shared utilities including `.train_gender_model()`
  (with mclust NULL fallback for zero-variance Y-fractions in all-female
  cohorts), `.gender_correct()`, `.get_mask()`, `.normalize_and_mask()`,
  `.train_pca()`, `.project_pc()`, `.predict_gender()`.

* `R/rwisecondorx_cbs.R` — `.exec_cbs()` wrapper around DNAcopy's `segment()`
  with the upstream WisecondorX conventions (0→NA, 0 weights→1e-99, split
  segments at large NA gaps, recalculate weighted means). Supports ParDNAcopy
  for parallel segmentation.

## Rcpp acceleration for KNN reference building

* New `src/knn_reference.cpp` implementing `knn_reference_cpp()` and
  `null_ratios_cpp()` in C++ with OpenMP parallelization. Replaces the
  O(n_bins² × n_samples) pure-R double for-loop with compiled code. The KNN
  reference-finding step that previously took minutes now completes in seconds.

* `Rcpp` added to `Imports` and `LinkingTo` in DESCRIPTION. OpenMP flags in
  `src/Makevars` and `src/Makevars.win`.

## Fix: KNN index semantics in predict normalization

* Fixed a correctness bug in `.normalize_once()` where global KNN indexes
  (stored by `knn_reference_cpp()`) were used directly to index into
  `chr_data`, the leave-one-chromosome-out subset. Global index `g` is now
  correctly mapped to local space: `g` if before the excluded chromosome,
  `g - n_chr` if after. This bug caused reference bins to be looked up at
  wrong positions during within-sample normalization.

* Note: upstream WisecondorX has the inverse issue — it stores LOCAL indexes
  during `newref` but uses them as GLOBAL indexes in the `null_ratio`
  computation. Our implementation stores GLOBAL indexes and now correctly
  handles both predict (local conversion) and null_ratio (global direct use).

## Fix: null ratio column count mismatch

* Fixed `rwisecondorx_predict()` crash when the autosomal sub-reference ("A")
  and gonosomal sub-reference ("F") have different numbers of null-ratio
  columns (due to different sample counts per partition). Both are now truncated
  to `min(ncol(aut), ncol(gon))` before `rbind()`.

## Fix: mclust GMM on zero-variance Y-fractions

* `.train_gender_model()` now handles the case where all female samples have
  exactly 0.0 Y-fraction (no Y reads). `mclust::Mclust(..., G=2)` returns
  `NULL` in this scenario; the fix falls back to a gap-based cutoff:
  `min(nonzero_y_fractions) / 2`.

## Synthetic cohort generator

* New `generate_cohort()` creates synthetic BAM files for testing using
  "compressed" chromosome lengths (100bp per bin instead of 100kb). Produces
  identical bin COUNT structure to GRCh37 at 100kb resolution with ~435KB
  BAMs. Supports injecting trisomy signal (extra reads on target chromosome).

* New `COMPRESSED_BINSIZE` constant (100L) exported for use with
  `bam_convert(binsize = COMPRESSED_BINSIZE)`.

* `inst/scripts/make_cohort.R` — CLI wrapper for cohort generation.

## Test infrastructure cleanup

* All test files now use `library(RWisecondorX)` instead of `source()` hacks
  to load R code directly. This is required because functions call compiled
  C++ code (`knn_reference_cpp`, `null_ratios_cpp`) which is only available
  in an installed package.

* Fixed `tinytest` assertion: `expect_message(..., pattern = NA)` does NOT
  mean "expect no messages" — it fails if no message is emitted. Replaced
  with `expect_silent()`.

* New `inst/tinytest/test_rwisecondorx.R` — 76 assertions covering reference
  building (gender model, masking, PCA, KNN indexes, null ratios), prediction
  (normalization, CBS, aberration calling), and trisomy detection sensitivity.

* New `inst/tinytest/test_cohort_pipeline.R` — 31 assertions for end-to-end
  pipeline: cohort generation → binning → reference building → prediction →
  trisomy detection. Tests that trisomy 21, 18, and 13 are detected as gains
  and euploid samples produce no aberrations.

* Total test count: 345 assertions, all passing.

## Y-unique region ratio for sex prediction

* New `nipter_y_unique_ratio()` counts reads overlapping 7 Y-chromosome
  unique gene regions (HSFY1, BPY2, BPY2B, BPY2C, XKRY, PRY, PRY2) and
  computes the ratio to total nuclear genome reads. Uses DuckDB/duckhts
  index-based region queries (`read_bam(region := ...)`) for efficient
  BAM access. The bundled GRCh37 regions file can be replaced with a custom
  file for other assemblies.

* New `nipter_sex_model_y_unique()` fits a 2-component GMM on Y-unique
  ratios (one per BAM in a cohort), producing a `NIPTeRSexModel` compatible
  with `nipter_predict_sex()`.

* `nipter_predict_sex()` gains a `y_unique_ratio` parameter for passing
  a pre-computed Y-unique ratio when a `"y_unique"` model is included in
  the consensus vote. This enables the full 3-model majority-vote sex
  prediction pipeline (Y-unique ratio + Y fraction + XY fractions).

* Bundled `inst/extdata/grch37_Y_UniqueRegions.txt` — TSV of 7 GRCh37
  Y-unique regions used for the Y-unique ratio calculation.

## Sex prediction via Gaussian mixture models

* New `nipter_sex_model()` fits a 2-component GMM on sex chromosome fractions
  from a `NIPTeRControlGroup` using `mclust::Mclust()`. Supports `"y_fraction"`
  (univariate on Y-chromosome fraction) and `"xy_fraction"` (bivariate on X + Y
  fractions) methods. The male cluster is identified as the component with
  higher median Y fraction.

* New `nipter_predict_sex()` classifies a `NIPTeRSample` as male or female
  given one or more `NIPTeRSexModel` objects. Multiple models use majority vote
  consensus (tie defaults to "female" — conservative for NIPT).

* `mclust` added to `Suggests` in DESCRIPTION.

## `nipter_bin_bam_bed()` SeparatedStrands output

* `nipter_bin_bam_bed()` gains a `separate_strands` parameter. When `TRUE`,
  outputs a 7-column BED (`chrom`, `start`, `end`, `count`, `count_fwd`,
  `count_rev`, `corrected_count`) where `count = count_fwd + count_rev`. When
  `FALSE` (default), the 5-column BED format is unchanged.

## SeparatedStrands support

* `bam_convert()` gains a `separate_strands` parameter. When `TRUE`, returns a
  list with `fwd` (forward strand) and `rev` (reverse strand) data frames.

* `nipter_bin_bam()` gains `separate_strands = TRUE` support, producing
  `NIPTeRSample` objects with class `c("NIPTeRSample", "SeparatedStrands")`.
  Autosomal chromosome reads are stored as a list of two matrices (forward and
  reverse) with rownames `"1F".."22F"` and `"1R".."22R"`.

* All NIPTeR statistical functions now support SeparatedStrands samples:
  - `nipter_gc_correct()`: LOESS/bin-weight fitted on summed strand counts,
    corrections applied independently to each strand matrix.
  - `nipter_chi_correct()`: chi-squared computed on summed strands, correction
    applied per-strand.
  - `nipter_z_score()` and `nipter_ncv_score()`: use collapsed (F+R) fractions.
  - `nipter_regression()`: doubles the predictor pool (44 candidates:
    `"1F".."22F"`, `"1R".."22R"`) with complementary exclusion within each
    model (selecting `"5F"` excludes both `"5F"` and `"5R"`).

## WisecondorX CLI wrapper fixes

* `wisecondorx_predict()`: removed hallucinated `ref_binsize` / `--binsize`
  parameter that does not exist in the upstream `predict` subcommand.

* `wisecondorx_newref()`: fixed `ref_binsize` default from `50000L` to
  `100000L` to match upstream default.

* `wisecondorx_predict()`: added `add_plot_title` parameter mapping to upstream
  `--add-plot-title`.
