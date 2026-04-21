# RWisecondorX 0.0.0.9001 (development version)

## Contract and workflow hardening

* Native RWisecondorX paths now fail more explicitly on invalid analytical
  states instead of silently degrading. This includes zero/non-finite sample
  totals during gender prediction, zero/non-finite reference coverage during
  masking/normalization, NA-gap CBS split mismatches, and inconsistent
  gonosomal remapping/null-ratio contracts during prediction.

* `WisecondorXReference` / `WisecondorXPrediction` validation is stricter.
  Reference branch metadata now checks integer/non-negative bin counts,
  cumulative masked-bin consistency, and mask cardinality. Prediction objects
  now enforce finite non-negative per-chromosome bin counts aligned with the
  returned chromosome result lists.

* `rwisecondorx_ref_qc()` now returns richer readiness metadata in addition to
  the existing PASS/WARN/FAIL verdict: per-branch summaries, explicit missing
  mean-distance counts, branch-level prediction readiness flags, and a compact
  top-level readiness summary.

## QC and cohort reporting

* Preprocess and prediction cohort scripts now emit more machine-readable QC
  surfaces. Sample-level outputs include metric readiness/missingness flags,
  GC loess valid/invalid bin summaries, optional GC-curve data TSVs, and new
  `sample_readiness.tsv`, `qc_flag_summary.tsv`, and
  `prediction_readiness.tsv` artifacts alongside the existing summary tables.

* `sample_metrics` and NIPTeR QC/report helpers now propagate explicit
  degraded-state signals instead of burying them in sparse side fields. This
  includes missing metric counts, readiness flags, control-group outlier
  counts, sex-model readiness counts, and clearer plotting subtitles for
  sample QC views.

## Script-path and test determinism

* Shipped scripts `convert_sample.R`, `build_reference.R`, and
  `precompute_gc.R` now use `this.path` for script-location discovery instead
  of relying on working-directory assumptions.

* Native RWisecondorX tinytests were tightened to match the stricter contract
  behavior. The full synthetic reference-build test now uses an explicit
  manual sex cutoff for deterministic execution, while coverage/gender edge
  cases now assert explicit failures rather than permissive fallback behavior.

* The package and script-level reference-building CLIs are now aligned more
  strictly. `build_reference.R` now defaults WisecondorX references to a
  100kb target binsize, keeps NIPTeR on its 50kb binsize path, and errors on
  NIPTeR BED-mode flags (`--ref-binsize`, `--gc-table`, `--fasta`,
  `--separate-strands`) that were previously accepted but ignored.

* New `nipter_control_group_qc()` returns a typed `NIPTControlGroupQC`
  object for control-cohort diagnostics. It exposes per-chromosome
  mean/SD/CV/Shapiro summaries, per-sample matching/outlier summaries,
  sex-stratified X/Y readiness metrics, and optional bin-level chi-profile
  summaries including expected counts, scaled mean/SD/CV, chi scores, and
  correction factors.

* The chi-correction implementation now reuses an explicit internal
  chi-profile builder instead of recomputing and discarding bin-level
  statistics. This keeps the correction path and the new QC/report path on
  the same calculation.

* Gaunosome model building and scoring now fail explicitly on invalid
  reference subsets and missing inputs. Sex-specific NCV/regression builders
  error when the filtered female or male subset is too small, NCV scoring now
  checks that the requested numerator and denominator chromosome columns are
  actually present, regression sex scores are centred on the learned
  `mean_ratio`, and tied sex-model votes now warn instead of silently picking
  a class.

* The NIPTeR S7 layer has been tightened up around the control-group
  abstraction. `NIPTControlGroup` objects now carry optional explicit
  `sample_sex` annotations plus a `sex_source`, and the public
  `nipter_as_control_group()` / `nipter_control_group_from_beds()` builders
  accept those labels directly.

* New `nipter_reference_frame()` builds the chromosome-level training table
  (`NChrReads_*` / `FrChrReads_*`) directly from a control group, with an
  optional `SampleSex` column. This is the package-level foundation for
  future gaunosome model building instead of pushing production-specific data
  munging into the core sample classes.

* The NIPTeR modelling layer now has typed S7 reference/model objects:
  `NIPTReferenceFrame`, `NIPTeRSexModel`, `NIPTeRSexPrediction`, and
  `NIPTReferenceModel`. New `nipter_build_reference()` packages a control
  group, reference frame, and built sex models into one validated reference
  object, and `nipter_predict_sex()` now accepts that reference object
  directly.

* `nipter_build_reference()` now enriches its embedded `NIPTReferenceFrame`
  with production-style sex-scoring columns (`ConsensusGender`, `RR_X`,
  `RR_Y`, sex-class leave-one-out z-scores, and `IsRefSexOutlier`) instead of
  leaving gaunosome bookkeeping to downstream application code. New
  `NIPTSexScore` / `nipter_sex_score()` score X/Y directly against the
  sex-matched non-outlier reference subset owned by that typed reference
  object.

* New `nipter_build_sex_ncv_models()` /
  `nipter_build_sex_regression_models()` attach typed gaunosome NCV and
  regression models directly to `NIPTReferenceModel`, and
  `nipter_ncv_sex_score()` / `nipter_regression_sex_score()` now score X/Y
  from that package-owned reference state instead of requiring loose
  production-script data frames and model lists.

* New `nipter_build_gaunosome_models()` and `nipter_gaunosome_score()`
  provide the package-level gaunosome API on top of `NIPTReferenceModel`:
  one precomputation step for X/Y NCV and regression models, and one typed
  aggregate score object bundling sex z-scores, NCV scores, regression
  scores, and a compact per-chromosome summary table.

* New `nipter_gaunosome_report()` and `write_nipter_gaunosome_output()`
  provide the cohort/reporting layer on top of the typed gaunosome API:
  batch sample scoring now returns one validated report object with a
  flattened summary table, and that summary can be written directly as a
  tab-separated output file without rebuilding downstream script logic.

* The S7 implementation now exposes compatibility accessors for the legacy
  `autosomal_chromosome_reads`, `sex_chromosome_reads`, `sample_name`, and
  correction-status fields, and the control-group fraction cache is now a
  real environment-backed cache instead of a dead placeholder.

* The NIPTeR S7 compatibility layer no longer registers bare `$` methods for
  upstream `NIPTeR` class names. This fixes a recursive dispatch bug that
  could overflow the C stack during cross-package conformance tests, while
  keeping compatibility for `RWisecondorX`'s own namespaced S7 objects.

* Raw NIPTeR bin-count samples now preserve integer storage for uncorrected
  autosomal and sex matrices. Only corrected matrices are promoted to double,
  which restores the original count semantics in round-trip and binning tests.

* WisecondorX objects now use a real S7 hierarchy: `WisecondorXSample`,
  `WisecondorXReference`, `WisecondorXPrediction`, and
  `WisecondorXReferenceQC` are validated S7 classes that intentionally
  subclass `list`, so existing list-style code keeps working while the object
  boundary is now typed and explicit.

* New `simulate_trisomy_bam()` creates a BAM-level trisomy simulation from a
  donor BAM/CRAM by thinning non-target chromosomes with an htslib-backed
  native implementation. It preserves the donor read layout and sort order,
  writes a real BAM, and indexes it with `Rduckhts`. New
  `inst/scripts/simulate_trisomy_bam.R` exposes the same primitive as a CLI.

* New `simulate_trisomy_cohort()` batches that donor-BAM simulator across a
  donor set and simulation grid, writes a cohort manifest, and reuses one
  `Rduckhts` indexing session for the whole run. New
  `inst/scripts/simulate_trisomy_cohort.R` exposes the cohort generator as a
  CLI.

* New `rwisecondorx_ref_qc()` mirrors the upstream Python `ref_qc.py`
  heuristics for native `WisecondorXReference` objects. It returns a
  structured PASS/WARN/FAIL report in R and can optionally write the report as
  JSON for later inspection, without requiring the Python `wisecondorx`
  package.

* `wisecondorx_newref()` now defaults to `cpus = 4L`, matching the native
  `rwisecondorx_newref()` implementation and the package CLI defaults. This
  removes an unnecessary single-thread slow path in the public Python wrapper.

* The `wisecondorx_*()` CLI wrappers now isolate `condathis`/libmamba cache
  state from the shared `~/.cache/mamba` tree by setting a session-local
  `XDG_CACHE_HOME` during environment creation and command execution. This
  avoids stale `proc.lock` failures from unrelated mamba sessions.

* The Python `wisecondorx_*()` wrappers now cache successful `condathis`
  environment setup for the rest of the R session and force `condathis` to run
  from a known-existing working directory. This avoids redundant `create_env()`
  calls and fixes the libmamba `cannot get current path` failure when a prior
  subprocess or test left R with an invalid or deleted current directory.

* Python conformance tests no longer convert live `wisecondorx` runtime
  failures into skips. Missing optional dependencies (`condathis`,
  `reticulate`, `numpy`) still skip, but environment creation, upstream
  `wisecondorx convert/newref/predict`, and Python-side NPZ loading now fail
  loudly with contextual errors.

* `bam_convert()` now keeps the native dense-bin contract from `Rduckhts`:
  chromosomes present in the BAM header come back as trailing-zero-padded
  vectors across the full header span, while `bam_convert_npz()` adds only the
  extra upstream WisecondorX `int(length / binsize + 1)` NPZ padding quirk at
  serialisation time. The NPZ tests now assert exact padded lengths instead of
  comparing only shared prefixes.

* The bundled `nipter_conformance_fixture.bam` is now treated as a real package
  invariant rather than an optional convenience. The generator validates that
  every chromosome 1-22/X/Y has both forward and reverse reads, fixture tests
  assert that property. Live cross-package checks against the installed
  `NIPTeR` package now skip with an explicit message when upstream
  `NIPTeR::bin_bam_sample()` is incompatible with the current
  `Rsamtools`/Bioconductor stack; the always-on inline formula and fixture
  tests remain strict.

## `bam_convert()` now delegates to native `bam_bin_counts(...)`

* `bam_convert()` and the downstream WisecondorX/NIPTeR binning layer now use
  the native `bam_bin_counts(...)` kernel bundled in `Rduckhts` instead of the
  previous SQL/window-function implementation over `read_bam()` and
  `FILE_OFFSET`.

* Public behavior is preserved: `rmdup = "streaming"` still matches the
  upstream WisecondorX `larp` / `larp2` semantics, while `rmdup = "flag"` and
  `rmdup = "none"` continue to expose duplicate-flag and no-dedup modes.

* This validates the new native counting primitive in real downstream code and
  removes RWisecondorX's direct dependency on the old SQL/file-order dedup
  workaround.

* `bam_convert_npz()` now accepts the same native filter knobs as
  `bam_convert()` and `bam_convert_bed()` (`mapq`, `require_flags`,
  `exclude_flags`). `inst/scripts/convert_sample.R` now forwards those filters
  in NPZ mode instead of silently ignoring them.

* `nipter_bin_bam_bed()` is now a strict bin-and-write convenience wrapper.
  GC-corrected BED export goes through `nipter_sample_to_bed()` instead of an
  awkward `corrected =` side path. This removes the double-binning trap from
  `inst/scripts/convert_sample.R` and `inst/scripts/build_reference.R`.

* The CLI scripts under `inst/scripts/` no longer depend on `%||%` before the
  package namespace is loaded. `make_cohort.R` also uses a simpler, more
  robust dev-tree fallback path.

* The package now bundles a whole-genome `nipter_conformance_fixture.bam` and
  installed-package NIPTeR conformance tests use it by default. The
  `NIPTER_CONFORMANCE_BAM` environment variable is now only an override for a
  custom pre-filtered BAM, not a hard gate.

## SeparatedStrands BED expanded to 9 columns; BED readers switched to `read_tabix()`

* `nipter_bin_bam_bed(separate_strands = TRUE)` now writes **9-column** BED
  files: `chrom`, `start`, `end`, `count`, `count_fwd`, `count_rev`,
  `corrected_count`, `corrected_fwd`, `corrected_rev`. The previous 7-column
  format lacked per-strand corrected values, so GC-corrected SeparatedStrands
  samples lost their per-strand correction information on the BED round-trip.

* `bed_to_sample()` and `bed_to_nipter_sample()` now use `read_tabix()`
  instead of `read_bed()` for reading BED.gz files. `read_tabix()` returns
  all columns as VARCHAR with no BED-schema type coercion, avoiding the
  problem where `read_bed()` maps columns 7-8 to `thick_start`/`thick_end`
  (INTEGER) and rejects double-precision corrected values.

* `bed_to_nipter_sample()` format detection now uses `tryCatch()` around the
  `column5` probe query, so it correctly falls back to CombinedStrands when
  the file has only 5 columns (where `column5` does not exist in `read_tabix`).

* New tests for corrected per-strand round-trip: CombinedStrands corrected
  values and SeparatedStrands per-strand corrected values (with independent
  forward/reverse multipliers) survive the write→read cycle within
  floating-point tolerance. Independence of forward and reverse corrections
  is explicitly verified.

* Fixed R CMD check NOTEs: `utils::globalVariables()` for mclust symbols in
  `R/aaa.R`; `utils::write.table()` qualified in `R/convert.R` and
  `R/nipter_bin.R`.

* Total test count: 391 assertions, all passing.

## BED.gz reader functions — close the round-trip

* New `bed_to_sample()` reads a 4-column BED.gz file (written by
  `bam_convert_bed()`) back into the named-list-of-integer-vectors format
  expected by `rwisecondorx_newref()` and `rwisecondorx_predict()`. This
  closes the WisecondorX round-trip: bin once with `bam_convert_bed()`, store
  the BED.gz, and reload for analysis without re-reading the BAM.

* New `bed_to_nipter_sample()` reads a 5-column (CombinedStrands) or 9-column
  (SeparatedStrands) BED.gz file (written by `nipter_bin_bam_bed()`) into a
  `NIPTeRSample` object compatible with all NIPTeR statistical functions. Column
  count is auto-detected. Handles literal `"NA"` strings in the `corrected_count`
  field via `TRY_CAST`. Sample name is inferred from the filename or set
  explicitly.

* Both functions accept an optional DuckDB connection for reuse across multiple
  files, creating one internally (with `allow_unsigned_extensions = "true"`)
  when none is supplied.

* New `inst/tinytest/test_bed_reader.R` — 46 assertions covering WisecondorX
  4-column round-trip, NIPTeR CombinedStrands 5-column round-trip,
  SeparatedStrands 9-column round-trip (all four matrices), corrected per-strand
  round-trip, sample name inference, and integration with `scale_sample()`.

* Total test count: 391 assertions (46 in test_bed_reader.R), all passing.

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
  outputs a 9-column BED (`chrom`, `start`, `end`, `count`, `count_fwd`,
  `count_rev`, `corrected_count`, `corrected_fwd`, `corrected_rev`) where
  `count = count_fwd + count_rev`. When `FALSE` (default), the 5-column BED
  format is unchanged.

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
