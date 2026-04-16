# RWisecondorX Agent Guidelines

This document provides guidance for AI agents working on `RWisecondorX`.

## Project Goal

Build an R toolkit for copy number analysis and trisomy prediction in
non-invasive prenatal testing (NIPT), comparable in scope to `NIPTUtils`
but built entirely on `Rduckhts`/`DuckDB` rather than
`Rsamtools`/Python. The package ports two upstream implementations:

- **WisecondorX** — copy number variation detection for cfDNA
  (Python/pysam upstream)
- **NIPTeR** — fast trisomy prediction via chromosomal fraction
  Z-scores, NCV scores, regression, and chi-squared correction
  (R/Rsamtools upstream)

Design priorities:

- **Exact conformance** with upstream implementations where documented
  and testable.
- **Performance** — DuckDB SQL replaces per-read Python loops and
  Rsamtools scans; large BAMs run in seconds, not minutes.
- **Interoperable file formats** — bgzipped, tabix-indexed BED files as
  the intermediate layer between binning and analysis, consumable by
  DuckDB, R, Python, or any tabix-aware tool.
- **No Python runtime dependency** — `Rduckhts` replaces `pysam` for all
  HTS operations; `reticulate` and `condathis` are optional and used
  only for conformance testing.

------------------------------------------------------------------------

## Implementation Status

### Completed

**Shared binning engine — `R/convert.R`** -
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md):
DuckDB/SQL read-counting core. Supports `binsize`, `mapq`,
`require_flags` (samtools `-f`), `exclude_flags` (samtools `-F`),
`rmdup` (`"streaming"` / `"flag"` / `"none"`), `separate_strands`
(returns list of forward/reverse data frames), CRAM via `reference`.
Achieves exact bin-for-bin conformance with WisecondorX on
HG00106.chrom11 (25,115 non-zero bins, 0 mismatches). -
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md):
bgzipped + tabix-indexed 4-column BED output (`chrom`, `start`, `end`,
`count`). Uses
[`Rduckhts::rduckhts_bgzip()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_bgzip.html)
and `rduckhts_tabix_index()` — no external tools.

**WisecondorX layer** - `R/wisecondorx_cli.R`:
[`wisecondorx_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_convert.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md),
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)
— thin `condathis` wrappers delegating to the official bioconda
package. - `R/wisecondorx_npz.R`:
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md)
— WisecondorX-compatible NPZ output via `reticulate`.

**NIPTeR binning layer — `R/nipter_bin.R`** -
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md):
produces `NIPTeRSample` objects. With `separate_strands = FALSE`
(default), class is `c("NIPTeRSample", "CombinedStrands")`; with
`separate_strands = TRUE`, class is
`c("NIPTeRSample", "SeparatedStrands")` with
`autosomal_chromosome_reads` as a list of two matrices (forward/reverse,
rownames `"1F".."22F"` / `"1R".."22R"`). Exposes `mapq`,
`require_flags`, `exclude_flags`, `rmdup` for pre-filtering matching
real-world NIPT pipelines (e.g. `mapq=40L, exclude_flags=1024L`). -
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md):
BED output with `separate_strands` support. When `FALSE` (default):
5-column BED (`chrom`, `start`, `end`, `count`, `corrected_count`). When
`TRUE`: 9-column BED (`chrom`, `start`, `end`, `count`, `count_fwd`,
`count_rev`, `corrected_count`, `corrected_fwd`, `corrected_rev`).

**Tests** - `inst/tinytest/test_fixtures.R` (41 assertions): synthetic
BAM/CRAM fixtures, all three `rmdup` modes, CRAM reference round-trip. -
`inst/tinytest/test_nipter.R` (47 assertions): `NIPTeRSample` structure,
matrix dimensions and rownames, MAPQ filter, `exclude_flags=1024` dedup,
BED.gz output shape, and bin-for-bin conformance against
[`NIPTeR::bin_bam_sample()`](https://rdrr.io/pkg/NIPTeR/man/bin_bam_sample.html)
using the bundled whole-genome fixture by default
(`NIPTER_CONFORMANCE_BAM` remains an override). -
`inst/tinytest/test_cli_args.R`: CLI argument builder unit tests. -
`inst/tinytest/test_npz.R`: NPZ round-trip tests. -
`inst/tinytest/test_integration.R`: WisecondorX bin-for-bin conformance
via condathis. - `inst/tinytest/test_sra_metadata.R`: SRA URL helpers.

**Infrastructure** - `R/aaa.R`: SRA metadata helpers. -
`scripts/make_fixtures.sh`: reproducible synthetic BAM/CRAM fixture
generation. - `DESCRIPTION`: NIPTeR authors (Dirk de Weerd, Lennart
Johansson) listed as `cph`. NIPTeR added to `Suggests`.

**NIPTeR statistical layer**

Each file is strictly separate; never mix NIPTeR and WisecondorX code.
All statistical functions support both CombinedStrands and
SeparatedStrands samples.

- `R/nipter_control.R` —
  [`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md):
  constructs a `NIPTeRControlGroup` from a list of `NIPTeRSample`
  objects with validation and dedup.
  [`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md):
  per-chromosome Z-scores and Shapiro-Wilk normality tests across the
  control group.
  [`nipter_match_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_control_group.md):
  selects the best-fit control samples for a test sample; distance
  computation is an Rcpp+OpenMP `nipter_ssd_scores_cpp()` kernel (one
  query vs N columns of the pre-extracted 22×N fractions matrix).
  [`nipter_match_matrix()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_matrix.md):
  computes the full N×N pairwise SSD matrix via
  `nipter_ssd_matrix_cpp()` — replaces the production
  `lapply(1:N, match_control_group(..., mode="report"))` hot loop with a
  single vectorized call.
  [`nipter_control_group_from_beds()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_control_group_from_beds.md):
  loads a directory of TSV.bgz files into a control group in one DuckDB
  pass via `rduckhts_tabix_multi()`. Internal helpers
  `.sample_chr_fractions()` (44-element for SeparatedStrands, 22 for
  Combined) and `.sample_chr_fractions_collapsed()` (always 22-element)
  handle strand dispatch.
- `R/nipter_gc.R` —
  [`nipter_gc_precompute()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_precompute.md):
  runs `rduckhts_fasta_nuc()` once and writes a 5-column TSV.bgz+tbi
  (`chrom`, `start`, `end`, `pct_gc`, `seq_len`); use this to pay the
  FASTA scan cost once per reference build rather than once per sample.
  [`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md):
  LOESS and bin-weight GC correction accepting a `gc_table` parameter —
  either a path to the precomputed TSV.bgz, an in-memory list, or `NULL`
  (falls back to the `fasta` path). When correcting a
  `NIPTeRControlGroup`, the GC table is resolved once and reused for
  every sample. SeparatedStrands: LOESS fitted on
  `Reduce("+", auto_list)` summed counts, corrections applied
  independently to each strand matrix via
  [`lapply()`](https://rdrr.io/r/base/lapply.html). Sex chromosome
  correction via nearest-neighbour (LOESS) or same-bucket weights
  (bin-weight).
- `R/nipter_chi.R` —
  [`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md):
  chi-squared overdispersion correction applied simultaneously to sample
  and control group. SeparatedStrands: chi computed on summed strand
  counts, correction applied per-strand via
  [`lapply()`](https://rdrr.io/r/base/lapply.html).
- `R/nipter_score.R` —
  [`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md):
  chromosomal fraction Z-score with Shapiro-Wilk normality test. Uses
  collapsed fractions for SeparatedStrands.
  [`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md):
  normalised chromosome value with brute-force denominator search using
  [`utils::combn()`](https://rdrr.io/r/utils/combn.html) (replaces
  [`sets::set_combn()`](https://rdrr.io/pkg/sets/man/set.html)
  dependency).
- `R/nipter_regression.R` —
  [`nipter_regression()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md):
  forward stepwise regression Z-score with train/test split, practical
  vs theoretical CV selection. Supports both CombinedStrands and
  SeparatedStrands. SeparatedStrands doubles the predictor pool (44
  candidates: `"1F".."22F","1R".."22R"`) with complementary exclusion
  (selecting `"5F"` excludes both `"5F"` and `"5R"` from the same
  model).
- `R/nipter_sex.R` —
  [`nipter_sex_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md):
  2-component GMM sex prediction via
  [`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html),
  supporting `"y_fraction"` and `"xy_fraction"` methods.
  [`nipter_sex_model_y_unique()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model_y_unique.md):
  2-component GMM on Y-unique region read ratios (BAM-level, uses
  [`nipter_y_unique_ratio()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_y_unique_ratio.md)).
  [`nipter_y_unique_ratio()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_y_unique_ratio.md):
  counts reads overlapping 7 Y-chromosome unique gene regions via
  DuckDB/duckhts index-based region queries and returns ratio to total
  nuclear reads.
  [`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md):
  classifies a sample as male/female using one or more models with
  majority-vote consensus; `y_unique_ratio` parameter enables the
  Y-unique model in the vote. Internal `.mclust_fit()` wrapper evaluates
  in mclust namespace to work around `mclustBIC()` not being
  namespace-qualified in upstream mclust.
- `inst/tinytest/test_nipter_stats.R` — 105 assertions covering all
  statistical functions including SeparatedStrands variants.
- `inst/tinytest/test_nipter_sex.R` — 45 assertions covering sex model
  building, Y-unique model, 3-model consensus, classification accuracy,
  prediction, and edge cases.

**Native WisecondorX implementation (newref + predict)**

Pure R/Rcpp port of the WisecondorX `newref` and `predict` pipelines,
with performance-critical KNN reference-finding in compiled C++. No
Python runtime dependency.

- `R/rwisecondorx_utils.R` — Shared utilities:
  [`.train_gender_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-train_gender_model.md)
  (2-component GMM on Y-fractions with zero-variance fallback),
  [`.gender_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-gender_correct.md),
  [`.get_mask()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-get_mask.md)
  (5% median coverage threshold),
  [`.normalize_and_mask()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize_and_mask.md),
  [`.train_pca()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-train_pca.md)
  (5 components, ratio correction),
  [`.project_pc()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-project_pc.md),
  [`.predict_gender()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-predict_gender.md).
  Also exports
  [`scale_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/scale_sample.md)
  for rescaling bin sizes.
- `R/rwisecondorx_newref.R` —
  [`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md):
  reference building pipeline. Gender model training → global bin mask →
  per-partition (A/F/M) normalize/PCA/distance-filter/KNN/null-ratios.
  [`.build_sub_reference()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-build_sub_reference.md)
  orchestrates each partition.
  [`.get_reference()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-get_reference.md)
  delegates to Rcpp.
- `R/rwisecondorx_predict.R` —
  [`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md):
  prediction pipeline. Rescale → predict gender → normalize autosomes
  (coverage + PCA + weights + optimal cutoff + 3-pass within-sample
  normalization with aberration masking) → normalize gonosomes → combine
  → inflate → log-transform → blacklist → CBS → segment Z-scores →
  aberration calling.
  [`.normalize()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize.md),
  [`.normalize_repeat()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize_repeat.md),
  [`.normalize_once()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize_once.md)
  implement the multi-pass normalization. Global-to-local index
  conversion in
  [`.normalize_once()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize_once.md)
  handles the index-space mapping correctly.
- `R/rwisecondorx_cbs.R` —
  [`.exec_cbs()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-exec_cbs.md):
  CBS wrapper; `parallel = TRUE` (now the default) uses
  `ParDNAcopy::parSegment(num.cores = cpus)` with an explicit thread
  count; falls back to
  [`DNAcopy::segment()`](https://rdrr.io/pkg/DNAcopy/man/segment.html)
  with a message if ParDNAcopy is absent. Matches upstream conventions
  (0→NA, 0 weights→1e-99, split segments at large NA gaps, recalculate
  weighted means).
  [`.get_segment_zscores()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-get_segment_zscores.md):
  segment-level Z-scores from null ratio distributions.
- `R/rwisecondorx_output.R` — BED and statistics output writers for
  `WisecondorXPrediction` objects.
- `src/knn_reference.cpp` — `knn_reference_cpp()`: OpenMP-parallelized
  KNN reference bin finding (leave-one-chromosome-out squared Euclidean
  distance, partial sort for K nearest). Stores **global** 1-based
  indexes. `null_ratios_cpp()`: OpenMP-parallelized null ratio
  computation using global indexes directly against full sample vectors.
- `src/nipter_matching.cpp` — `nipter_ssd_scores_cpp()`: one query vs N
  columns of a fractions matrix, returns N SSD scores (OpenMP).
  `nipter_ssd_matrix_cpp()`: symmetric N×N pairwise SSD matrix (OpenMP,
  symmetric schedule). Both accept `cpus` for thread count.
- `src/RcppExports.cpp`, `R/RcppExports.R` — auto-generated Rcpp glue.
- `R/zzz_rcpp.R` — `@useDynLib RWisecondorX, .registration = TRUE` and
  `@importFrom Rcpp sourceCpp` roxygen directives.
- `src/Makevars`, `src/Makevars.win` — OpenMP compilation flags.

**Synthetic cohort generator**

- `R/synthetic_cohort.R` —
  [`generate_cohort()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/generate_cohort.md):
  creates synthetic BAMs with “compressed” chromosome lengths (100bp per
  bin) for testing. Produces ~435KB BAMs per sample. Trisomy simulation
  uses the Nguyen et al. 2023 fetal fraction model
  (<doi:10.1101/2023.11.24.568620>): reads are stochastically removed
  from non-target chromosomes with probability `p = k*f/(1+k*f)` to
  simulate relative enrichment of the trisomy chromosome, where
  `k = 0.5` (non-mosaic) or `k = 0.25` (mosaic). Default fetal fraction
  is 10%. Exports `COMPRESSED_BINSIZE` constant (100L).
- `inst/scripts/make_cohort.R` — CLI wrapper for batch cohort
  generation.

**Tests — native WisecondorX and pipeline**

- `inst/tinytest/test_rwisecondorx.R` — 76 assertions: reference
  building (gender model, masking, PCA, KNN indexes and distances, null
  ratios), prediction (normalization, CBS, Z-scores, aberration
  calling), trisomy detection sensitivity (T21/T18/T13 detected as
  gains, euploid negative control clean).
- `inst/tinytest/test_cohort_pipeline.R` — 31 assertions: end-to-end
  pipeline (cohort generation → binning → newref → predict → aberration
  detection). All test files use
  [`library(RWisecondorX)`](https://sounkou-bioinfo.github.io/RWisecondorX/)
  instead of [`source()`](https://rdrr.io/r/base/source.html) hacks.

**BED.gz reader functions — round-trip from stored BED files**

- `R/bed_reader.R` —
  [`bed_to_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_sample.md):
  reads a 4-column BED.gz (from
  [`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md))
  into the named-list-of-integer-vectors format for
  [`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md)
  /
  [`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md).
  [`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md):
  reads a 5-column (CombinedStrands) or 9-column (SeparatedStrands)
  BED.gz (from
  [`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md))
  into a `NIPTeRSample` object. Auto-detects column count via
  `read_tabix()` probe; handles literal “NA” in `corrected_count` via
  `TRY_CAST`. Per-strand corrected values (`corrected_fwd`,
  `corrected_rev`) are read independently for SeparatedStrands.
- `inst/tinytest/test_bed_reader.R` — 46 assertions covering WisecondorX
  and NIPTeR round-trips, SeparatedStrands 9-column BED, corrected
  per-strand round-trip, sample name inference, and
  [`scale_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/scale_sample.md)
  integration.

### Open architectural questions

**NIPTeR layer**

- **Sex-stratified NCV for X/Y chromosomes**: The user’s clinical
  pipeline computes sex-stratified NCV denominators for X and Y
  (separate models for males vs females). Not yet implemented.
- **Sex-stratified regression for X/Y**: Forward stepwise
  [`lm()`](https://rdrr.io/r/stats/lm.html) models for X and Y
  fractions, stratified by predicted sex. Not yet implemented.
- **Multi-chromosome NIPTeR conformance fixture**: bundled as
  `inst/extdata/nipter_conformance_fixture.bam` and regenerated by
  `make fixtures` / `make nipter-fixture`. It contains all 24
  chromosomes, no unmapped reads, and no same-position collisions, so
  installed-package conformance tests do not require an external BAM.

**Native WisecondorX layer — what remains**

The pipeline (`newref` + `predict`) is functionally complete and all 76
unit tests + 31 cohort pipeline tests pass. The remaining work is
validation and extension:

- **End-to-end conformance vs upstream Python**: Bin-for-bin comparison
  of
  [`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md) +
  [`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md)
  output against the official Python `wisecondorx newref` + `predict` on
  real multi-sample BAMs. The reference object (PCA components, KNN
  indexes, null ratios) and final aberration calls both need to match.
  Requires a controlled multi-sample BAM set and the bioconda
  `wisecondorx` package via `condathis`. Not yet automated.
- **WisecondorX reference as interoperable file**: The
  `WisecondorXReference` object is currently serialized as an RDS file.
  The matrices (PCA components, mask, indexes, distances, null_ratios)
  are dense numeric arrays that don’t map naturally to a flat BED/TSV. A
  structured HDF5 or multiple-TSV.bgz serialization would allow Python
  tooling to read the reference. Not a priority unless the user has a
  Python consumer.
- **[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md)
  multi-file input via tabix_multi**: ~~Currently takes a list of
  in-memory sample lists. Could accept a bed_dir argument using
  `rduckhts_tabix_multi()` to load all WisecondorX 4-column BEDs at
  once, analogous to
  [`nipter_control_group_from_beds()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_control_group_from_beds.md).~~
  **DONE** —
  [`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md)
  now accepts `bed_dir` and `bed_pattern` parameters.
- **Threading budget propagation**:
  [`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md)
  accepts `cpus` for KNN finding.
  [`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md)
  now accepts `cpus` and passes it to `parSegment()`. Both should be
  wired together coherently when called from a pipeline script — use the
  same `cpus` value throughout.
- **Beta-mode aberration calling**:
  [`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md)
  supports `zscore` and `beta` (purity-based) cutoffs. The `beta` path
  (`ratio > 2^beta`) is implemented but not yet tested with real data.

------------------------------------------------------------------------

## Source Layout Expectations

- `R/convert.R` — shared BAM/CRAM binning engine; used by both
  WisecondorX and NIPTeR layers.
- `R/nipter_bin.R` — NIPTeR binning layer.
- `R/nipter_gc.R` — GC correction (LOESS and bin-weight).
- `R/nipter_score.R` — Z-score and NCV scoring.
- `R/nipter_chi.R` — chi-squared overdispersion correction.
- `R/nipter_control.R` — control group construction, diagnostics, and
  matching.
- `R/nipter_regression.R` — forward stepwise regression Z-score.
- `R/nipter_sex.R` — sex prediction via Gaussian mixture models
  (mclust).
- `R/rwisecondorx_utils.R` — shared utilities for native WisecondorX
  implementation.
- `R/rwisecondorx_newref.R` — native WisecondorX reference building.
- `R/rwisecondorx_predict.R` — native WisecondorX prediction pipeline.
- `R/rwisecondorx_cbs.R` — CBS segmentation wrapper
  (DNAcopy/ParDNAcopy).
- `R/rwisecondorx_output.R` — BED/statistics output generation.
- `R/bed_reader.R` — BED.gz reader functions
  ([`bed_to_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_sample.md),
  [`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md))
  for round-tripping from stored BED files.
- `R/synthetic_cohort.R` — synthetic cohort generator for testing.
- `R/simulate_trisomy_cohort.R` — donor BAM/CRAM cohort simulator for
  trisomy positives.
- `R/wisecondorx_cli.R` — CLI wrappers (condathis-based conformance
  tools).
- `R/wisecondorx_npz.R` — WisecondorX NPZ output.
- `R/aaa.R` — SRA metadata helpers.
- `R/RcppExports.R` — auto-generated Rcpp R wrappers.
- `R/zzz_rcpp.R` — roxygen directives for useDynLib + importFrom Rcpp.
- `src/knn_reference.cpp` — Rcpp + OpenMP KNN reference finding and null
  ratios.
- `src/RcppExports.cpp` — auto-generated Rcpp C++ glue.
- `src/Makevars`, `src/Makevars.win` — OpenMP compilation flags.
- `inst/tinytest/` — unit tests (one file per feature family).
- `inst/extdata/` — synthetic BAM/CRAM fixtures (including
  `nipter_conformance_fixture.bam`) and bundled reference data
  (`grch37_Y_UniqueRegions.txt`).
- `inst/scripts/make_cohort.R` — CLI wrapper for cohort generation.
- `inst/scripts/build_reference.R` — optparse CLI for building
  WisecondorX references or NIPTeR control groups from BAM/CRAM or
  BED.gz files. Supports `--mode wisecondorx|nipter`,
  `--bam-dir`/`--bam-list` (bins then builds), `--bed-dir`/`--bed-list`
  (builds from pre-binned BEDs), `--gc-table` or `--fasta` (NIPTeR GC
  correction during binning, mutually exclusive), `--bed-out-dir`, and
  all key parameters (binsize, mapq, rmdup, refsize, cpus, etc.).
- `inst/scripts/convert_sample.R` — optparse CLI for converting a single
  BAM/CRAM to a binned BED.gz or NPZ file. Supports
  `--mode wisecondorx|nipter`, `--npz` (WisecondorX NPZ output for
  Python interop), `--gc-table` or `--fasta` (NIPTeR GC correction,
  mutually exclusive), `--separate-strands`, and all binning parameters
  (binsize, mapq, rmdup, exclude-flags, require-flags, reference).
- `inst/scripts/precompute_gc.R` — optparse CLI for precomputing per-bin
  GC content from a reference FASTA. Wraps
  [`nipter_gc_precompute()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_precompute.md).
  Accepts `--fasta`, `--out`, `--binsize`. Output is a TSV.bgz+tbi that
  can be passed to `--gc-table` in the other scripts or to
  `nipter_gc_correct(gc_table = ...)`.
- The WisecondorX upstream algorithm reference is `.sync/WisecondorX/`.
- The NIPTeR upstream algorithm reference is `.sync/NIPTeR/`.
- The WisecondorX conformance script is
  `../../duckhts/scripts/wisecondorx_convert_conformance.py`.
- There is NO `inst/python/` directory.

------------------------------------------------------------------------

## Interoperable File Formats

Our “BED.gz” files are **BGZF-compressed, tabix-indexed TSV files** that
use BED coordinate conventions (0-based half-open intervals, no `chr`
prefix). They are not strict BED files — they have more than 3 columns
and the extra columns are not standard BED fields. Always read them with
`read_tabix()` (DuckDB SQL) or
[`Rduckhts::rduckhts_tabix_multi()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_tabix_multi.html),
never with [`read.table()`](https://rdrr.io/r/utils/read.table.html) or
BED-aware parsers that assume 3-column BED semantics.

- **4-column WisecondorX TSV.bgz**: `chrom`, `start`, `end`, `count`
  (written by
  [`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md)).
- **5-column NIPTeR TSV.bgz (CombinedStrands)**: `chrom`, `start`,
  `end`, `count`, `corrected_count` (written by
  [`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md)).
- **9-column NIPTeR TSV.bgz (SeparatedStrands)**: `chrom`, `start`,
  `end`, `count`, `count_fwd`, `count_rev`, `corrected_count`,
  `corrected_fwd`, `corrected_rev` (written by
  `nipter_bin_bam_bed(separate_strands = TRUE)`).
  `count = count_fwd + count_rev`; `corrected_*` columns are `NA` until
  a GC-corrected sample is supplied.
- **5-column GC table TSV.bgz**: `chrom`, `start`, `end`, `pct_gc`,
  `seq_len` (written by
  [`nipter_gc_precompute()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_precompute.md)).
  Precompute once per reference genome + bin size, pass path to
  `nipter_gc_correct(gc_table = ...)`.
- All files are bgzipped (BGZF) and tabix-indexed via
  [`Rduckhts::rduckhts_bgzip()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_bgzip.html)
  and
  [`Rduckhts::rduckhts_tabix_index()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_tabix_index.html).
  Do not use [`gzfile()`](https://rdrr.io/r/base/connections.html) or
  external tools.
- [`bed_to_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_sample.md)
  reads 4-column TSV.bgz back into the WisecondorX in-memory format.
  [`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md)
  reads 5- or 9-column TSV.bgz back into a `NIPTeRSample`.
  [`nipter_control_group_from_beds()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_control_group_from_beds.md)
  loads a whole directory via `rduckhts_tabix_multi()` in one DuckDB
  pass. All reading uses the SQL `read_tabix()` table function — never
  `read_bed()`.
- The WisecondorX reference (`WisecondorXReference`) is still RDS-only.
  Its matrices (PCA, KNN indexes, null ratios) are dense numeric arrays
  that cannot sensibly be represented in a flat TSV. Use
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html).

------------------------------------------------------------------------

## Convert Step: Exact SQL Replication of pysam larp/larp2

The WisecondorX convert step uses a streaming dedup state machine
(`larp`/`larp2`). The exact SQL equivalent uses DuckDB window functions
over `FILE_OFFSET`:

``` sql
LAG(pos) OVER (ORDER BY file_offset)                                          -- prev_pos (larp)
LAST_VALUE(CASE WHEN is_paired != 0 THEN pnext END IGNORE NULLS)              -- prev_pnext (larp2)
    OVER (ORDER BY file_offset ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING)
```

Key subtleties reproduced exactly:

- Improper pairs (`is_paired AND NOT is_proper`) are invisible to both
  larp and larp2. This filter is internal to `rmdup = "streaming"` — it
  is intrinsic to the WisecondorX algorithm, not a user flag option.
- Unpaired reads update larp but NOT larp2 — hence the `IGNORE NULLS`
  trick.
- `larp` is never reset between chromosomes — hence no `PARTITION BY`.
- Bin assignment uses integer division `pos // binsize` matching
  Python’s `int(pos / binsize)`.
- `pos` is 0-based: subtract 1 from duckhts 1-based POS before dividing.

------------------------------------------------------------------------

## Flag Filtering API

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md),
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md),
and
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md)
expose samtools-style `-f`/`-F` flag filtering:

- `require_flags` (integer bitmask, default `0L`): only reads where
  `(flag & require_flags) == require_flags` are kept. Equivalent to
  `samtools view -f`.
- `exclude_flags` (integer bitmask, default `0L`): reads where
  `(flag & exclude_flags) != 0` are dropped. Equivalent to
  `samtools view -F`.

Common real-world NIPT pre-filter pattern (matches
`samtools view --min-MQ 40 -F 1024`):

``` r
nipter_bin_bam("sample.dm.bam", mapq = 40L, exclude_flags = 1024L)
```

------------------------------------------------------------------------

## NIPTeR Layer: Known Issues and Design Decisions

### NIPTeR’s scanBam approach is “concat and pretend single-ended”

NIPTeR’s `bin_bam_sample()` calls `scanBam()` without filtering
parameters, then splits all read records by strand (`+` / `-`),
discarding only reads with strand `"*"` (unmapped). It then applies
[`unique()`](https://rdrr.io/r/base/unique.html) on positions per
chromosome per strand and sums forward and reverse counts into a single
matrix.

The consequence is that **paired-end data is treated as if you
concatenated both mates into a single-ended read pool**: R1 (forward)
and R2 (reverse) are each counted independently at their respective
start positions. A paired-end library produces roughly 2× the bin counts
of a single-end library for the same sequencing depth, because both
mates of every fragment are counted.

This creates a fundamental incompatibility: **a NIPTeR reference
(control group) built from single-end samples cannot be used to test a
paired-end sample, and vice versa.** In practice this goes unnoticed
when all samples in a cohort use the same library preparation, but
mixing single-end and paired-end data produces wildly incorrect Z-scores
and NCV values.

Our
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
replicates this behaviour (both mates counted) as the default
(`rmdup = "none"`, no flag filters), keeping it compatible with NIPTeR
control groups built the same way. Users who need fragment-level counts
(one count per fragment, not per read) should pre-filter to R1 only
using `require_flags = 0x41L` (paired + read1) or count only unique
fragments — this is a future enhancement.

### NIPTeR’s implicit positional dedup is similar to WisecondorX streaming dedup, but different

NIPTeR’s `bin_reads()` applies
[`unique()`](https://rdrr.io/r/base/unique.html) to read positions per
chromosome per strand before binning. The effect is that two reads with
the same start position on the same chromosome and strand are counted
once, not twice. This is position-based deduplication.

WisecondorX’s `larp` also collapses consecutive reads at the same
position. The critical differences:

|                          | NIPTeR [`unique()`](https://rdrr.io/r/base/unique.html) | WisecondorX larp/larp2                     |
|--------------------------|---------------------------------------------------------|--------------------------------------------|
| Scope                    | Per-chromosome, per-strand                              | Global file order, no reset between chroms |
| Pairing awareness        | None — operates on positions only                       | larp2 uses pnext for paired reads          |
| Improper pair handling   | None                                                    | Improper pairs excluded entirely           |
| Equivalent samtools flag | Not equivalent to any flag                              | Not expressible as a single flag           |

Both collapse same-position reads, but NIPTeR’s
[`unique()`](https://rdrr.io/r/base/unique.html) also silently merges
any two reads that happen to share a start position, including distinct
reads that are not PCR duplicates. On the `hg00106_chr11_fixture.bam`
fixture (1,018 paired-end reads) NIPTeR’s approach would count 1,003 vs
1,018 from `rmdup = "none"` — a difference of 15 reads, all at positions
where the duplicate flag (`0x400`) is set on one mate.

### Two additional bugs in NIPTeR’s scanBam call

**Bug 1 — Split length mismatch on BAMs with unmapped reads.**

``` r
# NIPTeR source — bin_sample.R
splitted_reads <- split(x = reads, f = droplevels(strands[strands != "*"]))
```

`strands[strands != "*"]` has length `N − U` (U = unmapped reads) but
`reads` has length `N`. R’s
[`split()`](https://rdrr.io/r/base/split.html) requires equal lengths,
so this errors or silently misbehaves on any BAM containing unmapped
reads. NIPTeR implicitly requires `samtools view -F 4` pre-filtering.

**Bug 2 — Positional dedup via
[`unique()`](https://rdrr.io/r/base/unique.html) is flag-ignorant.**

``` r
# NIPTeR source — bin_sample.R / bin_reads()
reads <- sapply(X = unique(reads_data_frame[min_read:max_read]),
                FUN = getbins, bin_size = bin_size)
bins <- tabulate(reads, nbins = n_bins)
```

[`tabulate()`](https://rdrr.io/r/base/tabulate.html) here counts unique
*positions per bin*, not reads. A PCR duplicate at a unique position
survives; two distinct reads at the same position are silently merged.
This is not equivalent to `-F 1024`.

**Consequence for conformance testing.** Exact bin-for-bin match against
[`NIPTeR::bin_bam_sample()`](https://rdrr.io/pkg/NIPTeR/man/bin_bam_sample.html)
requires a BAM with no unmapped reads and no two reads sharing a start
position per strand. The package now bundles such a fixture as
`inst/extdata/nipter_conformance_fixture.bam`; `NIPTER_CONFORMANCE_BAM`
is only an override for a custom fixture with the same constraints. Our
`rmdup = "none"` is the closest mode but counts each read independently
rather than unique positions.

------------------------------------------------------------------------

## Agent Working Instructions

Always read the existing implementation before changing it:

1.  Read the relevant R files, package metadata, and tests before
    proposing changes.
2.  When upstream WisecondorX behaviour is unclear, consult
    `../../duckhts/.sync/WisecondorX/`. Do NOT look for `inst/python/` —
    it has been removed.
3.  When upstream NIPTeR behaviour is unclear, consult
    `../../duckhts/.sync/NIPTeR/`.
4.  For conformance testing, use the official `wisecondorx` bioconda
    package via `condathis`. Do NOT vendor Python source.
5.  Preserve CRAN-friendly package behaviour. Do not make package
    loading depend on a preconfigured external Python installation.
6.  Keep package runtime logic and conformance tooling separate.
7.  Keep WisecondorX and NIPTeR layers in **separate files**
    (`nipter_*.R` vs `wisecondorx_*.R`). The shared engine is
    `convert.R`.

------------------------------------------------------------------------

## Documentation Conventions

- `README.Rmd` is the primary editable README source. Run `make readme`
  to regenerate `README.md`.
- R documentation and `NAMESPACE` are generated from roxygen2 tags. Run
  `make rd` after changing `@export`, `@param`, or `@seealso` tags.
- Do not manually edit `NAMESPACE` or `README.md`.

## Build and Workflow Rules

- Use make targets: `make rd`, `make test`, `make readme`,
  `make fixtures`.
- Keep build steps deterministic and non-interactive.

## Conformance Testing Rules

- Conformance tests skip cleanly when `condathis`, `reticulate`, or
  required external tools/packages are unavailable.
- `WISECONDORX_TEST_BAM` — a human BAM for WisecondorX integration
  tests.
- `NIPTER_CONFORMANCE_BAM` — optional override for the bundled
  coordinate-sorted, pre-filtered (no unmapped reads, no same-position
  strand collisions) whole-genome BAM used for NIPTeR conformance.

## Testing Rules

- Use `tinytest`. One file per feature family under `inst/tinytest/`.
- Update tests when R-facing behaviour changes.
- Tests requiring Python or `reticulate` must skip cleanly when
  unavailable.

## Copyright and Attribution

- NIPTeR authors Dirk de Weerd and Lennart Johansson are listed as `cph`
  in `DESCRIPTION`.
- WisecondorX authors are listed as `ctb`.
- Any new upstream algorithm ported must have its original authors
  credited with an appropriate role.

------------------------------------------------------------------------

## KNN Index Semantics — Critical Design Note

Our Rcpp implementation (`knn_reference_cpp()`) stores **global**
1-based indexes into the full masked-bin array. Upstream Python
(`get_ref_for_bins()`) stores **local** indexes relative to `chr_data`
(the leave-one-chromosome-out subset).

The consequence:

| Context                | Upstream Python                                               | Our R/Rcpp                                                     |
|------------------------|---------------------------------------------------------------|----------------------------------------------------------------|
| KNN index storage      | Local (into `chr_data`)                                       | Global (into full array)                                       |
| Predict normalization  | Correct: rebuilds `chr_data`, uses local indexes              | Correct: rebuilds `chr_data`, converts global→local via offset |
| Null ratio computation | **Bug**: uses local indexes as global (no `chr_data` rebuild) | Correct: uses global indexes against full sample vector        |

The global-to-local conversion in
[`.normalize_once()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize_once.md)
is: - Global index `g < chr_start` → local `g` (before excised
chromosome, unchanged) - Global index `g > chr_cum` → local `g - n_chr`
(after excised chromosome, shifted down) - Global index within
`[chr_start, chr_cum]` should never occur (same-chromosome bins excluded
during KNN building)

Do not change the index storage convention without updating both
[`.normalize_once()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize_once.md)
in predict and `null_ratios_cpp()` in newref.
