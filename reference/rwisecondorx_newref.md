# Build a WisecondorX reference from binned samples

Native R implementation of the WisecondorX `newref` pipeline. Takes a
list of binned samples (as returned by
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
or loaded from NPZ via reticulate) and builds a PCA-based reference
suitable for
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md).

## Usage

``` r
rwisecondorx_newref(
  samples = NULL,
  binsize = 100000L,
  sample_binsizes = NULL,
  nipt = FALSE,
  refsize = 300L,
  yfrac = NULL,
  gender_model_names = "V",
  female_partition_min_samples = 5L,
  male_partition_min_samples = 5L,
  mask_min_median_coverage_fraction = 0.05,
  gender_model_grid_min = 0,
  gender_model_grid_max = 0.02,
  gender_model_grid_length = 5000L,
  pca_components = 5L,
  pca_distance_min_class_bins = 10L,
  pca_distance_autosome_mad_multiplier = 20,
  pca_distance_autosome_floor = 10,
  pca_distance_chrX_mad_multiplier = 20,
  pca_distance_chrX_floor = 10,
  pca_distance_chrY_mad_multiplier = 50,
  pca_distance_chrY_floor = 15,
  null_ratio_max_samples = 100L,
  cpus = 4L,
  bed_dir = NULL,
  bed_pattern = "*.bed.gz",
  con = NULL
)
```

## Arguments

- samples:

  List of sample objects, each a named list of integer vectors keyed by
  chromosome (`"1"`–`"24"`), as returned by
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).
  At least 10 samples are required. Mutually exclusive with `bed_dir`.

- binsize:

  Integer; the target bin size in base pairs. All samples are rescaled
  to this size. Default `100000L`. Inferred from BED files when
  `bed_dir` is supplied and `sample_binsizes` is `NULL`.

- sample_binsizes:

  Optional integer vector of per-sample bin sizes. If `NULL` (default),
  all samples are assumed to already be at `binsize`.

- nipt:

  Logical; if `TRUE`, NIPT mode (no gender correction, no male gonosomal
  reference). Default `FALSE`.

- refsize:

  Integer; number of reference bin locations per target bin. Default
  `300L`.

- yfrac:

  Optional numeric; manual Y-fraction cutoff for gender classification.
  If `NULL` (default), the cutoff is derived from a GMM.

- gender_model_names:

  Character scalar; mclust covariance model name used for the Y-fraction
  Gaussian mixture during gender-model fitting. Default `"V"`, matching
  the current upstream WisecondorX implementation.

- female_partition_min_samples:

  Integer; minimum female samples required before applying
  female-specific masking and building the female gonosomal partition.
  In NIPT mode this is also the minimum female count required to
  proceed. Default `5L`.

- male_partition_min_samples:

  Integer; minimum male samples required before applying male-specific
  masking and building the male gonosomal partition in non-NIPT mode.
  Default `5L`.

- mask_min_median_coverage_fraction:

  Numeric scalar; bins with summed normalized coverage below this
  fraction of the global median are masked before model fitting. Default
  `0.05`.

- gender_model_grid_min:

  Numeric scalar; lower bound of the density grid used to locate the
  Y-fraction cutoff when `yfrac` is not provided. Default `0`.

- gender_model_grid_max:

  Numeric scalar; upper bound of that density grid. Default `0.02`.

- gender_model_grid_length:

  Integer; number of grid points used to locate the Y-fraction cutoff.
  Default `5000L`.

- pca_components:

  Integer; PCA components retained during between-sample normalization.
  Default `5L`.

- pca_distance_min_class_bins:

  Integer; minimum masked bins required before applying class-specific
  PCA-distance pruning. Default `10L`.

- pca_distance_autosome_mad_multiplier:

  Numeric; MAD multiplier used for autosomal PCA-distance pruning.
  Default `20`.

- pca_distance_autosome_floor:

  Numeric; minimum autosomal PCA-distance cutoff. Default `10`.

- pca_distance_chrX_mad_multiplier:

  Numeric; MAD multiplier used for chrX PCA-distance pruning. Default
  `20`.

- pca_distance_chrX_floor:

  Numeric; minimum chrX PCA-distance cutoff. Default `10`.

- pca_distance_chrY_mad_multiplier:

  Numeric; MAD multiplier used for chrY PCA-distance pruning. Default
  `50`.

- pca_distance_chrY_floor:

  Numeric; minimum chrY PCA-distance cutoff. Default `15`.

- null_ratio_max_samples:

  Integer; maximum number of cohort samples used when sampling columns
  for null-ratio estimation. Default `100L`, matching the current
  upstream implementation.

- cpus:

  Integer; number of threads for reference bin finding. Default `4L`.

- bed_dir:

  Optional character; path to a directory of 4-column bgzipped BED files
  (as written by
  [`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md)).
  All files matching `bed_pattern` are loaded in a single DuckDB pass
  via `rduckhts_tabix_multi()`. Mutually exclusive with `samples`.

- bed_pattern:

  Glob pattern for matching BED files inside `bed_dir`. Default
  `"*.bed.gz"`.

- con:

  Optional existing DuckDB connection. Used only when `bed_dir` is
  supplied. A temporary connection is created (and closed) if `NULL`.

## Value

A list-like `WisecondorXReference` S7 object containing autosomal and
(optionally) gonosomal sub-references. See Details for the full
structure.

## Details

The pipeline trains a gender model (2-component GMM on Y-fractions),
optionally applies gender correction for non-NIPT workflows, computes a
global bin mask, then builds three sub-references: autosomes (A), female
gonosomes (F), and male gonosomes (M). Each sub-reference includes PCA
components, within-sample reference bin indices and distances, and null
ratios for between-sample Z-scoring.

This is a faithful port of the upstream Python `wisecondorx newref`,
crediting the original WisecondorX authors.

The returned reference object contains:

- binsize:

  Integer; the reference bin size.

- is_nipt:

  Logical; whether NIPT mode was used.

- trained_cutoff:

  Numeric; Y-fraction gender cutoff.

- has_female:

  Logical; whether a female gonosomal reference exists.

- has_male:

  Logical; whether a male gonosomal reference exists.

- algorithm_params:

  Named list of native reference-building parameters used to create the
  object.

- mask, bins_per_chr, masked_bins_per_chr, masked_bins_per_chr_cum,
  pca_components, pca_mean, indexes, distances, null_ratios:

  Autosomal reference components.

- mask.F, ..., null_ratios.F:

  Female gonosomal reference (if present).

- mask.M, ..., null_ratios.M:

  Male gonosomal reference (if present).

## See also

[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md),
[`scale_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/scale_sample.md),
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
