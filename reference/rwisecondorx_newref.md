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

A list (the reference object) with class `"WisecondorXReference"`,
containing autosomal and (optionally) gonosomal sub-references. See
Details for the full structure.

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
