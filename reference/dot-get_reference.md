# Find within-sample reference bins (KNN)

For each target bin, finds the `refsize` most similar bins from other
chromosomes, measured by Euclidean distance in PCA-corrected space. Also
computes null ratios for between-sample Z-scoring.

## Usage

``` r
.get_reference(
  pca_corrected,
  masked_bins_per_chr,
  masked_bins_per_chr_cum,
  refsize,
  gender,
  cpus
)
```

## Arguments

- pca_corrected:

  Numeric matrix `(n_masked_bins, n_samples)`.

- masked_bins_per_chr:

  Integer vector.

- masked_bins_per_chr_cum:

  Integer vector (cumulative).

- refsize:

  Number of reference bins per target.

- gender:

  `"A"`, `"F"`, or `"M"`.

- cpus:

  Number of threads for parallel computation. Default `1L`.

## Value

List with `indexes` (integer matrix), `distances` (numeric matrix),
`null_ratios` (numeric matrix).

## Details

Mirrors `newref_tools.get_reference()` and `get_ref_for_bins()`. Uses
Rcpp + OpenMP for performance-critical distance computation and null
ratio calculation.
