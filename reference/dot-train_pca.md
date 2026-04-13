# Train PCA and return ratio-corrected data

Fits a PCA model with `n_comp` components to the masked reference data,
then corrects the data by dividing by the PCA reconstruction (ratio
correction, not subtractive). Mirrors `newref_tools.train_pca()`.

## Usage

``` r
.train_pca(ref_data, n_comp = 5L)
```

## Arguments

- ref_data:

  Numeric matrix `(n_masked_bins, n_samples)`.

- n_comp:

  Number of PCA components (default 5).

## Value

A list with:

- corrected:

  Numeric matrix `(n_masked_bins, n_samples)` — ratio- corrected data.

- components:

  Numeric matrix `(n_comp, n_masked_bins)` — PCA rotation (loadings).

- center:

  Numeric vector of length `n_masked_bins` — PCA column means.
