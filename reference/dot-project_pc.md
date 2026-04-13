# Project a single sample through stored PCA and apply ratio correction

Mirrors `predict_tools.project_pc()` in upstream WisecondorX.

## Usage

``` r
.project_pc(sample_data, pca_components, pca_mean)
```

## Arguments

- sample_data:

  Numeric vector of length `n_masked_bins` (coverage- normalized and
  masked).

- pca_components:

  Matrix `(n_comp, n_masked_bins)`.

- pca_mean:

  Numeric vector of length `n_masked_bins`.

## Value

Numeric vector of length `n_masked_bins` — ratio-corrected.
