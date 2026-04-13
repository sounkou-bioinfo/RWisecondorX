# Stack samples into a matrix for a subset of chromosomes and apply a mask

Stack samples into a matrix for a subset of chromosomes and apply a mask

## Usage

``` r
.normalize_and_mask(samples, chr_range, mask)
```

## Arguments

- samples:

  List of sample objects.

- chr_range:

  Integer vector of chromosome indices (1-based).

- mask:

  Logical vector covering all bins in `chr_range`.

## Value

Numeric matrix of shape `(n_masked_bins, n_samples)`.
