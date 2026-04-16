# Extract the autosomal read-count matrix from a NIPTSample

Returns a 22 x n_bins numeric matrix (rows = chromosomes "1"-"22", cols
= bin indices). For `SeparatedStrandsSample`, forward and reverse counts
are summed.

## Usage

``` r
autosomal_matrix(x, ...)
```

## Arguments

- x:

  A `NIPTSample` object.

## Value

A 22-row numeric matrix.
