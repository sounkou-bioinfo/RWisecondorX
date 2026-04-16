# Extract the sex-chromosome read-count matrix from a NIPTSample

Returns a 2 x n_bins numeric matrix (row 1 = X, row 2 = Y). For
`SeparatedStrandsSample`, forward and reverse counts are summed.

## Usage

``` r
sex_matrix(x, ...)
```

## Arguments

- x:

  A `NIPTSample` object.

- ...:

  Reserved for S7 method dispatch; currently unused.

## Value

A 2-row numeric matrix.
