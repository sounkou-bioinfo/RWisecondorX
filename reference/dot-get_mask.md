# Compute the global bin mask from a set of reference samples

Identifies bins with sufficient coverage across all reference samples.
Bins where the summed normalized coverage is less than 5% of the median
are masked out. Mirrors `newref_tools.get_mask()` in upstream
WisecondorX.

## Usage

``` r
.get_mask(samples, ref_bins_per_chr = NULL)
```

## Arguments

- samples:

  List of sample objects (each a named list of integer vectors keyed by
  chromosome `"1"`–`"24"`).

- ref_bins_per_chr:

  Optional integer vector of length 24 giving the minimum number of bins
  per chromosome. When supplied (e.g. from the global mask), each
  chromosome is zero-padded to at least this many bins so that the
  returned mask has the same length as the global mask. This avoids
  length mismatches when combining gender-specific masks with `&`.

## Value

A list with two elements:

- mask:

  Logical vector of length `sum(bins_per_chr)`; `TRUE` for bins with
  sufficient coverage.

- bins_per_chr:

  Integer vector of length 24 giving the number of bins per chromosome.
