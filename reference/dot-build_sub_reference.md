# Build a sub-reference for one gender partition

Build a sub-reference for one gender partition

## Usage

``` r
.build_sub_reference(samples, gender, total_mask, bins_per_chr, refsize, cpus)
```

## Arguments

- samples:

  List of sample objects.

- gender:

  `"A"` (autosomes), `"F"` (female gonosomes), `"M"` (male gonosomes).

- total_mask:

  Logical mask over all 24 chromosomes.

- bins_per_chr:

  Integer vector of length 24.

- refsize:

  Number of reference bins per target.

- cpus:

  Number of threads.

## Value

Named list with mask, bins_per_chr, masked_bins_per_chr, etc.
