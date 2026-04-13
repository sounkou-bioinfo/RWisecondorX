# Run the full normalization chain for one gender partition

Run the full normalization chain for one gender partition

## Usage

``` r
.normalize(sample, ref, ref_gender, maskrepeats)
```

## Arguments

- sample:

  Named list of bin counts.

- ref:

  Reference list.

- ref_gender:

  `"A"`, `"F"`, or `"M"`.

- maskrepeats:

  Number of optimal cutoff iterations.

## Value

List with results_r, results_z, results_w, ref_sizes, m_lr, m_z.
