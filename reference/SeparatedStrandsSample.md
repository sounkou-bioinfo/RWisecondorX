# NIPTeR sample with separate forward and reverse read counts

NIPTeR sample with separate forward and reverse read counts

## Usage

``` r
SeparatedStrandsSample(
  sample_name = character(0),
  binsize = integer(0),
  correction = NIPTCorrectionRecord(),
  auto_fwd = NULL,
  auto_rev = NULL,
  sex_fwd = NULL,
  sex_rev = NULL
)
```

## Slots

- `auto_fwd`:

  22 x n_bins numeric matrix (forward strand).

- `auto_rev`:

  22 x n_bins numeric matrix (reverse strand).

- `sex_fwd`:

  2 x n_bins numeric matrix (forward strand; row1=X, row2=Y).

- `sex_rev`:

  2 x n_bins numeric matrix (reverse strand).
