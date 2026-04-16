# NIPTeR sample with combined forward+reverse read counts

NIPTeR sample with combined forward+reverse read counts

## Usage

``` r
CombinedStrandsSample(
  sample_name = character(0),
  binsize = integer(0),
  correction = NIPTCorrectionRecord(),
  auto_matrix = NULL,
  sex_matrix_ = NULL
)
```

## Arguments

- sample_name:

  Sample identifier.

- binsize:

  Positive integer bin width in base pairs.

- correction:

  A `NIPTCorrectionRecord` describing the applied corrections.

- auto_matrix:

  22 x n_bins numeric matrix for autosomes.

- sex_matrix\_:

  2 x n_bins numeric matrix for X/Y.

## Slots

- `auto_matrix`:

  22 x n_bins numeric matrix (rows = chr 1-22).

- `sex_matrix_`:

  2 x n_bins numeric matrix (row 1 = X, row 2 = Y).
