# Correction record for a NIPTeR sample

Tracks the sequence of corrections applied to the autosomal and sex
chromosome matrices of a `NIPTSample` object.

## Usage

``` r
NIPTCorrectionRecord(
  autosomal = .nipt_raw_correction_path(),
  sex = .nipt_raw_correction_path()
)
```

## Arguments

- autosomal:

  A `NIPTCorrectionPath` describing the autosomal correction history.
  Defaults to a single `"raw"` step.

- sex:

  A `NIPTCorrectionPath` describing the sex-chromosome correction
  history. Defaults to a single `"raw"` step.
