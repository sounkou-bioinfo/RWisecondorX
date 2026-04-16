# Write gaunosome summary output

Writes the flattened gaunosome summary table to a tab-separated text
file. Accepts either a single `NIPTGaunosomeScore` or a batch
`NIPTGaunosomeReport`; single scores are written as one-sample reports.

## Usage

``` r
write_nipter_gaunosome_output(x, outprefix)
```

## Arguments

- x:

  A `NIPTGaunosomeScore` or `NIPTGaunosomeReport`.

- outprefix:

  Path prefix for the output file.

## Value

The written summary-table path, invisibly.
