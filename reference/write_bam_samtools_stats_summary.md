# Write a samtools-stats-compatible `SN` summary section

Writes the implemented summary-number subset from
[`bam_samtools_stats_summary()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_samtools_stats_summary.md)
in a `samtools stats`-compatible `SN` text format.

## Usage

``` r
write_bam_samtools_stats_summary(x, file)
```

## Arguments

- x:

  One-row `data.frame` returned by
  [`bam_samtools_stats_summary()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_samtools_stats_summary.md).

- file:

  Output path for the text file.

## Value

`file`, invisibly.

## Details

This writer emits the `SN` section only. It does not attempt to
reproduce the full `samtools stats` output surface such as insert-size
histograms, GC-depth, or cycle-level tables.
