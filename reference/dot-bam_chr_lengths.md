# Query BAM header for chromosome lengths, keyed "1"-"24"

Returns a named list: keys are "1"-"24" (matching bam_convert()
convention), values are integer chromosome lengths. Only chromosomes
1-22/X/Y present in the BAM header are included; contigs, decoys, etc.
are ignored.

## Usage

``` r
.bam_chr_lengths(con, bam)
```

## Arguments

- con:

  Open DBI connection with duckhts loaded.

- bam:

  Path to BAM/CRAM.

## Value

Named list of integer lengths.
