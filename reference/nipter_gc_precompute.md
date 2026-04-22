# Pre-compute and save per-bin GC content to a TSV.bgz file

Runs `rduckhts_fasta_nuc()` once and writes the GC percentage table to a
bgzipped, tabix-indexed TSV file. Pass the resulting path to
`nipter_gc_correct(gc_table = ...)` to avoid recomputing GC content for
every sample in a large cohort.

## Usage

``` r
nipter_gc_precompute(fasta, binsize = 50000L, out, con = NULL)
```

## Arguments

- fasta:

  Path to an indexed reference FASTA file (.fa/.fasta with .fai).

- binsize:

  Bin size in base pairs (default `50000L`).

- out:

  Path for the output file. The tabix index is written alongside as
  `<out>.tbi`.

- con:

  Optional open DBI connection with duckhts loaded.

## Value

`out` invisibly.

## Details

The output is a 5-column, tab-delimited TSV.bgz: `chrom`, `start`,
`end`, `pct_gc`, `seq_len`. Coordinates are 0-based half-open intervals
(BED convention). Chromosomes use no `chr` prefix (`1`–`22`, `X`, `Y`).
GC is computed on callable `A/C/G/T` sequence only, excluding
`N`/ambiguous bases from the denominator. Bins with no callable sequence
are written with `pct_gc = NA`.

## See also

[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)

## Examples

``` r
if (FALSE) { # \dontrun{
nipter_gc_precompute("hg38.fa", binsize = 50000L, out = "hg38_gc_50k.tsv.bgz")
cg <- nipter_gc_correct(cg, gc_table = "hg38_gc_50k.tsv.bgz")
} # }
```
