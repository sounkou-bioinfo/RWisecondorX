# Write NIPTeR-style bin counts to a bgzipped BED file

Bins a BAM/CRAM file with NIPTeR defaults and writes the result to a
bgzipped, tabix-indexed BED file with five columns: `chrom`, `start`,
`end`, `count`, `corrected_count`.

## Usage

``` r
nipter_bin_bam_bed(
  bam,
  bed,
  binsize = 50000L,
  mapq = 0L,
  require_flags = 0L,
  exclude_flags = 0L,
  rmdup = c("none", "flag"),
  corrected = NULL,
  con = NULL,
  reference = NULL,
  index = TRUE
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- bed:

  Path for the output `.bed.gz` file. The tabix index is written
  alongside as `paste0(bed, ".tbi")`.

- binsize:

  Bin size in base pairs (default 50000).

- mapq:

  Minimum mapping quality (default `0L`). Set to `40L` to match common
  NIPT pre-filtering (`samtools view --min-MQ 40`).

- require_flags:

  Integer bitmask; only reads with all bits set are kept (samtools
  `-f`). Default `0L`.

- exclude_flags:

  Integer bitmask; reads with any bit set are dropped (samtools `-F`).
  Default `0L`. Set to `1024L` to exclude duplicate-flagged reads
  (`samtools view -F 1024`).

- rmdup:

  Duplicate removal strategy: `"none"` (default) or `"flag"`.

- corrected:

  Optional `NIPTeRSample` already processed by
  [`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md).
  When supplied, its corrected counts populate the `corrected_count`
  column; otherwise the column is `NA`.

- con:

  Optional open DBI connection with duckhts already loaded.

- reference:

  Optional FASTA reference path for CRAM inputs.

- index:

  Logical; create a tabix index (default `TRUE`).

## Value

`bed` (invisibly).

## Details

`corrected_count` is filled with `NA` until
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)
is available. Once GC correction is ported, the typical workflow will
be:
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
→
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)
→ `nipter_bin_bam_bed()` with the corrected sample passed as
`corrected`.

## See also

[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md),
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md),
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md)

## Examples

``` r
if (FALSE) { # \dontrun{
nipter_bin_bam_bed("sample.bam", "sample.nipter.bed.gz")

# With pre-filtering matching a typical NIPT pipeline
nipter_bin_bam_bed("sample.dm.bam", "sample.nipter.bed.gz",
                   mapq = 40L, exclude_flags = 1024L)
} # }
```
