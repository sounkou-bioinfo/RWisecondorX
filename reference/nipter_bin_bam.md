# Bin a BAM/CRAM file — NIPTeR style

Replicates NIPTeR's `bin_bam_sample()`: counts reads in fixed-width bins
across chromosomes 1-22, X, and Y. The original NIPTeR counts all mapped
reads with no MAPQ filter; real-world NIPT pipelines typically
pre-filter with `samtools view --min-MQ 40 -F 1024` before binning. Both
modes are supported through the `mapq`, `exclude_flags`, and
`require_flags` parameters.

## Usage

``` r
nipter_bin_bam(
  bam,
  binsize = 50000L,
  mapq = 0L,
  require_flags = 0L,
  exclude_flags = 0L,
  rmdup = c("none", "flag"),
  separate_strands = FALSE,
  con = NULL,
  reference = NULL
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- binsize:

  Bin size in base pairs. Default 50000 (NIPTeR's fixed bin size). Must
  match the binsize used when building any control group.

- mapq:

  Minimum mapping quality. Default `0L` (NIPTeR original: all mapped
  reads). Set to `40L` to match common NIPT pipeline pre-filtering
  (`samtools view --min-MQ 40`).

- require_flags:

  Integer bitmask; only reads with all bits set are kept (samtools
  `-f`). Default `0L` (no requirement).

- exclude_flags:

  Integer bitmask; reads with any bit set are dropped (samtools `-F`).
  Default `0L`. Set to `1024L` (`0x400`) to exclude reads marked as
  duplicates (`samtools view -F 1024`).

- rmdup:

  Duplicate removal strategy. `"none"` (default, NIPTeR standard) or
  `"flag"` (equivalent to `exclude_flags = 1024L`; for BAMs already
  processed by Picard / sambamba). NIPTeR does not perform streaming
  deduplication.

- separate_strands:

  Not yet implemented. Set `FALSE` (default).

- con:

  Optional open DBI connection with duckhts already loaded.

- reference:

  Optional FASTA reference path for CRAM inputs.

## Value

An object of class `NIPTeRSample`: a named list with
`autosomal_chromosome_reads` (a list of one 22-row integer matrix, rows
named `"1"`–`"22"`, columns are bins), `sex_chromosome_reads` (a list of
one 2-row integer matrix, rows named `"X"` and `"Y"`),
`correction_status_autosomal` (`"Uncorrected"`), `correction_status_sex`
(`"Uncorrected"`), and `sample_name`.

## Details

The result is reshaped into a `NIPTeRSample` object whose structure
parallels NIPTeR's `NIPTSample`: autosomal reads as a chromosome-by-bin
matrix and sex chromosome reads as a separate two-row matrix.

Strand separation (`separate_strands = TRUE`) is not yet implemented; it
requires a forward/reverse split in the SQL layer and will be added when
the NIPTeR regression layer is ported.

## See also

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md),
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# NIPTeR original defaults: all mapped reads, no dedup
sample <- nipter_bin_bam("sample.bam", binsize = 50000L)

# Common NIPT pipeline: MAPQ >= 40, exclude duplicate-flagged reads
sample <- nipter_bin_bam("sample.dm.bam", binsize = 50000L,
                         mapq = 40L, exclude_flags = 1024L)

sample$autosomal_chromosome_reads[[1]]["21", ]   # chr21 bin counts
} # }
```
