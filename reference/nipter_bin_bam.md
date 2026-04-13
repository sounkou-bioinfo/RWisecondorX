# Bin a BAM/CRAM file — NIPTeR style

Replicates NIPTeR's `bin_bam_sample()`: counts reads in fixed-width bins
across chromosomes 1-22, X, and Y, with no MAPQ filter and no
pair-status filter (all mapped reads are counted), matching the
Rsamtools `scanBam` defaults used by the original.

## Usage

``` r
nipter_bin_bam(
  bam,
  binsize = 50000L,
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

- rmdup:

  Duplicate removal strategy. `"none"` (default, NIPTeR standard) or
  `"flag"` (for BAMs already processed by Picard / sambamba). NIPTeR
  does not perform streaming deduplication.

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

Internally calls
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
with `mapq = 0L`, `filter_improper_pairs = FALSE`, and the
caller-supplied `rmdup`. The result is reshaped into a `NIPTeRSample`
object whose structure parallels NIPTeR's `NIPTSample`: autosomal reads
as a chromosome-by-bin matrix and sex chromosome reads as a separate
two-row matrix.

Strand separation (`separate_strands = TRUE`) is not yet implemented; it
requires a forward/reverse split in the SQL layer and will be added when
the NIPTeR regression layer is ported.

## See also

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md),
`nipter_gc_correct()`

## Examples

``` r
if (FALSE) { # \dontrun{
sample <- nipter_bin_bam("sample.bam", binsize = 50000L)
sample$autosomal_chromosome_reads[[1]]["21", ]   # chr21 bin counts
} # }
```
