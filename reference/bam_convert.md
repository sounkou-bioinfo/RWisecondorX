# Count reads per bin from a BAM or CRAM file

Core read-counting engine shared by the WisecondorX and NIPTeR layers.
Reads aligned reads from a BAM or CRAM file and returns per-bin read
counts for chromosomes 1-22 and the sex chromosomes (X mapped to key
"23", Y to "24").

## Usage

``` r
bam_convert(
  bam,
  binsize = 5000L,
  mapq = 1L,
  rmdup = c("streaming", "flag", "none"),
  filter_improper_pairs = TRUE,
  con = NULL,
  reference = NULL
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- binsize:

  Bin size in base pairs. Default 5000 (WisecondorX); use 50000 for
  NIPTeR-style workflows.

- mapq:

  Minimum mapping quality to retain a read. Default `1L` (WisecondorX).
  Set to `0L` to disable MAPQ filtering (NIPTeR).

- rmdup:

  Duplicate removal strategy: `"streaming"` — WisecondorX larp/larp2
  consecutive-position dedup (default; not meaningful when
  `filter_improper_pairs = FALSE`); `"flag"` — exclude reads with SAM
  flag `0x400` (pre-marked by Picard / sambamba); `"none"` — no
  duplicate removal.

- filter_improper_pairs:

  When `TRUE` (default, WisecondorX behaviour) paired reads that are not
  properly paired (`FLAG & 0x2 == 0`) are excluded. Set to `FALSE` to
  include all mapped reads regardless of pair status (NIPTeR behaviour).

- con:

  An optional open DBI connection with the duckhts extension already
  loaded. If `NULL` (default), a temporary in-memory DuckDB connection
  is created for this call.

- reference:

  Optional FASTA reference path for CRAM inputs.

## Value

A named list with one integer vector per chromosome key (`"1"`–`"22"`,
`"23"` for X, `"24"` for Y). Each vector contains per-bin read counts
(bin 0 = positions 0 to binsize-1). Chromosomes absent from the BAM are
`NULL`.

## Details

Filter parameters let callers replicate the exact behaviour of both
tools. WisecondorX defaults: `mapq = 1`, `rmdup = "streaming"`,
`filter_improper_pairs = TRUE`. NIPTeR defaults: `mapq = 0`,
`rmdup = "none"`, `filter_improper_pairs = FALSE`.

When `rmdup = "streaming"` the WisecondorX larp/larp2 state machine is
reproduced exactly: improper pairs are invisible to the dedup logic
(they do not update larp2), unpaired reads update larp but not larp2,
and larp is never reset between chromosomes. Bin assignment uses
truncating integer division matching Python's `int(pos / binsize)`.

## See also

[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md),
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# WisecondorX defaults
bins <- bam_convert("sample.bam", binsize = 5000L, rmdup = "streaming")

# NIPTeR defaults
bins <- bam_convert("sample.bam", binsize = 50000L, mapq = 0L,
                    rmdup = "none", filter_improper_pairs = FALSE)
} # }
```
