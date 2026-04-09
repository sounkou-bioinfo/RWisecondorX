# Convert BAM/CRAM to WisecondorX bin counts

Replicates the WisecondorX `convert` step: reads aligned reads from a
BAM or CRAM file, applies the same filtering and duplicate-removal
strategy as the upstream Python implementation, and returns per-bin read
counts for chromosomes 1-22, X (as 23), and Y (as 24).

## Usage

``` r
bam_convert(
  bam,
  binsize = 5000L,
  rmdup = c("streaming", "none", "flag"),
  con = NULL,
  reference = NULL
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- binsize:

  Bin size in base pairs (default 5000, matching WisecondorX convert
  default). The reference bin size should be a multiple of this value.

- rmdup:

  Duplicate removal strategy. One of:

  `"streaming"`

  :   Consecutive-position dedup matching WisecondorX's larp/larp2
      streaming state machine (default). Recommended when the BAM has
      not been pre-processed with a dedup tool.

  `"none"`

  :   No duplicate removal. Use for NIPT where read depth is low and
      duplicate removal is not recommended (corresponds to WisecondorX
      `--normdup` flag).

  `"flag"`

  :   Use the SAM 0x400 duplicate flag. Use when the BAM has been
      processed with Picard, sambamba, or similar tools.

- con:

  An optional open DBI connection with the duckhts extension already
  loaded. If `NULL` (default), a temporary in-memory DuckDB connection
  is created for this call.

- reference:

  Optional FASTA reference path for CRAM inputs. Passed to
  `read_bam(..., reference := ...)`. Leave `NULL` for BAM inputs.

## Value

A named list with one integer vector per chromosome key ("1"-"22", "23"
for X, "24" for Y). Each vector has length
`floor(chr_length / binsize) + 1` and contains per-bin read counts.
Chromosomes with no reads present in the BAM are `NULL`.

## Details

The default `rmdup = "streaming"` exactly reproduces the pysam
larp/larp2 deduplication used by WisecondorX. Key subtleties preserved:

- Improper pairs are invisible to the dedup state machine (they never
  update larp/larp2 in pysam's `continue` branch).

- Unpaired reads update larp but NOT larp2.

- larp is never reset between chromosomes.

- Bin assignment matches Python's `int(pos / binsize)` (truncating
  division).

## See also

WisecondorX paper: Huijsdens-van Amsterdam et al. (2018). Conformance
reference: `wisecondorx_convert_conformance.py` in the duckhts
repository, which validates exact bin-for-bin agreement on real NIPT
data.

## Examples

``` r
if (FALSE) { # \dontrun{
bins <- bam_convert("sample.bam", binsize = 5000, rmdup = "streaming")
bins[["11"]]  # bin counts for chromosome 11
} # }
```
