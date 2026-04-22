# Compute samtools-stats-style BAM summary numbers

Runs a native htslib-backed summary pass over a BAM or CRAM file and
returns a subset of the summary-number (`SN`) fields produced by
`samtools stats`. The core mismatch, duplicate-flag, and
CIGAR-denominator logic is ported from upstream `samtools` `stats.c`.

## Usage

``` r
bam_samtools_stats_summary(
  path,
  reference = NULL,
  min_mapq = 0L,
  require_flags = 0L,
  exclude_flags = 0L,
  decompression_threads = 0L,
  report_filtered_stream_bookkeeping = FALSE
)
```

## Arguments

- path:

  Path to the input BAM or CRAM file.

- reference:

  Optional reference FASTA path for CRAM input when required.

- min_mapq:

  Integer; additional minimum MAPQ filter applied before the
  samtools-style summary logic. This is a package-level extension for
  preprocessing workflows; native `samtools stats` does not expose a
  direct MAPQ threshold option.

- require_flags:

  Integer; required SAM flag mask, equivalent in spirit to samtools
  `-f`.

- exclude_flags:

  Integer; excluded SAM flag mask, equivalent in spirit to samtools
  `-F`.

- decompression_threads:

  Integer; number of htslib decompression worker threads to use. `0`
  disables worker threads.

- report_filtered_stream_bookkeeping:

  Logical. When `FALSE` (default), `raw_total_sequences` stays anchored
  to the original BAM/CRAM input and `filtered_sequences` reports the
  number excluded by `min_mapq` / flag-based filtering. When `TRUE`,
  these two fields instead follow the semantics of
  `samtools view ... | samtools stats -`, where the retained stream is
  treated as the input and `filtered_sequences` is therefore `0`.

## Value

A one-row `data.frame`.

## Details

`error_rate` follows upstream `samtools stats` semantics exactly:
`mismatches_from_nm / bases_mapped_cigar`, where mismatches are
accumulated from the `NM` auxiliary tag and `bases_mapped_cigar` counts
`M`, `I`, `=`, and `X` CIGAR operations.

`reads_duplicated` and the derived duplicate fractions are
duplicate-flag summaries. They are **not** Picard `PERCENT_DUPLICATION`.

For exact CLI comparison when using `min_mapq`, compare against a
pipeline of the form `samtools view -u -q ... | samtools stats -`.
