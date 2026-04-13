# Count reads per bin from a BAM or CRAM file

Core read-counting engine shared by the WisecondorX and NIPTeR layers.
Reads aligned reads from a BAM or CRAM file and returns per-bin read
counts for chromosomes 1-22 and the sex chromosomes (X mapped to key
`"23"`, Y to `"24"`).

## Usage

``` r
bam_convert(
  bam,
  binsize = 5000L,
  mapq = 1L,
  require_flags = 0L,
  exclude_flags = 0L,
  rmdup = c("streaming", "flag", "none"),
  con = NULL,
  reference = NULL
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- binsize:

  Bin size in base pairs. Default `5000L` (WisecondorX); use `50000L`
  for NIPTeR-style workflows.

- mapq:

  Minimum mapping quality. Default `1L` (WisecondorX / samtools
  default). Set to `0L` to retain all reads regardless of MAPQ (NIPTeR).

- require_flags:

  Integer bitmask. Only reads for which
  `(FLAG & require_flags) == require_flags` are retained. `0L` (default)
  imposes no requirement. Example: `require_flags = 0x2L` keeps only
  properly paired reads.

- exclude_flags:

  Integer bitmask. Reads for which `(FLAG & exclude_flags) != 0` are
  dropped. `0L` (default) drops nothing. Example:
  `exclude_flags = 0xF04L` drops unmapped, secondary, QC-fail and
  supplementary reads (common samtools pre-filter).

- rmdup:

  Duplicate removal strategy. `"streaming"` (default) applies the
  WisecondorX larp/larp2 algorithm (also excludes improper pairs).
  `"flag"` drops reads with SAM flag `0x400`. `"none"` keeps all reads
  that pass the other filters.

- con:

  Optional open DBI connection with duckhts already loaded. If `NULL`
  (default) a temporary in-memory DuckDB connection is created.

- reference:

  Optional FASTA reference path for CRAM inputs.

## Value

A named list with one integer vector per chromosome key (`"1"`– `"22"`,
`"23"` for X, `"24"` for Y). Each vector contains per-bin read counts
(bin 0 = positions 0 to `binsize - 1`). Chromosomes absent from the BAM
are `NULL`.

## Details

Read filtering mirrors the `samtools view` convention: `mapq` sets the
minimum mapping quality; `require_flags` is a bitmask of flags that must
**all** be set (equivalent to `samtools view -f`); `exclude_flags` is a
bitmask of flags that must **all** be clear (equivalent to
`samtools view -F`). Use the duckhts UDF `sam_flag_bits(flag)` to
inspect named flag fields, or `sam_flag_has(flag, bit)` to test
individual bits.

`rmdup` controls duplicate removal independently of the flag filters:
`"streaming"` applies the WisecondorX larp/larp2 consecutive-position
state machine and also enforces the WisecondorX improper-pair rule
(paired reads that are not properly paired are excluded from both
counting and the dedup state — this is intrinsic to the algorithm, not a
flag option); `"flag"` additionally excludes reads with SAM flag `0x400`
set (Picard / sambamba pre-marked duplicates); `"none"` applies no
deduplication.

## See also

[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md),
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# WisecondorX defaults
bins <- bam_convert("sample.bam")

# NIPTeR defaults — all mapped reads, no dedup, 50 kb bins
bins <- bam_convert("sample.bam", binsize = 50000L, mapq = 0L,
                    rmdup = "none")

# Pre-filtered BAM: skip unmapped + secondary + supplementary, flag dedup
bins <- bam_convert("sample.bam",
                    exclude_flags = bitwOr(0x4L, bitwOr(0x100L, 0x800L)),
                    rmdup = "flag")
} # }
```
