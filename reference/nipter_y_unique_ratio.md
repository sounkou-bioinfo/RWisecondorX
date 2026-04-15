# Compute Y-unique region read ratio from a BAM file

Counts reads overlapping the 7 Y-chromosome unique gene regions (HSFY1,
BPY2, BPY2B, BPY2C, XKRY, PRY, PRY2) and divides by total nuclear genome
reads (chromosomes 1–22, X, Y). This ratio is a strong univariate sex
predictor used in clinical NIPT pipelines.

## Usage

``` r
nipter_y_unique_ratio(
  bam,
  mapq = 1L,
  require_flags = 0L,
  exclude_flags = 0L,
  regions_file = NULL,
  con = NULL,
  reference = NULL
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- mapq:

  Minimum mapping quality. Default `1L`.

- require_flags:

  Integer bitmask; only reads with all bits set are kept (samtools
  `-f`). Default `0L`.

- exclude_flags:

  Integer bitmask; reads with any of these flags set are dropped
  (samtools `-F`). Default `0L`.

- regions_file:

  Path to a TSV file of Y-unique regions with columns `Chromosome`,
  `Start`, `End`, `GeneName`. Defaults to the bundled GRCh37 file.
  Supply a custom file for GRCh38 or other assemblies.

- con:

  Optional open DBI connection with duckhts loaded. If `NULL` (default)
  a temporary in-memory DuckDB connection is created.

- reference:

  Optional FASTA reference path for CRAM inputs.

## Value

A list with elements:

- ratio:

  Numeric; Y-unique reads / total nuclear reads. `0` if total nuclear
  reads is zero.

- y_unique_reads:

  Integer; reads overlapping any of the 7 Y-unique regions.

- total_nuclear_reads:

  Integer; total reads on chromosomes 1–22, X, Y.

- regions:

  Data frame of the regions used (columns: `Chromosome`, `Start`, `End`,
  `GeneName`).

## Details

The regions are defined in the bundled file
`extdata/grch37_Y_UniqueRegions.txt` (GRCh37 coordinates). The BAM must
be indexed (`.bai` or `.csi`).

Read counting uses DuckDB/duckhts with index-based region queries
(`read_bam(region := ...)`) for the Y-unique intervals, and a separate
full-genome scan for the nuclear total. Both queries apply the same MAPQ
and flag filters.

## See also

[`nipter_sex_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md),
[`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
yr <- nipter_y_unique_ratio("sample.bam", mapq = 1L)
yr$ratio
yr$y_unique_reads
yr$total_nuclear_reads
} # }
```
