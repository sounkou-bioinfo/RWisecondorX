# Compute Y-unique region read ratio from a BAM file

Counts reads overlapping a Y-chromosome interval set and divides by
total nuclear genome reads (chromosomes 1–22, X, Y). This ratio is a
strong univariate sex predictor used in the clinical NIPT pipeline
mirrored here.

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

  Path to a Y-unique regions file. Headered TSV files with columns
  `Chromosome`, `Start`, `End`, `GeneName` are accepted, as are
  headerless BED-like files whose first 3 columns are
  chromosome/start/end and optional 4th column is a gene or interval
  label. Defaults to the bundled GRCh37 file. Supply a custom file for
  GRCh38 or other assemblies.

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

  Integer; reads overlapping the supplied Y interval set.

- total_nuclear_reads:

  Integer; total reads on chromosomes 1–22, X, Y.

- regions:

  Data frame of the intervals used (columns: `Chromosome`, `Start`,
  `End`, `GeneName`).

- regions_file:

  Normalized path to the regions file used.

## Details

For GRCh37, the default interval set is the bundled legacy
`extdata/grch37_Y_chrom_blacklist.bed`, which matches the historical
`samtools stats -t human_g1k_v37.fasta_Y_chrom_blacklist.bed` path used
by the production shell pipeline. A smaller 7-gene interval file can
still be supplied explicitly through `regions_file`, but it is no longer
the default because it does not reproduce the legacy
`YUniqueRatioFiltered` semantics.

The BAM must be indexed (`.bai` or `.csi`).

Read counting uses DuckDB/duckhts with index-based region queries
(`read_bam(region := ...)`) issued once per Y-unique interval, and a
separate full-genome scan for the nuclear total. Interval chromosome
names are matched to the BAM header so both `Y` and `chrY`-style contigs
work. Both query paths apply the same MAPQ and flag filters.

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
