# Convert BAM/CRAM to a bgzipped BED bin-count file

Runs
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
and writes the per-bin read counts to a bgzipped, tab-delimited BED file
(four columns: `chrom`, `start`, `end`, `count`), then creates a tabix
index (`.tbi`) alongside it via Rduckhts. This is the language-agnostic
intermediate format for the RWisecondorX native pipeline; the file can
be queried directly with duckhts/DuckDB or any tabix-aware tool.

## Usage

``` r
bam_convert_bed(
  bam,
  bed,
  binsize = 5000L,
  mapq = 1L,
  require_flags = 0L,
  exclude_flags = 0L,
  rmdup = c("streaming", "flag", "none"),
  con = NULL,
  reference = NULL,
  index = TRUE
)
```

## Arguments

- bam:

  Path to an indexed BAM or CRAM file.

- bed:

  Path for the output `.bed.gz` file (created or overwritten). The tabix
  index is written to `paste0(bed, ".tbi")`.

- binsize:

  Bin size in base pairs (default 5000).

- mapq:

  Minimum mapping quality (default `1L`).

- require_flags:

  Integer bitmask; only reads with all bits set are kept (samtools
  `-f`). Default `0L` (no requirement).

- exclude_flags:

  Integer bitmask; reads with any bit set are dropped (samtools `-F`).
  Default `0L` (nothing dropped).

- rmdup:

  Duplicate-removal strategy passed to
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

- con:

  Optional open DBI connection with duckhts already loaded. If `NULL`
  (default), a temporary in-memory DuckDB connection is created.

- reference:

  Optional FASTA reference path for CRAM inputs.

- index:

  Logical; when `TRUE` (default) a tabix index is created alongside the
  bgzipped output.

## Value

`bed` (invisibly).

## Details

Coordinates are 0-based half-open intervals matching the BED convention
(`start = bin_index * binsize`, `end = start + binsize`). Chromosomes
are written as `1`–`22`, `X`, `Y` (no `chr` prefix) in numeric order.
All header-defined bins are written, including those with a count of
zero.

bgzip and tabix indexing are performed via
[`Rduckhts::rduckhts_bgzip()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_bgzip.html)
and
[`Rduckhts::rduckhts_tabix_index()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_tabix_index.html),
so no external tools are required.

## See also

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md),
[`wisecondorx_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_convert.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bam_convert_bed("sample.bam", "sample.bed.gz", binsize = 5000, rmdup = "streaming")
} # }
```
