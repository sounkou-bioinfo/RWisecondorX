# Bin a BAM/CRAM and write NIPTeR-style counts to a bgzipped BED

Convenience wrapper that bins a BAM/CRAM with
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
and immediately writes the result via
[`nipter_sample_to_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sample_to_bed.md).
The BAM is read exactly once.

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
  separate_strands = FALSE,
  con = NULL,
  reference = NULL,
  index = TRUE,
  metadata = NULL
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

- separate_strands:

  Logical; when `TRUE`, includes separate `count_fwd` and `count_rev`
  columns. Default `FALSE`.

- con:

  Optional open DBI connection with duckhts already loaded.

- reference:

  Optional FASTA reference path for CRAM inputs.

- index:

  Logical; create a tabix index (default `TRUE`).

- metadata:

  Optional named list or named atomic vector of provenance metadata to
  write as leading `##RWX_<key>=<value>` lines before the BED body. The
  data rows remain headerless. Default `NULL` writes no metadata.

## Value

`bed` (invisibly).

## Details

If you need to GC-correct the sample before writing, call
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
and
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)
yourself, then pass both to
[`nipter_sample_to_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sample_to_bed.md).
That avoids reading the BAM a second time.

Column layout:

- **CombinedStrands** (default): 5 columns — `chrom`, `start`, `end`,
  `count`, `corrected_count`.

- **SeparatedStrands** (`separate_strands = TRUE`): 9 columns — `chrom`,
  `start`, `end`, `count`, `count_fwd`, `count_rev`, `corrected_count`,
  `corrected_fwd`, `corrected_rev`.

## See also

[`nipter_sample_to_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sample_to_bed.md),
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

# Strand-separated output (9 columns)
nipter_bin_bam_bed("sample.bam", "sample.stranded.bed.gz",
                   separate_strands = TRUE)
} # }
```
