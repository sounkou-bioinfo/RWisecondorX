# Write an in-memory NIPTeRSample to a bgzipped BED file

Serialises an already-binned `NIPTeRSample` — optionally GC-corrected —
to a bgzipped, tabix-indexed BED file without re-reading any BAM.

## Usage

``` r
nipter_sample_to_bed(
  sample,
  bed,
  binsize,
  corrected = NULL,
  con = NULL,
  index = TRUE
)
```

## Arguments

- sample:

  A `NIPTeRSample` object (`CombinedStrands` or `SeparatedStrands`) from
  [`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md).

- bed:

  Path for the output `.bed.gz` file.

- binsize:

  Bin size in base pairs used when `sample` was created. Required to
  reconstruct genomic coordinates.

- corrected:

  Optional GC-corrected `NIPTeRSample` from
  [`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md).
  Populated corrected-count columns when supplied; those columns are
  `NA` otherwise.

- con:

  Optional open DBI connection with duckhts already loaded.

- index:

  Logical; write a tabix index alongside the BED (default `TRUE`).

## Value

`bed` (invisibly).

## Details

Use this when you have already called
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
and optionally
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md):
the BAM is read exactly once and the result is written in one step. For
a single one-shot BAM → BED conversion with no intermediate sample
object, see
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md).

Column layout mirrors
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md):

- **CombinedStrands** (default): 5 columns — `chrom`, `start`, `end`,
  `count`, `corrected_count`.

- **SeparatedStrands** (`separate_strands = TRUE`): 9 columns — `chrom`,
  `start`, `end`, `count`, `count_fwd`, `count_rev`, `corrected_count`,
  `corrected_fwd`, `corrected_rev`.

## See also

[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md),
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md),
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md),
[`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Bin once, write — no GC correction
s <- nipter_bin_bam("sample.bam", binsize = 50000L)
nipter_sample_to_bed(s, "sample.bed.gz", binsize = 50000L)

# Bin once, GC-correct, write — BAM read only once
s    <- nipter_bin_bam("sample.bam", binsize = 50000L)
corr <- nipter_gc_correct(s, gc_table = "hg38_gc_50k.tsv.bgz")
nipter_sample_to_bed(s, "sample.bed.gz", binsize = 50000L, corrected = corr)
} # }
```
