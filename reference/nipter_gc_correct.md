# GC-correct a NIPTeR sample or control group

Adjusts bin counts for GC-content bias using either LOESS regression or
bin-weight normalisation.

## Usage

``` r
nipter_gc_correct(
  object,
  fasta = NULL,
  method = c("loess", "bin"),
  span = 0.75,
  include_sex = FALSE,
  binsize = 50000L,
  gc_table = NULL,
  con = NULL
)
```

## Arguments

- object:

  A `NIPTeRSample` or `NIPTeRControlGroup`.

- fasta:

  Path to an indexed reference FASTA file (.fa/.fasta with .fai index).
  Ignored when `gc_table` is supplied.

- method:

  GC correction method: `"loess"` (default) or `"bin"` (bin-weight).

- span:

  LOESS smoothing parameter (only used when `method = "loess"`). Default
  `0.75`.

- include_sex:

  Logical; correct sex chromosomes (X, Y) as well? Default `FALSE`.

- binsize:

  Bin size used when binning the sample (default 50000L). Must match the
  binsize of the sample. Ignored when `gc_table` is a list (bin size is
  already encoded in the vector lengths).

- gc_table:

  Pre-computed GC table. Either a path to a TSV.bgz file (from
  [`nipter_gc_precompute`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_precompute.md))
  or the in-memory named list returned by a previous `.get_gc_table()`
  call. When `NULL` (default), GC content is computed from `fasta`.

- con:

  Optional open DBI connection with duckhts loaded. If `NULL`, a
  temporary connection is created.

## Value

A corrected copy of `object` with the same class. Correction status is
updated from `"Uncorrected"` to `"GC Corrected"`.

## Details

GC content can be supplied in three ways via the `gc_table` parameter:

- Pre-computed file:

  Path to a TSV.bgz produced by
  [`nipter_gc_precompute`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_precompute.md).
  Fastest for large cohorts — compute once, reuse for every sample.

- In-memory list:

  Named list of numeric vectors (one per chromosome) as returned by
  `.get_gc_table()`. Useful when chaining corrections within a session.

- FASTA path via `fasta`:

  Compute GC on-the-fly for every call. Convenient for single-sample
  use; slow for many samples.

**LOESS method** (default): Fits a LOESS curve of read counts vs GC
percentage across all autosomal bins with known GC and non-zero reads.
Each bin is then scaled by `median(counts) / fitted(loess)`, so that all
bins are normalised to the genome-wide median. This is the NIPTeR
default method.

**Bin-weight method**: Groups bins by GC percentage (0.1\\ computes the
mean read count per GC bucket, then scales each bin by
`global_mean / bucket_mean`. Faster than LOESS but less smooth.

Sex chromosome correction (when `include_sex = TRUE`) uses a
nearest-neighbour lookup against the autosomal LOESS curve (LOESS
method) or the same GC bucket weights (bin-weight method).

## See also

[`nipter_gc_precompute()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_precompute.md),
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md),
[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# One-shot: compute GC and correct in one call
corrected <- nipter_gc_correct(sample, fasta = "hg38.fa")

# Recommended for cohorts: precompute once, reuse
nipter_gc_precompute("hg38.fa", binsize = 50000L, out = "hg38_gc_50k.tsv.bgz")
cg <- nipter_gc_correct(cg, gc_table = "hg38_gc_50k.tsv.bgz")
test_sample <- nipter_gc_correct(test_sample, gc_table = "hg38_gc_50k.tsv.bgz")
} # }
```
