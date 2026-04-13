# GC-correct a NIPTeR sample or control group

Adjusts bin counts for GC-content bias using either LOESS regression or
bin-weight normalisation. GC content per bin is computed on-the-fly from
the reference FASTA via
[`rduckhts_fasta_nuc`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_fasta_nuc.html)
rather than from bundled precomputed tables (as the original NIPTeR does
with `sysdata.rda`).

## Usage

``` r
nipter_gc_correct(
  object,
  fasta,
  method = c("loess", "bin"),
  span = 0.75,
  include_sex = FALSE,
  binsize = 50000L,
  con = NULL
)
```

## Arguments

- object:

  A `NIPTeRSample` or `NIPTeRControlGroup`.

- fasta:

  Path to an indexed reference FASTA file (.fa/.fasta with .fai index).

- method:

  GC correction method: `"loess"` (default) or `"bin"` (bin-weight).

- span:

  LOESS smoothing parameter (only used when `method = "loess"`). Default
  `0.75`.

- include_sex:

  Logical; correct sex chromosomes (X, Y) as well? Default `FALSE`.

- binsize:

  Bin size used when binning the sample (default 50000L). Must match the
  binsize of the sample.

- con:

  Optional open DBI connection with duckhts loaded. If `NULL`, a
  temporary connection is created.

## Value

A corrected copy of `object` with the same class. Correction status is
updated from `"Uncorrected"` to `"GC corrected"`.

## Details

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

[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md),
[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sample <- nipter_bin_bam("sample.bam")
corrected <- nipter_gc_correct(sample, fasta = "hg38.fa")

# Correct an entire control group
cg <- nipter_as_control_group(samples)
cg_corrected <- nipter_gc_correct(cg, fasta = "hg38.fa")
} # }
```
