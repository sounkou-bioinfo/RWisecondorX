# Read a NIPTeR-format BED file into a NIPTeRSample

Reads a 5-column or 7-column bgzipped BED file (as written by
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md))
and returns a `NIPTeRSample` object suitable for all NIPTeR statistical
functions:
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md),
[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md),
[`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md),
[`nipter_regression()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md),
and
[`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md).

## Usage

``` r
bed_to_nipter_sample(bed, name = NULL, binsize = NULL, con = NULL)
```

## Arguments

- bed:

  Path to a bgzipped (or plain) BED file.

- name:

  Optional sample name. If `NULL` (default), derived from the BED file
  basename (e.g. `"sample"` from `"sample.bed.gz"`).

- binsize:

  Optional integer; bin size in base pairs. If `NULL` (default),
  inferred from the first row of the BED file.

- con:

  Optional open DBI connection with duckhts already loaded.

## Value

An object of class `c("NIPTeRSample", <strand_type>)`:

**`CombinedStrands`** (from 5-column BED): same structure as
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
with `separate_strands = FALSE`.

**`SeparatedStrands`** (from 7-column BED): same structure as
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
with `separate_strands = TRUE`.

## Details

A 5-column BED (`chrom`, `start`, `end`, `count`, `corrected_count`)
produces a `CombinedStrands` sample. A 7-column BED (`chrom`, `start`,
`end`, `count`, `count_fwd`, `count_rev`, `corrected_count`) produces a
`SeparatedStrands` sample with independent forward/reverse count
matrices. The number of columns is detected automatically.

## See also

[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md),
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md),
[`bed_to_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_sample.md),
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md),
[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Write NIPTeR-style bin counts to BED, then read back
nipter_bin_bam_bed("sample.bam", "sample.bed.gz")
sample <- bed_to_nipter_sample("sample.bed.gz")

# Use with NIPTeR scoring
samples <- lapply(bed_files, bed_to_nipter_sample)
cg <- nipter_as_control_group(samples)
z21 <- nipter_z_score(samples[[1]], cg, chromo_focus = 21)
} # }
```
