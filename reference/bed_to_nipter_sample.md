# Read a NIPTeR-format BED file into a NIPTeRSample

Reads a 5-column or 9-column bgzipped BED file (as written by
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
bed_to_nipter_sample(
  bed,
  name = NULL,
  binsize = NULL,
  autosomal_source = c("auto", "raw", "corrected"),
  sex_source = c("match", "raw", "corrected"),
  con = NULL
)
```

## Arguments

- bed:

  Path to a bgzipped (or plain) BED file with a `.tbi` index.

- name:

  Optional sample name. If `NULL` (default), derived from the BED file
  basename (e.g. `"sample"` from `"sample.bed.gz"`).

- binsize:

  Optional integer; bin size in base pairs. If `NULL` (default),
  inferred from the first row of the BED file.

- autosomal_source:

  Which autosomal counts to realize in the returned sample. `"auto"`
  (default) uses corrected autosomal columns when present, otherwise raw
  counts. `"raw"` always uses raw count columns. `"corrected"` requires
  corrected columns to be present.

- sex_source:

  Which sex-chromosome counts to realize in the returned sample.
  `"match"` (default) follows `autosomal_source`. `"raw"` always uses
  raw count columns. `"corrected"` requires corrected sex columns to be
  present.

- con:

  Optional open DBI connection with duckhts already loaded.

## Value

An object of class `c("NIPTeRSample", <strand_type>)`:

**`CombinedStrands`** (from 5-column BED): same structure as
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
with `separate_strands = FALSE`.

**`SeparatedStrands`** (from 9-column BED): same structure as
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
with `separate_strands = TRUE`.

## Details

A 5-column BED (`chrom`, `start`, `end`, `count`, `corrected_count`)
produces a `CombinedStrands` sample. A 9-column BED (`chrom`, `start`,
`end`, `count`, `count_fwd`, `count_rev`, `corrected_count`,
`corrected_fwd`, `corrected_rev`) produces a `SeparatedStrands` sample
with independent forward/reverse count matrices. Optional leading
`##RWX_<key>=<value>` metadata lines are ignored automatically. The
number of columns is detected automatically.

When corrected columns contain non-NA values (i.e. the BED was written
with a GC-corrected `corrected` argument), `autosomal_source = "auto"`
(default) realizes corrected autosomal counts and `sex_source = "match"`
follows that same choice for X/Y. Callers can override this and request
raw or corrected sex-chromosome values independently of the autosomes.
For `SeparatedStrands`, the per-strand corrected values
(`corrected_fwd`, `corrected_rev`) are used when corrected sex counts
are selected.

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
