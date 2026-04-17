# Load a NIPTeR control group from a directory of TSV.bgz files

Reads all `.bed.gz` (or `.tsv.bgz`) files in `bed_dir` using
`rduckhts_tabix_multi()` — a single multi-file DuckDB scan — and
constructs a `NIPTeRControlGroup` from the results. This is much faster
than `lapply(files, bed_to_nipter_sample)` for large cohorts because all
files are read in one pass.

## Usage

``` r
nipter_control_group_from_beds(
  bed_dir,
  pattern = "*.bed.gz",
  binsize = NULL,
  autosomal_source = c("auto", "raw", "corrected"),
  sex_counts = c("match", "raw", "corrected"),
  description = "General control group",
  sample_sex = NULL,
  sex_source = NULL,
  con = NULL
)
```

## Arguments

- bed_dir:

  Character; directory containing one `.bed.gz` or `.tsv.bgz` file per
  control sample, each produced by
  [`nipter_bin_bam_bed`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md).

- pattern:

  Glob pattern for file discovery (default `"*.bed.gz"`).

- binsize:

  Optional integer; bin size in base pairs. If `NULL` (default),
  inferred from the first row of the first file.

- autosomal_source:

  Which autosomal counts to realize in each imported sample. `"auto"`
  (default) uses corrected autosomal columns when present, otherwise raw
  counts. `"raw"` always uses raw count columns. `"corrected"` requires
  corrected columns.

- sex_counts:

  Which sex-chromosome counts to realize in each imported sample.
  `"match"` (default) follows `autosomal_source`. `"raw"` always uses
  raw count columns. `"corrected"` requires corrected sex columns.

- description:

  Label for the resulting control group (default
  `"General control group"`).

- sample_sex:

  Optional character vector of known sex labels for the samples in
  `bed_dir`. Names must match the inferred sample names when supplied as
  a named vector.

- sex_source:

  Optional string describing the provenance of `sample_sex`.

- con:

  Optional open DBI connection with duckhts loaded.

## Value

A `NIPTeRControlGroup`.

## Details

The column count (5 or 9) is detected automatically from the first file:
5-column BEDs produce a `CombinedStrands` control group; 9-column BEDs
(written by `nipter_bin_bam_bed(separate_strands = TRUE)`) produce a
`SeparatedStrands` control group. All files in the directory must have
the same column count. When corrected BED columns are present,
`autosomal_source` and `sex_counts` control whether the imported samples
use raw or corrected values for each compartment before constructing the
control group.

## See also

[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md),
[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md),
[`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Bin all controls to BED once
for (bam in bam_files) {
  nipter_bin_bam_bed(bam, file.path("controls/", sub(".bam$", ".bed.gz", basename(bam))))
}
# Load them all at once
cg <- nipter_control_group_from_beds("controls/")
} # }
```
