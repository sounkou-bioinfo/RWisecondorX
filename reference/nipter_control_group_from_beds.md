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
  description = "General control group",
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

- description:

  Label for the resulting control group (default
  `"General control group"`).

- con:

  Optional open DBI connection with duckhts loaded.

## Value

A `NIPTeRControlGroup`.

## Details

The column count (5 or 9) is detected automatically from the first file:
5-column BEDs produce a `CombinedStrands` control group; 9-column BEDs
(written by `nipter_bin_bam_bed(separate_strands = TRUE)`) produce a
`SeparatedStrands` control group. All files in the directory must have
the same column count.

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
