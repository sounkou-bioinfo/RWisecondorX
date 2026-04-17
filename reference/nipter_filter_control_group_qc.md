# Filter control samples by read-depth and GC QC

Applies hard sample-level QC gates before reference/model building. This
is intended for cohort curation steps such as excluding controls with
too few unique reads or aberrant GC content.

## Usage

``` r
nipter_filter_control_group_qc(
  control_group,
  sample_qc,
  sample_col = NULL,
  total_unique_reads_col = NULL,
  gc_col = NULL,
  min_total_unique_reads = NULL,
  max_total_unique_reads = NULL,
  gc_mad_cutoff = NULL
)
```

## Arguments

- control_group:

  A `NIPTeRControlGroup`.

- sample_qc:

  Data frame containing at least a sample-name column and, if the
  corresponding filters are enabled, QC columns for total unique reads
  and/or GC content.

- sample_col:

  Optional sample-name column in `sample_qc`. When `NULL`, common names
  such as `sample_name` and `Sample` are inferred.

- total_unique_reads_col:

  Optional total-unique-reads column in `sample_qc`. Required when
  either read-depth threshold is supplied. If `NULL`, common names such
  as `TotalUniqueReads` are inferred.

- gc_col:

  Optional GC column in `sample_qc`. Required when `gc_mad_cutoff` is
  supplied. If `NULL`, common names such as `GCPCTAfterFiltering` are
  inferred.

- min_total_unique_reads:

  Optional minimum allowed total unique reads. Samples below this
  threshold are dropped.

- max_total_unique_reads:

  Optional maximum allowed total unique reads. Samples above this
  threshold are dropped.

- gc_mad_cutoff:

  Optional robust MAD cutoff for the GC column. Samples more than this
  many MADs from the cohort median are dropped.

## Value

A list with:

- control_group:

  The filtered control group.

- retained_samples:

  Character vector of retained sample names.

- excluded_samples:

  Character vector of excluded sample names.

- exclusion_table:

  Data frame describing dropped samples and reasons.

- matched_qc:

  The QC rows matched to the control-group sample order.

- settings:

  Named list of applied QC filter settings.

## See also

[`nipter_drop_control_group_samples()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_drop_control_group_samples.md),
[`nipter_build_reference()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md)
