# Iteratively prune aberrant controls using chi-correction diagnostics

Mirrors the legacy production cleaning loop for NIPTeR control cohorts:
chi-correct the current control group, diagnose chromosomal-fraction
outliers, drop flagged samples, and repeat until no further outliers
remain or the minimum control size would be violated.

## Usage

``` r
nipter_prune_control_group_outliers(
  control_group,
  sample = NULL,
  chi_cutoff = 3.5,
  collapse_strands = FALSE,
  z_cutoff = 3,
  max_aberrant_chromosomes = 2L,
  outlier_rule = c("any_aberrant_score", "bidirectional_or_multichromosome"),
  min_controls = 10L,
  max_iterations = 100L,
  verbose = FALSE
)
```

## Arguments

- control_group:

  A `NIPTeRControlGroup`.

- sample:

  Optional `NIPTeRSample` used as the dummy sample passed into
  [`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md).
  Defaults to the first control sample.

- chi_cutoff:

  Numeric chi-correction cutoff passed to
  [`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md).
  Default `3.5`.

- collapse_strands:

  Logical; keep 44 strand-resolved diagnostics for `SeparatedStrands`
  (`FALSE`, default) or collapse to 22 autosomal fractions.

- z_cutoff:

  Absolute Z-score cutoff passed to
  [`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md).

- max_aberrant_chromosomes:

  Maximum number of distinct aberrant chromosomes allowed before a
  sample is dropped. Default `2L`.

- outlier_rule:

  Character scalar controlling how aberrant samples are dropped.
  `"any_aberrant_score"` (default) removes any sample appearing in
  `abberant_scores`, matching the upstream NIPTeR vignette.
  `"bidirectional_or_multichromosome"` drops a sample only when both
  strands of one chromosome are aberrant or when more than
  `max_aberrant_chromosomes` distinct chromosomes are aberrant.

- min_controls:

  Minimum allowed retained control count. The iteration stops before
  dropping samples if doing so would leave fewer controls than this
  threshold.

- max_iterations:

  Maximum number of pruning iterations. Default `100L`.

- verbose:

  Logical; emit per-iteration messages.

## Value

A named list with:

- control_group:

  The final uncorrected control group after dropping flagged samples.

- chi_corrected_control_group:

  The final chi-corrected control group used for the terminal diagnostic
  pass.

- diagnostics:

  The terminal diagnostic list from
  [`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md).

- dropped_samples:

  Character vector of unique dropped sample names in drop order.

- iteration_log:

  Data frame recording flagged and retained control counts per
  iteration.

- converged:

  Logical; `TRUE` when the loop ended without flagged samples.

- stop_reason:

  Text reason for the stopping condition.

## Details

For `SeparatedStrands` control groups with `collapse_strands = FALSE`
(default), the pruning rule matches the legacy script: drop a sample if
both strands of any chromosome are aberrant, or if more than
`max_aberrant_chromosomes` distinct chromosomes are aberrant.

## See also

[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md),
[`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md),
[`nipter_drop_control_group_samples()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_drop_control_group_samples.md)
