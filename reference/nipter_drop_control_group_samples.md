# Drop samples from a NIPTeR control group by sample name

Convenience helper for quickly pruning a control cohort before
rebuilding downstream reference artifacts.

## Usage

``` r
nipter_drop_control_group_samples(control_group, sample_names)
```

## Arguments

- control_group:

  A `NIPTeRControlGroup` object.

- sample_names:

  Character vector of sample names to drop.

## Value

A `NIPTeRControlGroup` of the same strand type with the requested
samples removed.

## See also

[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md),
[`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md)
