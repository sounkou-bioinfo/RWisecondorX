# Write NIPTeR reference QC plots

Writes PNG plots for the existing NIPTeR QC bundle and reference
sex-model spaces.

## Usage

``` r
write_nipter_reference_plots(qc, reference, outprefix, dpi = 150L)
```

## Arguments

- qc:

  A `NIPTControlGroupQC` from
  [`nipter_control_group_qc`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_control_group_qc.md).

- reference:

  A `NIPTReferenceModel` from
  [`nipter_build_reference`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md).

- outprefix:

  Output file prefix.

- dpi:

  Output PNG DPI. Default `150L`.

## Value

A named character vector of written file paths, invisibly.
