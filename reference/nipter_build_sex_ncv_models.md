# Build sex-chromosome NCV models on a typed reference

Searches denominator sets for X and Y separately on the female and male
non-outlier subsets of a
[`NIPTReferenceModel`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/NIPTReferenceModel.md).
The resulting `NIPTSexNCVModel` objects are attached to the reference as
`$sex_ncv_models`.

## Usage

``` r
nipter_build_sex_ncv_models(
  reference,
  candidate_chromosomes = c(1:12, 14:16, 20, 22),
  min_elements = 6L,
  max_elements = 9L,
  focus_chromosomes = c("X", "Y")
)
```

## Arguments

- reference:

  A `NIPTReferenceModel` from
  [`nipter_build_reference`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_build_reference.md).

- candidate_chromosomes:

  Integer vector of autosomes allowed in the denominator pool. Defaults
  to the production-style set `c(1:12, 14:16, 20, 22)`.

- min_elements:

  Minimum denominator-set size. Default `6L`.

- max_elements:

  Maximum denominator-set size. Default `9L`.

- focus_chromosomes:

  Character vector; any subset of `c("X", "Y")`.

## Value

The input `NIPTReferenceModel`, with a typed `$sex_ncv_models` field
containing nested `female`/`male` and `X`/`Y` model entries.
