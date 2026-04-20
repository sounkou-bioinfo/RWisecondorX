# Build a typed QC report for a NIPTeR control group

Summarises the stability of a `NIPTControlGroup` at the chromosome,
sample, and optional autosomal-bin level. The chromosome summary exposes
the mean, standard deviation, coefficient of variation, and Shapiro-Wilk
normality statistic for each autosome, across multiple scoring spaces in
one long table: post-chi chromosomal fractions, NCV, and the four RBZ
predictor sets. The sample summary adds SSD-based matching scores plus
strand-aware chromosomal alert counts from the post-chi diagnostic pass.
When sex labels are available, the sample summary is enriched with the
same per-sample sex-cluster metrics written into the typed reference
frame (`ConsensusGender`, sex-outlier flag, and `Z_X_XX`/`Z_X_XY`/
`Z_Y_XX`/`Z_Y_XY`), and the report also includes sex-stratified
`RR_X`/`RR_Y` spread metrics relevant for gaunosome model readiness plus
a dedicated long-form sex-model summary for XX/XY sex-chromosome Z, NCV,
and regression models. Optional bin-level output exposes the
scaled-count CV and chi-correction profile that underlies
[`nipter_chi_correct`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md).

## Usage

``` r
nipter_control_group_qc(
  control_group,
  sample_sex = NULL,
  reference_model = NULL,
  chi_cutoff = 3.5,
  z_cutoff = 3,
  collapse_strands = FALSE,
  max_aberrant_chromosomes = 2L,
  outlier_rule = c("any_aberrant_score", "bidirectional_or_multichromosome"),
  include_bins = FALSE
)
```

## Arguments

- control_group:

  A `NIPTControlGroup`.

- sample_sex:

  Optional character vector overriding the sex labels stored on
  `control_group`. Accepted values are `"female"`, `"male"`,
  `"ambiguous"`, and `"unknown"`.

- reference_model:

  Optional `NIPTReferenceModel` built from the same controls. When
  supplied, the QC bundle also extracts the explicit XX/XY
  sex-chromosome Z, NCV, and regression model summaries from it.

- chi_cutoff:

  Numeric scalar used to flag overdispersed bins in the optional chi
  profile. Default `3.5`.

- z_cutoff:

  Absolute post-chi z-score threshold used when counting aberrant
  chromosome rows in the sample summary. Default `3`.

- collapse_strands:

  Logical; when `FALSE` (default), separated strands retain their legacy
  F/R-specific diagnostic rows and the sample summary reports both
  distinct chromosomes and bidirectional aberrations. When `TRUE`,
  sample-level aberration counts collapse to 22 autosomes.

- max_aberrant_chromosomes:

  Maximum number of distinct aberrant chromosomes allowed before
  `is_chromosomal_outlier` is set in the sample summary. Default `2L`.

- outlier_rule:

  Character scalar describing how to convert `abberant_scores` into
  sample-level outlier flags. `"any_aberrant_score"` (default) removes
  any sample with at least one post-chi chromosome row beyond
  `z_cutoff`. `"bidirectional_or_multichromosome"` drops a sample only
  when both strands of one chromosome are aberrant or when too many
  distinct chromosomes are aberrant.

- include_bins:

  Logical; when `TRUE`, compute the autosomal bin-level scaled-count CV
  and chi profile. Default `FALSE`.

## Value

A typed `NIPTControlGroupQC` object.

## See also

[`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md),
[`nipter_match_matrix()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_matrix.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)
