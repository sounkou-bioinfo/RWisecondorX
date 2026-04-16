# Predict Fetal Fraction with SeqFF

Predict Fetal Fraction with SeqFF

## Usage

``` r
seqff_predict(
  input,
  input_type = c("bam", "sam", "counts"),
  mapq = NULL,
  require_flags = NULL,
  exclude_flags = NULL,
  con = NULL,
  reference = NULL
)
```

## Arguments

- input:

  BAM path, SAM path, counts path, or a data frame with `binName` and
  `counts` columns.

- input_type:

  One of `bam`, `sam`, or `counts`.

- mapq:

  Optional integer mapping-quality threshold for BAM input.

- require_flags:

  Optional integer bitmask; only reads with all bits set are kept for
  BAM input.

- exclude_flags:

  Optional integer bitmask; reads with any bit set are dropped for BAM
  input.

- con:

  Optional open DBI connection with duckhts already loaded.

- reference:

  Optional FASTA reference path for CRAM input.

## Value

Named numeric vector with `SeqFF`, `Enet`, and `WRSC`.
