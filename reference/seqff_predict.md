# Predict Fetal Fraction with SeqFF

Predict Fetal Fraction with SeqFF

## Usage

``` r
seqff_predict(
  input,
  input_type = c("bam", "sam", "counts"),
  samtools_bin = "samtools",
  samtools_exclude_flags = NULL,
  samtools_min_mapq = NULL
)
```

## Arguments

- input:

  BAM path, SAM path, counts path, or a data frame with `binName` and
  `counts` columns.

- input_type:

  One of `bam`, `sam`, or `counts`.

- samtools_bin:

  Samtools executable to use for BAM input.

- samtools_exclude_flags:

  Optional integer flag passed to `samtools view -F`.

- samtools_min_mapq:

  Optional integer mapping-quality threshold passed to
  `samtools view -q`.

## Value

Named numeric vector with `SeqFF`, `Enet`, and `WRSC`.
