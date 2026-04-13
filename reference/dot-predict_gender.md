# Predict gender from Y-chromosome read fraction

Classifies a sample as male or female based on the fraction of reads
mapping to chrY, using a cutoff from the reference GMM model. Mirrors
`predict_tools.predict_gender()` in upstream WisecondorX.

## Usage

``` r
.predict_gender(sample, trained_cutoff)
```

## Arguments

- sample:

  Named list of integer/numeric vectors keyed by chromosome.

- trained_cutoff:

  Numeric; Y-fraction threshold from reference.

## Value

`"M"` or `"F"`.
