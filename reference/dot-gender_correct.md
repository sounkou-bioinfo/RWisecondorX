# Correct gonosomal read counts for male samples

Doubles chrX (key `"23"`) and chrY (key `"24"`) counts for male samples
to level with autosomal diploid coverage. This is a no-op for female
samples. Mirrors `overall_tools.gender_correct()` in upstream
WisecondorX.

## Usage

``` r
.gender_correct(sample, gender)
```

## Arguments

- sample:

  Named list of integer/numeric vectors keyed by chromosome.

- gender:

  Character; `"M"` for male, `"F"` for female.

## Value

Modified `sample` with doubled chrX/Y if male.
