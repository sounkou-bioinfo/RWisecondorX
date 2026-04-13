# Train a gender model using GMM on Y-fractions

Fits a 2-component Gaussian mixture model to the Y-chromosome read
fractions of a set of reference samples. Uses
[`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html)
as the R-native replacement for `sklearn.GaussianMixture`. The local
minimum of the fitted density is used as the male/female cutoff.

## Usage

``` r
.train_gender_model(samples, yfrac = NULL)
```

## Arguments

- samples:

  List of sample objects.

- yfrac:

  Optional numeric; if given, used as the manual cutoff instead of the
  GMM-derived one.

## Value

A list with:

- genders:

  Character vector of `"M"` / `"F"` per sample.

- cutoff:

  Numeric; the Y-fraction cutoff used.

## Details

Mirrors `newref_tools.train_gender_model()` in upstream WisecondorX.
