# Scale a binned sample to a new bin size

Rescales per-chromosome bin-count vectors by aggregating bins. The new
bin size must be a positive integer multiple of the original. This is
the R equivalent of `overall_tools.scale_sample()` in the upstream
WisecondorX Python package.

## Usage

``` r
scale_sample(sample, from_size, to_size)
```

## Arguments

- sample:

  Named list of integer vectors keyed by chromosome (`"1"`–`"24"`), as
  returned by
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md).

- from_size:

  Current bin size in base pairs.

- to_size:

  Target bin size in base pairs. Must be an integer multiple of
  `from_size`.

## Value

Named list structured like `sample` but with bins aggregated to the new
size.

## Examples

``` r
if (FALSE) { # \dontrun{
bins_5k  <- bam_convert("sample.bam", binsize = 5000L)
bins_100k <- scale_sample(bins_5k, from_size = 5000, to_size = 100000)
} # }
```
