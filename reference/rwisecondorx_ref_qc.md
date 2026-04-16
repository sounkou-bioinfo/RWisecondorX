# QC a native WisecondorX reference

Native R implementation of the upstream WisecondorX `ref_qc.py`
heuristics for `WisecondorXReference` objects. Inspects within-sample
reference-bin distances, compatibility between autosomal and gonosomal
mask prefixes, and (for male references) chrY usability metrics.

## Usage

``` r
rwisecondorx_ref_qc(reference, min_ref_bins = 150L, output_json = NULL)
```

## Arguments

- reference:

  A `WisecondorXReference` object from
  [`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md).

- min_ref_bins:

  Integer; minimum number of usable reference bins per target bin before
  the report flags a warning. Default `150L`, matching the upstream
  Python QC heuristic.

- output_json:

  Optional path for writing the QC report as JSON. When supplied,
  requires the `jsonlite` package.

## Value

A named list with class `"WisecondorXReferenceQC"` containing:

- overall_verdict:

  `"PASS"`, `"WARN"`, or `"FAIL"`.

- worst_severity:

  Integer severity code: `0L` pass, `1L` warn, `2L` fail.

- compat_issues:

  Character vector of autosomal-prefix compatibility issues between
  autosomal and gonosomal sub-references.

- metrics:

  Named list of per-branch metrics (`A`, `F`, `M`) including verdict,
  message, mean/std of per-bin mean distances, outlier counts, and
  low-reference counts. Male reports also include chrY metrics.

## Details

The function returns a structured QC report and can optionally write the
report as JSON for later inspection.

## See also

[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md)
