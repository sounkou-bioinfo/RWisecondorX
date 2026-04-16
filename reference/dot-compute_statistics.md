# Compute per-chromosome statistics

Mirrors `predict_output._generate_chr_statistics_file()`.

## Usage

``` r
.compute_statistics(
  results_r,
  results_w,
  results_c,
  results_nr,
  bins_per_chr,
  binsize,
  ref_gender,
  gender,
  n_reads
)
```

## Note

**Upstream bug replicated for conformance (last-bin drop).** Upstream
Python (`predict_output.py:210`) sets `end = bins_per_chr - 1` when
constructing whole-chromosome segments. Combined with `get_z_score`'s
`array[s:e]` slicing (Python half-open), this drops the last bin of
every chromosome from the whole-chromosome Z-score calculation. We
replicate this exactly (`end = bins_per_chr - 1L`) for conformance with
Python WisecondorX output.
