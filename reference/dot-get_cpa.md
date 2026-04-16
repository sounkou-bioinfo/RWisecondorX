# Copy number profile abnormality (CPA)

Copy number profile abnormality (CPA)

## Usage

``` r
.get_cpa(results_c, binsize)
```

## Note

**Upstream bug replicated for conformance (CPA +1 overcount).** Upstream
Python (`overall_tools.py:146`) computes segment length as
`segment[2] - segment[1] + 1`. Given that CBS returns half-open
`[start, end)` intervals, the correct length is `end - start`, so the
`+ 1` overcounts by one bin per segment. We replicate this exactly
(`end - start + 1L`) for conformance with Python WisecondorX output.
