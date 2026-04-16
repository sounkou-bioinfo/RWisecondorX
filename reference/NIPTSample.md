# Abstract base class for NIPTeR bin-count samples

Never instantiated directly. Subclasses are `CombinedStrandsSample`
(reads are summed across strands) and `SeparatedStrandsSample` (forward
and reverse counts stored independently).

## Usage

``` r
NIPTSample(
  sample_name = character(0),
  binsize = integer(0),
  correction = NIPTCorrectionRecord()
)
```
