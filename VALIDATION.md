# RWisecondorX Validation Record

This document satisfies [rewrites.bio](https://rewrites.bio) principles
**4.1** (test and benchmark with real data) and **4.3** (pin versions
and document). It records what has been validated, the exact version
each claim applies to, how to reproduce the result, and what remains
unvalidated.

------------------------------------------------------------------------

## Upstream Versions

| Tool                  | Version | Source                                |
|-----------------------|---------|---------------------------------------|
| WisecondorX           | 1.2.x   | `bioconda::wisecondorx` via condathis |
| NIPTeR                | 1.3.x   | CRAN `NIPTeR`                         |
| htslib (via Rduckhts) | 1.21    | bundled in Rduckhts                   |

> **Policy**: When a new upstream version is released, re-run the
> conformance tests listed below and update this document. Breaking
> changes (format changes, algorithm updates) trigger a version bump in
> RWisecondorX; bug-fix updates are re-validated silently.

------------------------------------------------------------------------

## What Has Been Validated

### 1. `bam_convert()` — bin-for-bin agreement with Python WisecondorX

**Status**: Validated ✓  
**Claim**:
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
with `rmdup="streaming"` produces bin counts identical (bin-for-bin,
zero mismatches) to `wisecondorx convert` on the same BAM. The current
implementation delegates counting to the native
[`Rduckhts::rduckhts_bam_bin_counts()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_bam_bin_counts.html)
kernel and reshapes that output into the WisecondorX chromosome-keyed
format.

**Data**: `HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam`  
**Validated bins**: 25,115 non-zero bins on chromosome 11  
**Mismatches**: 0  
**Test file**: `inst/tinytest/test_integration.R`

**Command used** (reference):

``` bash
# Python (via condathis)
wisecondorx convert HG00106.chrom11.bam HG00106.npz --binsize 5000

# R
bam_convert("HG00106.chrom11.bam", binsize = 5000L, rmdup = "streaming")
```

**How to reproduce**:

``` r
# Requires: condathis, reticulate, numpy, wisecondorx bioconda env
tinytest::run_test_file("inst/tinytest/test_integration.R")
```

------------------------------------------------------------------------

### 2. NIPTeR statistical layer — formula verification

**Status**: Validated against inline reference ✓  
**Claim**:
[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md),
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)
reproduce the NIPTeR formulas exactly (to floating-point precision).

**Method**: Inline reference computations in `test_nipter_conformance.R`
(Arm A) replicate the NIPTeR source formulas directly in R, then compare
against our function output. No external package dependency.

**Key assertions**: - Z-score:
`(sample_frac - mean(ctrl_fracs)) / sd(ctrl_fracs)` matches to
`< 1e-10` - Chi correction: non-overdispersed bins unchanged;
overdispersed bins not inflated - NCV: denominator set does not include
focus chromosome - Regression: returns numeric z-score with predictor
names

**Test file**: `inst/tinytest/test_nipter_conformance.R` (Arm A, always
runs)

**How to reproduce**:

``` r
tinytest::run_test_file("inst/tinytest/test_nipter_conformance.R")
```

------------------------------------------------------------------------

### 3. NIPTeR statistical layer — cross-check against NIPTeR R package

**Status**: Conditional ⚠ (requires NIPTeR package + conformance BAM)  
**Claim**:
[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md)
and
[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)
agree with `NIPTeR::chromosomal.zscore()` and `NIPTeR::chi.correct()`
within tolerance 0.01 on the same binned data.

**Constraint**: The conformance BAM must have: 1. No unmapped reads
(pre-filter with `samtools view -F 4`) 2. No two reads sharing the same
start position on the same chromosome per strand

These constraints avoid two known bugs in NIPTeR (see AGENTS.md for
details).

**How to reproduce**:

``` bash
# Build multi-chromosome synthetic fixture (satisfies both constraints)
make nipter-fixture
# Then run conformance tests with the fixture
NIPTER_CONFORMANCE_BAM=inst/extdata/nipter_conformance_fixture.bam make test
# Or use a real pre-filtered whole-genome BAM:
NIPTER_CONFORMANCE_BAM=/path/to/prefiltered.bam make test
```

------------------------------------------------------------------------

### 4. Trisomy detection — synthetic smoke test cohort

**Status**: Smoke-tested ✓  
**Claim**: The native R pipeline (`rwisecondorx_newref` +
`rwisecondorx_predict`) detects T21, T18, and T13 on the package’s
synthetic smoke-test cohort and produces zero aberration calls on the
chosen euploid negative control in that same synthetic setup.

**Data**: 50-sample synthetic cohort from
[`generate_cohort()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/generate_cohort.md)
(compressed coordinates, `COMPRESSED_BINSIZE = 100`)

**Key assertions**: - T21 sample: `pred$aberrations` contains at least
one `chr=="21"` + `type=="gain"` row - T18 sample: at least one
`chr=="18"` + `type=="gain"` row - T13 sample: at least one
`chr=="13"` + `type=="gain"` row - First euploid sample:
`nrow(pred$aberrations) == 0`

**Important limitation**:
[`generate_cohort()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/generate_cohort.md)
is a fast synthetic generator for package testing. It is not a
calibrated clinical simulator and does not replace real-data conformance
or performance validation.

**Test files**: - `inst/tinytest/test_cohort_pipeline.R` (Z-score level
assertions) - `inst/tinytest/test_wisecondorx_e2e.R` (aberration call
level assertions)

**How to reproduce**:

``` r
make conformance
# or
tinytest::run_test_file("inst/tinytest/test_wisecondorx_e2e.R")
```

------------------------------------------------------------------------

## What Is NOT Yet Validated On Real Data

### A. `rwisecondorx_newref/predict` segment and aberration agreement vs Python

**Status**: Not validated ❌  
**Gap**: We have not run both the native R pipeline and the Python
WisecondorX pipeline on the same set of real NIPT BAMs and compared
their segment boundaries or aberration calls.

**Blocking**: Requires a cohort of ≥10 real NIPT BAMs (whole genome, not
exome) with at least one sample of known trisomy status.

**How to run** (once BAMs are available):

``` bash
export RWXCONF_CONTROL_BAMS="ctrl1.bam:ctrl2.bam:..."
export RWXCONF_TEST_BAM="test_t21.bam"
export RWXCONF_CPUS=8
Rscript inst/scripts/real_data_conformance.R
```

**Acceptance criteria**: - Bin correlation (R vs Python on convert
step): Pearson r \> 0.999 - Chr-level aberration call Jaccard (R vs
Python on same reference): ≥ 0.8 - Segment Jaccard on trisomy
chromosome: ≥ 0.7 (CBS is stochastic; use same seed)

Update this file with the results after running.

------------------------------------------------------------------------

### B. Beta mode (purity-based calling)

**Status**: Implemented, not validated ❌  
**Gap**: The `beta` parameter in
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md)
implements ratio-based aberration cutoffs
(`log2((ploidy ± beta/2) / ploidy)`) for somatic CNV analysis. The
formula is correct (matches upstream WisecondorX) but has not been
tested on real tumour samples with known purity.

**Not appropriate for NIPT**: NIPT pipelines should always use Z-score
calling (the default).

------------------------------------------------------------------------

### C. GC correction numeric agreement with NIPTeR

**Status**: Known divergence, documented ❌  
**Gap**:
[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)
computes GC content on-the-fly from a reference FASTA via
`rduckhts_fasta_nuc()`. `NIPTeR::gc.correct()` uses bundled GC tables
from `sysdata.rda`. These tables differ in precision and GC-binning
strategy. Numeric equality is **not expected or required**. Both
implementations produce the same qualitative correction (GC-biased bins
are downweighted), but the corrected values will differ.

------------------------------------------------------------------------

### D. Sex-stratified NCV and regression for X/Y chromosomes

**Status**: Not implemented ❌  
**Gap**: The clinical pipeline computes sex-stratified NCV denominators
and regression models for X and Y chromosomes (separate models for males
vs females). Not yet implemented.

------------------------------------------------------------------------

## NIPTeR Known Bugs (Upstream)

Two bugs in NIPTeR `bin_bam_sample()` constrain the conformance BAM
requirements:

1.  **Split-length mismatch**:
    `split(reads, droplevels(strands[strands != "*"]))` drops
    unmapped-read entries from the factor but not from `reads`, causing
    silent misbehaviour on BAMs with unmapped reads.

2.  **Implicit positional dedup**: `bin_reads()` calls
    [`unique()`](https://rdrr.io/r/base/unique.html) on positions per
    chromosome per strand before binning, silently treating two reads at
    the same start position as one. This is NOT equivalent to `-F 1024`
    (duplicate flag).

These are bugs in the upstream NIPTeR package, not in RWisecondorX. Our
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
is documented to match the upstream behaviour when the BAM satisfies the
above constraints (see AGENTS.md for full details).

------------------------------------------------------------------------

## How To Run All Conformance Tests

``` bash
# Synthetic tests (always run; no external data needed)
make install
make test

# WisecondorX E2E conformance (synthetic; Python arm conditional on condathis)
make conformance

# Build multi-chr NIPTeR fixture for Arm B tests
make nipter-fixture
NIPTER_CONFORMANCE_BAM=inst/extdata/nipter_conformance_fixture.bam make test

# Full real-data conformance (server; requires real BAMs)
export RWXCONF_CONTROL_BAMS="bam1.bam:bam2.bam:..."
export RWXCONF_TEST_BAM="test.bam"
Rscript inst/scripts/real_data_conformance.R
```

------------------------------------------------------------------------

## AI Assistance Disclosure

This package was written with the assistance of AI coding agents (Claude
Sonnet 4.6 / Claude Code). Correctness is validated by comparing output
against the upstream tools on a suite of synthetic and real sequencing
datasets — not by manual code review alone. The AI generated the
implementation; humans defined the validation criteria and verified the
results.

Where validation is incomplete (see “Not Yet Validated” section above),
this is explicitly documented per [rewrites.bio](https://rewrites.bio)
principle 2.3.
