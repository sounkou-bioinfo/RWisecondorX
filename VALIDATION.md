# RWisecondorX Validation Record

This document satisfies [rewrites.bio](https://rewrites.bio) principles
**4.1** (test and benchmark with real data) and **4.3** (pin versions
and document). It records what has been validated, how to reproduce it,
and what remains unvalidated.

------------------------------------------------------------------------

## Upstream Versions

| Tool                  | Version | Source                                  |
|-----------------------|---------|-----------------------------------------|
| WisecondorX           | 1.2.x   | `bioconda::wisecondorx` via `condathis` |
| NIPTeR                | 1.3.x   | CRAN `NIPTeR`                           |
| htslib (via Rduckhts) | 1.21    | bundled in `Rduckhts`                   |

> **Policy**: when upstream versions change, rerun the relevant
> conformance checks below and update this file.

------------------------------------------------------------------------

## What Has Been Validated

### 1. WisecondorX convert semantics in the native R path

**Status**: Validated ✓

**Claim**: the native convert layer reproduces the upstream WisecondorX
convert semantics that matter for bin counts:

- streaming duplicate handling (`larp` / `larp2`)
- internal `MAPQ >= 1` behavior
- chromosome-keyed sample structure used by the native `rwisecondorx`
  path

**Scope**:

- exact convert semantics are implemented in
  [`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
  [`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md),
  and
  [`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md)
- upstream-facing CLI behavior is preserved in
  [`wisecondorx_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_convert.md)
- current real-data conformance workflow compares prepared native BED.gz
  test cases against prepared upstream NPZ test cases through
  `inst/scripts/wisecondorx_conformance.R`

**Reproduction**:

``` bash
Rscript inst/scripts/wisecondorx_conformance.R \
  --rwisecondorx-ref /path/to/rwisecondorx_reference.rds \
  --wisecondorx-ref /path/to/wisecondorx_reference.npz \
  --case-manifest /path/to/cases.tsv \
  --sample-binsize 100000 \
  --out-dir /path/to/report_dir
```

`cases.tsv` must contain:

``` tsv
sample  bed npz
sample1 /path/to/native/sample1.bed.gz  /path/to/upstream/sample1.npz
sample2 /path/to/native/sample2.bed.gz  /path/to/upstream/sample2.npz
```

The script runs both `predict` paths and compares the resulting per-bin
`_bins.bed` files.

### 2. NIPTeR statistical layer formula verification

**Status**: Validated ✓

**Claim**: the core NIPTeR statistical functions match their intended
formulas:

- [`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md)
- [`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)
- [`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md)
- [`nipter_regression()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md)

**Method**:

- inline reference computations in
  `inst/tinytest/test_nipter_conformance.R`
- broader functional coverage in `inst/tinytest/test_nipter_stats.R`
- sex-model and gaunosome coverage in `inst/tinytest/test_nipter_sex.R`

These tests are mostly synthetic/unit-level by design; they verify
formulas and control-flow, not cohort realism.

**Reproduction**:

``` bash
make test-nipter
```

### 3. Real-BAM structural validation for the NIPTeR path

**Status**: Validated on internal real data ✓

**Claim**: the NIPTeR preprocessing and object-building path works on
real BAMs when driven from a manifest:

- binning
- BED round-trip
- control-group construction
- basic cohort structural checks

**Test file**:

- `inst/tinytest/test_nipter_real_manifest.R`

**Default internal manifest**:

- `/mnt/data/BixCTF/NiptSeqNeo/all_bam_list_sample_500.txt`
- default limit: `50`

**Reproduction**:

``` bash
make test-real-nipter
```

or explicitly:

``` bash
THREADS=20 \
NIPTER_REAL_BAM_ENABLE=1 \
NIPTER_REAL_BAM_LIST=/mnt/data/BixCTF/NiptSeqNeo/all_bam_list_sample_500.txt \
NIPTER_REAL_BAM_LIMIT=50 \
R -e "tinytest::run_test_file('inst/tinytest/test_nipter_real_manifest.R')"
```

### 4. SeqFF support data and wrapper integration

**Status**: Validated ✓

**Claim**: the package-local SeqFF support data loads correctly and
[`seqff_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/seqff_predict.md)
is wired into the package as a usable helper.

**Test file**:

- `inst/tinytest/test_seqff.R`

**Reproduction**:

``` bash
make test-seqff
```

### 5. Reference/control-group build CLIs match the current artifact boundaries

**Status**: Validated structurally ✓

**Claim**:

- `inst/scripts/convert_sample.R` differentiates:
  - `rwisecondorx`
  - `wisecondorx`
  - `nipter`
- `inst/scripts/build_reference.R` differentiates:
  - `--mode rwisecondorx` from BED.gz
  - `--mode wisecondorx` from NPZ
  - `--mode nipter` from BED.gz
- `inst/scripts/preprocess_cohort.R` stages real cohorts into reusable
  artifact directories instead of mixing conversion, reference building,
  and scoring in one script
- the same preprocessing script now emits SeqFF fetal-fraction estimates
  as a cohort-side preprocessing artifact

This is structural validation of the workflow boundary, not numeric
cross-implementation equivalence by itself.

------------------------------------------------------------------------

## What Is Not Yet Fully Validated

### A. End-to-end native vs upstream WisecondorX agreement on a full real NIPT cohort

**Status**: In progress ❌

The package now has the correct artifact-based conformance workflow
(`preprocess_cohort.R` → `build_reference.R` →
`wisecondorx_conformance.R`), but this document does not yet record a
completed full-cohort result with:

- exact prepared sample set
- exact native and upstream references
- exact summary output from `wisecondorx_conformance_summary.tsv`

What remains to record:

- number of test cases
- exact match rate on `_bins.bed`
- any persistent per-bin or per-chromosome discrepancies

### B. Clinical positive-case validation on real trisomy samples

**Status**: Not yet validated ❌

The package is currently set up for real-data preprocessing and
control-group building, but this file does not yet contain validation on
confirmed positive NIPT cases for:

- T21
- T18
- T13
- sex chromosome abnormalities

Until those are run and recorded, the package should be treated as
structurally validated on real data, not clinically validated.

### C. Native vs upstream WisecondorX segment and aberration agreement

**Status**: Not yet recorded ❌

The current conformance CLI compares per-bin prediction outputs. This is
the right first acceptance gate, but this file does not yet record:

- segment-level agreement
- aberration-call agreement
- disagreement patterns on noisy real samples

Those should be added once the positive-case cohort is available.

### D. GC correction numeric equality with upstream NIPTeR

**Status**: Expected divergence; not a target ❌

[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)
computes GC content on the fly or from precomputed package-side GC
tables. Upstream `NIPTeR::gc.correct()` uses bundled GC data and
different implementation details. Exact numeric equality is not expected
and is not currently a validation target.

### E. Sex-stratified X/Y scoring in the NIPTeR clinical path

**Status**: Not yet complete ❌

The package still needs broader real-data validation, and in some cases
further implementation work, for sex-aware X/Y scoring and related
gaunosome reporting in the clinical NIPT setting.

------------------------------------------------------------------------

## What Was Removed From The Validation Story

The following older validation claims are no longer part of the current
package story and should not be cited:

- bundled BAM-fixture conformance as the main validation path
- `make nipter-fixture`
- synthetic cohort / `generate_cohort()` smoke validation
- `inst/scripts/real_data_conformance.R`
- `inst/tinytest/test_wisecondorx_e2e.R`

The package has been refocused on real BAM cohorts and reusable prepared
artifacts instead.

------------------------------------------------------------------------

## How To Run The Current Validation Slices

``` bash
# fast package-local checks
make test-fast

# BED / NPZ / I/O checks
make test-io

# NIPTeR statistical and control-group checks
make test-nipter

# native WisecondorX unit checks
make test-rwisecondorx

# SeqFF support
make test-seqff

# internal real-BAM structural validation for NIPTeR
make test-real-nipter
```

For real cohort preprocessing and conformance:

``` bash
# 1. preprocess BAMs into native BED.gz, optional upstream NPZ, and NIPTeR BED.gz
Rscript inst/scripts/preprocess_cohort.R \
  --bam-list /path/to/bams.txt \
  --out-root /path/to/workdir \
  --fasta /path/to/reference.fa \
  --wcx-write-npz

# 2. build the native and upstream references from prepared artifacts
Rscript inst/scripts/build_reference.R \
  --mode rwisecondorx \
  --bed-dir /path/to/workdir/wcx_beds \
  --out /path/to/workdir/refs/rwisecondorx_reference.rds

Rscript inst/scripts/build_reference.R \
  --mode wisecondorx \
  --npz-dir /path/to/workdir/wcx_npz \
  --out /path/to/workdir/refs/wisecondorx_reference.npz

# 3. compare both predict paths on prepared test cases
Rscript inst/scripts/wisecondorx_conformance.R \
  --rwisecondorx-ref /path/to/workdir/refs/rwisecondorx_reference.rds \
  --wisecondorx-ref /path/to/workdir/refs/wisecondorx_reference.npz \
  --case-manifest /path/to/cases.tsv \
  --out-dir /path/to/workdir/wisecondorx_conformance
```

------------------------------------------------------------------------

## AI Assistance Disclosure

This package was written with the assistance of AI coding agents.
Correctness is validated by comparison against upstream tools, formula
checks, and real-data structural tests, not by generated code alone.

Where validation is incomplete, this file states that explicitly.
