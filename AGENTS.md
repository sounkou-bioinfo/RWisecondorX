# RWisecondorX Agent Guide

This guide is for AI agents working in `RWisecondorX`.

## 1. Project Summary

`RWisecondorX` is an R toolkit for copy-number analysis and trisomy
prediction in non-invasive prenatal testing (NIPT). It aims to provide:

- **Upstream conformance** where behavior is documented and testable
- **High performance** via `Rduckhts` + DuckDB SQL instead of per-read
  Python/Rsamtools loops
- **Interoperable intermediate files** via bgzipped, tabix-indexed
  TSV/BED-like files
- **No required Python runtime** for the native R implementation

The package ports two upstream systems:

- **WisecondorX** — cfDNA copy-number aberration detection
- **NIPTeR** — trisomy prediction via chromosome-fraction statistics

## 2. Core Working Principles

1.  **Read existing code and tests first.**
2.  **Keep WisecondorX and NIPTeR code separated by file family**:
    - `nipter_*.R`
    - `rwisecondorx_*.R` / `wisecondorx_*.R`
    - shared binning only in `R/convert.R`
3.  **Prefer explicit contracts over defensive magic.** This project is
    still pre-1.0; do not add awkward internal backwards-compatibility
    layers just to preserve accidental behavior.
4.  **Fail loudly on invalid inputs, missing artifacts, impossible
    parameter combinations, or degenerate statistical states.** Do not
    silently substitute defaults, alternate branches, uniform weights,
    or fallback outputs unless that behavior is intentionally documented
    as upstream-conformance logic.
5.  **Keep runtime logic and conformance tooling separate.**
6.  **Prefer upstream conformance unless an intentional divergence is
    already documented.**
7.  **Design steps to be as stateless and deterministic as possible.**
    Functions and scripts should clearly define inputs, outputs,
    assumptions, and artifact contracts.
8.  **Rich QC is mandatory.** Prefer exposing metrics, intermediate
    summaries, and visualizations rather than hiding weak or degenerate
    results behind permissive fallbacks.
9.  **Preserve CRAN-friendly behavior.** Do not make package loading
    depend on an external Python install.
10. **When behavior is unclear, inspect upstream references before
    changing code.**
11. **If filesystem access or permissions are unclear, stop and ask
    before assuming.**

## 3. Where To Look When Behavior Is Unclear

### Upstream references

- WisecondorX upstream: `.sync/WisecondorX/`
- NIPTeR upstream: `.sync/NIPTeR/`

### Conformance tools

- Use the official bioconda `wisecondorx` package via `condathis`
- Do **not** vendor Python source into this repo
- Do **not** create or rely on an `inst/python/` directory

## 4. Current Architecture

## 4.1 Shared binning layer

### `R/convert.R`

Shared BAM/CRAM binning engine used by both WisecondorX and NIPTeR.

Implemented features: -
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md) -
DuckDB/SQL read counting - supports `binsize`, `mapq`, `require_flags`,
`exclude_flags` - duplicate modes: `"streaming"`, `"flag"`, `"none"` -
optional strand separation - CRAM support via `reference` -
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md) -
writes bgzipped + tabix-indexed 4-column TSV/BED-like files

## 4.2 WisecondorX layers

### CLI interoperability

- `R/wisecondorx_cli.R`
  - [`wisecondorx_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_convert.md)
  - [`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md)
  - [`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)
- thin `condathis` wrappers around official upstream CLI

### NPZ compatibility

- `R/wisecondorx_npz.R`
  - [`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md)
- kept for upstream CLI conformance and interoperability, not native
  runtime

### Native R implementation

Pure R/Rcpp implementation of WisecondorX `newref` + `predict`.

Key files: - `R/rwisecondorx_utils.R` - shared utilities - gender model
training - mask construction - PCA helpers -
[`scale_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/scale_sample.md) -
`R/rwisecondorx_newref.R` - native reference building - A/F/M partition
construction - KNN/null-ratio orchestration -
`R/rwisecondorx_predict.R` - native prediction pipeline - rescaling,
gender prediction, autosome/gonosome normalization, CBS, aberration
calling - `R/rwisecondorx_cbs.R` - DNAcopy / ParDNAcopy wrapper -
segment z-scores and summary statistics - `R/rwisecondorx_output.R` -
BED/statistics output writers

Compiled code: - `src/knn_reference.cpp` - `knn_reference_cpp()` -
`null_ratios_cpp()` - `src/RcppExports.cpp`, `R/RcppExports.R` -
`src/Makevars`, `src/Makevars.win` - `R/zzz_rcpp.R`

## 4.3 NIPTeR layers

### Binning

- `R/nipter_bin.R`
  - [`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
  - [`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md)

### Statistics and reference operations

- `R/nipter_control.R`
  - control groups, diagnostics, matching, matrix SSD helpers
- `R/nipter_gc.R`
  - GC precompute and correction
- `R/nipter_chi.R`
  - chi-squared correction
- `R/nipter_score.R`
  - z-scores and NCV
- `R/nipter_regression.R`
  - regression models
- `R/nipter_sex.R`
  - mclust-based sex prediction and Y-unique model

### Additional sex / gaunosome reporting layers

- `R/nipter_reference.R`
- `R/nipter_gaunosome_models.R`
- `R/nipter_gaunosome_score.R`
- `R/nipter_gaunosome_report.R`
- `R/nipter_plot.R`
- `R/nipter_control_qc.R`
- `R/nipter_sex_score.R`

## 5. File Format Contracts

Our so-called “BED.gz” files are actually **BGZF-compressed,
tabix-indexed TSV files** using BED-like genomic coordinates.

Important rules: - 0-based, half-open coordinates - chromosome names
have **no `chr` prefix** - these are **not strict 3-column BED files** -
always read them with: - `read_tabix()` -
[`Rduckhts::rduckhts_tabix_multi()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_tabix_multi.html) -
never use generic BED parsers that assume strict BED semantics

### Supported TSV.bgz layouts

- **WisecondorX 4-column**: `chrom`, `start`, `end`, `count`
- **NIPTeR 5-column**: `chrom`, `start`, `end`, `count`,
  `corrected_count`
- **NIPTeR 9-column**: `chrom`, `start`, `end`, `count`, `count_fwd`,
  `count_rev`, `corrected_count`, `corrected_fwd`, `corrected_rev`
- **GC table 5-column**: `chrom`, `start`, `end`, `pct_gc`, `seq_len`

### Readers

- [`bed_to_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_sample.md)
  reads WisecondorX 4-column files
- [`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md)
  reads NIPTeR 5- or 9-column files
- [`nipter_control_group_from_beds()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_control_group_from_beds.md)
  loads a whole BED directory in one pass

### Reference persistence

- `WisecondorXReference` is currently **RDS-only**
- do not try to flatten PCA / KNN / null-ratio matrices into a single
  TSV

## 6. Source Layout Expectations

### Core package files

- `R/convert.R`
- `R/nipter_bin.R`
- `R/nipter_gc.R`
- `R/nipter_score.R`
- `R/nipter_chi.R`
- `R/nipter_control.R`
- `R/nipter_regression.R`
- `R/nipter_sex.R`
- `R/rwisecondorx_utils.R`
- `R/rwisecondorx_newref.R`
- `R/rwisecondorx_predict.R`
- `R/rwisecondorx_cbs.R`
- `R/rwisecondorx_output.R`
- `R/bed_reader.R`
- `R/synthetic_cohort.R`
- `R/simulate_trisomy_cohort.R`
- `R/wisecondorx_cli.R`
- `R/wisecondorx_npz.R`
- `R/aaa.R`

### Scripts

- `inst/scripts/convert_sample.R`
- `inst/scripts/build_reference.R`
- `inst/scripts/precompute_gc.R`
- `inst/scripts/make_cohort.R`

### Tests and data

- `inst/tinytest/`
- `inst/extdata/`
- `inst/extdata/write_compat_npz.py`

## 7. Testing, Build, and Documentation Rules

## 7.1 Testing

- use **tinytest**
- keep one test file per feature family where practical
- update tests when R-facing behavior changes
- tests requiring Python / `reticulate` / external tools must skip
  cleanly

### Relevant environment variables

- `WISECONDORX_TEST_BAM`
- `NIPTER_CONFORMANCE_BAM`

## 7.2 Build workflow

Use make targets: - `make rd` - `make test` - `make readme` -
`make fixtures`

Keep builds deterministic and non-interactive.

## 7.3 Documentation

- `README.Rmd` is the editable source
- regenerate `README.md` with `make readme`
- roxygen generates Rd + `NAMESPACE`
- do **not** edit `README.md` or `NAMESPACE` by hand

## 8. Dependency Guidance

### 8.1 `Rduckhts`

`RWisecondorX` tracks the **edge** of `Rduckhts`, not just CRAN.

If you need to install/update it for development, prefer r-universe:

``` r
install.packages(
  "Rduckhts",
  repos = c(
    "https://rgenomicsetl.r-universe.dev",
    "https://cloud.r-project.org"
  )
)
```

Do not “fix” behavior by lowering expectations to an older CRAN build.

### 8.2 `this.path`

Use [`this.path`](https://CRAN.R-project.org/package=this.path) as the
preferred solution for locating the executing script when relative
script paths matter.

Project rule: - prefer `this.path` over ad hoc path heuristics tied to
[`commandArgs()`](https://rdrr.io/r/base/commandArgs.html), IDE-specific
globals, working-directory assumptions, or hand-rolled source-path
detection - especially for `inst/scripts/*.R`, use explicit path
discovery and then derive sibling resources relative to the script
location - avoid brittle
[`getwd()`](https://rdrr.io/r/base/getwd.html)-relative behavior in
scripts

Design intent: - define inputs, outputs, and resource roots explicitly -
make path resolution reproducible across `Rscript`, IDEs, notebooks, and
sourced execution - remove hidden environment coupling where practical

Because this project is still pre-1.0, do not carry awkward internal
path-compatibility shims just to preserve old incidental behavior.

## 9. Convert-Step Conformance Notes

The WisecondorX convert step replicates upstream `larp` / `larp2`
streaming dedup behavior in SQL.

Key SQL pattern:

``` sql
LAG(pos) OVER (ORDER BY file_offset)
LAST_VALUE(CASE WHEN is_paired != 0 THEN pnext END IGNORE NULLS)
  OVER (ORDER BY file_offset ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING)
```

Important semantics: - improper pairs are invisible to both `larp` and
`larp2` - unpaired reads update `larp` but not `larp2` - `larp` is never
reset by chromosome - binning uses integer division matching Python
`pos // binsize` - `duckhts` POS is 1-based; convert to 0-based before
binning

## 10. Flag Filtering API

These functions expose samtools-style `-f` / `-F` filters: -
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md) -
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md) -
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md) -
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md)

### Parameters

- `require_flags`
  - keep reads where `(flag & require_flags) == require_flags`
- `exclude_flags`
  - drop reads where `(flag & exclude_flags) != 0`

Example:

``` r
nipter_bin_bam("sample.dm.bam", mapq = 40L, exclude_flags = 1024L)
```

## 11. NIPTeR-Specific Caveats

## 11.1 Upstream NIPTeR treats paired-end data like concatenated single-end reads

`bin_bam_sample()` counts both mates independently by strand and then
bins them. This means: - paired-end libraries produce roughly ~2x counts
relative to fragment counting - single-end and paired-end references are
not safely interchangeable

Our default behavior preserves that compatibility.

## 11.2 Upstream NIPTeR positional dedup is not duplicate-flag filtering

Upstream [`unique()`](https://rdrr.io/r/base/unique.html) on
positions: - collapses same-position reads per chromosome + strand - is
not equivalent to `-F 1024` - can merge biologically distinct reads
sharing a start site

## 11.3 Upstream NIPTeR scanBam bugs

### Bug 1: unmapped-read split mismatch

Upstream code can break when BAMs include unmapped reads because the
split factor length does not match the read vector length.

### Bug 2: `unique()` is flag-ignorant

A duplicate at a unique position survives; two distinct reads at the
same position collapse.

### Conformance implication

Exact NIPTeR conformance requires a BAM with: - no unmapped reads - no
same-position strand collisions

That is why this repo bundles
`inst/extdata/nipter_conformance_fixture.bam`.

## 12. Native WisecondorX Critical Notes

## 12.1 KNN index semantics

Our Rcpp KNN implementation stores **global 1-based indexes into the
full masked array**.

Upstream Python stores **local indexes into the leave-one-chromosome-out
array**.

### Consequences

| Context                | Upstream Python                          | Our R/Rcpp                    |
|------------------------|------------------------------------------|-------------------------------|
| Stored indexes         | Local                                    | Global                        |
| Predict normalization  | Uses local indexes                       | Converts global to local      |
| Null-ratio computation | Uses local indexes incorrectly as global | Uses global indexes correctly |

### Global-to-local conversion in `.normalize_once()`

- if `g < chr_start`, local index is `g`
- if `g > chr_cum`, local index is `g - n_chr`
- indexes inside `[chr_start, chr_cum]` should never exist

Do **not** change index semantics unless you update both: -
[`.normalize_once()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-normalize_once.md)
in `R/rwisecondorx_predict.R` - `null_ratios_cpp()` in
`src/knn_reference.cpp`

## 12.2 Known upstream WisecondorX deviations intentionally replicated

These are deliberate conformance quirks.

- **CPA +1 overcount**
  - upstream computes segment length as `end - start + 1` for half-open
    segments
  - replicated in
    [`.get_cpa()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-get_cpa.md)
- **Whole-chromosome last-bin drop**
  - upstream computes whole-chromosome z-scores with
    `end = bins_per_chr - 1`
  - replicated in
    [`.compute_statistics()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/dot-compute_statistics.md)

## 13. Implementation Status Snapshot

## 13.1 Implemented and stable

### Shared / infrastructure

- shared DuckDB binning
- CRAM support
- tabix-indexed BED/TSV outputs
- fixture generation
- SRA metadata helpers

### WisecondorX

- CLI wrappers via `condathis`
- NPZ compatibility path
- native `newref` + `predict`
- CBS integration
- native output writers
- Rcpp/OpenMP KNN + null-ratio code

### NIPTeR

- binning
- control groups
- GC correction
- chi correction
- z-score / NCV
- regression
- sex models including Y-unique ratios
- BED round-trips

### Cohort simulation

- synthetic cohort generation
- cohort pipeline tests

## 13.2 Open / remaining work

### WisecondorX

- end-to-end conformance vs upstream Python on real multi-sample
  datasets
- interoperable non-RDS WisecondorX reference format, if ever needed
- consistent threading-budget propagation through pipeline scripts
- real-data validation of beta-mode aberration calling
- continued hardening of explicit artifact contracts, QC summaries, and
  visualization outputs for cohort workflows

### NIPTeR

- additional clinical sex-stratified X/Y NCV and regression work as
  needed
- continued conformance validation against real-world fixtures
- continued expansion of QC/reporting surfaces rather than burying edge
  cases in permissive behavior

## 14. Agent Policy on Fallbacks, Contracts, and QC

When changing package code or scripts:

- **Do not hide bugs with fallbacks.**
  - Bad: silently switching to another reference branch, inventing
    weights, swallowing missing files, or returning cosmetically valid
    output from broken internals.
  - Good: error clearly, report diagnostics, and emit QC artifacts that
    explain what failed or degraded.
- **Make assumptions explicit.**
  - expected bin size
  - chromosome naming assumptions
  - required metadata columns
  - expected sample/reference compatibility
  - thread and memory assumptions when relevant
- **Define outputs explicitly.**
  - primary result object
  - sidecar QC tables
  - plots
  - serialized artifacts
  - metadata needed for downstream steps
- **Keep steps stateless where possible.**
  - inputs should be passed in or declared
  - outputs should be written to declared locations
  - avoid hidden dependence on working directory, IDE state, or mutable
    global options
- **Prefer rich QC over permissive success.**
  - expose low ref-bin counts
  - expose masked-bin fractions
  - expose sex-model confidence / ambiguity
  - expose normalization and segmentation warnings
  - write plots and summary tables when a workflow already has a QC
    surface

## 15. Copyright and Attribution

- NIPTeR authors Dirk de Weerd and Lennart Johansson are listed as `cph`
  in `DESCRIPTION`
- WisecondorX authors are listed as `ctb`
- preserve attribution for any additional upstream ports
