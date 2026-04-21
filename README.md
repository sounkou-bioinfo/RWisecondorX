
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWisecondorX

<!-- badges: start -->

[![R-CMD-check](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml)
[![rewrites.bio](https://rewrites.bio/badges/rewrites-bio.svg)](https://rewrites.bio)
<!-- badges: end -->

`RWisecondorX` is an R toolkit for copy number analysis and trisomy
prediction in non-invasive prenatal testing (NIPT), built on top of
`Rduckhts` and DuckDB. It covers two upstream method families:

  - [WisecondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX)
    for CNV detection
  - [NIPTeR](https://github.com/molgenis/NIPTeR) for
    chromosomal-fraction based trisomy prediction

The package is structured around three workflows:

  - **Native R WisecondorX**: `bam_convert()`, `bam_convert_bed()`,
    `rwisecondorx_newref()`, `rwisecondorx_predict()`
  - **Native R NIPTeR**: `nipter_bin_bam()`, GC correction, chi
    correction, Z-score, NCV, regression, sex prediction, reference QC,
    and plotting
  - **Optional upstream conformance**: `bam_convert_npz()` plus
    `wisecondorx_convert()`, `wisecondorx_newref()`, and
    `wisecondorx_predict()`

The native convert step is backed by the `bam_bin_counts(...)` kernel in
`Rduckhts`. It reproduces the upstream WisecondorX `larp` / `larp2`
streaming deduplication behaviour, while keeping all core processing
inside R and DuckDB. Both native analysis paths also use Rcpp/OpenMP in
their hot loops: WisecondorX for KNN reference finding and null-ratio
construction, and NIPTeR for control matching, NCV search, and
regression model search.

## Installation

`RWisecondorX` is developed against the current edge of `Rduckhts`, not
an older CRAN snapshot. Install `Rduckhts` from r-universe first, then
install `RWisecondorX`:

``` r
install.packages(
  "Rduckhts",
  repos = c(
    "https://rgenomicsetl.r-universe.dev",
    "https://cloud.r-project.org"
  )
)

# install.packages("pak")
pak::pak("sounkou-bioinfo/RWisecondorX")
```

`RWisecondorX` imports `Rduckhts`, `DBI`, and `duckdb`. `reticulate` is
needed only to write `.npz` files. `condathis` is needed only for the
upstream CLI wrappers.

## Workflow Overview

| Goal                                         | Primary functions                                             | Main output                           |
| -------------------------------------------- | ------------------------------------------------------------- | ------------------------------------- |
| Convert BAM/CRAM to native WisecondorX bins  | `bam_convert()`, `bam_convert_bed()`                          | in-memory sample or 4-column BED.gz   |
| Convert BAM/CRAM to native NIPTeR bins       | `nipter_bin_bam()`, `nipter_bin_bam_bed()`                    | `NIPTeRSample` or 5/9-column BED.gz   |
| Build a native WisecondorX reference         | `rwisecondorx_newref()`                                       | `WisecondorXReference`                |
| Predict native WisecondorX CNVs              | `rwisecondorx_predict()`                                      | `WisecondorXPrediction`               |
| Build a NIPTeR-style reference and QC bundle | `nipter_build_reference()`, `nipter_build_gaunosome_models()` | typed reference object                |
| Run upstream WisecondorX for conformance     | `bam_convert_npz()`, `wisecondorx_*()`                        | upstream `.npz` files and CLI outputs |

If you only want the native R workflow, you can ignore `reticulate`,
`condathis`, and the `wisecondorx_*()` wrappers.

## Convert BAM To Bin Counts

### WisecondorX-style bins

``` r
bins <- bam_convert(
  "sample.dm.bam",
  binsize = 100000L,
  mapq = 1L,
  rmdup = "streaming"
)
```

The `rmdup` argument controls duplicate handling:

  - `"streaming"` reproduces upstream WisecondorX `larp` / `larp2`
  - `"none"` skips deduplication entirely
  - `"flag"` respects the SAM duplicate flag (`0x400`)

### BED.gz output and round-trip

`bam_convert_bed()` writes a bgzipped, tabix-indexed BED-like TSV with
columns `chrom`, `start`, `end`, and `count`. This is the main file
format used by the native R WisecondorX workflow.

``` r
bed_file <- "sample.wisecondorx.bed.gz"

bam_convert_bed(
  bam     = "sample.dm.bam",
  bed     = bed_file,
  binsize = 100000L,
  rmdup   = "streaming"
)

sample_reloaded <- bed_to_sample(bed_file, binsize = 100000L)
```

For the NIPTeR side, `bed_to_nipter_sample()` reads a 5-column or
9-column BED.gz back into a `NIPTeRSample` object.

## Native WisecondorX

`rwisecondorx_newref()` and `rwisecondorx_predict()` are pure R/Rcpp
implementations of the reference-building and prediction stages. This is
the main production path when you want WisecondorX-style CNV calling
without a Python runtime dependency.

``` r
control_bams <- c("ctrl_01.dm.bam", "ctrl_02.dm.bam", "ctrl_03.dm.bam")
control_samples <- lapply(
  control_bams,
  function(bam) bam_convert(bam, binsize = 100000L, mapq = 1L, rmdup = "streaming")
)
names(control_samples) <- sub("\\.bam$", "", basename(control_bams))

ref <- rwisecondorx_newref(
  samples = control_samples,
  binsize = 100000L,
  nipt = TRUE,
  cpus = 4L
)

test_sample <- bam_convert("test.dm.bam", binsize = 100000L, mapq = 1L, rmdup = "streaming")
pred <- rwisecondorx_predict(test_sample, ref, zscore = 5, seed = 42L, cpus = 4L)
pred$aberrations
```

### Key prediction controls

The native predictor exposes the controls most commonly needed in
practice:

  - `beta` for purity-based ratio calling
  - `gender = "F"` / `"M"` to override the Y-fraction classifier
  - `parallel = TRUE` with `ParDNAcopy`, or `parallel = FALSE` for
    serial `DNAcopy::segment()`

`gender` does **not** have the same operational meaning in every
reference mode. Current native behaviour matches upstream WisecondorX:

  - with a non-NIPT reference, `gender` changes the operational
    gonosomal branch and applies male sex-chromosome leveling when
    forced to `"M"`
  - with an NIPT-mode reference, `gender` only overrides the
    reported/predicted gender label; gonosomal normalization still uses
    the female/NIPT branch internally

This asymmetry is upstream-compatible rather than ideal. A planned
native extension is to support an explicit, opt-in operational branch
override even in NIPT mode, to better align WisecondorX-side
sex-chromosome handling with the native NIPTeR flow.

``` r
pred_forced <- rwisecondorx_predict(
  test_sample,
  ref,
  gender = "F",
  beta = NULL,
  seed = 42L
)
```

### Outputs and QC

`write_wisecondorx_output()` writes:

  - `<outprefix>_bins.bed`
  - `<outprefix>_segments.bed`
  - `<outprefix>_aberrations.bed`
  - `<outprefix>_statistics.txt`

`rwisecondorx_ref_qc()` computes the upstream-style reference QC
heuristics and returns a structured report, with optional JSON output.

``` r
qc <- rwisecondorx_ref_qc(ref, output_json = "reference_qc.json")
write_wisecondorx_output(pred, outprefix = "results/sample_01")
```

The package does **not yet** provide native ggplot helpers for
`WisecondorXPrediction` objects. If you need the upstream PNG/PDF
plotting behaviour, use `wisecondorx_predict(..., plot = TRUE)` through
the optional CLI wrapper path.

### Building a reference from BED files

`rwisecondorx_newref()` can also load a directory of 4-column BED.gz
files produced by `bam_convert_bed()`:

``` r
ref <- rwisecondorx_newref(
  bed_dir = "bed_counts/",
  binsize = 100000L,
  nipt    = TRUE,
  cpus    = 4L
)
```

### Somatic beta mode

`rwisecondorx_predict()` supports a `beta` parameter for purity-based
aberration calling:

``` r
pred_somatic <- rwisecondorx_predict(
  sample    = tumour_sample,
  reference = ref,
  beta      = 0.4,
  seed      = 1L
)
```

When `beta` is supplied, ratio cutoffs replace Z-score thresholds:

``` text
gain cutoff: log2((ploidy + beta/2) / ploidy)
loss cutoff: log2((ploidy - beta/2) / ploidy)
```

This mode is for somatic CNV analysis, not routine NIPT trisomy calling.

## NIPTeR Workflow

The native NIPTeR layer covers binning, control-group construction, GC
correction, chi correction, autosomal scoring, sex prediction,
sex-chromosome models, and reference QC plotting. Like the native
WisecondorX path, it uses Rcpp-accelerated kernels in the expensive
search/matching steps rather than pure R loops throughout.

### Binning

`nipter_bin_bam()` produces `NIPTeRSample` objects compatible with all
downstream NIPTeR-style functions.

``` r
samples <- lapply(bam_files, nipter_bin_bam, binsize = 50000L)
```

### Control group construction

`nipter_as_control_group()` validates and deduplicates the reference
cohort. `nipter_diagnose_control_group()` reports chromosome-level
diagnostics. `nipter_match_control_group()` and `nipter_match_matrix()`
support sample matching and outlier QC.

``` r
cg <- nipter_as_control_group(samples)
diag <- nipter_diagnose_control_group(cg)
matched_cg <- nipter_match_control_group(test_sample, cg, n = 50L)
mm <- nipter_match_matrix(cg, cpus = 4L)
```

### GC correction

`nipter_gc_correct()` adjusts bin counts for GC-content bias. For
repeated workflows, `nipter_gc_precompute()` stores the GC table once as
a TSV.bgz file.

``` r
nipter_gc_precompute(
  fasta = "hg38.fa",
  binsize = 50000L,
  out = "gc_table.tsv.bgz"
)

cg <- nipter_gc_correct(cg, gc_table = "gc_table.tsv.bgz")
test_sample <- nipter_gc_correct(test_sample, gc_table = "gc_table.tsv.bgz")
```

### Chi-squared correction

`nipter_chi_correct()` identifies overdispersed bins across the control
group and applies the correction to both the sample and the control
cohort.

``` r
corrected   <- nipter_chi_correct(test_sample, cg, chi_cutoff = 3.5)
test_sample <- corrected$sample
cg          <- corrected$control_group
```

### Scoring

`nipter_z_score()` computes chromosomal fraction Z-scores:

``` r
z21 <- nipter_z_score(test_sample, cg, chromo_focus = 21)
z21$sample_z_score
```

`nipter_ncv_score()` computes normalized chromosome values:

``` r
ncv21 <- nipter_ncv_score(test_sample, cg, chromo_focus = 21)
ncv21$sample_ncv_score
ncv21$denominators
```

`nipter_regression()` uses forward stepwise regression:

``` r
reg21 <- nipter_regression(test_sample, cg, chromo_focus = 21)
reg21$models[[1]]$z_score
reg21$models[[1]]$predictors
```

### Reference building, QC, and plots

For production use, the package exposes a typed reference-building and
QC layer:

  - `nipter_build_reference()` packages the control group, sex models,
    and reference frame
  - `nipter_build_gaunosome_models()` adds sex-stratified X/Y NCV and
    regression models
  - `nipter_control_group_qc()` computes chromosome-, sample-,
    sex-model-, and optional bin-level QC summaries
  - `write_nipter_reference_plots()` writes the current plot bundle

<!-- end list -->

``` r
reference <- nipter_build_reference(
  cg,
  sample_sex = sample_sex,
  y_unique_ratios = y_unique_ratios
)
reference <- nipter_build_gaunosome_models(reference)

qc <- nipter_control_group_qc(
  reference$control_group,
  reference_model = reference,
  include_bins = TRUE
)

plot_paths <- write_nipter_reference_plots(
  qc,
  reference,
  outprefix = "nipter_reference_qc"
)
```

The current plotting helpers are:

  - `nipter_plot_qc_chromosomes()`
  - `nipter_plot_qc_samples()`
  - `nipter_plot_qc_bins(..., metric = "cv_scaled" | "chi_z")`
  - `nipter_plot_reference_sex_boxplots()`
  - `nipter_plot_reference_sex_scatter(..., space = "fraction" |
    "ratio")`

### Sex prediction

`nipter_sex_model()` trains Gaussian-mixture models for sex prediction
from Y-fraction, X/Y fraction pairs, or Y-unique ratios.
`nipter_predict_sex()` combines one or more models by majority vote.

``` r
cg <- nipter_as_control_group(samples)
cg <- nipter_gc_correct(cg, fasta = "hg38.fa")
corrected <- nipter_chi_correct(samples[[1]], cg)
cg <- corrected$control_group

model_y  <- nipter_sex_model(cg, method = "y_fraction")
model_xy <- nipter_sex_model(cg, method = "xy_fraction")

nipter_predict_sex(test_sample, models = list(model_y, model_xy))
```

## Optional Upstream WisecondorX CLI

The `wisecondorx_*()` wrappers are for upstream interoperability and
conformance, not for the main native R workflow.

Use them when you specifically need:

  - upstream `.npz` files
  - upstream Python plotting output
  - direct comparison against the official `wisecondorx` implementation

### Writing `.npz` files

``` r
fixture_npz <- "sample.npz"
np <- reticulate::import("numpy", convert = FALSE)

bam_convert_npz(
  bam     = "sample.dm.bam",
  npz     = fixture_npz,
  binsize = 100000L,
  rmdup   = "streaming",
  np      = np
)
```

### Wrapping the upstream CLI

The wrappers expose the main upstream knobs explicitly:

  - `wisecondorx_newref(..., yfrac = ..., plotyfrac = ...)`
  - `wisecondorx_predict(..., beta = ..., gender = ..., plot = TRUE,
    add_plot_title = TRUE, ylim = c(-0.5, 0.5), cairo = TRUE)`

<!-- end list -->

``` r
# Requires: condathis installed; wisecondorx bioconda package auto-installed

npz_files <- list.files("controls/", "\\.npz$", full.names = TRUE)

wisecondorx_newref(
  npz_files   = npz_files,
  output      = "reference.npz",
  ref_binsize = 100000L,
  nipt        = TRUE,
  refsize     = 300L,
  cpus        = 4L
)

wisecondorx_predict(
  npz           = "sample.npz",
  ref           = "reference.npz",
  output_prefix = "results/sample",
  bed           = TRUE,
  seed          = 1L
)
```

For cohort workflows, the package also provides:

  - `inst/scripts/preprocess_cohort.R`
  - `inst/scripts/build_reference.R`
  - `inst/scripts/predict_cohort.R`
  - `inst/scripts/wisecondorx_conformance.R`

## SRA Metadata Helpers

`sra_runinfo_url()`, `download_sra_runinfo()`, and `read_sra_runinfo()`
are small helpers for staging NCBI SRA run metadata.

``` r
sra_runinfo_url("PRJNA400134")
#> [1] "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA400134"
```

## Validation

Per [rewrites.bio](https://rewrites.bio) principle 4.3, see
[`VALIDATION.md`](VALIDATION.md) for a full account of what has been
validated, on what data, against which upstream versions, and what
remains unvalidated.

Current high-level status:

  - `bam_convert()` is bin-for-bin identical to Python WisecondorX on
    the current convert fixture
  - the NIPTeR statistical layer is covered by dedicated native tests
  - native and upstream WisecondorX conformance is exercised on prepared
    real-cohort artifacts via `inst/scripts/wisecondorx_conformance.R`
  - beta mode is implemented but still needs more clinical validation

## Development

`README.Rmd` is the editable source for this document. Run `make readme`
to rerender `README.md` and `make rd` to regenerate documentation from
roxygen2 comments.
