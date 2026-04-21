
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWisecondorX

<!-- badges: start -->

[![R-CMD-check](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml)
[![rewrites.bio](https://rewrites.bio/badges/rewrites-bio.svg)](https://rewrites.bio)
<!-- badges: end -->

`RWisecondorX` is an R toolkit for copy number analysis and trisomy
prediction in non-invasive prenatal testing (NIPT), built on top of
`Rduckhts` and DuckDB. It ports both
[WisecondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX)
(CNV detection) and [NIPTeR](https://github.com/molgenis/NIPTeR)
(trisomy prediction via Z-scores, NCV, regression, and chi-squared
correction), following the [rewrites.bio](https://rewrites.bio)
principles: exact emulation of upstream outputs, full credit to the
original authors, and transparency about AI assistance.

The core convert step now runs through the native `bam_bin_counts(...)`
kernel bundled in `Rduckhts`, with no Python dependency. `bam_convert()`
returns per-bin read counts in memory; `bam_convert_bed()` writes them
to a bgzipped, tabix-indexed BED file readable by DuckDB or any
tabix-aware tool; and `bam_convert_npz()` serialises them to a
WisecondorX-compatible `.npz` via `reticulate` for Python CLI
conformance and interop. The convert implementation exactly replicates
the upstream `larp` / `larp2` streaming deduplication behaviour,
achieving bin-for-bin agreement with the Python implementation.

The NIPTeR statistical layer provides GC correction (via on-the-fly
FASTA computation, not bundled tables), chi-squared overdispersion
correction, chromosomal fraction Z-scores, normalised chromosome values
(NCV), and forward stepwise regression — all operating on `NIPTeRSample`
objects produced by `nipter_bin_bam()`.

For pipelines that need the full upstream toolchain, thin `condathis`
wrappers cover every stage: `wisecondorx_convert()`,
`wisecondorx_newref()`, and `wisecondorx_predict()` each delegate to the
official bioconda package without requiring a manual conda setup.

## Installation

`RWisecondorX` is developed against the current edge of `Rduckhts`, not
an older CRAN snapshot. Install `Rduckhts` from the r-universe first,
then install `RWisecondorX`:

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
needed only to write `.npz` files; `condathis` is needed only for the
upstream CLI wrappers. Neither is required for the core convert step.

## Convert A BAM To Bin Counts

``` r
bins <- bam_convert(
  "sample.dm.bam",
  binsize = 100000L,
  mapq = 1L,
  rmdup = "streaming"
)
```

The `rmdup` argument controls duplicate handling. `"streaming"` exactly
replicates the upstream WisecondorX `larp` / `larp2` state machine and
is the default. `"none"` skips deduplication entirely, matching
`wisecondorx convert --normdup` for NIPT data where read depth is low.
`"flag"` uses the SAM duplicate flag (`0x400`) for BAMs already
processed by Picard or sambamba.

## BED.gz Output And Round-Trip

`bam_convert_bed()` is the language-agnostic output path. It produces a
bgzipped BED file (`chrom`, `start`, `end`, `count`; 0-based
coordinates) and a tabix index, using `rduckhts_bgzip` and
`rduckhts_tabix_index` from `Rduckhts` — no external tools or Python
required.

``` r
bed_file <- "sample.wisecondorx.bed.gz"

bam_convert_bed(
  bam     = "sample.dm.bam",
  bed     = bed_file,
  binsize = 100000L,
  rmdup   = "streaming"
)

# Round-trip: reload into the named-list format for rwisecondorx_newref()
sample_reloaded <- bed_to_sample(bed_file, binsize = 100000L)
```

`bed_to_nipter_sample()` reads a 5-column NIPTeR BED.gz into a
`NIPTeRSample` object compatible with all NIPTeR statistical functions.
Column count is auto-detected (5 = CombinedStrands, 9 =
SeparatedStrands).

## Native WisecondorX Pipeline (No Python)

`rwisecondorx_newref()` and `rwisecondorx_predict()` are pure R/Rcpp
implementations of the full WisecondorX algorithm — reference building,
PCA correction, within-sample normalisation, CBS segmentation, and
aberration calling — with no Python dependency. Performance-critical KNN
reference-bin finding is compiled via Rcpp with OpenMP parallelisation.
CBS segmentation uses DNAcopy directly, or `ParDNAcopy` when `parallel =
TRUE`.

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

`rwisecondorx_predict(parallel = TRUE)` requires `ParDNAcopy`. If that
package is not installed, call `rwisecondorx_predict(..., parallel =
FALSE)` to use serial `DNAcopy::segment()`.

The native predictor also exposes the two upstream predict controls most
often needed in practice:

  - `beta` for purity-based ratio calling
  - `gender = "F"` / `"M"` to force the gonosomal branch instead of
    using the Y-fraction classifier

<!-- end list -->

``` r
pred_forced <- rwisecondorx_predict(
  test_sample,
  ref,
  gender = "F",
  beta = NULL,
  seed = 42L
)
```

### Native WisecondorX QC And Outputs

On the native R side, prediction output is currently tabular rather than
graphical:

  - `write_wisecondorx_output()` writes
      - `<outprefix>_bins.bed`
      - `<outprefix>_segments.bed`
      - `<outprefix>_aberrations.bed`
      - `<outprefix>_statistics.txt`
  - `rwisecondorx_ref_qc()` runs the upstream-style reference QC
    heuristics and returns a structured report, with optional JSON
    output

<!-- end list -->

``` r
qc <- rwisecondorx_ref_qc(ref, output_json = "reference_qc.json")
write_wisecondorx_output(pred, outprefix = "results/sample_01")
```

The package does **not yet** provide native ggplot helpers for
`WisecondorXPrediction` objects. If you need the upstream PNG/PDF
plotting behaviour, use the CLI wrapper `wisecondorx_predict(..., plot =
TRUE)` instead.

### Multi-File BED Input For Reference Building

`rwisecondorx_newref()` also accepts a directory of 4-column BED.gz
files produced by `bam_convert_bed()` via the `bed_dir` parameter. All
matching files are loaded in a single DuckDB pass — no in-memory sample
list needed:

``` r
ref <- rwisecondorx_newref(
  bed_dir = "bed_counts/",     # directory of *.bed.gz files
  binsize = 100000L,
  nipt    = TRUE,
  cpus    = 4L
)
```

## Optional NPZ And CLI Workflow

`bam_convert_npz()` and the `wisecondorx_*()` CLI wrappers are the
conformance and interoperability path. The production native path is
`bam_convert_bed()`/`bam_convert()` with `rwisecondorx_newref()` and
`rwisecondorx_predict()`. When the full upstream Python WisecondorX
pipeline is needed, `wisecondorx_convert()` wraps `wisecondorx convert`
via `condathis` and produces an `.npz` through the official
implementation. All three wrappers handle conda environment creation
automatically on first use.

``` r
wisecondorx_convert(
  bam     = "sample.bam",
  npz     = "sample.npz",
  binsize = 100000L,
  normdup = FALSE          # TRUE → --normdup (NIPT)
)
```

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

The CLI wrappers expose the upstream arguments documented in the
WisecondorX README and delegate to the official `wisecondorx` bioconda
package.

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

The wrappers keep the upstream plotting and manual override controls
explicit:

  - `wisecondorx_newref(..., yfrac = ..., plotyfrac = ...)`
  - `wisecondorx_predict(..., beta = ..., gender = ..., plot = TRUE,
    add_plot_title = TRUE, ylim = c(-0.5, 0.5), cairo = TRUE)`

For cohort preprocessing and reference building, the package now keeps
the native and upstream WisecondorX paths explicit:

  - `inst/scripts/preprocess_cohort.R` stages a cohort manifest and
    writes:
      - SeqFF fetal-fraction estimates
      - native `rwisecondorx` BED.gz files in `rwcx_beds/`
      - optional upstream `wisecondorx` NPZ files in `wisecondorx_npz/`,
        produced by the Python CLI wrapper
      - NIPTeR BED.gz files
  - `inst/scripts/build_reference.R --mode rwisecondorx` builds a native
    RDS reference from BED.gz files
  - `inst/scripts/build_reference.R --mode wisecondorx` builds an
    upstream NPZ reference from NPZ files
  - `inst/scripts/predict_cohort.R` scores a preprocessed cohort against
    native RWisecondorX and/or NIPTeR references
  - `inst/scripts/wisecondorx_conformance.R` takes those prepared
    artifacts (`bed` + `npz` test cases, plus both references), runs
    both `predict` paths, and compares the resulting per-bin `_bins.bed`
    outputs

## Beta Mode (Somatic CNV)

`rwisecondorx_predict()` supports a `beta` parameter for purity-based
aberration calling, appropriate for somatic CNV analysis where tumour
purity is known. When `beta` is supplied, ratio cutoffs replace Z-score
thresholds:

    gain cutoff: log2((ploidy + beta/2) / ploidy)
    loss cutoff: log2((ploidy - beta/2) / ploidy)

`beta` is the tumour purity estimate (0–1; 0 = most liberal calling, 1 =
most conservative; optimally close to true purity). This mode is **not
appropriate for NIPT**, which uses the default Z-score calling.

``` r
# Somatic CNV: 40% purity tumour sample
pred_somatic <- rwisecondorx_predict(
  sample    = tumour_sample,
  reference = ref,
  beta      = 0.4,   # ~tumour purity
  seed      = 1L
)
```

## NIPTeR-Style Trisomy Prediction

`RWisecondorX` ports the statistical analysis layer from
[NIPTeR](https://github.com/molgenis/NIPTeR) (de Weerd et al.),
replacing `Rsamtools` with `Rduckhts` and bundled `sysdata.rda` tables
with on-the-fly computation via `rduckhts_fasta_nuc()`. The analysis
pipeline follows NIPTeR’s workflow: bin samples, build a control group,
correct for GC bias and overdispersion, then compute Z-scores or NCV
values for trisomy prediction.

### Binning

`nipter_bin_bam()` produces `NIPTeRSample` objects compatible with all
downstream functions. Pre-filtering flags mirror `samtools view -f`/`-F`
conventions, so real-world NIPT pipelines that mark duplicates with
Picard can pass `exclude_flags = 1024L` directly.

``` r
samples <- lapply(bam_files, nipter_bin_bam, binsize = 50000L)
```

### Control Group Construction

A control group is the reference cohort against which test samples are
scored. `nipter_as_control_group()` validates and deduplicates the input
samples. `nipter_diagnose_control_group()` reports per-chromosome
Z-scores and Shapiro-Wilk normality p-values across the group, useful
for identifying outlier samples before scoring. When working with a
large cohort, `nipter_match_control_group()` selects the
closest-matching subset for a given test sample based on chromosomal
fraction distance (Rcpp+OpenMP accelerated). `nipter_match_matrix()`
computes the full N×N pairwise SSD matrix — useful for QC outlier
detection before building the control group.

``` r
cg <- nipter_as_control_group(samples)
diag <- nipter_diagnose_control_group(cg)
matched_cg <- nipter_match_control_group(test_sample, cg, n = 50L)
mm <- nipter_match_matrix(cg, cpus = 4L)
```

### GC Correction

`nipter_gc_correct()` adjusts bin counts for GC-content bias. GC
percentages are computed per bin from the reference FASTA via
`rduckhts_fasta_nuc()` rather than from static bundled tables. Two
methods are available: LOESS regression (default, matching NIPTeR) and
bin-weight normalisation. For cohort workflows, `nipter_gc_precompute()`
pre-computes and saves GC fractions to a TSV.bgz file so the FASTA is
only scanned once.

``` r
# Pre-compute GC table (one-time per genome / binsize)
nipter_gc_precompute(fasta = "hg38.fa", binsize = 50000L,
                     out = "gc_table.tsv.bgz")

# Apply GC correction using pre-computed table (fast)
cg          <- nipter_gc_correct(cg, gc_table = "gc_table.tsv.bgz")
test_sample <- nipter_gc_correct(test_sample, gc_table = "gc_table.tsv.bgz")
```

### Chi-Squared Correction

`nipter_chi_correct()` identifies overdispersed bins across the control
group using a normalised chi-squared test and downweights them. The
correction is applied simultaneously to both the test sample and control
group, maintaining consistency for downstream scoring.

``` r
corrected   <- nipter_chi_correct(test_sample, cg, chi_cutoff = 3.5)
test_sample <- corrected$sample
cg          <- corrected$control_group
```

### Scoring

`nipter_z_score()` computes the chromosomal fraction Z-score for a focus
chromosome. A Z-score above 3 is the conventional threshold for trisomy.

``` r
z21 <- nipter_z_score(test_sample, cg, chromo_focus = 21)
z21$sample_z_score
z21$control_statistics
```

`nipter_ncv_score()` computes the normalised chromosome value, searching
for the denominator chromosome set that minimises control-group variance
before computing the Z-score. More robust than the standard Z-score when
the control group is small.

``` r
ncv21 <- nipter_ncv_score(test_sample, cg, chromo_focus = 21)
ncv21$sample_ncv_score
ncv21$denominators
```

`nipter_regression()` uses forward stepwise linear regression to predict
the focus chromosome’s fraction from other chromosomes, then scores the
residual.

``` r
reg21 <- nipter_regression(test_sample, cg, chromo_focus = 21)
reg21$models[[1]]$z_score
reg21$models[[1]]$predictors
```

### Reference Building, QC, And Plots

For production NIPTeR workflows, the package now exposes a typed
reference-building and QC layer rather than only the primitive scoring
functions:

  - `nipter_build_reference()` packages the control group, sex models,
    and reference frame
  - `nipter_build_gaunosome_models()` adds sex-stratified X/Y NCV and
    regression models
  - `nipter_control_group_qc()` computes chromosome-, sample-,
    sex-model-, and optional bin-level QC summaries
  - `write_nipter_reference_plots()` writes the current QC/sex-model
    plot bundle

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

The current NIPTeR plotting helpers are:

  - `nipter_plot_qc_chromosomes()` for chromosome CV summaries across Z,
    NCV, and RBZ, including XX/XY sex-model rows
  - `nipter_plot_qc_samples()` for retained-vs-outlier control
    diagnostics
  - `nipter_plot_qc_bins(..., metric = "cv_scaled" | "chi_z")` for
    autosomal bin-level QC
  - `nipter_plot_reference_sex_boxplots()` for
    fraction/Y-unique/sex-cluster-Z distributions
  - `nipter_plot_reference_sex_scatter(..., space = "fraction" |
    "ratio")` for X/Y and RR spaces

## Sex Prediction

`nipter_sex_model()` trains a 2-component Gaussian mixture model on
chromosome fraction data from a control group for sex prediction. Three
model types are available: Y-fraction alone, bivariate X+Y fractions,
and Y-unique region read ratios from specific Y-chromosome genes.
`nipter_predict_sex()` classifies a test sample using one or more models
with majority-vote consensus (ties → female, the conservative choice for
NIPT).

``` r
# Requires: mclust package; whole-genome BAMs; FASTA for GC correction
cg <- nipter_as_control_group(samples)
cg <- nipter_gc_correct(cg, fasta = "hg38.fa")
corrected <- nipter_chi_correct(samples[[1]], cg)
cg <- corrected$control_group

model_y  <- nipter_sex_model(cg, method = "y_fraction")
model_xy <- nipter_sex_model(cg, method = "xy_fraction")

nipter_predict_sex(test_sample, models = list(model_y, model_xy))
```

## SRA Metadata Helpers

`sra_runinfo_url()`, `download_sra_runinfo()`, and `read_sra_runinfo()`
are small helpers for staging NCBI SRA run metadata.

``` r
sra_runinfo_url("PRJNA400134")
#> [1] "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA400134"
```

## Validation And Conformance

Per [rewrites.bio](https://rewrites.bio) principle 4.3, see
[`VALIDATION.md`](VALIDATION.md) for a full account of what has been
validated, on what data, against which upstream versions, and what
remains unvalidated.

**Summary of current status:**

  - `bam_convert()`: bin-for-bin identical to Python WisecondorX on
    `HG00106.chrom11` (0 mismatches, 25,115 bins).
  - NIPTeR statistical layer (Z-score, NCV, chi, regression): formulas
    verified against inline reference implementations in
    `test_nipter_conformance.R`.
  - Native and upstream WisecondorX conformance is run on prepared
    real-cohort artifacts via `inst/scripts/wisecondorx_conformance.R`,
    comparing the resulting per-bin `_bins.bed` outputs from both
    `predict` implementations.
  - Beta mode (purity-based calling): implemented, not yet tested with
    clinical samples.

## Development

`README.Rmd` is the editable source for this document. Run `make readme`
to rerender `README.md` and `make rd` to regenerate documentation from
roxygen2 comments.
