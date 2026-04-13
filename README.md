
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWisecondorX

<!-- badges: start -->

[![R-CMD-check](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sounkou-bioinfo/RWisecondorX/actions/workflows/R-CMD-check.yaml)
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

The core convert step runs entirely in R/SQL via `Rduckhts` with no
Python dependency. `bam_convert()` returns per-bin read counts in
memory; `bam_convert_bed()` writes them to a bgzipped, tabix-indexed BED
file readable by DuckDB or any tabix-aware tool; and `bam_convert_npz()`
serialises them to a WisecondorX-compatible `.npz` via `reticulate`. The
convert implementation exactly replicates the upstream `larp` / `larp2`
streaming deduplication behaviour, achieving bin-for-bin agreement with
the Python implementation.

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

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("sounkou-bioinfo/RWisecondorX")
```

`RWisecondorX` imports `Rduckhts`, `DBI`, and `duckdb`. `reticulate` is
needed only to write `.npz` files; `condathis` is needed only for the
upstream CLI wrappers. Neither is required for the core convert step.

## Convert A BAM To Bin Counts

The package ships with a small chromosome 11 BAM fixture for examples
and tests.

``` r
library(RWisecondorX)

fixture_bam <- system.file(
  "extdata",
  "hg00106_chr11_fixture.bam",
  package = "RWisecondorX"
)

bins <- bam_convert(
  fixture_bam,
  binsize = 5000L,
  rmdup = "streaming"
)

names(Filter(Negate(is.null), bins))
#> [1] "11"

length(bins[["11"]])
#> [1] 1198

sum(bins[["11"]])
#> [1] 982
```

The `rmdup` argument controls duplicate handling. `"streaming"` exactly
replicates the upstream WisecondorX `larp` / `larp2` state machine and
is the default. `"none"` skips deduplication entirely, matching
`wisecondorx convert --normdup` for NIPT data where read depth is low.
`"flag"` uses the SAM duplicate flag (`0x400`) for BAMs already
processed by Picard or sambamba.

## Optional BED.gz Output

`bam_convert_bed()` is the language-agnostic output path. It produces a
bgzipped BED file (`chrom`, `start`, `end`, `count`; 0-based
coordinates) and a tabix index alongside it, using `rduckhts_bgzip` and
`rduckhts_tabix_index` from `Rduckhts` — no external tools or Python
required. The resulting file can be queried by region directly from
DuckDB.

``` r
bam_convert_bed(
  bam     = "sample.bam",
  bed     = "sample.bed.gz",   # → sample.bed.gz + sample.bed.gz.tbi
  binsize = 5000L,
  rmdup   = "streaming"
)
```

## Optional NPZ And CLI Workflow

When the full upstream WisecondorX pipeline is needed,
`wisecondorx_convert()` wraps `wisecondorx convert` via `condathis` and
produces an `.npz` through the official Python implementation.
`wisecondorx_newref()` and `wisecondorx_predict()` cover the remaining
stages. All three handle conda environment creation automatically on
first use.

``` r
wisecondorx_convert(
  bam     = "sample.bam",
  npz     = "sample.npz",
  binsize = 5000L,
  normdup = FALSE          # TRUE → --normdup (NIPT)
)
```

``` r
library(RWisecondorX)

fixture_npz <- tempfile(fileext = ".npz")
np <- reticulate::import("numpy", convert = FALSE)

bam_convert_npz(
  bam     = fixture_bam,
  npz     = fixture_npz,
  binsize = 5000L,
  rmdup   = "streaming",
  np      = np
)

c(
  npz_created = file.exists(fixture_npz),
  npz_size    = file.info(fixture_npz)$size
)
#> npz_created    npz_size 
#>           1         392
```

The CLI wrappers expose the upstream arguments documented in the
WisecondorX README. They delegate to the official `wisecondorx` bioconda
package installed automatically via `condathis`.

``` r
library(RWisecondorX)

npz_files <- list.files("controls/", "\\.npz$", full.names = TRUE)

wisecondorx_newref(
  npz_files   = npz_files,
  output      = "reference.npz",
  ref_binsize = 100000L,
  nipt        = TRUE,
  refsize     = 300L,
  plotyfrac   = "yfrac_plot.png",
  cpus        = 4L
)

wisecondorx_predict(
  npz             = "sample.npz",
  ref             = "reference.npz",
  output_prefix   = "results/sample",
  bed             = TRUE,
  plot            = TRUE,
  add_plot_title  = TRUE,
  seed            = 1L
)
```

## Native WisecondorX Pipeline (No Python)

`rwisecondorx_newref()` and `rwisecondorx_predict()` are pure R/Rcpp
implementations of the full WisecondorX algorithm — reference building,
PCA correction, within-sample normalisation, CBS segmentation, and
aberration calling — with no Python dependency. Performance-critical KNN
reference-bin finding is compiled via Rcpp with OpenMP parallelisation.
CBS segmentation uses DNAcopy (or ParDNAcopy for parallel operation)
directly instead of calling the upstream `CBS.R` subprocess.

The typical workflow is: bin samples with `bam_convert()`, build a
reference with `rwisecondorx_newref()`, then call aberrations with
`rwisecondorx_predict()`. The reference is a plain R list that can be
serialised with `saveRDS()` and shared between sessions.

``` r
library(RWisecondorX)

bam_files <- list.files("controls/", "\\.bam$", full.names = TRUE)

samples <- lapply(bam_files, bam_convert, binsize = 5000L, rmdup = "streaming")

ref <- rwisecondorx_newref(
  samples  = samples,
  binsize  = 100000L,
  nipt     = TRUE,
  refsize  = 300L,
  cpus     = 4L
)

saveRDS(ref, "reference.rds")

test_sample <- bam_convert("test.bam", binsize = 5000L, rmdup = "streaming")

prediction <- rwisecondorx_predict(
  sample         = test_sample,
  reference      = ref,
  sample_binsize = 5000L,
  outprefix      = "results/test",
  zscore         = 5,
  alpha          = 1e-4,
  seed           = 1L
)

prediction$aberrations
prediction$statistics
```

### Synthetic Cohort For Testing

`generate_cohort()` creates a synthetic BAM cohort suitable for
end-to-end pipeline testing. It produces 50 samples using “compressed”
chromosome lengths (100bp per 100kb genomic bin), so each BAM is ~435KB.
The cohort includes euploid females and males plus three trisomy samples
(T21, T18, T13) for sensitivity validation.

``` r
library(RWisecondorX)

out_dir <- tempfile("cohort_")
generate_cohort(out_dir)

bams <- list.files(out_dir, "\\.bam$", full.names = TRUE)
length(bams)  # 50

samples <- lapply(bams, bam_convert,
                  binsize = COMPRESSED_BINSIZE, rmdup = "none")
ref <- rwisecondorx_newref(samples, binsize = COMPRESSED_BINSIZE,
                           nipt = TRUE, cpus = 4L)

t21_bam <- file.path(out_dir, "sample_048_T21.bam")
t21 <- bam_convert(t21_bam, binsize = COMPRESSED_BINSIZE, rmdup = "none")
pred <- rwisecondorx_predict(t21, ref,
                             sample_binsize = COMPRESSED_BINSIZE,
                             zscore = 5, seed = 42L)

pred$aberrations[pred$aberrations$type == "gain", c("chr", "start", "end", "zscore")]
```

## Sex Prediction

`nipter_sex_model()` trains a 2-component Gaussian mixture model on
chromosome fraction data from a control group for sex prediction. Three
model types are available: Y-fraction alone, bivariate X+Y fractions,
and Y-unique region read ratios from specific Y-chromosome genes.
`nipter_predict_sex()` classifies a test sample using one or more models
with majority-vote consensus.

``` r
library(RWisecondorX)

cg <- nipter_as_control_group(samples)
cg <- nipter_gc_correct(cg, fasta = "hg38.fa")
corrected <- nipter_chi_correct(samples[[1]], cg)
cg <- corrected$control_group

model_y  <- nipter_sex_model(cg, method = "y_fraction")
model_xy <- nipter_sex_model(cg, method = "xy_fraction")

nipter_predict_sex(test_sample, models = list(model_y, model_xy))
```

## NIPTeR-Style Trisomy Prediction

`RWisecondorX` also ports the statistical analysis layer from
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
fraction distance.

``` r
cg <- nipter_as_control_group(samples)
diag <- nipter_diagnose_control_group(cg)
matched_cg <- nipter_match_control_group(test_sample, cg, n_controls = 50L)
```

### GC Correction

`nipter_gc_correct()` adjusts bin counts for GC-content bias. GC
percentages are computed per bin from the reference FASTA via
`rduckhts_fasta_nuc()` rather than from static bundled tables. Two
methods are available: LOESS regression (default, matching NIPTeR) and
bin-weight normalisation. The function accepts either a single sample or
an entire control group.

``` r
cg <- nipter_gc_correct(cg, fasta = "hg38.fa", method = "loess")
test_sample <- nipter_gc_correct(test_sample, fasta = "hg38.fa")
```

### Chi-Squared Correction

`nipter_chi_correct()` identifies overdispersed bins across the control
group using a normalised chi-squared test and downweights them. The
correction is applied simultaneously to both the test sample and control
group, maintaining consistency for downstream scoring.

``` r
corrected <- nipter_chi_correct(test_sample, cg, chi_cutoff = 3.5)
test_sample <- corrected$sample
cg <- corrected$control_group
```

### Scoring

`nipter_z_score()` computes the chromosomal fraction Z-score for a focus
chromosome, testing how far the sample deviates from the control
distribution. A Z-score above 3 is the conventional threshold for
trisomy. Each result includes a Shapiro-Wilk p-value for the control
distribution so that normality violations are visible.

``` r
z21 <- nipter_z_score(test_sample, cg, chromo_focus = 21)
z21$sample_z_score
z21$control_statistics
```

`nipter_ncv_score()` computes the normalised chromosome value, which
searches for the denominator chromosome set that minimises control-group
variance before computing the Z-score. This is more robust than the
standard Z-score when the control group is small or the chromosomal
fraction has high variance.

``` r
ncv21 <- nipter_ncv_score(test_sample, cg, chromo_focus = 21)
ncv21$sample_ncv_score
ncv21$denominators
```

`nipter_regression()` uses forward stepwise linear regression to predict
the focus chromosome’s fraction from other chromosomes, then scores the
residual. It builds multiple models with a train/test split and selects
the CV (practical vs theoretical) that gives the more conservative
estimate.

``` r
reg21 <- nipter_regression(test_sample, cg, chromo_focus = 21)
reg21$models[[1]]$z_score
reg21$models[[1]]$predictors
```

## SRA Metadata Helpers

`sra_runinfo_url()`, `download_sra_runinfo()`, and `read_sra_runinfo()`
are small helpers for staging NCBI SRA run metadata used by conformance
fixtures.

``` r
library(RWisecondorX)

sra_runinfo_url("PRJNA400134")
#> [1] "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA400134"
```

## Development

`README.Rmd` is the editable source for this document. Run `make readme`
to rerender `README.md`, `make rd` to regenerate documentation from
roxygen2 comments, and `make fixtures` to rebuild the bundled test data.
