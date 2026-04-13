# RWisecondorX

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
Python dependency.
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md)
returns per-bin read counts in memory;
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md)
writes them to a bgzipped, tabix-indexed BED file readable by DuckDB or
any tabix-aware tool; and
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md)
serialises them to a WisecondorX-compatible `.npz` via `reticulate`. The
convert implementation exactly replicates the upstream `larp` / `larp2`
streaming deduplication behaviour, achieving bin-for-bin agreement with
the Python implementation.

The NIPTeR statistical layer provides GC correction (via on-the-fly
FASTA computation, not bundled tables), chi-squared overdispersion
correction, chromosomal fraction Z-scores, normalised chromosome values
(NCV), and forward stepwise regression — all operating on `NIPTeRSample`
objects produced by
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md).

For pipelines that need the full upstream toolchain, thin `condathis`
wrappers cover every stage:
[`wisecondorx_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_convert.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md),
and
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)
each delegate to the official bioconda package without requiring a
manual conda setup.

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

[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md)
is the language-agnostic output path. It produces a bgzipped BED file
(`chrom`, `start`, `end`, `count`; 0-based coordinates) and a tabix
index alongside it, using `rduckhts_bgzip` and `rduckhts_tabix_index`
from `Rduckhts` — no external tools or Python required. The resulting
file can be queried by region directly from DuckDB.

``` r
bam_convert_bed(
  bam     = "sample.bam",
  bed     = "sample.bed.gz",   # → sample.bed.gz + sample.bed.gz.tbi
  binsize = 5000L,
  rmdup   = "streaming"
)
```

## BED.gz Round-Trip

Once bin counts have been written to BED.gz,
[`bed_to_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_sample.md)
and
[`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md)
read them back into the in-memory formats expected by the analysis
pipelines. This means you can bin once, store the compact BED.gz files,
and reload them for any number of downstream analyses without touching
the original BAM again.

[`bed_to_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_sample.md)
reads a 4-column WisecondorX BED.gz into the named list of integer
vectors consumed by
[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md)
and
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md):

``` r
sample <- bed_to_sample("sample.bed.gz", binsize = 5000L)
```

[`bed_to_nipter_sample()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bed_to_nipter_sample.md)
reads a 5-column (CombinedStrands) or 9-column (SeparatedStrands) NIPTeR
BED.gz into a `NIPTeRSample` object compatible with all NIPTeR
statistical functions. Column count is auto-detected:

``` r
ns <- bed_to_nipter_sample("sample_nipter.bed.gz", binsize = 50000L)
class(ns)  # "NIPTeRSample" "CombinedStrands"
```

A typical production workflow bins BAMs once and stores the results,
then loads them for analysis across multiple sessions:

``` r
library(RWisecondorX)

bam_files <- list.files("bams/", "\\.bam$", full.names = TRUE)
bed_dir   <- "bed_counts/"
dir.create(bed_dir, showWarnings = FALSE)

for (bam in bam_files) {
  bed <- file.path(bed_dir, sub("\\.bam$", ".bed.gz", basename(bam)))
  bam_convert_bed(bam, bed, binsize = 5000L, rmdup = "streaming")
}

bed_files <- list.files(bed_dir, "\\.bed\\.gz$", full.names = TRUE)
samples   <- lapply(bed_files, bed_to_sample, binsize = 5000L)

ref <- rwisecondorx_newref(samples, binsize = 100000L, nipt = TRUE, cpus = 4L)
saveRDS(ref, "reference.rds")
```

## Optional NPZ And CLI Workflow

When the full upstream WisecondorX pipeline is needed,
[`wisecondorx_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_convert.md)
wraps `wisecondorx convert` via `condathis` and produces an `.npz`
through the official Python implementation.
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md)
and
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)
cover the remaining stages. All three handle conda environment creation
automatically on first use.

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

[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md)
and
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md)
are pure R/Rcpp implementations of the full WisecondorX algorithm —
reference building, PCA correction, within-sample normalisation, CBS
segmentation, and aberration calling — with no Python dependency.
Performance-critical KNN reference-bin finding is compiled via Rcpp with
OpenMP parallelisation. CBS segmentation uses DNAcopy (or ParDNAcopy for
parallel operation) directly instead of calling the upstream `CBS.R`
subprocess.

The typical workflow is: bin samples with
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
build a reference with
[`rwisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_newref.md),
then call aberrations with
[`rwisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/rwisecondorx_predict.md).
The reference is a plain R list that can be serialised with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and shared between
sessions.

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

[`generate_cohort()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/generate_cohort.md)
creates a synthetic BAM cohort suitable for end-to-end pipeline testing.
It produces 50 samples using “compressed” chromosome lengths (100bp per
100kb genomic bin), so each BAM is ~435KB. The cohort includes euploid
females and males plus three trisomy samples (T21, T18, T13) for
sensitivity validation.

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

[`nipter_sex_model()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_sex_model.md)
trains a 2-component Gaussian mixture model on chromosome fraction data
from a control group for sex prediction. Three model types are
available: Y-fraction alone, bivariate X+Y fractions, and Y-unique
region read ratios from specific Y-chromosome genes.
[`nipter_predict_sex()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_predict_sex.md)
classifies a test sample using one or more models with majority-vote
consensus.

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

[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
produces `NIPTeRSample` objects compatible with all downstream
functions. Pre-filtering flags mirror `samtools view -f`/`-F`
conventions, so real-world NIPT pipelines that mark duplicates with
Picard can pass `exclude_flags = 1024L` directly.

``` r
samples <- lapply(bam_files, nipter_bin_bam, binsize = 50000L)
```

### Control Group Construction

A control group is the reference cohort against which test samples are
scored.
[`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md)
validates and deduplicates the input samples.
[`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md)
reports per-chromosome Z-scores and Shapiro-Wilk normality p-values
across the group, useful for identifying outlier samples before scoring.
When working with a large cohort,
[`nipter_match_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_control_group.md)
selects the closest-matching subset for a given test sample based on
chromosomal fraction distance.

``` r
cg <- nipter_as_control_group(samples)
diag <- nipter_diagnose_control_group(cg)
matched_cg <- nipter_match_control_group(test_sample, cg, n_controls = 50L)
```

### GC Correction

[`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md)
adjusts bin counts for GC-content bias. GC percentages are computed per
bin from the reference FASTA via `rduckhts_fasta_nuc()` rather than from
static bundled tables. Two methods are available: LOESS regression
(default, matching NIPTeR) and bin-weight normalisation. The function
accepts either a single sample or an entire control group.

``` r
cg <- nipter_gc_correct(cg, fasta = "hg38.fa", method = "loess")
test_sample <- nipter_gc_correct(test_sample, fasta = "hg38.fa")
```

### Chi-Squared Correction

[`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md)
identifies overdispersed bins across the control group using a
normalised chi-squared test and downweights them. The correction is
applied simultaneously to both the test sample and control group,
maintaining consistency for downstream scoring.

``` r
corrected <- nipter_chi_correct(test_sample, cg, chi_cutoff = 3.5)
test_sample <- corrected$sample
cg <- corrected$control_group
```

### Scoring

[`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md)
computes the chromosomal fraction Z-score for a focus chromosome,
testing how far the sample deviates from the control distribution. A
Z-score above 3 is the conventional threshold for trisomy. Each result
includes a Shapiro-Wilk p-value for the control distribution so that
normality violations are visible.

``` r
z21 <- nipter_z_score(test_sample, cg, chromo_focus = 21)
z21$sample_z_score
z21$control_statistics
```

[`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md)
computes the normalised chromosome value, which searches for the
denominator chromosome set that minimises control-group variance before
computing the Z-score. This is more robust than the standard Z-score
when the control group is small or the chromosomal fraction has high
variance.

``` r
ncv21 <- nipter_ncv_score(test_sample, cg, chromo_focus = 21)
ncv21$sample_ncv_score
ncv21$denominators
```

[`nipter_regression()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md)
uses forward stepwise linear regression to predict the focus
chromosome’s fraction from other chromosomes, then scores the residual.
It builds multiple models with a train/test split and selects the CV
(practical vs theoretical) that gives the more conservative estimate.

``` r
reg21 <- nipter_regression(test_sample, cg, chromo_focus = 21)
reg21$models[[1]]$z_score
reg21$models[[1]]$predictors
```

## SRA Metadata Helpers

[`sra_runinfo_url()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/sra_runinfo.md),
[`download_sra_runinfo()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/sra_runinfo.md),
and
[`read_sra_runinfo()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/sra_runinfo.md)
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
