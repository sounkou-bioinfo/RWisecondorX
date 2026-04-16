# RWisecondorX — S7 Class Hierarchy Design
*2026-04-16*

## 0. Corrections to the previous code review

### P2-3 was partially wrong — `include_XY` for GC correction IS correct

`BinReads.R` (the upstream preprocessing script) performs LOESS GC correction on sex
chromosome bins using nearest-neighbor interpolation from the autosomal LOESS fit.
The output is saved as `*.dm.filtered.gc.losses.includeXY.rds` — confirming that
`include_XY = TRUE` IS the production GC correction setting.

Our current `nipter_gc_correct(include_sex = TRUE)` path implements this correctly:
it interpolates the correction factor from the nearest autosomal GC value in the
fitted LOESS curve, which exactly matches `BinReads.R` lines 172–184.

What IS confirmed non-conformant: `include_XY = FALSE` is the invariant for chi
correction in both `BuildCoverageModels2025.R` (line 163) and
`CoverageProjectionSCA_Reports.R` (line 252). The `include_sex = TRUE` path in
`nipter_chi_correct()` is never used in the production pipeline and should be removed.

### Sex chromosome scoring uses a completely separate pipeline

The production pipeline (`CoverageProjectionSCA_Reports.R`) does NOT use NIPTeR's
`calculate_z_score()`, `prepare_ncv()`, or `perform_regression()` for sex
chromosomes. Sex chromosome analysis is entirely separate:

1. **Sex classification**: 3-model GMM (mclust) majority vote
   - Y unique ratio model
   - XY fraction model (2D GMM on log(X frac) vs log(Y frac))
   - Y fraction model (1D GMM)
   - plus k-mer model from niptmer tool

2. **Sex chromosome Z-scores**: computed from sex-matched flat chromosome fraction
   data frames (not NIPTeR objects), using `ZScoreMinusSelf()` (leave-one-out)

3. **Sex chromosome NCV**: custom brute-force denominator search on raw count columns
   (`NChrReads_X`, `NChrReads_Y`) against sex-matched reference sets, completely
   separate from NIPTeR's `prepare_ncv()` which works on NIPTeRSample objects

4. **Sex chromosome regression**: forward stepwise `lm` + `step()` on sex-matched
   flat data frames, NOT NIPTeR's `perform_regression()`

This means `nipter_sex.R` in its current form is targeting the wrong abstraction.
Sex chromosome analysis operates on chromosome-fraction data frames derived from the
chi-corrected controls, not on NIPTeRSample objects directly.

---

## 1. Architecture overview

The current S3 design has three problems that S7 fixes:

1. **No accessor abstraction**: `CombinedStrands` and `SeparatedStrands` dispatch
   is copy-pasted across ≥5 files. S7 generics make this a single method definition.

2. **No property validation**: correction state is an unenforced character vector.
   S7 validators enforce the allowed correction pipeline ordering.

3. **No lazy caching**: chromosome fractions are recomputed on every scoring call.
   S7 `new_property()` with a read-only computed property handles this cleanly.

---

## 2. S7 class definitions

### 2.1 Read-count leaf classes

```r
library(S7)

# ---- Correction record -------------------------------------------------------

NIPTCorrectionRecord <- new_class(
  "NIPTCorrectionRecord",
  properties = list(
    autosomal = new_property(class_character, default = "Uncorrected"),
    sex       = new_property(class_character, default = "Uncorrected")
  ),
  validator = function(self) {
    valid_auto <- c("Uncorrected", "GC_LOESS", "GC_bin", "Chi_corrected")
    if (!all(self@autosomal %in% valid_auto))
      sprintf("autosomal must be in {%s}", paste(valid_auto, collapse = ", "))
  }
)

# ---- Abstract NIPTSample -----------------------------------------------------
#
# Never instantiated directly. Provides:
#   - shared properties: sample_name, binsize, correction
#   - declared generics that subclasses must implement

NIPTSample <- new_class(
  "NIPTSample",
  abstract = TRUE,
  properties = list(
    sample_name = class_character,
    binsize     = new_property(class_integer,
                               validator = function(value)
                                 if (value < 1L) "binsize must be positive"),
    correction  = new_property(NIPTCorrectionRecord,
                               default = quote(NIPTCorrectionRecord()))
  )
)

# ---- Generics ----------------------------------------------------------------
#
# autosomal_matrix(): always returns 22 x n_bins double matrix
# sex_matrix()      : always returns  2 x n_bins double matrix (rows: X, Y)
# strand_type()     : "combined" or "separated"
# n_bins()          : ncol(autosomal_matrix(x))

autosomal_matrix <- new_generic("autosomal_matrix", "x")
sex_matrix       <- new_generic("sex_matrix",       "x")
strand_type      <- new_generic("strand_type",      "x")
n_bins           <- new_generic("n_bins",            "x")

method(n_bins, NIPTSample) <- function(x) ncol(autosomal_matrix(x))

# ---- CombinedStrandsSample ---------------------------------------------------

CombinedStrandsSample <- new_class(
  "CombinedStrandsSample",
  parent = NIPTSample,
  properties = list(
    auto_matrix = new_property(
      class_double,           # stored as double matrix (22 x n_bins)
      validator = function(value) {
        if (!is.matrix(value) || nrow(value) != 22L)
          "auto_matrix must be a 22-row matrix"
      }
    ),
    sex_matrix_ = new_property(
      class_double,           # 2 x n_bins (row 1 = X, row 2 = Y)
      validator = function(value) {
        if (!is.matrix(value) || nrow(value) != 2L)
          "sex_matrix_ must be a 2-row matrix (X, Y)"
      }
    )
  )
)

method(autosomal_matrix, CombinedStrandsSample) <- function(x) x@auto_matrix
method(sex_matrix,       CombinedStrandsSample) <- function(x) x@sex_matrix_
method(strand_type,      CombinedStrandsSample) <- function(x) "combined"

# ---- SeparatedStrandsSample --------------------------------------------------

SeparatedStrandsSample <- new_class(
  "SeparatedStrandsSample",
  parent = NIPTSample,
  properties = list(
    auto_fwd = new_property(class_double,
                            validator = function(v)
                              if (!is.matrix(v) || nrow(v) != 22L)
                                "auto_fwd must be 22-row matrix"),
    auto_rev = new_property(class_double,
                            validator = function(v)
                              if (!is.matrix(v) || nrow(v) != 22L)
                                "auto_rev must be 22-row matrix"),
    sex_fwd  = new_property(class_double,
                            validator = function(v)
                              if (!is.matrix(v) || nrow(v) != 2L)
                                "sex_fwd must be 2-row matrix"),
    sex_rev  = new_property(class_double,
                            validator = function(v)
                              if (!is.matrix(v) || nrow(v) != 2L)
                                "sex_rev must be 2-row matrix")
  )
)

method(autosomal_matrix, SeparatedStrandsSample) <- function(x) x@auto_fwd + x@auto_rev
method(sex_matrix,       SeparatedStrandsSample) <- function(x) x@sex_fwd + x@sex_rev
method(strand_type,      SeparatedStrandsSample) <- function(x) "separated"
```

### 2.2 Control group with cached fractions

```r
# ---- NIPTControlGroup --------------------------------------------------------
#
# Lazy fractions cache: computed on first access, invalidated on sample mutation.
# fractions_auto(): 22 x N matrix of autosomal chromosome fractions
# fractions_sep() : 44 x N matrix (1F..22F, 1R..22R) for SeparatedStrands

NIPTControlGroup <- new_class(
  "NIPTControlGroup",
  abstract = TRUE,
  properties = list(
    samples     = class_list,        # list of NIPTSample (same strand type)
    description = new_property(class_character, default = "General control group"),
    # Internal lazy cache — not user-facing, prefixed with dot
    .frac_cache = new_property(
      new_union(class_NULL, class_double),
      default = NULL
    )
  ),
  validator = function(self) {
    if (length(self@samples) < 2L) "need >= 2 samples"
    # all samples same strand type (enforced in constructor)
  }
)

fractions_auto <- new_generic("fractions_auto", "cg")
fractions_for_regression <- new_generic("fractions_for_regression", "cg")
n_samples <- new_generic("n_samples", "cg")
control_sample_names <- new_generic("control_sample_names", "cg")

method(n_samples, NIPTControlGroup) <- function(cg) length(cg@samples)
method(control_sample_names, NIPTControlGroup) <- function(cg)
  vapply(cg@samples, function(s) s@sample_name, character(1L))

# ---- CombinedControlGroup / SeparatedControlGroup ---------------------------

CombinedControlGroup <- new_class(
  "CombinedControlGroup",
  parent = NIPTControlGroup
)

SeparatedControlGroup <- new_class(
  "SeparatedControlGroup",
  parent = NIPTControlGroup
)

# fractions_auto: 22 x N collapsed fractions (same for both subtypes)
method(fractions_auto, CombinedControlGroup) <- function(cg) {
  if (!is.null(cg@.frac_cache)) return(cg@.frac_cache)
  mats <- lapply(cg@samples, autosomal_matrix)
  totals <- vapply(mats, sum, numeric(1L))
  fracs <- mapply(function(m, t) rowSums(m) / t, mats, totals)  # 22 x N
  cg@.frac_cache <- fracs
  fracs
}

# SeparatedControlGroup needs 44-row fractions for regression
method(fractions_auto, SeparatedControlGroup) <- function(cg) {
  if (!is.null(cg@.frac_cache)) return(cg@.frac_cache)
  mats <- lapply(cg@samples, autosomal_matrix)  # each 22 x n_bins (sum F+R)
  totals <- vapply(mats, sum, numeric(1L))
  fracs <- mapply(function(m, t) rowSums(m) / t, mats, totals)  # 22 x N
  cg@.frac_cache <- fracs
  fracs
}

method(fractions_for_regression, CombinedControlGroup) <- function(cg) {
  fractions_auto(cg)  # 22 x N
}

method(fractions_for_regression, SeparatedControlGroup) <- function(cg) {
  # 44 x N: rows are "1F".."22F","1R".."22R"
  mats_fwd <- lapply(cg@samples, function(s) s@auto_fwd)
  mats_rev <- lapply(cg@samples, function(s) s@auto_rev)
  totals <- vapply(cg@samples, function(s) sum(autosomal_matrix(s)), numeric(1L))
  fwd <- mapply(function(m, t) rowSums(m) / t, mats_fwd, totals)  # 22 x N
  rev <- mapply(function(m, t) rowSums(m) / t, mats_rev, totals)  # 22 x N
  rownames(fwd) <- paste0(1:22, "F")
  rownames(rev) <- paste0(1:22, "R")
  rbind(fwd, rev)  # 44 x N
}
```

### 2.3 NCV template

```r
NCVTemplate <- new_class(
  "NCVTemplate",
  properties = list(
    focus_chromosome  = class_integer,
    denominators      = class_character,
    ctrl_mean         = class_double,
    ctrl_sd           = class_double,
    ctrl_cv           = class_double,
    shapiro_p         = class_double,
    train_z_scores    = class_double,
    test_z_scores     = class_double,
    train_sample_names = class_character,
    test_sample_names  = class_character
  )
)
```

### 2.4 Sex models (new, not in current package)

```r
# ---- SexGMMModels ------------------------------------------------------------
# Three mclust models for sex classification (majority vote)

SexGMMModels <- new_class(
  "SexGMMModels",
  properties = list(
    yuniq_model        = class_list,   # Mclust object (1D, Y unique ratio)
    xy_fractions_model = class_list,   # Mclust object (2D, X + Y fractions)
    y_fraction_model   = class_list,   # Mclust object (1D, Y fraction)
    male_class_yuniq   = class_integer,
    male_class_xy      = class_integer,
    male_class_y       = class_integer
  )
)

sex_predict <- new_generic("sex_predict", "model")

# ---- SexNCVModel / SexNCVModelSet --------------------------------------------
# NCV denominators for sex chromosome scoring (custom, NOT NIPTeR prepare_ncv)
# Operates on flat chromosome-count data frames (NChrReads_* columns), not
# NIPTeRSample objects.

SexNCVModel <- new_class(
  "SexNCVModel",
  properties = list(
    chr_focus        = class_character,    # "NChrReads_X" or "NChrReads_Y"
    denominators     = class_character,    # e.g. c("NChrReads_1", "NChrReads_3", ...)
    reference_sex    = class_character,    # "female" or "male"
    ctrl_mean        = class_double,
    ctrl_sd          = class_double,
    ctrl_cv          = class_double,
    shapiro_p        = class_double
  )
)

SexNCVModelSet <- new_class(
  "SexNCVModelSet",
  properties = list(
    x_female = SexNCVModel,
    y_female = SexNCVModel,
    x_male   = SexNCVModel,
    y_male   = SexNCVModel
  )
)

# ---- SexRegressionModels -----------------------------------------------------
# Forward stepwise lm models for sex chromosome fraction prediction
# Fit separately per sex per target chromosome

SexRegressionModels <- new_class(
  "SexRegressionModels",
  properties = list(
    x_female = class_list,   # lm object after step()
    y_female = class_list,
    x_male   = class_list,
    y_male   = class_list
  )
)
```

### 2.5 Reference model (top-level pre-built object)

```r
# ---- NIPTReferenceModel ------------------------------------------------------
# The serialised reference object produced by nipter_build_reference().
# Replaces the untyped ControlObject list currently used in the production
# pipeline. Validates its structure on construction; provides typed accessors.

NIPTReferenceModel <- new_class(
  "NIPTReferenceModel",
  properties = list(
    # NIPTeR controls (GC-corrected, chi-corrected)
    nipter_controls         = NIPTControlGroup,
    chi_corrected_controls  = NIPTControlGroup,

    # NCV templates for autosomes (22-element list, index = chromosome number)
    ncv_templates_auto = class_list,

    # Sex-specific models
    sex_gmm_models      = SexGMMModels,
    sex_ncv_models      = SexNCVModelSet,
    sex_regression      = SexRegressionModels,

    # Reference QC data frame (one row per reference sample)
    # Columns: Sample_name, FrChrReads_1..22, FrChrReads_X/Y,
    #          NChrReads_1..22, NChrReads_X/Y, ConsensusGender,
    #          YUniqueRatioFiltered, IsRefSexOutlier, etc.
    reference_qc_df = class_data.frame,

    # Optional: WisecondorX reference (list or NULL)
    wisecondorx_ref = new_union(class_NULL, class_list),

    # Provenance
    build_date   = class_character,
    build_params = class_list
  ),
  validator = function(self) {
    if (length(self@ncv_templates_auto) != 22L)
      "ncv_templates_auto must have exactly 22 entries"
  }
)
```

---

## 3. Generics and the correction pipeline

```r
# ---- Correction generics -----------------------------------------------------

gc_correct_loess <- new_generic("gc_correct_loess", "x")
gc_correct_bin   <- new_generic("gc_correct_bin",   "x")
chi_correct      <- new_generic("chi_correct",       "x")

# x can be NIPTSample or NIPTControlGroup
# For NIPTSample methods, returns a new NIPTSample with updated correction record
# For NIPTControlGroup, returns a new NIPTControlGroup

# ---- Scoring generics --------------------------------------------------------

z_score     <- new_generic("z_score",     "sample")     # → ZScoreResult
ncv_score   <- new_generic("ncv_score",   "sample")     # → NCVResult
rbz_score   <- new_generic("rbz_score",   "sample")     # → RBZResult
sex_score   <- new_generic("sex_score",   "sample")     # → SexResult
```

All scoring generics take `sample` (NIPTSample) and `ref` (NIPTReferenceModel) as
arguments — not a raw NIPTControlGroup. This enforces that callers always work from a
typed, pre-validated reference object rather than assembling inputs ad-hoc.

---

## 4. Strand-type compatibility enforcement

The current `nipter_regression()` strand-type mismatch bug (P2-1) is structurally
impossible in the S7 design because:

1. `CombinedControlGroup` and `SeparatedControlGroup` are distinct S7 classes
2. `fractions_for_regression()` dispatches to the correct implementation
3. The scoring generics take `NIPTReferenceModel`, not raw `NIPTControlGroup` — the
   reference model is built from samples of a single strand type at construction time
4. If a mismatch is somehow attempted, the validator on `NIPTControlGroup` rejects it

No explicit strand-type guard is needed in scoring functions — the type system enforces it.

---

## 5. Removal of `nipter_sex.R` sex chromosome via-NIPTeR path

The production pipeline analyses sex chromosomes via flat chromosome-count data frames
and sex-matched reference sets (`ReferenceSexPredictionsDF`), NOT via NIPTeR's
`calculate_z_score()` or `prepare_ncv()`.

The current `nipter_sex.R` (which uses `nipter_z_score(sample, cg, chromo_focus = 23)`
style calls on sex chromosome rows) is incorrect for the production workflow. It should
be replaced with:

```r
# New: sex_score() generic dispatches on NIPTReferenceModel
method(sex_score, CombinedStrandsSample) <- function(sample, ref) {
  # 1. Classify sex via GMM majority vote (uses ref@sex_gmm_models + sample QC)
  # 2. Select sex-matched reference subset from ref@reference_qc_df
  # 3. Compute X/Y Z-scores from FrChrReads_X/Y columns (ZScoreMinusSelf)
  # 4. Compute X/Y NCV scores using ref@sex_ncv_models
  # 5. Compute X/Y regression scores using ref@sex_regression
  # Returns SexResult object
}
```

---

## 6. Correction-state validation

The `NIPTCorrectionRecord` validator enforces allowed states:
- `"Uncorrected"`
- `"GC_LOESS"` — set after `gc_correct_loess()`
- `"GC_bin"` — set after `gc_correct_bin()`
- `"Chi_corrected"` — set after `chi_correct()`

Downstream scoring functions check the correction record via `validate()`:

```r
.require_corrected <- function(sample, what = "autosomal") {
  rec <- sample@correction
  state <- if (what == "autosomal") rec@autosomal else rec@sex
  if (!any(grepl("^GC_|^Chi_", state))) {
    stop(sprintf(
      "Sample '%s' has not been GC or chi corrected (status: '%s'). ",
      sample@sample_name, state
    ), call. = FALSE)
  }
}
```

This replaces the current silent pass-through where a non-corrected sample produces
wrong Z-scores without any warning.

---

## 7. CRAN requirements

The S7 design is CRAN-compatible because:

1. `S7` is on CRAN (v0.2.1 as of 2025-11-14) — add to `Imports`
2. No use of `methods::setClass()` or `Rcpp` S4 integration required
3. `methods_register()` in `.onLoad()` is all that's required for method dispatch
4. Abstract classes (`abstract = TRUE`) are pure R — no compilation
5. Properties with validators are pure R closures
6. Minimum R version can stay at R >= 4.1.0 (S7 requires >= 3.5.0)

Add to `DESCRIPTION`:
```
Imports: S7 (>= 0.2.0), ...
```

---

## 8. Migration path (S3 → S7)

The migration should be done in phases, not in one commit:

**Phase 1 (P1)**: Define S7 classes alongside existing S3 classes. Add `as_niptsample()`
converters from the old S3 structure. Add `autosomal_matrix()`, `sex_matrix()`,
`strand_type()`, `n_bins()` generics with methods for both old S3 and new S7 classes.
Tests continue to pass.

**Phase 2 (P2)**: Port `nipter_bin_bam()` and `bed_to_nipter_sample()` to return S7
objects. Port `nipter_as_control_group()` to return S7 control group. All internal
scoring functions accept both old and new types via generics.

**Phase 3 (P3)**: Port scoring functions to use `NIPTReferenceModel`. Add
`nipter_build_reference()` that produces a validated `NIPTReferenceModel`. Remove
`include_sex = TRUE` from `nipter_chi_correct()`.

**Phase 4 (P4)**: Remove old S3 class definitions and all `inherits(x, "CombinedStrands")`
dispatch branches. Port `nipter_sex.R` to the flat-data-frame approach.
Port `nipter_gc_correct()` to S7. Full test coverage on S7 objects.

---

## 9. What NOT to change from the current implementation

- `rwisecondorx_*` (WisecondorX layer) — completely separate from NIPTeR; already has
  its own `WisecondorXReference` S3 class that works well
- `src/knn_reference.cpp` and all other Rcpp code — interface through existing R wrappers;
  S7 objects passed to Rcpp need their matrices extracted before the C++ call
- The `bam_convert()` / `nipter_bin_bam()` DuckDB/BAM reading logic — unchanged
- `generate_cohort()` — unchanged (synthetic data utility)
- Test fixtures in `inst/tinytest/` — updated to use S7 constructors but logic unchanged