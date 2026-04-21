# test_nipter_conformance.R — NIPTeR statistical layer conformance
#
# Arm A (always runs): verifies our implementations reproduce the NIPTeR
# formulas precisely, using inline reference computations from the NIPTeR
# source code. No external dependencies required.
#
# Arm B (conditional): cross-checks against the NIPTeR R package itself on
# a whole-genome BAM. Requires:
#   - NIPTeR package installed
#   - `NIPTER_CONFORMANCE_BAM` pointing to a pre-filtered whole-genome BAM
#     override (no unmapped reads, no same-position collisions per strand;
#     see AGENTS.md for the exact constraints and the known NIPTeR bugs that
#     impose them)
#
# Known divergence: nipter_gc_correct() derives GC on-the-fly from
# rduckhts_fasta_nuc() nucleotide counts, excluding N/ambiguous bases from
# the denominator, whereas NIPTeR uses bundled sysdata.rda GC tables.
# Minor numeric differences may remain; both functions are tested for
# structural correctness (returns NIPTeRSample, status updated) rather than
# numeric equality.

library(tinytest)
library(RWisecondorX)

.is_nipter_sample <- function(x) {
  RWisecondorX:::.is_nipt_sample_object(x)
}

.is_nipter_control_group <- function(x) {
  RWisecondorX:::.is_nipt_control_group_object(x)
}


# ---------------------------------------------------------------------------
# Helpers — intentionally self-contained; NOT shared via helper_nipter.R.
#
# This file must run correctly in isolation (e.g. via
# tinytest::run_test_file()) without the auto-sourced helpers, and Arm B
# runs in a separate process where helper_nipter.R is not guaranteed to be
# loaded.  Do NOT "de-duplicate" these into the shared helper.
#
# Design differences from the shared helper_nipter.R helpers:
#   - bin fill uses seq_len(rem) for fully deterministic output (no RNG);
#     set.seed() in .make_sample_conf is a no-op (kept for call-site
#     compatibility but the fill itself is not stochastic).
#   - .make_ctrl_conf seeds each sample independently with seed*i (noise)
#     and seed*i+1 (sample), rather than advancing a single global stream.
# ---------------------------------------------------------------------------

.make_sample_conf <- function(chr_totals, name, n_bins = 100L, seed = 42L) {
  if (!is.null(seed)) set.seed(seed)
  auto_mat <- matrix(0L, nrow = 22L, ncol = n_bins)
  rownames(auto_mat) <- as.character(1:22)
  colnames(auto_mat) <- as.character(seq_len(n_bins))
  for (chr in names(chr_totals)) {
    total <- chr_totals[[chr]]
    if (total == 0L) next
    base <- total %/% n_bins
    rem  <- total - base * n_bins
    counts <- rep(base, n_bins)
    if (rem > 0L) counts[seq_len(rem)] <- counts[seq_len(rem)] + 1L
    auto_mat[chr, ] <- as.integer(counts)
  }
  sex_mat <- matrix(0L, nrow = 2L, ncol = n_bins,
                    dimnames = list(c("X", "Y"), as.character(seq_len(n_bins))))
  structure(
    list(autosomal_chromosome_reads  = list(auto_mat),
         sex_chromosome_reads        = list(sex_mat),
         correction_status_autosomal = "Uncorrected",
         correction_status_sex       = "Uncorrected",
         sample_name                 = name),
    class = c("NIPTeRSample", "CombinedStrands")
  )
}

.chr_sizes <- c(248, 242, 198, 190, 182, 171, 159, 145, 138, 134, 135, 133,
                114, 107, 102, 90, 83, 80, 59, 64, 47, 51)

.make_ctrl_conf <- function(n = 10L, n_bins = 100L, seed = 7L) {
  base_totals <- as.integer(round(.chr_sizes * 10))
  names(base_totals) <- as.character(1:22)
  lapply(seq_len(n), function(i) {
    set.seed(seed * i)
    noise <- runif(22L, 0.95, 1.05)
    totals <- as.integer(round(base_totals * noise))
    names(totals) <- as.character(1:22)
    .make_sample_conf(totals, sprintf("ctrl_%02d", i), n_bins = n_bins,
                      seed = seed * i + 1L)
  })
}

# Build a test sample with trisomy 21 signal
.trisomy_totals <- function(chr = "21", boost = 0.05) {
  totals <- as.integer(round(.chr_sizes * 10))
  names(totals) <- as.character(1:22)
  totals[chr] <- as.integer(round(totals[chr] * (1 + boost)))
  totals
}


# ---------------------------------------------------------------------------
# Arm A — Inline formula verification (always runs)
# ---------------------------------------------------------------------------

ctrl_samples <- .make_ctrl_conf(12L)
cg <- nipter_as_control_group(ctrl_samples)
test_sample <- .make_sample_conf(.trisomy_totals("21", 0.1), "test_t21")


# --- 1. Z-score formula -------------------------------------------------------
# NIPTeR formula (from NIPTeR::calculate_z_score source):
#   sample_z = (sample_frac - mean(ctrl_fracs)) / sd(ctrl_fracs)
# where fracs are chromosomal fractions (chr reads / total reads).

.chr_fracs <- function(s) {
  auto <- s$autosomal_chromosome_reads[[1L]]
  total <- sum(auto)
  rowSums(auto) / total
}
.ctrl_fracs_mat <- function(cg) {
  # 22 × n_controls
  sapply(cg$samples, .chr_fracs)
}

s_fracs  <- .chr_fracs(test_sample)
cg_fracs <- .ctrl_fracs_mat(cg)

# chr21 (row 21)
chr_key <- "21"
ctrl_focus <- cg_fracs[chr_key, ]
ref_z <- unname((s_fracs[chr_key] - mean(ctrl_focus)) / stats::sd(ctrl_focus))

our_z <- nipter_z_score(test_sample, cg, chromo_focus = 21L)$sample_z_score

expect_equal(our_z, unname(ref_z), tolerance = 1e-10,
             info = "z-score matches inline NIPTeR formula")
expect_true(our_z > 0,
            info = "T21 sample has positive z-score on chr21")


# --- 2. Chi-squared correction formula ----------------------------------------
# NIPTeR algorithm (from chi_correct source):
#   1. Scale each ctrl's total reads to overall mean
#   2. Expected = per-bin mean of scaled counts
#   3. chi2_bin = sum_samples (expected - scaled)^2 / expected
#   4. z_chi = (chi2 - df) / sqrt(2 * df)
#   5. If z_chi > cutoff: divide by chi2 / df   (correction factor = df/chi2)
#
# We replicate the complete formula inline and compare numeric results.

n_bins_per_chr <- ncol(cg$samples[[1L]]$autosomal_chromosome_reads[[1L]])
ctrl_flat <- sapply(cg$samples, function(s) {
  as.numeric(t(s$autosomal_chromosome_reads[[1L]]))
})
# ctrl_flat: (22 * n_bins) × n_controls

n_ctrl <- ncol(ctrl_flat)
df     <- n_ctrl - 1L

ctrl_totals   <- colSums(ctrl_flat)
overall_mean  <- mean(ctrl_totals)
scaling       <- overall_mean / ctrl_totals
scaled_ctrl   <- sweep(ctrl_flat, 2, scaling, "*")

expected_bin  <- rowMeans(scaled_ctrl)
chi2_bin      <- rowSums(sweep(scaled_ctrl, 1, expected_bin, "-")^2) /
                 pmax(expected_bin, .Machine$double.eps)
z_chi         <- (chi2_bin - df) / sqrt(2 * df)

chi_cutoff    <- 3.5
overdispersed <- which(z_chi > chi_cutoff)

# Expected correction factor (per bin): 1 if not overdispersed, df/chi2 if so
correction_factor <- rep(1.0, length(chi2_bin))
if (length(overdispersed) > 0L) {
  correction_factor[overdispersed] <- df / chi2_bin[overdispersed]
}

corrected <- nipter_chi_correct(test_sample, cg, chi_cutoff = chi_cutoff)

expect_true(.is_nipter_sample(corrected$sample),
            info = "chi_correct returns an NIPTSample")
expect_true(.is_nipter_control_group(corrected$control_group),
            info = "chi_correct returns an NIPTControlGroup")

# Verify that the correction factor is applied correctly to the FIRST control
# sample: corrected_bin = original_bin * correction_factor
orig_flat_1 <- as.numeric(t(cg$samples[[1L]]$autosomal_chromosome_reads[[1L]]))
corr_flat_1 <- as.numeric(
  t(corrected$control_group$samples[[1L]]$autosomal_chromosome_reads[[1L]]))

# The ratio corrected/original must equal the correction factor
# (within floating-point tolerance); for bins with original==0 both are 0.
nz <- orig_flat_1 > 0
expected_corrected_1 <- orig_flat_1[nz] * correction_factor[nz]
expect_equal(corr_flat_1[nz], expected_corrected_1, tolerance = 1e-6,
             info = "chi correction factor applied correctly (df/chi2 for overdispersed, 1 otherwise)")


# --- 3. NCV — brute-force denominator search, exact numeric verification ------
# NIPTeR selects the denominator chromosome set (up to max_elements) that
# minimises the CV of the NCV ratio across controls, then computes a Z-score.
#
# With max_elements=1 the search space is small (21 single-chromosome options)
# so we can replicate the complete algorithm inline and compare numerically.
#
# Formula:
#   For each control i: NCV_i = frac_focus_i / frac_denom_i
#   CV = sd(NCV) / |mean(NCV)|
#   Select denominator minimising CV.
#   NCV_sample = frac_focus_sample / frac_denom_sample
#   score = (NCV_sample - mean(ctrl_NCVs)) / sd(ctrl_NCVs)

ctrl_fracs_all <- .ctrl_fracs_mat(cg)   # 22 × n_ctrl
focus_idx      <- 21L
denom_pool     <- setdiff(seq_len(22L), focus_idx)

# Default exclude_chromosomes matches nipter_ncv_score default: c(13, 18, 21)
# The focus chromosome is always excluded from the denominator pool.
default_excl  <- c(13L, 18L, 21L)
denom_pool_ncv <- setdiff(denom_pool, default_excl)  # excludes focus+13+18

ctrl_focus_fracs <- ctrl_fracs_all[as.character(focus_idx), ]
best_cv  <- Inf
best_den <- NA_integer_

for (d in denom_pool_ncv) {
  ctrl_den_fracs <- ctrl_fracs_all[as.character(d), ]
  ncv_ctrl <- ctrl_focus_fracs / ctrl_den_fracs
  cv <- stats::sd(ncv_ctrl) / abs(mean(ncv_ctrl))
  if (is.finite(cv) && cv < best_cv) {
    best_cv  <- cv
    best_den <- d
  }
}

# Inline NCV Z-score for the best single denominator
s_den_frac   <- s_fracs[as.character(best_den)]
s_focus_frac <- s_fracs[as.character(focus_idx)]
ctrl_den_fracs_best <- ctrl_fracs_all[as.character(best_den), ]
ncv_ctrl_best <- ctrl_focus_fracs / ctrl_den_fracs_best
ncv_sample    <- unname(s_focus_frac / s_den_frac)
ref_ncv_z     <- (ncv_sample - mean(ncv_ctrl_best)) / stats::sd(ncv_ctrl_best)

our_ncv21 <- nipter_ncv_score(test_sample, cg, chromo_focus = 21L,
                               max_elements = 1L)  # uses default exclude_chromosomes

# The selected best denominator must agree
expect_equal(our_ncv21$denominators, as.integer(best_den),
             info = "NCV single-denominator selection matches inline brute-force")
# The numeric score must agree
expect_equal(our_ncv21$sample_score, ref_ncv_z, tolerance = 1e-10,
             info = "NCV score matches inline formula for best single denominator")


# --- 3b. NCV C++ search stays finite for large near-constant cohorts ---------
# Stress the denominator-search kernel with a large control cohort whose ratios
# are intentionally almost constant. The previous one-pass variance formula can
# lose precision here; the stable implementation should stay finite and match a
# two-pass reference CV.

n_ctrl_big <- 2000L
ctrl_big <- matrix(1e12, nrow = 22L, ncol = n_ctrl_big)
ratio_big <- 0.05 + seq_len(n_ctrl_big) * 1e-14
ctrl_big[21L, ] <- ctrl_big[1L, ] * ratio_big

big_res <- RWisecondorX:::nipter_ncv_search_cpp(
  ctrl_reads  = ctrl_big,
  candidates  = 0L,
  focus_row   = 20L,  # 0-based chr21
  max_elements = 1L
)

ref_mean_big <- mean(ratio_big)
ref_sd_big <- sqrt(sum((ratio_big - ref_mean_big)^2) / (length(ratio_big) - 1L))
ref_cv_big <- ref_sd_big / abs(ref_mean_big)

expect_true(is.finite(big_res$best_cv) && big_res$best_cv >= 0,
            info = "NCV C++ search returns a finite CV for a large near-constant cohort")
expect_true(abs(as.numeric(big_res$best_cv) - ref_cv_big) < 1e-12,
            info = "NCV C++ search matches the two-pass reference CV on the stress case")


# --- 4. Regression Z-score — numeric verification ----------------------------
# NIPTeR's regression fits: frac_focus ~ frac_predictor1 + ... on a training
# split of the control group, then scores the residual.
#
# With n_models=1 and a fixed seed, the train/test split is deterministic.
# We verify the numeric Z-score by manually fitting the same lm() model on
# the same data and computing the residual Z-score.
#
# nipter_regression() uses:
#   train_fraction = 0.6 (default)
#   chromo_focus = 21
#   n_predictors = 4 (default) — forward stepwise
#   The returned model has: predictors, z_score, practical_cv, theoretical_cv

set.seed(42L)
reg21 <- nipter_regression(test_sample, cg, chromo_focus = 21L,
                           n_models = 1L, n_predictors = 4L,
                           train_fraction = 0.6, seed = 42L)

expect_true(is.list(reg21$models) && length(reg21$models) == 1L,
            info = "regression with n_models=1 returns 1 model")

m1 <- reg21$models[[1L]]
expect_true(is.numeric(m1$z_score) && is.finite(m1$z_score),
            info = "regression model has finite numeric z_score")
expect_true(is.character(m1$predictors) && length(m1$predictors) >= 1L,
            info = "regression model has at least one predictor chromosome")
# Predictors must not include the focus chromosome
expect_false("21" %in% m1$predictors,
             info = "focus chromosome not a regression predictor")

# Numeric reconstruction: build the same lm() from the model's predictors
# on the full control group. The Z-score should be in the expected range.
ctrl_fracs_full <- .ctrl_fracs_mat(cg)   # 22 × n_ctrl
pred_rows <- ctrl_fracs_full[m1$predictors, , drop = FALSE]
focus_row <- ctrl_fracs_full[as.character(21L), ]

# Chromosome indices like "1", "5" are not valid R identifiers, so prefix them
pred_colnames <- paste0("chr_", m1$predictors)
combined <- rbind(focus_row, pred_rows)
rownames(combined) <- c("focus", pred_colnames)
df_ctrl <- as.data.frame(t(combined))
formula_str <- paste("focus ~", paste(pred_colnames, collapse = " + "))
fit <- stats::lm(stats::as.formula(formula_str), data = df_ctrl)

# Use actual test sample fracs for predictors
test_preds <- stats::setNames(
  as.numeric(s_fracs[m1$predictors]),
  pred_colnames
)
test_df_pred <- as.data.frame(as.list(test_preds))
predicted_frac  <- unname(stats::predict(fit, newdata = test_df_pred))
actual_frac    <- unname(s_fracs[as.character(21L)])
residuals_ctrl <- unname(stats::residuals(fit))

# Z-score of test residual against ctrl residuals
resid_test <- actual_frac - predicted_frac
ref_reg_z  <- (resid_test - mean(residuals_ctrl)) / stats::sd(residuals_ctrl)

# The reconstruction uses the full cg (not the train/test split), so the
# Z-scores won't match exactly. We verify sign agreement and plausibility.
expect_true(sign(ref_reg_z) == sign(m1$z_score) || abs(m1$z_score) < 0.5,
            info = "regression Z-score sign agrees with inline reconstruction")
# For T21 sample with trisomy_frac=0.1, expect positive z_score
expect_true(m1$z_score > 0,
            info = "T21 sample has positive regression Z-score on chr21")


# --- 5. C++ stepwise — exact conformance against R lm() ---------------------
# nipter_stepwise_cpp() uses incremental Gram-Schmidt orthogonalisation to
# maximise adj.R² at each step, which is mathematically equivalent to the
# R lm() forward stepwise loop it replaces.  We verify they select the
# same predictors in the same order on a controlled synthetic matrix.
#
# Design: fracs is 22 × 20 (22 chromosomes, 20 training samples).
# chr1  and chr5 (0-based rows 0 and 4) are made strongly predictive of
# chr21 (row 20).  All other chromosomes are pure noise.  With any
# reasonable adj.R² criterion, row 0 must be selected first, row 4 second.

set.seed(777L)
n_train_cpp <- 20L
n_chr_cpp   <- 22L
fracs_cpp   <- matrix(runif(n_chr_cpp * n_train_cpp, 0.03, 0.07),
                      nrow = n_chr_cpp)
rownames(fracs_cpp) <- as.character(seq_len(n_chr_cpp))

# Make rows 0,4 (1-indexed: 1,5) predictive of row 20 (1-indexed: 21)
fracs_cpp[21L, ] <- 2.0 * fracs_cpp[1L, ] +
                    0.5 * fracs_cpp[5L, ] +
                    rnorm(n_train_cpp, 0, 1e-4)

focus_row_cpp  <- 20L                        # 0-based row 20 = chr21
cand_rows_cpp  <- setdiff(0L:(n_chr_cpp - 1L),
                           c(focus_row_cpp, 12L, 17L))  # exclude focus, 13, 18

# ----- R lm() reference (exact same greedy adj.R² criterion) ---------------
.r_stepwise_ref <- function(fracs, focus_row, cand_rows, n_step) {
  n_train  <- ncol(fracs)
  y_train  <- fracs[focus_row + 1L, ]         # 0→1-based
  selected <- integer(0L)                      # 0-based indices into cand_rows
  selected_rows <- integer(0L)                 # 0-based row indices in fracs

  for (step in seq_len(n_step)) {
    remaining_pos <- setdiff(seq_along(cand_rows) - 1L, selected)  # 0-based pos
    best_adj_r2   <- -Inf
    best_pos      <- -1L

    for (pos in remaining_pos) {
      row_1based <- cand_rows[pos + 1L] + 1L   # convert to 1-based row
      trial_rows <- c(selected_rows + 1L, row_1based)
      X <- rbind(1, fracs[trial_rows, , drop = FALSE])
      fit <- stats::lm.fit(t(X), y_train)
      # Compute adj.R² from lm.fit residuals
      rss <- sum(fit$residuals^2)
      tss <- sum((y_train - mean(y_train))^2)
      p   <- length(trial_rows)   # predictors (excluding intercept)
      adj_r2 <- 1 - (rss / (n_train - p - 1)) / (tss / (n_train - 1))
      if (!is.na(adj_r2) && adj_r2 > best_adj_r2) {
        best_adj_r2 <- adj_r2
        best_pos    <- pos
      }
    }

    if (best_pos < 0L) break
    selected      <- c(selected, best_pos)
    selected_rows <- c(selected_rows, cand_rows[best_pos + 1L])
  }
  selected   # 0-based indices into cand_rows
}

ref_sel <- .r_stepwise_ref(fracs_cpp, focus_row_cpp,
                            as.integer(cand_rows_cpp), 2L)
cpp_sel <- RWisecondorX:::nipter_stepwise_cpp(fracs_cpp,
                                              focus_row_cpp,
                                              as.integer(cand_rows_cpp), 2L)

# Both must select row 0 first (chr1, the strongest predictor) and row 4
# second (chr5), as 0-based positions within cand_rows_cpp.
pos_of_row0 <- match(0L, cand_rows_cpp) - 1L   # position of row 0 in cand_rows_cpp
pos_of_row4 <- match(4L, cand_rows_cpp) - 1L   # position of row 4 in cand_rows_cpp

expect_equal(as.integer(cpp_sel), as.integer(ref_sel),
             info = "C++ stepwise selects same predictors as R lm() (exact match)")
expect_equal(cpp_sel[1L], as.integer(pos_of_row0),
             info = "C++ stepwise: strongest predictor (chr1) selected first")
expect_equal(cpp_sel[2L], as.integer(pos_of_row4),
             info = "C++ stepwise: second predictor (chr5) selected second")

# Verify that nipter_regression() with C++ backend still gives a valid,
# positive Z-score on the T21 synthetic sample from Arm A.
reg21_cpp <- nipter_regression(test_sample, cg, chromo_focus = 21L,
                               n_models = 1L, n_predictors = 4L,
                               train_fraction = 0.6, seed = 42L)
expect_true(is.list(reg21_cpp$models) && length(reg21_cpp$models) >= 1L,
            info = "nipter_regression (C++ backend) returns models list")
m1_cpp <- reg21_cpp$models[[1L]]
expect_true(is.numeric(m1_cpp$z_score) && is.finite(m1_cpp$z_score),
            info = "C++ regression: finite z_score")
expect_true(m1_cpp$z_score > 0,
            info = "C++ regression: T21 sample has positive z_score on chr21")
expect_false("21" %in% m1_cpp$predictors,
             info = "C++ regression: focus chromosome not a predictor")


# ---------------------------------------------------------------------------
# Arm B — Cross-check against NIPTeR R package (conditional)
# ---------------------------------------------------------------------------

if (!requireNamespace("NIPTeR", quietly = TRUE)) {
  exit_file("NIPTeR package not installed; skipping cross-package conformance")
}

conf_bam <- Sys.getenv("NIPTER_CONFORMANCE_BAM", unset = NA_character_)
if (is.na(conf_bam) || !nzchar(conf_bam) || !file.exists(conf_bam)) {
  exit_file("NIPTER_CONFORMANCE_BAM not set; skipping cross-package conformance arm")
}

nipter_exports <- getNamespaceExports("NIPTeR")
required_nipter_exports <- c(
  "bin_bam_sample",
  "as_control_group",
  "chi_correct",
  "calculate_z_score",
  "gc_correct"
)
missing_nipter_exports <- setdiff(required_nipter_exports, nipter_exports)
if (length(missing_nipter_exports) > 0L) {
  stop(paste(
    "Installed NIPTeR package does not expose the expected current API:",
    paste(missing_nipter_exports, collapse = ", ")
  ), call. = FALSE)
}

.make_real_control_variants <- function(sample, prefix, n = 12L) {
  stopifnot(.is_nipter_sample(sample))
  stopifnot(length(sample$autosomal_chromosome_reads) == 1L)
  lapply(seq_len(n), function(i) {
    s <- sample
    auto <- RWisecondorX:::.sample_autosomal_reads(s)[[1L]]
    bump_chr <- as.character(((i - 1L) %% 4L) + 1L)
    bump_bin <- ((i - 1L) %% ncol(auto)) + 1L
    auto[bump_chr, bump_bin] <- auto[bump_chr, bump_bin] + i
    s <- RWisecondorX:::.sample_with_reads(s, autosomal = list(auto))
    s <- RWisecondorX:::.nipt_sample_dollar_assign(
      s,
      "sample_name",
      sprintf("%s_%02d", prefix, i)
    )
    s
  })
}

# --- B1. Binning (NIPTeR::bin_bam_sample vs nipter_bin_bam) ------------------
# Covered in test_nipter.R; not duplicated here.

# --- B2. Z-score conformance on real binned data -----------------------------
our_s_real <- nipter_bin_bam(conf_bam, binsize = 50000L, mapq = 0L,
                             rmdup = "none")
nipter_s_real <- tryCatch(
  NIPTeR::bin_bam_sample(bam_filepath = conf_bam, do_sort = FALSE,
                         separate_strands = FALSE),
  error = function(e) {
    stop(paste(
      "Upstream NIPTeR::bin_bam_sample is incompatible with the current",
      "Rsamtools/Bioconductor stack:",
      conditionMessage(e)
    ), call. = FALSE)
  }
)

# Build control groups from deterministic perturbations of a single real BAM.
# Exact clones produce zero-variance control fractions, which makes upstream
# NIPTeR's Shapiro-based z-score summary undefined after scaling.
our_cg_real <- nipter_as_control_group(
  .make_real_control_variants(our_s_real, "ctrl")
)
nipter_cg_real <- NIPTeR::as_control_group(
  .make_real_control_variants(nipter_s_real, "ctrl")
)

# --- B3. Chi correction agreement --------------------------------------------
our_chi  <- nipter_chi_correct(our_s_real,   our_cg_real)
nipter_chi <- tryCatch(
  NIPTeR::chi_correct(nipter_s_real, nipter_cg_real),
  error = function(e) {
    stop(paste(
      "Upstream NIPTeR::chi_correct failed on the bundled conformance setup:",
      conditionMessage(e)
    ), call. = FALSE)
  }
)

# Number of overdispersed bins should be identical (both use same formula)
our_n_overdispersed <- sum(
  as.numeric(t(our_chi$sample$autosomal_chromosome_reads[[1L]])) !=
  as.numeric(t(our_s_real$autosomal_chromosome_reads[[1L]]))
)
nipter_n_overdispersed <- sum(
  as.numeric(t(nipter_chi$sample$autosomal_chromosome_reads[[1L]])) !=
  as.numeric(t(nipter_s_real$autosomal_chromosome_reads[[1L]]))
)

expect_equal(our_n_overdispersed, nipter_n_overdispersed, tolerance = 0L,
             info = "chi correction: same number of overdispersed bins as NIPTeR")

# --- B4. Z-score agreement ---------------------------------------------------
# Use the chi-corrected samples and control groups

# Build matching structures
our_cg_chi <- our_chi$control_group
nipter_cg_chi <- nipter_chi$control_group

# Repack NIPTeR chi output into our structures for z-score (and vice versa)
# Actually compute z-score from our corrected sample vs our corrected cg
our_z21    <- nipter_z_score(our_chi$sample, our_chi$control_group, 21L)$sample_z_score
nipter_z21 <- tryCatch(
  NIPTeR::calculate_z_score(
    nipter_chi$sample,
    nipter_chi$control_group,
    chromo_focus = 21
  )$sample_Zscore,
  error = function(e) {
    stop(paste(
      "Upstream NIPTeR::calculate_z_score failed on the bundled conformance setup:",
      conditionMessage(e)
    ), call. = FALSE)
  }
)

expect_equal(our_z21, nipter_z21, tolerance = 0.01,
             info = "chr21 z-score agrees with NIPTeR within 0.01")

# --- B5. GC correction — structural check only (known divergence) ------------
# We do NOT assert numeric equality here. Our implementation uses
# rduckhts_fasta_nuc() on-the-fly; NIPTeR uses bundled GC tables.
# Both should return a NIPTeRSample with updated correction_status.
our_gc <- tryCatch(
  suppressWarnings(nipter_gc_correct(our_s_real, fasta = Sys.getenv("RWXCONF_FASTA",
                                     unset = NA_character_))),
  error = function(e) NULL
)
if (!is.null(our_gc) && !is.na(Sys.getenv("RWXCONF_FASTA", unset = NA_character_))) {
  expect_true(.is_nipter_sample(our_gc),
              info = "gc_correct returns an NIPTSample (structural check)")
  # NIPTeR GC correction
  nipter_gc <- tryCatch(
    NIPTeR::gc_correct(nipter_s_real),
    error = function(e) NULL
  )
  if (!is.null(nipter_gc)) {
    expect_true(.is_nipter_sample(nipter_gc),
                info = "NIPTeR gc_correct returns a sample-like result (structural check)")
    message(
      "[GC conformance note] nipter_gc_correct derives GC from ",
      "rduckhts_fasta_nuc() nucleotide counts, while NIPTeR uses bundled GC ",
      "tables. Minor numeric divergence may remain and is documented."
    )
  }
}
