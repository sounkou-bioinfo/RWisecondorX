library(tinytest)

# ---------------------------------------------------------------------------
# Tests for nipter_sex.R — sex prediction via mclust GMM
#
# These tests use synthetic NIPTeRSample objects with realistic sex
# chromosome signal (no BAM fixture required).
# ---------------------------------------------------------------------------

.source_candidate <- function(path) {
  candidates <- c(path, file.path("..", "..", path))
  candidates[file.exists(candidates)][1L]
}

src_files <- c("R/convert.R", "R/nipter_bin.R", "R/nipter_control.R",
               "R/nipter_gc.R", "R/nipter_chi.R", "R/nipter_score.R",
               "R/nipter_regression.R", "R/nipter_sex.R")
src_paths <- vapply(src_files, .source_candidate, character(1L))

if (all(!is.na(src_paths))) {
  for (p in src_paths) source(p)
} else if (requireNamespace("RWisecondorX", quietly = TRUE)) {
  for (fn in c("nipter_as_control_group", "nipter_sex_model",
               "nipter_predict_sex", "nipter_sex_model_y_unique",
               "nipter_y_unique_ratio")) {
    assign(fn, getExportedValue("RWisecondorX", fn))
  }
} else {
  stop("Unable to locate source files or installed RWisecondorX.", call. = FALSE)
}

# Skip if mclust is not available
if (!requireNamespace("mclust", quietly = TRUE)) {
  exit_file("mclust not installed — skipping sex prediction tests")
}

# ---------------------------------------------------------------------------
# Helpers: build synthetic samples with sex chromosome signal
# ---------------------------------------------------------------------------

# Approximate relative chromosome sizes (hg38, Mb, rounded)
.chr_sizes <- c(248, 242, 198, 190, 182, 171, 159, 145, 138, 134, 135, 133,
                114, 107, 102, 90, 83, 80, 59, 64, 47, 51)

.make_sex_sample <- function(name, sex = c("female", "male"),
                             n_bins = 100L, scale = 1, seed = NULL) {
  sex <- match.arg(sex)
  if (!is.null(seed)) set.seed(seed)

  # Autosomal reads (proportional to chromosome size)
  chr_totals <- as.integer(round(.chr_sizes * 10 * scale *
                                  runif(22, 0.95, 1.05)))
  auto_mat <- matrix(0L, nrow = 22L, ncol = n_bins)
  rownames(auto_mat) <- as.character(1:22)
  colnames(auto_mat) <- as.character(seq_len(n_bins))
  for (i in 1:22) {
    total <- chr_totals[i]
    base <- total %/% n_bins
    remainder <- total - base * n_bins
    counts <- rep(base, n_bins)
    if (remainder > 0L) counts[sample(n_bins, remainder)] <- counts[sample(n_bins, remainder)] + 1L
    auto_mat[i, ] <- as.integer(counts)
  }
  auto_total <- sum(auto_mat)

  # Sex chromosome reads — the signal that distinguishes male/female
  # Female (XX): X ~ 3-4% of autosomal, Y ~ 0.001%
  # Male (XY):   X ~ 1.5-2% of autosomal, Y ~ 0.2-0.5%
  if (sex == "female") {
    x_total <- as.integer(round(auto_total * runif(1, 0.030, 0.040)))
    y_total <- as.integer(round(auto_total * runif(1, 0.0001, 0.001)))
  } else {
    x_total <- as.integer(round(auto_total * runif(1, 0.015, 0.020)))
    y_total <- as.integer(round(auto_total * runif(1, 0.002, 0.005)))
  }

  sex_mat <- matrix(0L, nrow = 2L, ncol = n_bins)
  rownames(sex_mat) <- c("X", "Y")
  colnames(sex_mat) <- as.character(seq_len(n_bins))

  # Distribute X reads
  base_x <- x_total %/% n_bins
  rem_x  <- x_total - base_x * n_bins
  x_counts <- rep(base_x, n_bins)
  if (rem_x > 0L) x_counts[sample(n_bins, rem_x)] <- x_counts[sample(n_bins, rem_x)] + 1L
  sex_mat["X", ] <- as.integer(x_counts)

  # Distribute Y reads
  base_y <- y_total %/% n_bins
  rem_y  <- y_total - base_y * n_bins
  y_counts <- rep(base_y, n_bins)
  if (rem_y > 0L && y_total > 0L) y_counts[sample(n_bins, min(rem_y, n_bins))] <- y_counts[sample(n_bins, min(rem_y, n_bins))] + 1L
  sex_mat["Y", ] <- as.integer(y_counts)

  structure(
    list(
      autosomal_chromosome_reads  = list(auto_mat),
      sex_chromosome_reads        = list(sex_mat),
      correction_status_autosomal = "Uncorrected",
      correction_status_sex       = "Uncorrected",
      sample_name                 = name
    ),
    class = c("NIPTeRSample", "CombinedStrands")
  )
}

# Build a mixed-sex control group
.make_sex_control_set <- function(n_female = 15L, n_male = 15L,
                                  scale = 1, seed = 42L) {
  set.seed(seed)
  samples <- vector("list", n_female + n_male)
  for (i in seq_len(n_female)) {
    samples[[i]] <- .make_sex_sample(paste0("female_", i), "female",
                                     seed = seed + i)
  }
  for (i in seq_len(n_male)) {
    samples[[n_female + i]] <- .make_sex_sample(paste0("male_", i), "male",
                                                seed = seed + n_female + i)
  }
  samples
}

# ===================================================================
# BUILD SEX MODELS
# ===================================================================

sex_samples <- .make_sex_control_set(15L, 15L, seed = 123L)
sex_cg <- nipter_as_control_group(sex_samples)

# --- nipter_sex_model: Y-fraction method ---

model_y <- nipter_sex_model(sex_cg, method = "y_fraction")

expect_true(inherits(model_y, "NIPTeRSexModel"),
            info = "sex model has correct class")
expect_identical(model_y$method, "y_fraction",
                 info = "sex model method is y_fraction")
expect_true(model_y$male_cluster %in% c(1L, 2L),
            info = "male_cluster is 1 or 2")
expect_identical(length(model_y$classifications), 30L,
                 info = "classifications has 30 entries")
expect_true(all(model_y$classifications %in% c("male", "female")),
            info = "all classifications are male/female")

# The model should correctly classify most samples
true_sex <- c(rep("female", 15L), rep("male", 15L))
names(true_sex) <- names(model_y$classifications)
accuracy_y <- mean(model_y$classifications == true_sex)
expect_true(accuracy_y >= 0.80,
            info = paste0("Y-fraction model accuracy >= 80% (got ",
                         round(accuracy_y * 100), "%)"))

# Fractions matrix should be n_samples x 2
expect_identical(nrow(model_y$fractions), 30L,
                 info = "fractions matrix has 30 rows")
expect_identical(ncol(model_y$fractions), 2L,
                 info = "fractions matrix has 2 columns")
expect_identical(colnames(model_y$fractions),
                 c("x_fraction", "y_fraction"),
                 info = "fractions columns named correctly")

# Y fractions for males should be higher than females on average
male_y <- model_y$fractions[grepl("^male_", rownames(model_y$fractions)),
                            "y_fraction"]
female_y <- model_y$fractions[grepl("^female_", rownames(model_y$fractions)),
                              "y_fraction"]
expect_true(mean(male_y) > mean(female_y),
            info = "male Y-fractions higher than female Y-fractions")


# --- nipter_sex_model: XY-fraction method ---

model_xy <- nipter_sex_model(sex_cg, method = "xy_fraction")

expect_true(inherits(model_xy, "NIPTeRSexModel"),
            info = "XY model has correct class")
expect_identical(model_xy$method, "xy_fraction",
                 info = "sex model method is xy_fraction")

accuracy_xy <- mean(model_xy$classifications == true_sex)
expect_true(accuracy_xy >= 0.80,
            info = paste0("XY-fraction model accuracy >= 80% (got ",
                         round(accuracy_xy * 100), "%)"))


# ===================================================================
# PREDICT SEX
# ===================================================================

# Predict on a known female
test_female <- .make_sex_sample("test_female", "female", seed = 999L)
pred_f <- nipter_predict_sex(test_female, model_y)

expect_true(inherits(pred_f, "NIPTeRSexPrediction"),
            info = "prediction has correct class")
expect_true(pred_f$prediction %in% c("male", "female"),
            info = "prediction is male or female")
expect_identical(pred_f$prediction, "female",
                 info = "female sample predicted as female")
expect_identical(pred_f$sample_name, "test_female",
                 info = "sample_name preserved in prediction")
expect_true(is.numeric(pred_f$fractions),
            info = "fractions is numeric")
expect_identical(names(pred_f$fractions), c("x_fraction", "y_fraction"),
                 info = "fractions has correct names")

# Predict on a known male
test_male <- .make_sex_sample("test_male", "male", seed = 888L)
pred_m <- nipter_predict_sex(test_male, model_y)
expect_identical(pred_m$prediction, "male",
                 info = "male sample predicted as male")

# --- Consensus prediction from multiple models ---

pred_consensus <- nipter_predict_sex(test_female, model_y, model_xy)
expect_identical(length(pred_consensus$model_predictions), 2L,
                 info = "consensus has 2 model predictions")
expect_identical(pred_consensus$prediction, "female",
                 info = "consensus predicts female for female sample")

pred_consensus_m <- nipter_predict_sex(test_male, model_y, model_xy)
expect_identical(pred_consensus_m$prediction, "male",
                 info = "consensus predicts male for male sample")

# --- Predict with list of models ---

pred_list <- nipter_predict_sex(test_male, list(model_y, model_xy))
expect_identical(pred_list$prediction, "male",
                 info = "list-of-models prediction works")


# ===================================================================
# EDGE CASES
# ===================================================================

# Too few samples for sex model
small_cg <- nipter_as_control_group(sex_samples[1:3])
expect_error(nipter_sex_model(small_cg),
             info = "rejects control group with < 4 samples")

# No model provided to predict
expect_error(nipter_predict_sex(test_female),
             info = "rejects predict with no models")


# ===================================================================
# Y-UNIQUE MODEL (nipter_sex_model_y_unique)
# ===================================================================

# Simulate Y-unique ratios: females near 0, males near 0.003
set.seed(42L)
n_f <- 15L
n_m <- 15L
female_ratios <- abs(rnorm(n_f, mean = 0.00005, sd = 0.00002))
male_ratios   <- abs(rnorm(n_m, mean = 0.003,   sd = 0.0005))
all_ratios <- c(female_ratios, male_ratios)
names(all_ratios) <- c(paste0("female_", seq_len(n_f)),
                        paste0("male_",   seq_len(n_m)))
true_sex_yu <- c(rep("female", n_f), rep("male", n_m))
names(true_sex_yu) <- names(all_ratios)

model_yu <- nipter_sex_model_y_unique(all_ratios)

expect_true(inherits(model_yu, "NIPTeRSexModel"),
            info = "y_unique model has correct class")
expect_identical(model_yu$method, "y_unique",
                 info = "y_unique model method is correct")
expect_true(model_yu$male_cluster %in% c(1L, 2L),
            info = "y_unique male_cluster is 1 or 2")
expect_identical(length(model_yu$classifications), 30L,
                 info = "y_unique classifications has 30 entries")

# Classification accuracy
accuracy_yu <- mean(model_yu$classifications == true_sex_yu)
expect_true(accuracy_yu >= 0.90,
            info = paste0("Y-unique model accuracy >= 90% (got ",
                         round(accuracy_yu * 100), "%)"))

# Male cluster should have higher ratio values
male_idx <- model_yu$classifications == "male"
expect_true(median(all_ratios[male_idx]) > median(all_ratios[!male_idx]),
            info = "y_unique male cluster has higher ratios")

# Edge case: too few samples
expect_error(nipter_sex_model_y_unique(c(a = 0.001, b = 0.002, c = 0.003)),
             info = "y_unique rejects < 4 samples")

# Unnamed ratios get default names
unnamed_ratios <- unname(all_ratios)
model_unnamed <- nipter_sex_model_y_unique(unnamed_ratios)
expect_true(!is.null(names(model_unnamed$classifications)),
            info = "unnamed ratios get default sample names")


# ===================================================================
# PREDICT SEX WITH Y-UNIQUE MODEL
# ===================================================================

# Predict female sample using y_unique model via nipter_predict_sex
pred_yu_f <- nipter_predict_sex(test_female, model_yu,
                                y_unique_ratio = 0.00003)
expect_identical(pred_yu_f$prediction, "female",
                 info = "y_unique model predicts female correctly")
expect_identical(pred_yu_f$model_predictions[["y_unique"]], "female",
                 info = "y_unique model_predictions entry correct for female")

# Predict male sample using y_unique model
pred_yu_m <- nipter_predict_sex(test_male, model_yu,
                                y_unique_ratio = 0.004)
expect_identical(pred_yu_m$prediction, "male",
                 info = "y_unique model predicts male correctly")

# Three-model consensus with y_unique
pred_3 <- nipter_predict_sex(test_female, model_y, model_xy, model_yu,
                             y_unique_ratio = 0.00003)
expect_identical(length(pred_3$model_predictions), 3L,
                 info = "3-model consensus has 3 predictions")
expect_identical(pred_3$prediction, "female",
                 info = "3-model consensus predicts female correctly")

pred_3m <- nipter_predict_sex(test_male, model_y, model_xy, model_yu,
                              y_unique_ratio = 0.004)
expect_identical(pred_3m$prediction, "male",
                 info = "3-model consensus predicts male correctly")

# y_unique model without ratio → warning + NA + skipped in consensus
expect_warning(
  pred_skip <- nipter_predict_sex(test_female, model_yu),
  info = "warns when y_unique model present but ratio missing"
)
expect_true(is.na(pred_skip$model_predictions[["y_unique"]]),
            info = "y_unique prediction is NA when ratio missing")
# With only the y_unique model skipped, consensus defaults to female (tie)
expect_identical(pred_skip$prediction, "female",
                 info = "consensus defaults to female when only model skipped")


# ===================================================================
# nipter_y_unique_ratio — input validation
# ===================================================================

# Bad bam path
expect_error(nipter_y_unique_ratio("nonexistent.bam"),
             info = "y_unique_ratio rejects nonexistent BAM")

# Bad regions file — use a temp file as a fake BAM so we pass the file.exists
# check on the BAM path and hit the regions_file check instead.
.tmp_fake_bam <- tempfile(fileext = ".bam")
file.create(.tmp_fake_bam)
expect_error(nipter_y_unique_ratio(.tmp_fake_bam,
                                   regions_file = "nonexistent_regions.txt"),
             info = "y_unique_ratio rejects nonexistent regions file")
unlink(.tmp_fake_bam)
