library(tinytest)
library(RWisecondorX)

# ---------------------------------------------------------------------------
# Tests for nipter_sex.R — sex prediction via mclust GMM
#
# These tests use synthetic NIPTeRSample objects with realistic sex
# chromosome signal (no BAM fixture required).
# ---------------------------------------------------------------------------

library(mclust)

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

  CombinedStrandsSample(
    sample_name = name,
    binsize = as.integer(50000L),
    chrom_lengths = setNames(
      rep.int(as.integer(50000L) * n_bins, 24L),
      c(as.character(1:22), "X", "Y")
    ),
    auto_matrix = auto_mat,
    sex_matrix_ = sex_mat
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

expect_true(is.list(model_y),
            info = "sex model stays list-like")
expect_true(S7::S7_inherits(model_y, NIPTeRSexModel),
            info = "sex model uses the typed NIPTeRSexModel S7 class")
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

expect_true(S7::S7_inherits(model_xy, NIPTeRSexModel),
            info = "XY model uses the typed NIPTeRSexModel S7 class")
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

expect_true(is.list(pred_f),
            info = "prediction stays list-like")
expect_true(S7::S7_inherits(pred_f, NIPTeRSexPrediction),
            info = "prediction uses the typed NIPTeRSexPrediction S7 class")
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

model_y_tie <- RWisecondorX:::.as_nipter_sex_model(
  modifyList(
    as.list(model_y),
    list(male_cluster = if (model_y$male_cluster == 1L) 2L else 1L)
  )
)
expect_warning(
  pred_tie <- nipter_predict_sex(test_female, model_y_tie, model_xy),
  pattern = "tied",
  info = "sex prediction warns on a tied model vote"
)
expect_identical(pred_tie$prediction, "female",
                 info = "tied sex-model vote resolves conservatively to female")

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

imbalanced_ncv_samples <- c(
  list(.make_sex_sample("few_female_1", "female", seed = 2001L)),
  lapply(seq_len(4L), function(i) .make_sex_sample(paste0("few_male_", i), "male",
                                                   seed = 2010L + i))
)
imbalanced_ncv_labels <- c(
  few_female_1 = "female",
  few_male_1 = "male",
  few_male_2 = "male",
  few_male_3 = "male",
  few_male_4 = "male"
)
imbalanced_ncv_cg <- nipter_as_control_group(
  imbalanced_ncv_samples,
  sample_sex = imbalanced_ncv_labels,
  sex_source = "synthetic_truth"
)
imbalanced_ncv_ref <- nipter_build_reference(
  imbalanced_ncv_cg,
  sample_sex = imbalanced_ncv_labels,
  sex_source = "synthetic_truth",
  sex_methods = c("y_fraction", "xy_fraction")
)
expect_error(
  nipter_build_sex_ncv_models(
    imbalanced_ncv_ref,
    candidate_chromosomes = c(1L, 2L, 4L, 5L, 6L, 7L),
    min_elements = 2L,
    max_elements = 3L
  ),
  pattern = "female non-outlier subset has only 1 sample",
  info = "NCV builder errors clearly when one sex subset is too small"
)

imbalanced_reg_samples <- c(
  lapply(seq_len(3L), function(i) .make_sex_sample(paste0("reg_female_", i), "female",
                                                   seed = 2100L + i)),
  lapply(seq_len(4L), function(i) .make_sex_sample(paste0("reg_male_", i), "male",
                                                   seed = 2200L + i))
)
imbalanced_reg_samples <- unlist(imbalanced_reg_samples, recursive = FALSE)
imbalanced_reg_labels <- c(
  reg_female_1 = "female",
  reg_female_2 = "female",
  reg_female_3 = "female",
  reg_male_1 = "male",
  reg_male_2 = "male",
  reg_male_3 = "male",
  reg_male_4 = "male"
)
imbalanced_reg_cg <- nipter_as_control_group(
  imbalanced_reg_samples,
  sample_sex = imbalanced_reg_labels,
  sex_source = "synthetic_truth"
)
imbalanced_reg_ref <- nipter_build_reference(
  imbalanced_reg_cg,
  sample_sex = imbalanced_reg_labels,
  sex_source = "synthetic_truth",
  sex_methods = c("y_fraction", "xy_fraction")
)
expect_error(
  nipter_build_sex_regression_models(
    imbalanced_reg_ref,
    candidate_chromosomes = c(1L, 2L, 4L, 5L, 6L, 7L),
    n_models = 2L,
    n_predictors = 2L,
    extra_predictors = character()
  ),
  pattern = "female non-outlier subset has only .*need >= 4",
  info = "regression builder errors clearly when one sex subset is too small"
)

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

expect_true(S7::S7_inherits(model_yu, NIPTeRSexModel),
            info = "y_unique model uses the typed NIPTeRSexModel S7 class")
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

# y_unique model without ratio -> hard error
expect_error(
  nipter_predict_sex(test_female, model_yu),
  pattern = "y_unique_ratio",
  info = "predict errors when a y_unique model is supplied without an explicit ratio"
)


# ===================================================================
# TYPED REFERENCE MODEL
# ===================================================================

sex_labels <- c(
  setNames(rep("female", 15L), paste0("female_", seq_len(15L))),
  setNames(rep("male", 15L), paste0("male_", seq_len(15L)))
)

ref_model <- nipter_build_reference(
  sex_cg,
  sample_sex = sex_labels,
  sex_source = "synthetic_truth",
  sex_methods = c("y_fraction", "xy_fraction"),
  y_unique_ratios = all_ratios,
  build_params = list(source = "test_nipter_sex")
)

expect_true(is.list(ref_model),
            info = "reference model stays list-like")
expect_true(S7::S7_inherits(ref_model, NIPTReferenceModel),
            info = "reference model uses the typed NIPTReferenceModel S7 class")
expect_true(S7::S7_inherits(ref_model$reference_frame, NIPTReferenceFrame),
            info = "reference model contains a typed reference frame")
expect_identical(sort(names(ref_model$sex_models)),
                 c("xy_fraction", "y_fraction", "y_unique"),
                 info = "reference model stores all requested sex models")
expect_identical(ref_model$sample_sex_source, "synthetic_truth",
                 info = "reference model keeps sample-sex provenance")
expect_identical(ref_model$reference_frame$SampleSex,
                 unname(sex_labels[ref_model$reference_frame$Sample_name]),
                 info = "reference model frame carries aligned sample sex labels")
expect_identical(ref_model$reference_frame$ConsensusGender,
                 unname(sex_labels[ref_model$reference_frame$Sample_name]),
                 info = "reference model frame resolves consensus sex labels")
expect_true(all(c("RR_X", "RR_Y", "RR_X_SexClassMAD", "RR_Y_SexClassMAD",
                  "Z_X_XX", "Z_X_XY", "Z_Y_XX", "Z_Y_XY",
                  "IsRefSexOutlier", "YUniqueRatio") %in%
                  names(ref_model$reference_frame)),
            info = "reference model frame includes sex-scoring metadata columns")
expect_equal(ref_model$reference_frame$RR_X,
             ref_model$reference_frame$FrChrReads_X,
             info = "RR_X matches X fraction against autosomal totals")
expect_equal(ref_model$reference_frame$RR_Y,
             ref_model$reference_frame$FrChrReads_Y,
             info = "RR_Y matches Y fraction against autosomal totals")
expect_equal(ref_model$reference_frame$YUniqueRatio,
             unname(all_ratios[ref_model$reference_frame$Sample_name]),
             info = "reference model frame aligns Y-unique ratios to samples")
expect_true(is.logical(ref_model$reference_frame$IsRefSexOutlier),
            info = "reference model frame stores logical sex-outlier flags")
expect_true(all(is.finite(ref_model$reference_frame$Z_X_XX[
  ref_model$reference_frame$ConsensusGender == "female"
])),
info = "reference model frame stores finite XX-reference X z-scores for female controls")
expect_true(all(is.finite(ref_model$reference_frame$Z_Y_XY[
  ref_model$reference_frame$ConsensusGender == "male"
])),
info = "reference model frame stores finite XY-reference Y z-scores for male controls")

female_frame <- as.data.frame(ref_model$reference_frame, stringsAsFactors = FALSE)
female_non_outlier <- female_frame[
  female_frame$ConsensusGender == "female" & !female_frame$IsRefSexOutlier,
  ,
  drop = FALSE
]
female_target <- female_non_outlier$Sample_name[[1L]]
female_others <- female_non_outlier[
  female_non_outlier$Sample_name != female_target,
  ,
  drop = FALSE
]
manual_female_z_x_xx <- (female_non_outlier$FrChrReads_X[[1L]] -
  mean(female_others$FrChrReads_X)) / stats::sd(female_others$FrChrReads_X)
expect_equal(
  female_frame$Z_X_XX[female_frame$Sample_name == female_target],
  manual_female_z_x_xx,
  tolerance = 1e-10,
  info = "reference-frame Z_X_XX uses leave-one-out scoring within the female cluster"
)

male_target <- female_frame$Sample_name[female_frame$ConsensusGender == "male"][[1L]]
manual_male_z_x_xx <- (
  female_frame$FrChrReads_X[female_frame$Sample_name == male_target] -
    mean(female_non_outlier$FrChrReads_X)
) / stats::sd(female_non_outlier$FrChrReads_X)
expect_equal(
  female_frame$Z_X_XX[female_frame$Sample_name == male_target],
  manual_male_z_x_xx,
  tolerance = 1e-10,
  info = "reference-frame Z_X_XX scores male controls against the female cluster without leave-one-out"
)

nonbinary_labels <- sex_labels
nonbinary_labels[["female_1"]] <- "ambiguous"
nonbinary_labels[["male_1"]] <- "unknown"

ref_model_nonbinary <- nipter_build_reference(
  sex_cg,
  sample_sex = nonbinary_labels,
  sex_source = "synthetic_truth",
  sex_methods = c("y_fraction", "xy_fraction"),
  y_unique_ratios = all_ratios
)

expect_identical(
  ref_model_nonbinary$reference_frame$ConsensusGender[
    match(c("female_1", "male_1"), ref_model_nonbinary$reference_frame$Sample_name)
  ],
  c("ambiguous", "unknown"),
  info = "reference frame preserves explicit ambiguous and unknown consensus labels"
)
expect_true(all(c("ConsensusGender", "IsRefSexOutlier") %in%
                  names(ref_model_nonbinary$reference_frame)),
            info = "reference frame still carries sex metadata when non-binary labels are present")

ref_model_nonbinary_ncv <- nipter_build_sex_ncv_models(
  ref_model_nonbinary,
  min_elements = 2L,
  max_elements = 2L,
  candidate_chromosomes = 1:12
)
expect_true("sex_ncv_models" %in% names(ref_model_nonbinary_ncv),
            info = "gaunosome NCV models still build from remaining binary non-outlier subsets")

pred_from_ref <- nipter_predict_sex(test_male, ref_model, y_unique_ratio = 0.004)
expect_identical(pred_from_ref$prediction, "male",
                 info = "reference-model dispatch predicts male sample correctly")


# ===================================================================
# TYPED SEX SCORE
# ===================================================================

score_f <- nipter_sex_score(test_female, ref_model, y_unique_ratio = 0.00003)
score_m <- nipter_sex_score(test_male, ref_model, y_unique_ratio = 0.004)

expect_true(is.list(score_f),
            info = "sex score stays list-like")
expect_true(S7::S7_inherits(score_f, NIPTSexScore),
            info = "sex score uses the typed NIPTSexScore S7 class")
expect_identical(score_f$predicted_sex, "female",
                 info = "female sample scores against female controls")
expect_identical(score_m$predicted_sex, "male",
                 info = "male sample scores against male controls")
expect_identical(score_f$sample_name, "test_female",
                 info = "sex score preserves sample name")
expect_identical(score_m$sample_name, "test_male",
                 info = "sex score preserves sample name for males")
expect_identical(
  unname(score_f$z_scores[c("Z_FrChrReads_X", "Z_FrChrReads_Y")]),
  unname(score_f$z_scores[c("Z_FrChrReads_X_XX", "Z_FrChrReads_Y_XX")]),
  info = "female selected z-scores reuse the XX reference hypothesis"
)
expect_identical(
  unname(score_m$z_scores[c("Z_FrChrReads_X", "Z_FrChrReads_Y")]),
  unname(score_m$z_scores[c("Z_FrChrReads_X_XY", "Z_FrChrReads_Y_XY")]),
  info = "male selected z-scores reuse the XY reference hypothesis"
)

female_ref <- ref_model$reference_frame[
  ref_model$reference_frame$ConsensusGender == "female" &
    !ref_model$reference_frame$IsRefSexOutlier,
  ,
  drop = FALSE
]
male_ref <- ref_model$reference_frame[
  ref_model$reference_frame$ConsensusGender == "male" &
    !ref_model$reference_frame$IsRefSexOutlier,
  ,
  drop = FALSE
]

manual_female_y <- (score_f$sample_metrics[["FrChrReads_Y"]] -
                      mean(female_ref$FrChrReads_Y)) /
  stats::sd(female_ref$FrChrReads_Y)
manual_male_y <- (score_m$sample_metrics[["FrChrReads_Y"]] -
                    mean(male_ref$FrChrReads_Y)) /
  stats::sd(male_ref$FrChrReads_Y)
manual_female_cv_y <- 100 * stats::sd(female_ref$FrChrReads_Y) /
  mean(female_ref$FrChrReads_Y)
manual_male_cv_y <- 100 * stats::sd(male_ref$FrChrReads_Y) /
  mean(male_ref$FrChrReads_Y)

expect_equal(score_f$z_scores[["Z_FrChrReads_Y"]], manual_female_y,
             info = "female selected Y z-score matches manual reference calculation")
expect_equal(score_m$z_scores[["Z_FrChrReads_Y"]], manual_male_y,
             info = "male selected Y z-score matches manual reference calculation")
expect_equal(score_f$cv[["Z_FrChrReads_CV_Y"]], manual_female_cv_y,
             info = "female selected Y CV matches manual reference calculation")
expect_equal(score_m$cv[["Z_FrChrReads_CV_Y"]], manual_male_cv_y,
             info = "male selected Y CV matches manual reference calculation")
expect_identical(score_f$reference_sizes[["same_sex"]], nrow(female_ref),
                 info = "female score reports the filtered same-sex reference size")
expect_identical(score_m$reference_sizes[["same_sex"]], nrow(male_ref),
                 info = "male score reports the filtered same-sex reference size")
expect_true(all(score_f$reference_sample_names %in% female_ref$Sample_name),
            info = "female score uses only non-outlier female controls")
expect_true(all(score_m$reference_sample_names %in% male_ref$Sample_name),
            info = "male score uses only non-outlier male controls")


# ===================================================================
# TYPED SEX NCV / REGRESSION MODELS
# ===================================================================

small_candidate_pool <- c(1L, 2L, 4L, 5L, 6L, 7L)
ref_model_models <- nipter_build_sex_ncv_models(
  ref_model,
  candidate_chromosomes = small_candidate_pool,
  min_elements = 2L,
  max_elements = 3L
)
ref_model_models <- nipter_build_sex_regression_models(
  ref_model_models,
  candidate_chromosomes = small_candidate_pool,
  n_models = 2L,
  n_predictors = 2L,
  extra_predictors = character()
)

expect_true(S7::S7_inherits(ref_model_models, NIPTReferenceModel),
            info = "gaunosome builders keep the typed reference model class")
expect_true(S7::S7_inherits(ref_model_models$sex_ncv_models$female$X, NIPTSexNCVModel),
            info = "female X NCV model is typed")
expect_true(S7::S7_inherits(ref_model_models$sex_ncv_models$male$Y, NIPTSexNCVModel),
            info = "male Y NCV model is typed")
expect_true(all(vapply(ref_model_models$sex_regression_models$female$X,
                       function(x) S7::S7_inherits(x, NIPTSexRegressionModel),
                       logical(1L))),
            info = "female X regression models are typed")
expect_true(all(vapply(ref_model_models$sex_regression_models$male$Y,
                       function(x) S7::S7_inherits(x, NIPTSexRegressionModel),
                       logical(1L))),
            info = "male Y regression models are typed")

sample_rr <- function(sample) {
  auto_counts <- rowSums(sample$autosomal_chromosome_reads[[1L]])
  sex_counts <- rowSums(sample$sex_chromosome_reads[[1L]])
  auto_total <- sum(auto_counts)
  vals <- c(auto_counts, sex_counts)
  names(vals) <- c(paste0("RR_", as.character(1:22)), "RR_X", "RR_Y")
  vals / auto_total
}

ncv_x_f <- nipter_ncv_sex_score(test_female, ref_model_models,
                                focus_chromosome = "X",
                                y_unique_ratio = 0.00003)
ncv_y_m <- nipter_ncv_sex_score(test_male, ref_model_models,
                                focus_chromosome = "Y",
                                y_unique_ratio = 0.004)

expect_true(S7::S7_inherits(ncv_x_f, NIPTSexNCVScore),
            info = "female X NCV score is typed")
expect_true(S7::S7_inherits(ncv_y_m, NIPTSexNCVScore),
            info = "male Y NCV score is typed")
expect_identical(ncv_x_f$predicted_sex, "female",
                 info = "female X NCV score uses female reference set")
expect_identical(ncv_y_m$predicted_sex, "male",
                 info = "male Y NCV score uses male reference set")
expect_identical(ncv_x_f$sample_scores[["selected"]],
                 ncv_x_f$sample_scores[["female"]],
                 info = "selected NCV score reuses the female score for female samples")
expect_identical(ncv_y_m$sample_scores[["selected"]],
                 ncv_y_m$sample_scores[["male"]],
                 info = "selected NCV score reuses the male score for male samples")

female_x_ncv_model <- ref_model_models$sex_ncv_models$female$X
male_y_ncv_model <- ref_model_models$sex_ncv_models$male$Y
test_female_row <- c(rowSums(test_female$autosomal_chromosome_reads[[1L]]),
                     rowSums(test_female$sex_chromosome_reads[[1L]]))
names(test_female_row) <- c(paste0("NChrReads_", as.character(1:22)),
                            "NChrReads_X", "NChrReads_Y")
test_male_row <- c(rowSums(test_male$autosomal_chromosome_reads[[1L]]),
                   rowSums(test_male$sex_chromosome_reads[[1L]]))
names(test_male_row) <- c(paste0("NChrReads_", as.character(1:22)),
                          "NChrReads_X", "NChrReads_Y")

manual_ncv_f_x <- (test_female_row[["NChrReads_X"]] /
                     sum(test_female_row[female_x_ncv_model$denominators]) -
                     female_x_ncv_model$control_statistics[["mean"]]) /
  female_x_ncv_model$control_statistics[["sd"]]
manual_ncv_m_y <- (test_male_row[["NChrReads_Y"]] /
                     sum(test_male_row[male_y_ncv_model$denominators]) -
                     male_y_ncv_model$control_statistics[["mean"]]) /
  male_y_ncv_model$control_statistics[["sd"]]

expect_equal(ncv_x_f$sample_scores[["female"]], manual_ncv_f_x,
             info = "female X NCV score matches manual calculation")
expect_equal(ncv_y_m$sample_scores[["male"]], manual_ncv_m_y,
             info = "male Y NCV score matches manual calculation")

reg_x_f <- nipter_regression_sex_score(test_female, ref_model_models,
                                       focus_chromosome = "X",
                                       y_unique_ratio = 0.00003)
reg_y_m <- nipter_regression_sex_score(test_male, ref_model_models,
                                       focus_chromosome = "Y",
                                       y_unique_ratio = 0.004)

expect_true(S7::S7_inherits(reg_x_f, NIPTSexRegressionScore),
            info = "female X regression score is typed")
expect_true(S7::S7_inherits(reg_y_m, NIPTSexRegressionScore),
            info = "male Y regression score is typed")
expect_identical(reg_x_f$predicted_sex, "female",
                 info = "female regression score uses female models")
expect_identical(reg_y_m$predicted_sex, "male",
                 info = "male regression score uses male models")
expect_equal(reg_x_f$aggregate_scores[["selected_mean"]],
             reg_x_f$aggregate_scores[["female_mean"]],
             info = "female selected regression mean matches female model mean")
expect_equal(reg_y_m$aggregate_scores[["selected_mean"]],
             reg_y_m$aggregate_scores[["male_mean"]],
             info = "male selected regression mean matches male model mean")

female_x_reg_model <- ref_model_models$sex_regression_models$female$X[[1L]]
male_y_reg_model <- ref_model_models$sex_regression_models$male$Y[[1L]]
test_female_rr <- sample_rr(test_female)
test_male_rr <- sample_rr(test_male)

female_x_newdata <- as.data.frame(as.list(test_female_rr[c(female_x_reg_model$predictors,
                                                           female_x_reg_model$response_column)]))
male_y_newdata <- as.data.frame(as.list(test_male_rr[c(male_y_reg_model$predictors,
                                                       male_y_reg_model$response_column)]))
manual_reg_f_x_ratio <- test_female_rr[[female_x_reg_model$response_column]] /
  as.numeric(stats::predict(female_x_reg_model$fit, newdata = female_x_newdata))
manual_reg_m_y_ratio <- test_male_rr[[male_y_reg_model$response_column]] /
  as.numeric(stats::predict(male_y_reg_model$fit, newdata = male_y_newdata))
manual_reg_f_x <- (manual_reg_f_x_ratio -
                     female_x_reg_model$control_statistics[["mean_ratio"]]) /
  female_x_reg_model$control_statistics[["sd_ratio"]]
manual_reg_m_y <- (manual_reg_m_y_ratio -
                     male_y_reg_model$control_statistics[["mean_ratio"]]) /
  male_y_reg_model$control_statistics[["sd_ratio"]]

expect_equal(reg_x_f$scores$female[[1L]], manual_reg_f_x,
             info = "female X regression score matches manual first-model calculation")
expect_equal(reg_y_m$scores$male[[1L]], manual_reg_m_y,
             info = "male Y regression score matches manual first-model calculation")


# ===================================================================
# AGGREGATE GAUNOSOME API
# ===================================================================

ref_model_bundle <- nipter_build_gaunosome_models(
  ref_model,
  candidate_chromosomes = small_candidate_pool,
  ncv_min_elements = 2L,
  ncv_max_elements = 3L,
  regression_n_models = 2L,
  regression_n_predictors = 2L,
  regression_extra_predictors = character()
)

expect_true(S7::S7_inherits(ref_model_bundle, NIPTReferenceModel),
            info = "aggregate gaunosome builder keeps the typed reference class")
expect_true("sex_ncv_models" %in% names(ref_model_bundle),
            info = "aggregate gaunosome builder attaches NCV models")
expect_true("sex_regression_models" %in% names(ref_model_bundle),
            info = "aggregate gaunosome builder attaches regression models")

gauno_f <- nipter_gaunosome_score(test_female, ref_model_bundle,
                                  y_unique_ratio = 0.00003)
gauno_m <- nipter_gaunosome_score(test_male, ref_model_bundle,
                                  y_unique_ratio = 0.004)

expect_true(S7::S7_inherits(gauno_f, NIPTGaunosomeScore),
            info = "aggregate female gaunosome score is typed")
expect_true(S7::S7_inherits(gauno_m, NIPTGaunosomeScore),
            info = "aggregate male gaunosome score is typed")
expect_identical(gauno_f$predicted_sex, "female",
                 info = "aggregate female gaunosome score preserves predicted sex")
expect_identical(gauno_m$predicted_sex, "male",
                 info = "aggregate male gaunosome score preserves predicted sex")
expect_identical(sort(names(gauno_f$ncv_scores)), c("X", "Y"),
                 info = "aggregate gaunosome score includes X/Y NCV components")
expect_identical(sort(names(gauno_f$regression_scores)), c("X", "Y"),
                 info = "aggregate gaunosome score includes X/Y regression components")
expect_identical(gauno_f$summary$chromosome, c("X", "Y"),
                 info = "aggregate summary reports one row per requested gaunosome")
expect_identical(gauno_m$summary$chromosome, c("X", "Y"),
                 info = "aggregate male summary reports one row per requested gaunosome")

gauno_fx <- gauno_f$summary[gauno_f$summary$chromosome == "X", , drop = FALSE]
gauno_my <- gauno_m$summary[gauno_m$summary$chromosome == "Y", , drop = FALSE]

expect_equal(gauno_fx$z_score, gauno_f$sex_score$z_scores[["Z_FrChrReads_X"]],
             info = "aggregate summary X z-score matches the selected sex score")
expect_equal(gauno_my$z_score, gauno_m$sex_score$z_scores[["Z_FrChrReads_Y"]],
             info = "aggregate summary Y z-score matches the selected sex score")
expect_equal(gauno_fx$ncv_score_selected,
             gauno_f$ncv_scores$X$sample_scores[["selected"]],
             info = "aggregate summary X NCV score matches the component score")
expect_equal(gauno_my$ncv_score_selected,
             gauno_m$ncv_scores$Y$sample_scores[["selected"]],
             info = "aggregate summary Y NCV score matches the component score")
expect_equal(gauno_fx$regression_score_selected,
             gauno_f$regression_scores$X$aggregate_scores[["selected_mean"]],
             info = "aggregate summary X regression score matches the component score")
expect_equal(gauno_my$regression_score_selected,
             gauno_m$regression_scores$Y$aggregate_scores[["selected_mean"]],
             info = "aggregate summary Y regression score matches the component score")


# ===================================================================
# COHORT GAUNOSOME REPORT / OUTPUT
# ===================================================================

gauno_report <- nipter_gaunosome_report(
  list(test_female, test_male),
  ref_model_bundle,
  y_unique_ratios = c(test_female = 0.00003, test_male = 0.004)
)

expect_true(S7::S7_inherits(gauno_report, NIPTGaunosomeReport),
            info = "batch gaunosome report is typed")
expect_identical(sort(gauno_report$sample_names), c("test_female", "test_male"),
                 info = "batch gaunosome report records all sample names")
expect_identical(sort(names(gauno_report$scores)), c("test_female", "test_male"),
                 info = "batch gaunosome report names score entries by sample")
expect_identical(nrow(gauno_report$summary), 4L,
                 info = "batch gaunosome summary has one row per sample/chromosome")
expect_true(all(c("sample_name", "chromosome", "predicted_sex") %in%
                  names(gauno_report$summary)),
            info = "batch gaunosome summary contains the flattened reporting columns")

report_fx <- gauno_report$summary[
  gauno_report$summary$sample_name == "test_female" &
    gauno_report$summary$chromosome == "X",
  ,
  drop = FALSE
]
report_my <- gauno_report$summary[
  gauno_report$summary$sample_name == "test_male" &
    gauno_report$summary$chromosome == "Y",
  ,
  drop = FALSE
]

expect_equal(report_fx$ncv_score_selected,
             gauno_report$scores$test_female$ncv_scores$X$sample_scores[["selected"]],
             info = "batch report X NCV score matches the stored female component")
expect_equal(report_my$regression_score_selected,
             gauno_report$scores$test_male$regression_scores$Y$aggregate_scores[["selected_mean"]],
             info = "batch report Y regression score matches the stored male component")

gauno_outprefix <- tempfile("gaunosome_report_")
gauno_out <- write_nipter_gaunosome_output(gauno_report, gauno_outprefix)
expect_true(file.exists(gauno_out),
            info = "gaunosome output writer creates the TSV summary file")

gauno_written <- utils::read.delim(gauno_out, stringsAsFactors = FALSE)
expect_identical(nrow(gauno_written), nrow(gauno_report$summary),
                 info = "written gaunosome summary has the same row count as the report")
expect_identical(names(gauno_written), names(gauno_report$summary),
                 info = "written gaunosome summary preserves the report columns")
unlink(gauno_out)


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


# ===================================================================
# nipter_y_unique_ratio — BED coverage backend regression
# ===================================================================

.fixture_mixed_bam <- system.file("extdata", "fixture_mixed.bam", package = "Rduckhts")
.fixture_mixed_regions <- system.file("extdata", "fixture_mixed_regions.bed", package = "Rduckhts")

expect_true(file.exists(.fixture_mixed_bam),
            info = "Rduckhts mixed BAM fixture is available")
expect_true(file.exists(.fixture_mixed_regions),
            info = "Rduckhts mixed BED fixture is available")

fixture_yu <- nipter_y_unique_ratio(
  .fixture_mixed_bam,
  mapq = 0L,
  exclude_flags = 1796L,
  regions_file = .fixture_mixed_regions
)

expect_identical(fixture_yu$y_unique_reads, 4L,
                 info = "Y-unique backend sums post-filter BED coverage counts across intervals")
expect_identical(fixture_yu$total_nuclear_reads, 8L,
                 info = "Y-unique backend uses post-filter nuclear reads as the denominator")
expect_equal(fixture_yu$ratio, 0.5,
             info = "Y-unique ratio equals covered reads divided by total nuclear reads")
expect_identical(nrow(fixture_yu$regions), 3L,
                 info = "Y-unique ratio returns the intervals actually used")
