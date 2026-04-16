library(tinytest)
library(RWisecondorX)

.helper_nipter <- system.file("tinytest", "helper_nipter.R", package = "RWisecondorX")
if (nzchar(.helper_nipter)) {
  sys.source(.helper_nipter, envir = environment())
} else {
  sys.source("inst/tinytest/helper_nipter.R", envir = environment())
}

ctrl_samples <- .sim_nipter_control_set(4L, seed = 11L)
sex_labels <- c(
  ctrl_01 = "female",
  ctrl_02 = "female",
  ctrl_03 = "male",
  ctrl_04 = "male"
)

cg <- nipter_as_control_group(
  ctrl_samples,
  description = "Sex-annotated controls",
  sample_sex = sex_labels,
  sex_source = "explicit"
)

expect_identical(cg$sample_sex, sex_labels,
                 info = "control group stores explicit sex labels")
expect_identical(cg$sex_source, "explicit",
                 info = "control group stores sex label provenance")

ref <- nipter_reference_frame(cg)

expect_true(is.data.frame(ref),
            info = "reference frame returns a data.frame")
expect_true(S7::S7_inherits(ref, NIPTReferenceFrame),
            info = "reference frame returns the typed NIPTReferenceFrame S7 class")
expect_identical(nrow(ref), 4L,
                 info = "reference frame has one row per control")
expect_true(all(c("Sample_name", "SampleSex", "NChrReads_X", "NChrReads_Y",
                  "FrChrReads_X", "FrChrReads_Y") %in% names(ref)),
            info = "reference frame exposes sex-aware count and fraction columns")
expect_identical(ref$SampleSex, unname(sex_labels[ref$Sample_name]),
                 info = "reference frame carries aligned sample sex labels")
expect_true(all(vapply(ref[paste0("NChrReads_", c(1:22, "X", "Y"))],
                       is.numeric, logical(1L))),
            info = "reference frame count columns are numeric")
expect_true(all(vapply(ref[paste0("FrChrReads_", c(1:22, "X", "Y"))],
                       is.numeric, logical(1L))),
            info = "reference frame fraction columns are numeric")

subset_cg <- nipter_match_control_group(
  ctrl_samples[[1L]],
  cg,
  n = 2L,
  mode = "subset"
)
expect_identical(length(subset_cg$sample_sex), 2L,
                 info = "matched subset preserves sex labels")
expect_true(all(names(subset_cg$sample_sex) %in%
                  vapply(subset_cg$samples, function(s) s$sample_name, character(1L))),
            info = "matched subset sex labels stay aligned to kept samples")
