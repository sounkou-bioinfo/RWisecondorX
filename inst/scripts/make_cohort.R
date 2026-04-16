#!/usr/bin/env Rscript
# inst/scripts/make_cohort.R
#
# CLI wrapper for generate_cohort(). Generates 50 synthetic BAMs into a
# specified directory (defaults to tempdir()).
#
# Usage:
#   Rscript inst/scripts/make_cohort.R [output_dir]
#
# The real implementation lives in R/synthetic_cohort.R so that it is available as
# RWisecondorX::generate_cohort() when the package is installed.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[1] else tempdir()

if (requireNamespace("RWisecondorX", quietly = TRUE)) {
  RWisecondorX::generate_cohort(out_dir)
} else {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg) == 1L) sub("^--file=", "", script_arg) else "."
  source(file.path(dirname(script_path), "..", "..", "R", "synthetic_cohort.R"))
  generate_cohort(out_dir)
}
