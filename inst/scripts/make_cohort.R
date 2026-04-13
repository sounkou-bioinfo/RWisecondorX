#!/usr/bin/env Rscript
# inst/scripts/make_cohort.R
#
# CLI wrapper for generate_cohort(). Generates 50 synthetic BAMs into a
# specified directory (defaults to tempdir()).
#
# Usage:
#   Rscript inst/scripts/make_cohort.R [output_dir]
#
# The real implementation lives in R/cohort.R so that it is available as
# RWisecondorX::generate_cohort() when the package is installed.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[1] else tempdir()

if (requireNamespace("RWisecondorX", quietly = TRUE)) {
  RWisecondorX::generate_cohort(out_dir)
} else {
  # Fall back to sourcing directly (e.g. from a dev tree)
  source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "R", "cohort.R"))
  generate_cohort(out_dir)
}
