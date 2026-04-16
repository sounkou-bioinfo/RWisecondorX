#!/usr/bin/env Rscript
# inst/scripts/precompute_gc.R
#
# CLI script for precomputing per-bin GC content from a reference FASTA.
#
# Produces a bgzipped, tabix-indexed TSV file that can be passed to
# nipter_gc_correct(gc_table = ...) or to convert_sample.R / build_reference.R
# via --gc-table.
#
# Usage:
#   Rscript precompute_gc.R --fasta hg38.fa --out hg38_gc_50k.tsv.bgz
#   Rscript precompute_gc.R --fasta hg38.fa --out hg38_gc_5k.tsv.bgz --binsize 5000

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with: install.packages('optparse')",
       call. = FALSE)
}

library(optparse)

option_list <- list(
  make_option("--fasta", type = "character", default = NULL,
              help = "Indexed reference FASTA file (.fa/.fasta with .fai) [required]"),
  make_option("--out", type = "character", default = NULL,
              help = "Output path for the GC table (.tsv.bgz) [required]"),
  make_option("--binsize", type = "integer", default = 50000L,
              help = "Bin size in base pairs [default: %default]")
)

parser <- OptionParser(
  usage = "usage: %prog [options]",
  description = paste(
    "Precompute per-bin GC content from a reference FASTA.",
    "",
    "Output is a 5-column TSV.bgz (chrom, start, end, pct_gc, seq_len) with",
    "a tabix index (.tbi). Pass the output path to nipter_gc_correct(gc_table = ...)",
    "or to convert_sample.R / build_reference.R via --gc-table.",
    "",
    "Examples:",
    "  # Default 50kb bins (NIPTeR)",
    "  %prog --fasta hg38.fa --out hg38_gc_50k.tsv.bgz",
    "",
    "  # Custom bin size",
    "  %prog --fasta hg38.fa --out hg38_gc_5k.tsv.bgz --binsize 5000",
    sep = "\n"
  ),
  option_list = option_list
)

opts <- parse_args(parser)

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

if (is.null(opts$fasta)) {
  stop("--fasta is required (path to indexed reference FASTA file)", call. = FALSE)
}

if (!file.exists(opts$fasta)) {
  stop("FASTA file not found: ", opts$fasta, call. = FALSE)
}

if (is.null(opts$out)) {
  stop("--out is required (output path for GC table, e.g. hg38_gc_50k.tsv.bgz)", call. = FALSE)
}

if (opts$binsize < 1L) {
  stop("--binsize must be a positive integer, got: ", opts$binsize, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

library(RWisecondorX)

cat(sprintf("Precomputing GC table from %s (binsize=%d) ...\n",
            basename(opts$fasta), opts$binsize))

nipter_gc_precompute(
  fasta   = opts$fasta,
  binsize = opts$binsize,
  out     = opts$out
)

cat("Done.\n")
