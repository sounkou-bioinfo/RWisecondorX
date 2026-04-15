#!/usr/bin/env Rscript
# make_nipter_fixture.R — Generate multi-chromosome NIPTeR conformance fixture BAM
#
# Creates inst/extdata/nipter_conformance_fixture.bam — a synthetic whole-genome
# BAM suitable for cross-checking nipter_bin_bam() against NIPTeR::bin_bam_sample().
#
# Design constraints (from AGENTS.md / test_nipter.R):
#   1. No unmapped reads (avoids NIPTeR's split-length bug)
#   2. No two reads sharing the same start position on the same chromosome
#      (avoids NIPTeR's implicit unique() deduplication)
#   3. All 24 GRCh37 chromosomes present
#   4. Small file size suitable for inst/extdata (<2MB)
#
# Read placement: 1 read per 100kb bin at the bin start position + a small
# per-bin offset to avoid exact multiples (which could alias across binsize
# choices). This gives ~26,700 reads total across all chromosomes.
#
# Usage:
#   Rscript inst/scripts/make_nipter_fixture.R
#   # or via Makefile:
#   make nipter-fixture
#
# Output: inst/extdata/nipter_conformance_fixture.bam (+ .bai index)
# After building, set NIPTER_CONFORMANCE_BAM to this path to run Arm B tests.

stopifnot(nchar(Sys.which("samtools")) > 0L)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
pkg_root <- if (length(script_arg) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."))
} else {
  normalizePath(getwd())
}
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  stop("Could not locate package root from script path or working directory.", call. = FALSE)
}
out_bam  <- file.path(pkg_root, "inst", "extdata",
                       "nipter_conformance_fixture.bam")
out_sam  <- tempfile(fileext = ".sam")

# ---- GRCh37 chromosome sizes (bp) from BINS_PER_CHR * 100000 ---------------
# (Matches the BINS_PER_CHR internal constant in R/cohort.R)
BINS_PER_CHR <- c(
  "1"  = 2493L, "2"  = 2432L, "3"  = 1981L, "4"  = 1912L,
  "5"  = 1810L, "6"  = 1712L, "7"  = 1592L, "8"  = 1464L,
  "9"  = 1413L, "10" = 1356L, "11" = 1351L, "12" = 1339L,
  "13" = 1152L, "14" = 1074L, "15" = 1026L, "16" =  904L,
  "17" =  812L, "18" =  781L, "19" =  592L, "20" =  631L,
  "21" =  482L, "22" =  514L, "X"  = 1553L, "Y"  =  594L
)
CHR_LENGTHS <- BINS_PER_CHR * 100000L  # real GRCh37 lengths
CHR_NAMES   <- names(CHR_LENGTHS)

READ_LEN   <- 36L
READ_SEQ   <- paste(rep("A", READ_LEN), collapse = "")
READ_QUAL  <- paste(rep("F", READ_LEN), collapse = "")
READ_CIGAR <- paste0(READ_LEN, "M")
BIN_SIZE   <- 100000L  # one read per 100kb bin

set.seed(123L)

# ---- Build SAM ---------------------------------------------------------------
message("Building NIPTeR conformance fixture SAM ...")

hd  <- "@HD\tVN:1.6\tSO:coordinate"
sq  <- sprintf("@SQ\tSN:%s\tLN:%d", CHR_NAMES, CHR_LENGTHS)
rg  <- "@RG\tID:fixture\tSM:nipter_conformance\tPL:ILLUMINA"
header_lines <- c(hd, sq, rg)

records <- character(0L)
read_idx <- 0L

for (chr in CHR_NAMES) {
  chr_len <- CHR_LENGTHS[chr]
  n_bins  <- BINS_PER_CHR[chr]

  # One read per bin; position = bin_start + small deterministic offset
  # Offset ensures no read sits exactly on a bin boundary of a different binsize
  bin_starts <- (seq_len(n_bins) - 1L) * BIN_SIZE
  # Small offset: 5bp into the bin to avoid 0-based ambiguity
  positions <- bin_starts + 5L
  # Clip to valid range (must be < chr_len - READ_LEN + 1 for 1-based SAM)
  max_pos <- chr_len - READ_LEN
  positions <- pmin(positions, max_pos)
  positions <- positions + 1L  # convert to 1-based SAM coordinates
  flags <- ifelse(seq_len(n_bins) %% 2L == 0L, 16L, 0L)

  read_idx_start <- read_idx + 1L
  recs <- sprintf("r%d\t%d\t%s\t%d\t60\t%s\t*\t0\t0\t%s\t%s\tRG:Z:fixture",
                  read_idx_start + seq_len(n_bins) - 1L,
                  flags, chr, positions,
                  READ_CIGAR, READ_SEQ, READ_QUAL)
  records <- c(records, recs)
  read_idx <- read_idx + n_bins
}

message(sprintf("Writing %d reads across %d chromosomes ...", read_idx, length(CHR_NAMES)))
writeLines(c(header_lines, records), out_sam)

# ---- Convert to sorted, indexed BAM -----------------------------------------
message("Converting to BAM with samtools ...")
tmp_unsorted <- tempfile(fileext = ".bam")

ret_view <- system2("samtools", c("view", "-bS", out_sam, "-o", tmp_unsorted))
if (ret_view != 0L) stop("samtools view failed")

ret_sort <- system2("samtools", c("sort", tmp_unsorted, "-o", out_bam))
if (ret_sort != 0L) stop("samtools sort failed")

ret_idx  <- system2("samtools", c("index", out_bam))
if (ret_idx != 0L) stop("samtools index failed")

unlink(c(out_sam, tmp_unsorted))

# ---- Verify ------------------------------------------------------------------
stats_out <- system2("samtools", c("flagstat", out_bam), stdout = TRUE)
message("BAM flagstat:\n", paste(stats_out, collapse = "\n"))

message("\nNIPTeR conformance fixture written to:\n  ", out_bam)
message("\nInstalled-package tests now use this bundled fixture by default.")
message("Set NIPTER_CONFORMANCE_BAM only to override it with a custom BAM.")
