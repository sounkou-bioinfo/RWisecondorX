# test_nipter_matching.R — Tests for Rcpp SSD kernels, nipter_match_matrix(),
# nipter_control_group_from_beds(), nipter_gc_precompute(), and the gc_table
# parameter on nipter_gc_correct().
#
# All tests that need only synthetic data run unconditionally.
# Tests that need BAM/FASTA fixtures skip cleanly when fixtures are absent.
#
# Simulation helpers come from helper_nipter.R (auto-sourced by tinytest).

library(tinytest)
library(RWisecondorX)

.helper_nipter <- system.file("tinytest", "helper_nipter.R", package = "RWisecondorX")
if (nzchar(.helper_nipter)) {
  sys.source(.helper_nipter, envir = environment())
} else {
  sys.source("inst/tinytest/helper_nipter.R", envir = environment())
}

.is_nipter_sample <- function(x) {
  RWisecondorX:::.is_nipt_sample_object(x)
}

.is_nipter_control_group <- function(x) {
  RWisecondorX:::.is_nipt_control_group_object(x)
}

.is_combined_sample <- function(x) {
  inherits(x, "CombinedStrands") || RWisecondorX:::.is_s7_combined_sample(x)
}

.is_separated_sample <- function(x) {
  inherits(x, "SeparatedStrands") || RWisecondorX:::.is_s7_separated_sample(x)
}

.sample_name_of <- function(x) x$sample_name


# ---------------------------------------------------------------------------
# 1. nipter_ssd_scores_cpp — known inputs
# ---------------------------------------------------------------------------

# Simple 3-chromosome, 3-control case (all chromosomes compared)
# fracs matrix (3 × 3): each col is one control's fractions
fracs3 <- matrix(c(
  0.4, 0.3, 0.3,   # control 1
  0.5, 0.3, 0.2,   # control 2
  0.4, 0.4, 0.2    # control 3
), nrow = 3L, ncol = 3L)

query3 <- c(0.4, 0.3, 0.3)           # identical to control 1
cmp_idx <- as.integer(c(0L, 1L, 2L)) # all 3 rows (0-based)

scores <- RWisecondorX:::nipter_ssd_scores_cpp(fracs3, query3, cmp_idx, 1L)

expect_equal(length(scores), 3L,
             info = "ssd_scores_cpp returns one score per control")
expect_equal(scores[1L], 0.0,
             info = "SSD of query against identical control is 0")
expect_true(scores[2L] > 0.0,
            info = "SSD of query against different control is positive")
expect_true(scores[3L] > 0.0,
            info = "SSD of query against different control is positive")

# Manual check for control 2: (0.4-0.5)^2 + (0.3-0.3)^2 + (0.3-0.2)^2 = 0.02
expect_equal(scores[2L], 0.02, tolerance = 1e-10,
             info = "ssd_scores_cpp value matches manual calculation")

# Subset comparison: only rows 0 and 2 (exclude row 1)
scores_sub <- RWisecondorX:::nipter_ssd_scores_cpp(fracs3, query3,
                                                    as.integer(c(0L, 2L)), 1L)
expect_equal(scores_sub[2L], (0.4 - 0.5)^2 + (0.3 - 0.2)^2, tolerance = 1e-10,
             info = "ssd_scores_cpp subset excludes middle row correctly")


# ---------------------------------------------------------------------------
# 2. nipter_ssd_matrix_cpp — symmetry, diagonal, known values
# ---------------------------------------------------------------------------

ssd_mat <- RWisecondorX:::nipter_ssd_matrix_cpp(fracs3, cmp_idx, 1L)

expect_equal(dim(ssd_mat), c(3L, 3L),
             info = "ssd_matrix_cpp returns N×N matrix")
expect_equal(diag(ssd_mat), rep(0.0, 3L),
             info = "ssd_matrix_cpp diagonal is zero")
expect_equal(ssd_mat, t(ssd_mat),
             info = "ssd_matrix_cpp is symmetric")
expect_equal(ssd_mat[1L, 2L], ssd_mat[2L, 1L],
             info = "ssd_matrix_cpp[1,2] == ssd_matrix_cpp[2,1]")
# Known value: ctrl1 vs ctrl2 = (0.4-0.5)^2+(0.3-0.3)^2+(0.3-0.2)^2 = 0.02
expect_equal(ssd_mat[1L, 2L], 0.02, tolerance = 1e-10,
             info = "ssd_matrix_cpp off-diagonal matches manual calculation")

# Row means of the SSD matrix equal the score produced by match_control_group
# with mode="report" for each sample used as its own query (leave-self-in).
# Just verify the upper bound: row means are non-negative.
expect_true(all(rowMeans(ssd_mat) >= 0.0),
            info = "ssd_matrix_cpp row means are non-negative")


# ---------------------------------------------------------------------------
# 3. nipter_match_matrix() — structure and consistency with
#    nipter_match_control_group(mode = "report")
# ---------------------------------------------------------------------------

cg5 <- nipter_as_control_group(.sim_nipter_control_set(5L, n_bins = 50L))

mm <- nipter_match_matrix(cg5)

expect_equal(dim(mm), c(5L, 5L),
             info = "nipter_match_matrix returns 5×5 for 5-sample group")
expect_equal(rownames(mm), colnames(mm),
             info = "nipter_match_matrix rownames == colnames")
sample_names_cg5 <- vapply(cg5$samples, .sample_name_of, character(1L))
expect_equal(rownames(mm), sample_names_cg5,
             info = "nipter_match_matrix rownames match control sample names")
expect_equal(unname(diag(mm)), rep(0.0, 5L),
             info = "nipter_match_matrix diagonal is zero")
expect_equal(mm, t(mm),
             info = "nipter_match_matrix is symmetric")

# Compare nipter_match_matrix row 1 against nipter_match_control_group(mode="report")
# for the same query (sample 1 vs all 5 controls — note: includes itself, SSD=0).
report_1 <- nipter_match_control_group(cg5$samples[[1L]], cg5,
                                       n = 5L, mode = "report")
# report_1 is sorted ascending; mm[1,] is in order → extract by name
mm_row1 <- mm[1L, ]
expect_equal(sort(as.numeric(mm_row1)), sort(as.numeric(report_1)),
             tolerance = 1e-12,
             info = "nipter_match_matrix row agrees with match_control_group report")


# ---------------------------------------------------------------------------
# 4. nipter_match_control_group() — exclude_chromosomes respected in Rcpp path
# ---------------------------------------------------------------------------

# Build two samples that differ only on chromosome 13, 18, 21 (the excluded ones)
totals_a <- stats::setNames(rep(1000L, 22L), as.character(1:22))
totals_b <- totals_a
totals_b["13"] <- 2000L   # chr 13 is excluded by default
totals_b["18"] <- 2000L   # chr 18 is excluded by default
totals_b["21"] <- 2000L   # chr 21 is excluded by default

# Also add two "normal" samples as needed for a valid control group (min 2)
sample_a <- .sim_nipter_sample(totals_a, "a", n_bins = 50L, seed = 1L)
sample_b <- .sim_nipter_sample(totals_b, "b", n_bins = 50L, seed = 2L)
sample_c <- .sim_nipter_sample(totals_a, "c", n_bins = 50L, seed = 3L)  # same as a

cg_excl <- nipter_as_control_group(list(sample_a, sample_b, sample_c))

# Query is same as sample_a; with default exclusions (13,18,21), sample_a
# and sample_c should both score ~0, sample_b should score ~0 too (chr 13/18/21
# not included in distance). But without exclusions, sample_b scores higher.
scores_excl    <- nipter_match_control_group(sample_a, cg_excl, n = 3L,
                                              mode = "report")
scores_no_excl <- nipter_match_control_group(sample_a, cg_excl, n = 3L,
                                              mode = "report",
                                              exclude_chromosomes = integer(0L))

# With default exclusions: sample_b SSD should be very small (near 0)
# because chr 13/18/21 differences are excluded
expect_true(scores_excl["b"] < scores_no_excl["b"],
            info = "excluding chromosomes 13/18/21 reduces distance from sample_b")

# With no exclusions: sample_b should be the most distant
expect_true(scores_no_excl["b"] > scores_no_excl["a"],
            info = "without exclusions sample_b is more distant than sample_a")
expect_true(scores_no_excl["b"] > scores_no_excl["c"],
            info = "without exclusions sample_b is more distant than sample_c")


# ---------------------------------------------------------------------------
# 5. nipter_control_group_from_beds() — round-trip via BAM fixture
# ---------------------------------------------------------------------------

.fixture_path <- function(name) {
  f <- system.file("extdata", name, package = "RWisecondorX")
  if (nzchar(f)) f else NULL
}
test_bam <- .fixture_path("hg00106_chr11_fixture.bam")

if (is.null(test_bam)) {
  exit_file("hg00106_chr11_fixture.bam not available; skipping from_beds tests")
}

# Write 3 copies of the same BED (representing 3 "different" control samples
# by giving each a distinct filename).
bed_dir <- tempfile("nipter_cg_beds_")
dir.create(bed_dir)
# Note: no on.exit() here — tinytest evaluates expressions one-by-one, so
# on.exit() would fire immediately after this expression, deleting bed_dir
# before the for-loop below runs.  Temp dirs are cleaned up when R exits.

bed_paths <- file.path(bed_dir, paste0("ctrl_", 1:3, ".bed.gz"))
for (p in bed_paths) {
  nipter_bin_bam_bed(test_bam, p, binsize = 50000L)
}

expect_true(all(file.exists(bed_paths)),
            info = "BED.gz fixture files created")
expect_true(all(file.exists(paste0(bed_paths, ".tbi"))),
            info = "tabix indexes created")

# Load via multi-reader
cg_from_beds <- nipter_control_group_from_beds(bed_dir, binsize = 50000L)

expect_true(.is_nipter_control_group(cg_from_beds),
            info = "nipter_control_group_from_beds returns an NIPTControlGroup")
expect_equal(length(cg_from_beds$samples), 3L,
             info = "3 samples loaded (one per BED file)")

# Sample names derived from filenames (strip extension)
loaded_names <- sort(vapply(cg_from_beds$samples, .sample_name_of, character(1L)))
expect_equal(loaded_names, paste0("ctrl_", 1:3),
             info = "sample names derived from filenames without extension")

# Each sample is CombinedStrands NIPTeRSample
for (s in cg_from_beds$samples) {
  expect_true(.is_nipter_sample(s),
              info = paste("sample", s$sample_name, "inherits NIPTSample"))
  expect_true(.is_combined_sample(s),
              info = paste("sample", s$sample_name, "inherits CombinedStrandsSample"))
}

# Compare chr11 read counts against a direct bed_to_nipter_sample() load
single <- bed_to_nipter_sample(bed_paths[1L])
multi_s1 <- cg_from_beds$samples[[
  which(vapply(cg_from_beds$samples, .sample_name_of, character(1L)) == "ctrl_1")
]]

# Chr11 counts should match exactly between single and multi reads
chr11_single <- single$autosomal_chromosome_reads[[1L]]["11", ]
chr11_multi  <- multi_s1$autosomal_chromosome_reads[[1L]]["11", ]
# Trim to shorter (multi may have global n_bins >= single)
n_cmp <- min(length(chr11_single), length(chr11_multi))
expect_equal(chr11_single[seq_len(n_cmp)], chr11_multi[seq_len(n_cmp)],
             info = "multi-read chr11 counts match single-file load")

# All samples share the same matrix dimensions (global n_bins)
dims <- lapply(cg_from_beds$samples, function(s) {
  dim(s$autosomal_chromosome_reads[[1L]])
})
expect_true(length(unique(dims)) == 1L,
            info = "all control samples share identical matrix dimensions")

# nipter_match_matrix() works on a from_beds control group
mm_beds <- nipter_match_matrix(cg_from_beds)
expect_equal(dim(mm_beds), c(3L, 3L),
             info = "nipter_match_matrix works on from_beds control group")
expect_equal(unname(diag(mm_beds)), rep(0.0, 3L),
             info = "from_beds match matrix diagonal is zero")


# ---------------------------------------------------------------------------
# 5b. nipter_control_group_from_beds() — SeparatedStrands auto-detection
#
# Writes 9-column BED.gz files (separate_strands = TRUE) and verifies that
# nipter_control_group_from_beds() returns a SeparatedStrands control group.
# ---------------------------------------------------------------------------

ss_bed_dir <- tempfile("nipter_ss_beds_")
dir.create(ss_bed_dir)

ss_bed_paths <- file.path(ss_bed_dir, paste0("ss_ctrl_", 1:3, ".bed.gz"))
for (p in ss_bed_paths) {
  nipter_bin_bam_bed(test_bam, p, binsize = 50000L, separate_strands = TRUE)
}

expect_true(all(file.exists(ss_bed_paths)),
            info = "SS BED.gz fixture files created")
expect_true(all(file.exists(paste0(ss_bed_paths, ".tbi"))),
            info = "SS tabix indexes created")

# Verify 9 tab-delimited columns per row
ss_lines <- readLines(gzcon(file(ss_bed_paths[1L], "rb")), n = 3L)
ss_lines <- ss_lines[nzchar(ss_lines)]
if (length(ss_lines) > 0L) {
  expect_true(all(lengths(strsplit(ss_lines, "\t")) == 9L),
              info = "SeparatedStrands BED has 9 columns")
}

cg_ss_beds <- nipter_control_group_from_beds(ss_bed_dir, pattern = "*.bed.gz",
                                              binsize = 50000L)

expect_true(.is_nipter_control_group(cg_ss_beds),
            info = "SS from_beds returns an NIPTControlGroup")
expect_true(S7::S7_inherits(cg_ss_beds, SeparatedControlGroup),
            info = "SS from_beds returns a SeparatedControlGroup")
expect_equal(length(cg_ss_beds$samples), 3L,
             info = "3 SeparatedStrands samples loaded")

ss_names_loaded <- sort(vapply(cg_ss_beds$samples, .sample_name_of, character(1L)))
expect_equal(ss_names_loaded, paste0("ss_ctrl_", 1:3),
             info = "SS sample names derived from filenames")

for (s in cg_ss_beds$samples) {
  expect_true(.is_nipter_sample(s),
              info = paste("SS sample", s$sample_name, "inherits NIPTSample"))
  expect_true(.is_separated_sample(s),
              info = paste("SS sample", s$sample_name, "inherits SeparatedStrandsSample"))
  expect_identical(length(s$autosomal_chromosome_reads), 2L,
                   info = paste("SS sample", s$sample_name, "has 2 auto matrices"))
  expect_true(all(grepl("F$", rownames(s$autosomal_chromosome_reads[[1L]]))),
              info = paste("SS sample", s$sample_name, "fwd auto rownames end in F"))
  expect_true(all(grepl("R$", rownames(s$autosomal_chromosome_reads[[2L]]))),
              info = paste("SS sample", s$sample_name, "rev auto rownames end in R"))
  expect_identical(rownames(s$sex_chromosome_reads[[1L]]), c("XF", "YF"),
                   info = paste("SS sample", s$sample_name, "fwd sex rownames XF/YF"))
  expect_identical(rownames(s$sex_chromosome_reads[[2L]]), c("XR", "YR"),
                   info = paste("SS sample", s$sample_name, "rev sex rownames XR/YR"))
}

# Verify all samples share identical matrix dimensions (global n_bins sync)
ss_dims <- lapply(cg_ss_beds$samples, function(s) {
  lapply(s$autosomal_chromosome_reads, dim)
})
expect_true(length(unique(lapply(ss_dims, `[[`, 1L))) == 1L,
            info = "all SS samples share identical fwd auto matrix dimensions")

# fwd + rev reads per chr should match raw bin count from single-file load
ss_single <- bed_to_nipter_sample(ss_bed_paths[1L])
ss_multi_s1 <- cg_ss_beds$samples[[
  which(vapply(cg_ss_beds$samples, .sample_name_of, character(1L)) == "ss_ctrl_1")
]]
chr11_fwd_single <- ss_single$autosomal_chromosome_reads[[1L]]["11F", ]
chr11_fwd_multi  <- ss_multi_s1$autosomal_chromosome_reads[[1L]]["11F", ]
n_cmp_ss <- min(length(chr11_fwd_single), length(chr11_fwd_multi))
expect_equal(chr11_fwd_single[seq_len(n_cmp_ss)],
             chr11_fwd_multi[seq_len(n_cmp_ss)],
             info = "SS multi-read chr11 fwd counts match single-file load")

# nipter_match_matrix() must accept SeparatedStrands control groups
mm_ss <- nipter_match_matrix(cg_ss_beds)
expect_equal(dim(mm_ss), c(3L, 3L),
             info = "nipter_match_matrix works on SS from_beds control group")
expect_equal(unname(diag(mm_ss)), rep(0.0, 3L),
             info = "SS from_beds match matrix diagonal is zero")


# ---------------------------------------------------------------------------
# 6. nipter_gc_precompute() and gc_table= parameter on nipter_gc_correct()
# ---------------------------------------------------------------------------

fixture_ref <- .fixture_path("fixture_ref.fa")
if (is.null(fixture_ref)) {
  exit_file("fixture_ref.fa not available; run `make fixtures` to generate GC tests")
}

# Load one BED sample for GC correction
sample_for_gc <- bed_to_nipter_sample(bed_paths[1L], binsize = 50000L)

# 6a. nipter_gc_precompute() creates a TSV.bgz + tbi
gc_out <- tempfile(fileext = ".tsv.bgz")
# No on.exit() — see note above about tinytest per-expression eval.

nipter_gc_precompute(fixture_ref, binsize = 50000L, out = gc_out)

expect_true(file.exists(gc_out),
            info = "nipter_gc_precompute creates TSV.bgz output")
expect_true(file.exists(paste0(gc_out, ".tbi")),
            info = "nipter_gc_precompute creates tabix index")
expect_true(file.info(gc_out)$size > 0L,
            info = "nipter_gc_precompute output is non-empty")

# 6b. nipter_gc_correct() with gc_table = path produces a NIPTeRSample.
# The fixture_ref.fa is a tiny all-A sequence (~50 kb = 1 bin), which is below
# the 10-valid-bin minimum for LOESS fitting.  The function therefore returns
# the sample uncorrected with a warning.  We check that it returns a valid
# NIPTeRSample without error, and that correction_status is a known value.
corrected_path <- suppressWarnings(
  nipter_gc_correct(sample_for_gc, gc_table = gc_out, binsize = 50000L)
)

expect_true(.is_nipter_sample(corrected_path),
            info = "nipter_gc_correct(gc_table=path) returns an NIPTSample")
expect_true(corrected_path$correction_status_autosomal %in%
              c("Uncorrected", "GC Corrected"),
            info = "correction_status_autosomal is a known value after gc_table= call")

# 6c. gc_table = path and fasta = path give the same correction
# (both hit the <10 valid bins limit with the tiny fixture ref)
corrected_fasta <- suppressWarnings(
  nipter_gc_correct(sample_for_gc, fasta = fixture_ref, binsize = 50000L)
)

auto_path  <- corrected_path$autosomal_chromosome_reads[[1L]]
auto_fasta <- corrected_fasta$autosomal_chromosome_reads[[1L]]
expect_equal(auto_path, auto_fasta, tolerance = 1e-10,
             info = "gc_table=path correction matches fasta= correction")

# 6d. Passing the precomputed TSV.bgz path a second time gives identical result
corrected_path2 <- suppressWarnings(
  nipter_gc_correct(sample_for_gc, gc_table = gc_out, binsize = 50000L)
)
expect_equal(corrected_path$autosomal_chromosome_reads[[1L]],
             corrected_path2$autosomal_chromosome_reads[[1L]],
             tolerance = 1e-10,
             info = "repeated gc_table=path calls produce identical results")


# ---------------------------------------------------------------------------
# 7. Real GC correction test — synthetic sample with GC bias
#
# Construct a NIPTeRSample with 22 chromosomes × 50 bins where read counts
# are systematically biased by GC content: high-GC bins get more reads.
# Verify that nipter_gc_correct() reduces this bias for both LOESS and
# bin-weight methods.
# ---------------------------------------------------------------------------

n_gc_bins <- 50L
n_gc_chrs <- 22L
base_reads <- 100  # base read count per bin

set.seed(123)

# Create synthetic GC table: each chromosome gets GC values that vary
# systematically from ~30% to ~70% across its bins.
gc_table_syn <- list()
for (chr in as.character(1:n_gc_chrs)) {
  # Smooth gradient with slight per-chr offset
  offset <- (as.integer(chr) - 1L) * 0.005
  gc_table_syn[[chr]] <- seq(0.30 + offset, 0.70 + offset,
                              length.out = n_gc_bins)
}
# Add sex chromosomes (required by the gc_table format)
gc_table_syn[["X"]] <- seq(0.35, 0.65, length.out = n_gc_bins)
gc_table_syn[["Y"]] <- seq(0.35, 0.65, length.out = n_gc_bins)

# Build a NIPTeRSample with GC bias: count = base * (1 + 2*(gc - 0.5))
# This means bins with gc=0.3 get base*0.6=60, gc=0.5 get base*1.0=100,
# gc=0.7 get base*1.4=140.
col_names <- as.character(seq_len(n_gc_bins))
auto_mat_gc <- matrix(0L, nrow = n_gc_chrs, ncol = n_gc_bins,
                      dimnames = list(as.character(1:n_gc_chrs), col_names))
for (chr in as.character(1:n_gc_chrs)) {
  gc_vals <- gc_table_syn[[chr]]
  # GC bias: counts proportional to GC, plus small Poisson noise
  biased_counts <- as.integer(round(base_reads * (1 + 2 * (gc_vals - 0.5)) +
                                     rnorm(n_gc_bins, 0, 3)))
  biased_counts[biased_counts < 1L] <- 1L
  auto_mat_gc[chr, ] <- biased_counts
}

sex_mat_gc <- matrix(0L, nrow = 2L, ncol = n_gc_bins,
                     dimnames = list(c("X", "Y"), col_names))

gc_sample <- structure(
  list(autosomal_chromosome_reads  = list(auto_mat_gc),
       sex_chromosome_reads        = list(sex_mat_gc),
       correction_status_autosomal = "Uncorrected",
       correction_status_sex       = "Uncorrected",
       sample_name                 = "gc_biased_sample"),
  class = c("NIPTeRSample", "CombinedStrands")
)

# Verify the uncorrected sample has substantial GC bias.
# Compute correlation between GC and read count across all autosomal bins.
gc_flat <- unlist(gc_table_syn[as.character(1:n_gc_chrs)])
reads_flat_before <- as.numeric(t(gc_sample$autosomal_chromosome_reads[[1L]]))
cor_before <- cor(gc_flat, reads_flat_before)

expect_true(abs(cor_before) > 0.8,
            info = paste0("synthetic sample has strong GC bias before correction",
                          " (cor = ", round(cor_before, 3), ")"))

# 7a. LOESS GC correction
corrected_loess <- nipter_gc_correct(gc_sample, gc_table = gc_table_syn,
                                     method = "loess", span = 0.75)

expect_true(.is_nipter_sample(corrected_loess),
            info = "LOESS gc_correct returns an NIPTSample")
expect_equal(corrected_loess$correction_status_autosomal, "GC Corrected",
             info = "LOESS gc_correct sets correction_status to 'GC Corrected'")

# After LOESS correction, GC-count correlation should be substantially reduced
reads_flat_loess <- as.numeric(t(corrected_loess$autosomal_chromosome_reads[[1L]]))
cor_loess <- cor(gc_flat, reads_flat_loess)

expect_true(abs(cor_loess) < abs(cor_before) * 0.3,
            info = paste0("LOESS correction reduces GC-count correlation: ",
                          round(cor_before, 3), " -> ", round(cor_loess, 3)))

# Corrected counts should not all be identical to uncorrected (something changed)
expect_false(all(abs(reads_flat_loess - reads_flat_before) < 1e-10),
             info = "LOESS correction actually changes bin counts")

# Corrected counts should be positive (no negative counts)
expect_true(all(reads_flat_loess >= 0),
            info = "LOESS corrected counts are non-negative")

# 7b. Bin-weight GC correction
corrected_bin <- nipter_gc_correct(gc_sample, gc_table = gc_table_syn,
                                   method = "bin")

expect_true(.is_nipter_sample(corrected_bin),
            info = "bin-weight gc_correct returns an NIPTSample")
expect_equal(corrected_bin$correction_status_autosomal, "GC Corrected",
             info = "bin-weight gc_correct sets correction_status to 'GC Corrected'")

reads_flat_bin <- as.numeric(t(corrected_bin$autosomal_chromosome_reads[[1L]]))
cor_bin <- cor(gc_flat, reads_flat_bin)

expect_true(abs(cor_bin) < abs(cor_before) * 0.3,
            info = paste0("bin-weight correction reduces GC-count correlation: ",
                          round(cor_before, 3), " -> ", round(cor_bin, 3)))

# 7c. Both methods should produce different results from each other
# (LOESS is smooth, bin-weight is discrete)
expect_false(all(abs(reads_flat_loess - reads_flat_bin) < 1e-10),
             info = "LOESS and bin-weight produce different corrections")

# 7d. Verify coefficient of variation is reduced after correction.
# The CV of reads at similar GC should decrease after correction.
cv_before <- sd(reads_flat_before) / mean(reads_flat_before)
cv_loess  <- sd(reads_flat_loess) / mean(reads_flat_loess)
cv_bin    <- sd(reads_flat_bin) / mean(reads_flat_bin)

expect_true(cv_loess < cv_before,
            info = paste0("LOESS correction reduces CV: ",
                          round(cv_before, 4), " -> ", round(cv_loess, 4)))
expect_true(cv_bin < cv_before,
            info = paste0("bin-weight correction reduces CV: ",
                          round(cv_before, 4), " -> ", round(cv_bin, 4)))
