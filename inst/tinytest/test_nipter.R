library(tinytest)
library(RWisecondorX)

.fixture_path <- function(name) {
  f <- system.file("extdata", name, package = "RWisecondorX")
  if (nzchar(f)) f else NULL
}

test_bam <- .fixture_path("hg00106_chr11_fixture.bam")
if (is.null(test_bam)) {
  exit_file("hg00106_chr11_fixture.bam not available")
}

# ---------------------------------------------------------------------------
# Structure: NIPTeRSample object
# ---------------------------------------------------------------------------

sample <- nipter_bin_bam(test_bam, binsize = 50000L)

expect_true(inherits(sample, "NIPTeRSample"),    info = "result has class NIPTeRSample")
expect_true(inherits(sample, "CombinedStrands"), info = "result has class CombinedStrands")

expect_true(is.list(sample$autosomal_chromosome_reads), info = "autosomal_chromosome_reads is a list")
expect_true(is.list(sample$sex_chromosome_reads),       info = "sex_chromosome_reads is a list")
expect_identical(sample$correction_status_autosomal, "Uncorrected",
                 info = "autosomal correction status is Uncorrected")
expect_identical(sample$correction_status_sex, "Uncorrected",
                 info = "sex correction status is Uncorrected")
expect_true(is.character(sample$sample_name) && nzchar(sample$sample_name),
            info = "sample_name is a non-empty string")
expect_identical(sample$sample_name, "hg00106_chr11_fixture",
                 info = "sample_name is BAM basename without extension")

auto_mat <- sample$autosomal_chromosome_reads[[1L]]
sex_mat  <- sample$sex_chromosome_reads[[1L]]

expect_identical(nrow(auto_mat), 22L,           info = "autosomal matrix has 22 rows")
expect_identical(rownames(auto_mat), as.character(1:22),
                 info = "autosomal rownames are 1-22")
expect_identical(nrow(sex_mat), 2L,             info = "sex matrix has 2 rows")
expect_identical(rownames(sex_mat), c("X", "Y"), info = "sex rownames are X and Y")
expect_identical(ncol(auto_mat), ncol(sex_mat), info = "autosomal and sex matrices share bin width")
expect_true(is.integer(auto_mat), info = "autosomal matrix is integer")
expect_true(is.integer(sex_mat),  info = "sex matrix is integer")
expect_true(all(auto_mat >= 0L),  info = "autosomal counts are non-negative")
expect_true(all(sex_mat  >= 0L),  info = "sex counts are non-negative")

# The fixture BAM is chr11-only; all other autosomes should be zero.
expect_true(sum(auto_mat["11", ]) > 0L, info = "chr11 has reads in the fixture")
for (chr in setdiff(as.character(1:22), "11")) {
  expect_true(sum(auto_mat[chr, ]) == 0L,
              info = paste("chr", chr, "is zero (fixture is chr11-only)"))
}

# ---------------------------------------------------------------------------
# MAPQ filter reduces read count
# ---------------------------------------------------------------------------

sample_mapq40 <- nipter_bin_bam(test_bam, binsize = 50000L, mapq = 40L)

total_default <- sum(auto_mat) + sum(sex_mat)
total_mapq40  <- sum(sample_mapq40$autosomal_chromosome_reads[[1L]]) +
                 sum(sample_mapq40$sex_chromosome_reads[[1L]])

expect_true(total_mapq40 <= total_default, info = "mapq=40 keeps <= reads than mapq=0")
expect_true(total_mapq40 > 0L,            info = "mapq=40 retains some reads")

# ---------------------------------------------------------------------------
# exclude_flags=1024 removes duplicate-flagged reads
# ---------------------------------------------------------------------------

sample_nodup <- nipter_bin_bam(test_bam, binsize = 50000L, exclude_flags = 1024L)

total_nodup <- sum(sample_nodup$autosomal_chromosome_reads[[1L]]) +
               sum(sample_nodup$sex_chromosome_reads[[1L]])

expect_true(total_nodup <= total_default, info = "exclude_flags=1024 keeps <= reads than no filter")
# Fixture has 68 duplicates; total should drop
expect_true(total_nodup < total_default,  info = "exclude_flags=1024 removes some reads (fixture has duplicates)")

# ---------------------------------------------------------------------------
# BED output
# ---------------------------------------------------------------------------

bed_out <- tempfile(fileext = ".bed.gz")
on.exit(unlink(c(bed_out, paste0(bed_out, ".tbi"))), add = TRUE)

nipter_bin_bam_bed(test_bam, bed_out, binsize = 50000L)

expect_true(file.exists(bed_out),                    info = "BED.gz file created")
expect_true(file.exists(paste0(bed_out, ".tbi")),    info = "tabix index created")
expect_true(file.info(bed_out)$size > 0L,            info = "BED.gz is non-empty")

# Spot-check: 24 chromosomes × n_bins rows, 5 columns (tab-delimited)
lines <- readLines(gzcon(file(bed_out, "rb")), n = 5L)
expect_identical(lengths(strsplit(lines, "\t")), rep(5L, 5L),
                 info = "BED rows have 5 tab-delimited columns")

# ---------------------------------------------------------------------------
# Conformance against NIPTeR::bin_bam_sample()
#
# NIPTeR's bin_bam_sample() has two known bugs that constrain the conformance
# target:
#
#   1. Split length mismatch: split(reads, droplevels(strands[strands != "*"]))
#      drops unmapped-read entries from the factor but not from reads, so it
#      silently misbehaves on any BAM that contains unmapped reads.
#
#   2. Implicit positional dedup: bin_reads() calls unique() on positions per
#      chromosome per strand before binning.  Two reads at the same start
#      position are counted once, not twice — an implicit, flag-ignorant
#      deduplication that is NOT equivalent to -F 1024.
#
# For a meaningful bin-for-bin comparison the conformance BAM must therefore:
#   - Contain no unmapped reads  (avoids bug 1; pre-filter with -F 4)
#   - Have no two reads sharing the same start position per strand  (avoids
#     the unique() discrepancy when using rmdup = "none")
#
# The package ships a bundled whole-genome conformance fixture satisfying
# those constraints. Set NIPTER_CONFORMANCE_BAM only to override it with a
# custom pre-filtered BAM.
# ---------------------------------------------------------------------------

if (!requireNamespace("NIPTeR", quietly = TRUE)) {
  exit_file("NIPTeR not installed; skipping cross-package conformance")
}

conf_bam <- Sys.getenv("NIPTER_CONFORMANCE_BAM", unset = NA_character_)
if (is.na(conf_bam) || !nzchar(conf_bam) || !file.exists(conf_bam)) {
  conf_bam <- .fixture_path("nipter_conformance_fixture.bam")
}
if (is.null(conf_bam) || !file.exists(conf_bam)) {
  exit_file("NIPTeR conformance fixture not available; run `make fixtures`")
}

our_sample <- nipter_bin_bam(conf_bam, binsize = 50000L,
                             mapq = 0L, rmdup = "none")

ref_sample <- tryCatch(
  NIPTeR::bin_bam_sample(bam_filepath = conf_bam, do_sort = FALSE,
                         separate_strands = FALSE),
  error = function(e) {
    exit_file(paste("NIPTeR::bin_bam_sample failed:", conditionMessage(e)))
  }
)

ref_auto <- ref_sample$autosomal_chromosome_reads[[1L]]
our_auto <- our_sample$autosomal_chromosome_reads[[1L]]
ref_sex  <- ref_sample$sex_chromosome_reads[[1L]]
our_sex  <- our_sample$sex_chromosome_reads[[1L]]

# Total read counts must agree: NIPTeR drops strand="*" (unmapped),
# we drop rname IS NULL — equivalent filters.
expect_identical(sum(our_auto) + sum(our_sex),
                 sum(ref_auto) + sum(ref_sex),
                 info = "total read count matches NIPTeR")

# Bin-for-bin comparison; trim to the shorter column count (edge-bin tolerance).
n_auto <- min(ncol(ref_auto), ncol(our_auto))
n_sex  <- min(ncol(ref_sex),  ncol(our_sex))

for (chr in as.character(1:22)) {
  expect_identical(our_auto[chr, seq_len(n_auto)],
                   ref_auto[chr, seq_len(n_auto)],
                   info = paste("chr", chr, "autosomal bins match NIPTeR"))
}
expect_identical(our_sex["X", seq_len(n_sex)], ref_sex["X", seq_len(n_sex)],
                 info = "chrX bins match NIPTeR")
expect_identical(our_sex["Y", seq_len(n_sex)], ref_sex["Y", seq_len(n_sex)],
                 info = "chrY bins match NIPTeR")
