# test_bed_reader.R — Round-trip tests for bed_to_sample() and
# bed_to_nipter_sample()
#
# Tests write bin counts to bed.gz via bam_convert_bed() / nipter_bin_bam_bed(),
# then read them back via bed_to_sample() / bed_to_nipter_sample() and verify
# that the in-memory structures match.

library(tinytest)
library(RWisecondorX)


# ---------- Fixture setup -------------------------------------------------

.fixture_path <- function(name) {
  f <- system.file("extdata", name, package = "RWisecondorX")
  if (nzchar(f)) f else NULL
}

mixed_bam <- .fixture_path("fixture_mixed.bam")
if (is.null(mixed_bam)) {
  exit_file("Synthetic BAM fixtures not available; run `make fixtures`")
}


# ==========================================================================
# bed_to_sample() round-trip — WisecondorX format (4-column BED)
# ==========================================================================

bins_orig <- bam_convert(mixed_bam, binsize = 5000L, rmdup = "streaming")

bed_file <- tempfile(fileext = ".bed.gz")
bam_convert_bed(mixed_bam, bed_file, binsize = 5000L, rmdup = "streaming")

expect_true(file.exists(bed_file),
            info = "bed.gz file created by bam_convert_bed()")

# Read back
bins_rt <- bed_to_sample(bed_file)

expect_true(is.list(bins_rt),
            info = "bed_to_sample() returns a list")
expect_identical(length(bins_rt), 24L,
                 info = "bed_to_sample() returns 24 chromosome keys")
expect_identical(names(bins_rt), as.character(1:24),
                 info = "bed_to_sample() keys are '1'-'24'")

# Compare chr11 counts (the only chromosome with reads in this fixture)
expect_identical(bins_rt[["11"]], bins_orig[["11"]],
                 info = "bed_to_sample() chr11 round-trips exactly")

# Total read count preserved
total_orig <- sum(vapply(bins_orig, function(x) if (is.null(x)) 0L else sum(x),
                         integer(1L)))
total_rt   <- sum(vapply(bins_rt, function(x) if (is.null(x)) 0L else sum(x),
                         integer(1L)))
expect_identical(total_rt, total_orig,
                 info = "bed_to_sample() preserves total read count")

# Binsize inference
bins_rt2 <- bed_to_sample(bed_file, binsize = 5000L)
expect_identical(bins_rt2[["11"]], bins_orig[["11"]],
                 info = "bed_to_sample() with explicit binsize matches")


# ==========================================================================
# bed_to_nipter_sample() round-trip — CombinedStrands (5-column BED)
# ==========================================================================

nipter_orig <- nipter_bin_bam(mixed_bam, binsize = 5000L, rmdup = "none")

nipter_bed <- tempfile(fileext = ".bed.gz")
nipter_bin_bam_bed(mixed_bam, nipter_bed, binsize = 5000L, rmdup = "none")

expect_true(file.exists(nipter_bed),
            info = "NIPTeR bed.gz file created")

nipter_rt <- bed_to_nipter_sample(nipter_bed)

expect_true(inherits(nipter_rt, "NIPTeRSample"),
            info = "bed_to_nipter_sample() returns NIPTeRSample")
expect_true(inherits(nipter_rt, "CombinedStrands"),
            info = "5-column BED produces CombinedStrands")

# Matrix dimensions match
auto_orig <- nipter_orig$autosomal_chromosome_reads[[1L]]
auto_rt   <- nipter_rt$autosomal_chromosome_reads[[1L]]
expect_identical(nrow(auto_rt), 22L,
                 info = "autosomal matrix has 22 rows")
expect_identical(rownames(auto_rt), as.character(1:22),
                 info = "autosomal rownames are '1'-'22'")
expect_identical(ncol(auto_rt), ncol(auto_orig),
                 info = "autosomal column count matches original")

sex_orig <- nipter_orig$sex_chromosome_reads[[1L]]
sex_rt   <- nipter_rt$sex_chromosome_reads[[1L]]
expect_identical(nrow(sex_rt), 2L,
                 info = "sex matrix has 2 rows")
expect_identical(rownames(sex_rt), c("X", "Y"),
                 info = "sex rownames are 'X','Y'")
expect_identical(ncol(sex_rt), ncol(sex_orig),
                 info = "sex column count matches original")

# Bin counts round-trip exactly
expect_identical(as.integer(auto_rt), as.integer(auto_orig),
                 info = "autosomal bin counts round-trip exactly")
expect_identical(as.integer(sex_rt), as.integer(sex_orig),
                 info = "sex bin counts round-trip exactly")

# Correction status
expect_identical(nipter_rt$correction_status_autosomal, "Uncorrected",
                 info = "correction status is Uncorrected for raw BED")

# Sample name inference
expect_identical(nipter_rt$sample_name,
                 sub("\\.bed(\\.gz)?$", "", basename(nipter_bed)),
                 info = "sample name inferred from BED filename")

# Explicit name
nipter_rt_named <- bed_to_nipter_sample(nipter_bed, name = "my_sample")
expect_identical(nipter_rt_named$sample_name, "my_sample",
                 info = "explicit sample name honoured")


# ==========================================================================
# bed_to_nipter_sample() round-trip — SeparatedStrands (7-column BED)
# ==========================================================================

nipter_ss_orig <- nipter_bin_bam(mixed_bam, binsize = 5000L, rmdup = "none",
                                  separate_strands = TRUE)

ss_bed <- tempfile(fileext = ".bed.gz")
nipter_bin_bam_bed(mixed_bam, ss_bed, binsize = 5000L, rmdup = "none",
                   separate_strands = TRUE)

expect_true(file.exists(ss_bed),
            info = "SeparatedStrands bed.gz file created")

nipter_ss_rt <- bed_to_nipter_sample(ss_bed)

expect_true(inherits(nipter_ss_rt, "NIPTeRSample"),
            info = "SeparatedStrands bed_to_nipter_sample() returns NIPTeRSample")
expect_true(inherits(nipter_ss_rt, "SeparatedStrands"),
            info = "7-column BED produces SeparatedStrands")

# Forward autosomal
fwd_auto_orig <- nipter_ss_orig$autosomal_chromosome_reads[[1L]]
fwd_auto_rt   <- nipter_ss_rt$autosomal_chromosome_reads[[1L]]
expect_identical(nrow(fwd_auto_rt), 22L,
                 info = "SeparatedStrands fwd auto has 22 rows")
expect_identical(rownames(fwd_auto_rt), paste0(1:22, "F"),
                 info = "SeparatedStrands fwd auto rownames '1F'-'22F'")
expect_identical(as.integer(fwd_auto_rt), as.integer(fwd_auto_orig),
                 info = "SeparatedStrands fwd auto counts round-trip")

# Reverse autosomal
rev_auto_orig <- nipter_ss_orig$autosomal_chromosome_reads[[2L]]
rev_auto_rt   <- nipter_ss_rt$autosomal_chromosome_reads[[2L]]
expect_identical(rownames(rev_auto_rt), paste0(1:22, "R"),
                 info = "SeparatedStrands rev auto rownames '1R'-'22R'")
expect_identical(as.integer(rev_auto_rt), as.integer(rev_auto_orig),
                 info = "SeparatedStrands rev auto counts round-trip")

# Forward sex
fwd_sex_orig <- nipter_ss_orig$sex_chromosome_reads[[1L]]
fwd_sex_rt   <- nipter_ss_rt$sex_chromosome_reads[[1L]]
expect_identical(rownames(fwd_sex_rt), c("XF", "YF"),
                 info = "SeparatedStrands fwd sex rownames")
expect_identical(as.integer(fwd_sex_rt), as.integer(fwd_sex_orig),
                 info = "SeparatedStrands fwd sex counts round-trip")

# Reverse sex
rev_sex_orig <- nipter_ss_orig$sex_chromosome_reads[[2L]]
rev_sex_rt   <- nipter_ss_rt$sex_chromosome_reads[[2L]]
expect_identical(rownames(rev_sex_rt), c("XR", "YR"),
                 info = "SeparatedStrands rev sex rownames")
expect_identical(as.integer(rev_sex_rt), as.integer(rev_sex_orig),
                 info = "SeparatedStrands rev sex counts round-trip")


# ==========================================================================
# bed_to_sample() → rwisecondorx pipeline integration
# ==========================================================================

# Verify bed_to_sample() output can be passed to scale_sample()
bins_scaled <- scale_sample(bins_rt, from_size = 5000L, to_size = 10000L)
expect_true(is.list(bins_scaled),
            info = "bed_to_sample() output works with scale_sample()")
