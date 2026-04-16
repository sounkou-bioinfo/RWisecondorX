# helper_nipter.R — shared simulation helpers for NIPTeR tinytest suite
#
# Auto-sourced by tinytest before each test_*.R file that calls
# tinytest::test_package() or tinytest::test_dir(). Functions defined here are
# therefore available in every test file without duplication.
#
# Naming convention:
#   .sim_*   — pure simulation (no BAM / fixture dependency)
#   .fixture_path() — helper for inst/extdata fixture lookup
#
# Implementation note: reads-per-chromosome are distributed with the first
# `rem` bins getting one extra read (deterministic given a fixed total),
# not via sample() — so results are fully reproducible without set.seed()
# when using the deterministic chr_totals helpers.  Callers that need RNG-
# seeded totals should set.seed() before calling .sim_chr_totals().

# ---------------------------------------------------------------------------
# Chromosome size table (hg38, Mb rounded)
# ---------------------------------------------------------------------------

.SIM_CHR_SIZES <- c(248, 242, 198, 190, 182, 171, 159, 145, 138, 134, 135,
                    133, 114, 107, 102,  90,  83,  80,  59,  64,  47,  51)

# ---------------------------------------------------------------------------
# .sim_chr_totals(scale, trisomy_chr, trisomy_frac)
#
# Returns a named integer vector (keys "1"-"22") with per-chromosome total
# read counts proportional to chromosome size, optionally boosted on one
# chromosome to simulate a trisomy.
# ---------------------------------------------------------------------------

.sim_chr_totals <- function(scale = 1, trisomy_chr = NULL, trisomy_frac = 0.05) {
  totals <- as.integer(round(.SIM_CHR_SIZES * 10 * scale))
  names(totals) <- as.character(1:22)
  if (!is.null(trisomy_chr)) {
    k <- as.character(trisomy_chr)
    extra <- as.integer(round(totals[k] * trisomy_frac))
    totals[k] <- totals[k] + extra
  }
  totals
}

# ---------------------------------------------------------------------------
# .sim_fill_bins(total, n_bins)
#
# Distributes `total` reads across `n_bins` bins.  The first `rem` bins get
# base+1 reads; the rest get base.  Deterministic for a fixed total.
# ---------------------------------------------------------------------------

.sim_fill_bins <- function(total, n_bins) {
  if (total == 0L) return(integer(n_bins))
  base <- total %/% n_bins
  rem  <- total - base * n_bins
  counts <- rep(base, n_bins)
  if (rem > 0L) counts[seq_len(rem)] <- counts[seq_len(rem)] + 1L
  as.integer(counts)
}

# ---------------------------------------------------------------------------
# .sim_nipter_sample(chr_totals, name, n_bins, seed)
#
# Builds a CombinedStrands NIPTeRSample from a named integer vector of per-
# chromosome total read counts.
# ---------------------------------------------------------------------------

.sim_nipter_sample <- function(chr_totals, name, n_bins = 100L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  col_names <- as.character(seq_len(n_bins))

  auto_mat <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(as.character(1:22), col_names))
  for (chr in names(chr_totals)) {
    auto_mat[chr, ] <- .sim_fill_bins(chr_totals[chr], n_bins)
  }

  sex_mat <- matrix(0L, nrow = 2L, ncol = n_bins,
                    dimnames = list(c("X", "Y"), col_names))

  CombinedStrandsSample(
    sample_name = name,
    binsize = as.integer(50000L),
    auto_matrix = auto_mat,
    sex_matrix_ = sex_mat
  )
}

# ---------------------------------------------------------------------------
# .sim_nipter_control_set(n, scale, seed, n_bins)
#
# Returns a list of n CombinedStrands NIPTeRSamples with ±5% noise per
# chromosome.  Each sample's RNG seed is derived from `seed * i`.
# ---------------------------------------------------------------------------

.sim_nipter_control_set <- function(n = 10L, scale = 1, seed = 42L,
                                    n_bins = 100L) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    noise  <- runif(22L, 0.95, 1.05)
    totals <- as.integer(round(.sim_chr_totals(scale = scale) * noise))
    names(totals) <- as.character(1:22)
    .sim_nipter_sample(totals, sprintf("ctrl_%02d", i),
                       n_bins = n_bins, seed = seed * i)
  })
}

# ---------------------------------------------------------------------------
# .sim_nipter_ss_sample(chr_totals, name, n_bins, seed)
#
# Builds a SeparatedStrands NIPTeRSample.  Reads are split 50/50 between
# forward and reverse strands per chromosome.
# ---------------------------------------------------------------------------

.sim_nipter_ss_sample <- function(chr_totals, name, n_bins = 100L,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  col_names <- as.character(seq_len(n_bins))

  fwd_auto <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(paste0(1:22, "F"), col_names))
  rev_auto <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(paste0(1:22, "R"), col_names))

  for (chr in names(chr_totals)) {
    fwd_total <- as.integer(round(chr_totals[chr] * 0.5))
    rev_total <- chr_totals[chr] - fwd_total
    fwd_auto[paste0(chr, "F"), ] <- .sim_fill_bins(fwd_total, n_bins)
    rev_auto[paste0(chr, "R"), ] <- .sim_fill_bins(rev_total, n_bins)
  }

  fwd_sex <- matrix(0L, nrow = 2L, ncol = n_bins,
                    dimnames = list(c("XF", "YF"), col_names))
  rev_sex <- matrix(0L, nrow = 2L, ncol = n_bins,
                    dimnames = list(c("XR", "YR"), col_names))

  SeparatedStrandsSample(
    sample_name = name,
    binsize = as.integer(50000L),
    auto_fwd = fwd_auto,
    auto_rev = rev_auto,
    sex_fwd = fwd_sex,
    sex_rev = rev_sex
  )
}

# ---------------------------------------------------------------------------
# .sim_nipter_ss_control_set(n, scale, seed, n_bins)
#
# Returns a list of n SeparatedStrands NIPTeRSamples with ±5% noise.
# ---------------------------------------------------------------------------

.sim_nipter_ss_control_set <- function(n = 10L, scale = 1, seed = 42L,
                                       n_bins = 100L) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    noise  <- runif(22L, 0.95, 1.05)
    totals <- as.integer(round(.sim_chr_totals(scale = scale) * noise))
    names(totals) <- as.character(1:22)
    .sim_nipter_ss_sample(totals, sprintf("ss_ctrl_%02d", i),
                          n_bins = n_bins, seed = seed + i)
  })
}

# ---------------------------------------------------------------------------
# .fixture_path(name)
#
# Returns the path to an inst/extdata fixture file, or NULL if absent.
# ---------------------------------------------------------------------------

.fixture_path <- function(name) {
  f <- system.file("extdata", name, package = "RWisecondorX")
  if (nzchar(f)) f else NULL
}
