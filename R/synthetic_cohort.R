# synthetic_cohort.R — Synthetic BAM cohort generator for WisecondorX / NIPTeR smoke tests
#
# Generates coordinate-sorted, indexed BAM files using "compressed" chromosome
# lengths: each 100kb genomic bin maps to 100bp, so bin COUNT structure is
# identical to real GRCh37 at 100kb resolution, but BAMs are tiny (~435KB each).
#
# This generator is for fast package testing and algorithm smoke tests. It is
# intentionally not a calibrated clinical simulator, and should not be used to
# justify sensitivity/specificity claims on real NIPT cohorts.
#
# Cohort composition (default 50 samples, configurable via n_samples):
#   ~75% euploid female (XX, zero Y reads)
#   ~25% euploid male   (XY, ~50% X/Y relative to autosomes)
#    3 trisomy female (T21, T18, T13)
#
# GC bias simulation: each bin has a fixed pseudo-GC content (spatial random
# walk), and each sample has a random quadratic GC bias curve. This creates
# realistic per-bin coverage variance between samples, ensuring that reference
# null ratio distributions have proper spread for Z-score computation.
#
# Trisomy simulation follows Nguyen et al. 2023 (doi:10.1101/2023.11.24.568620):
# instead of naively multiplying reads on the trisomy chromosome by 1.5x, we
# model a fetal fraction f and REMOVE reads from non-target chromosomes so that
# the relative read coverage of the trisomy chromosome reflects an extra copy
# in the fetal genome.
#
# For non-mosaic trisomy at fetal fraction f:
#   - Reads on trisomy chromosome h are kept unchanged
#   - Reads on every other chromosome i are randomly removed with probability
#     p = 0.5*f / (1 + 0.5*f)
#   - This makes relative coverage C'h = Ch * (1 + 0.5*f), matching the
#     expected enrichment from an extra fetal copy
#
# For mosaic trisomy (50% placental mosaicism):
#   - Removal probability is p = 0.25*f / (1 + 0.25*f)
#   - Relative enrichment is C''h = Ch * (1 + 0.25*f)
#
# Default fetal fraction is 10% (f = 0.10), giving ~5% relative enrichment
# for non-mosaic trisomy — realistic for NIPT.
#
# Deterministic: sample i uses set.seed(42 + i).
# NIPTeR-compatible: no unmapped reads, unique positions per chromosome.
# Uses a native htslib-backed BAM writer plus Rduckhts BAM indexing — no
# external samtools dependency.


# ---- GRCh37 bin structure (package-internal constants) --------------------

# Number of 100kb bins per chromosome (ceiling(chr_length / 100000))
BINS_PER_CHR <- c(
  "1"  = 2493L, "2"  = 2432L, "3"  = 1981L, "4"  = 1912L,
  "5"  = 1810L, "6"  = 1712L, "7"  = 1592L, "8"  = 1464L,
  "9"  = 1413L, "10" = 1356L, "11" = 1351L, "12" = 1339L,
  "13" = 1152L, "14" = 1074L, "15" = 1026L, "16" =  904L,
  "17" =  812L, "18" =  781L, "19" =  592L, "20" =  631L,
  "21" =  482L, "22" =  514L, "X"  = 1553L, "Y"  =  594L
)

#' Compressed bin size for synthetic cohort BAMs
#'
#' Each 100kb genomic bin is compressed to 100bp, yielding identical bin
#' counts to GRCh37 at 100kb resolution.
#'
#' @export
COMPRESSED_BINSIZE <- 100L

# Compressed chromosome lengths: bins * 100bp per bin
COMPRESSED_CHR_LENGTHS <- BINS_PER_CHR * COMPRESSED_BINSIZE

# SAM chromosome names (no "chr" prefix, matching GRCh37/b37 convention)
CHR_NAMES <- names(BINS_PER_CHR)

# Read properties
READ_LEN <- 36L
# Default average reads per bin for euploid autosomes
READS_PER_BIN_DEFAULT <- 3.0


# ---- Sample profile constructors (internal) ------------------------------

.make_profile <- function(id, sex, chr_mult, trisomy_chr = NULL,
                          fetal_fraction = NULL, mosaic = FALSE) {
  list(id = id, sex = sex, chr_mult = chr_mult,
       trisomy_chr = trisomy_chr, fetal_fraction = fetal_fraction,
       mosaic = mosaic)
}

.euploid_female <- function(id) {
  mult <- stats::setNames(rep(1.0, 24), CHR_NAMES)
  mult["Y"] <- 0.0
  .make_profile(id, "F", mult)
}

.euploid_male <- function(id) {
  mult <- stats::setNames(rep(1.0, 24), CHR_NAMES)
  mult["X"] <- 0.5
  mult["Y"] <- 0.5
  .make_profile(id, "M", mult)
}

#' Build a trisomy profile using Nguyen et al. fetal fraction model.
#'
#' Euploid read counts are generated first (all chr_mult = 1.0 except Y=0),
#' then reads on non-target chromosomes are stochastically removed in
#' `.sample_positions()` to simulate the relative enrichment of the trisomy
#' chromosome.
#' @noRd
.trisomy_female <- function(id, trisomy_chr, fetal_fraction = 0.10,
                            mosaic = FALSE) {
  # chr_mult is 1.0 for all chromosomes — the trisomy signal is created

  # by REMOVING reads from non-target chromosomes in .sample_positions()
  mult <- stats::setNames(rep(1.0, 24), CHR_NAMES)
  mult["Y"] <- 0.0
  .make_profile(id, "F", mult, trisomy_chr = as.character(trisomy_chr),
                fetal_fraction = fetal_fraction, mosaic = mosaic)
}


# ---- Cohort definition ---------------------------------------------------

.build_cohort_profiles <- function(n_samples = 50L, fetal_fraction = 0.10,
                                   mosaic = FALSE) {
  n_samples <- as.integer(n_samples)
  stopifnot(n_samples >= 6L)  # need room for 3 trisomy + at least 1M + 2F

  n_trisomy <- 3L
  n_euploid <- n_samples - n_trisomy
  n_male    <- max(1L, round(n_euploid * 0.25))
  n_female  <- n_euploid - n_male

  profiles <- vector("list", n_samples)
  idx <- 1L

  for (i in seq_len(n_female)) {
    profiles[[idx]] <- .euploid_female(sprintf("euploid_f_%03d", i))
    idx <- idx + 1L
  }
  for (i in seq_len(n_male)) {
    profiles[[idx]] <- .euploid_male(sprintf("euploid_m_%03d", i))
    idx <- idx + 1L
  }

  profiles[[idx]] <- .trisomy_female("trisomy_21", "21",
                                     fetal_fraction, mosaic); idx <- idx + 1L
  profiles[[idx]] <- .trisomy_female("trisomy_18", "18",
                                     fetal_fraction, mosaic); idx <- idx + 1L
  profiles[[idx]] <- .trisomy_female("trisomy_13", "13",
                                     fetal_fraction, mosaic)

  profiles
}


# ---- GC-like coverage variance (internal) ---------------------------------
#
# Real NIPT samples exhibit per-bin coverage variance driven by GC content:
# some bins consistently get more/fewer reads depending on their nucleotide
# composition, and the strength/shape of this bias varies between samples
# (different library preps, extraction methods, sequencing runs).
#
# We simulate this with two components:
#   1. A fixed "GC landscape" — a per-bin pseudo-GC content in [0.3, 0.7],
#      generated once with a fixed seed. This creates spatial structure
#      (neighboring bins have correlated GC) that mimics real genomes.
#   2. A per-sample "GC bias curve" — a quadratic function mapping GC content
#      to a coverage multiplier, with coefficients drawn from a normal
#      distribution. Each sample's curve differs, creating realistic
#      between-sample variance while keeping the mean multiplier near 1.0.
#
# The result: bins with extreme GC (high or low) have more variable coverage
# across samples than bins near 50% GC, matching real sequencing behavior.

#' Generate a per-bin GC landscape for the compressed genome.
#'
#' Returns a named list of numeric vectors (one per chromosome), each of
#' length `BINS_PER_CHR[chr]`, with values between 0.3 and 0.7 representing
#' pseudo-GC content. The landscape is deterministic (fixed seed = 1337).
#'
#' Uses a smooth random walk per chromosome to create spatial correlation
#' between neighboring bins, matching real GC structure.
#' @noRd
.gc_landscape <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    saved_seed <- .Random.seed
    on.exit(assign(".Random.seed", saved_seed, envir = .GlobalEnv))
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }
  set.seed(1337L)

  gc <- vector("list", length(CHR_NAMES))
  names(gc) <- CHR_NAMES

  for (chr_name in CHR_NAMES) {
    n <- BINS_PER_CHR[chr_name]
    # Random walk with drift toward 0.5 (mean-reverting)
    walk <- numeric(n)
    walk[1] <- 0.5
    for (i in 2:n) {
      walk[i] <- walk[i - 1L] + stats::rnorm(1, sd = 0.02) +
        0.01 * (0.5 - walk[i - 1L])
    }
    # Clamp and rescale to [0.3, 0.7]
    walk <- pmin(pmax(walk, 0), 1)
    gc[[chr_name]] <- 0.3 + 0.4 * (walk - min(walk)) /
      max(max(walk) - min(walk), 1e-10)
  }
  gc
}

#' Generate per-bin coverage multipliers for one sample.
#'
#' Given a GC landscape and a sample-specific seed, generates a quadratic
#' GC bias curve and applies it to get per-bin multipliers. The multipliers
#' have mean ~1.0 (rescaled) so total read count is preserved.
#'
#' @param gc_landscape Named list from `.gc_landscape()`.
#' @param sample_seed Integer seed for this sample's bias curve.
#' @param gc_bias_strength Scaling factor for GC bias amplitude (default 0.3).
#'   Higher values produce more between-sample variance. 0 = no GC bias.
#' @return Named list of numeric vectors, one per chromosome.
#' @noRd
.gc_bias_multipliers <- function(gc_landscape, sample_seed,
                                  gc_bias_strength = 0.3) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    saved_seed <- .Random.seed
    on.exit(assign(".Random.seed", saved_seed, envir = .GlobalEnv))
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }
  set.seed(sample_seed)

  # Quadratic bias: mult(gc) = 1 + a*(gc - 0.5) + b*(gc - 0.5)^2
  # a controls linear tilt, b controls curvature
  a <- stats::rnorm(1, sd = gc_bias_strength)
  b <- stats::rnorm(1, sd = gc_bias_strength * 2)

  mults <- vector("list", length(gc_landscape))
  names(mults) <- names(gc_landscape)

  all_vals <- numeric(0)
  for (chr_name in names(gc_landscape)) {
    gc <- gc_landscape[[chr_name]]
    centered_gc <- gc - 0.5
    m <- 1.0 + a * centered_gc + b * centered_gc^2
    m <- pmax(m, 0.1)  # floor at 0.1 to avoid negative/zero
    mults[[chr_name]] <- m
    all_vals <- c(all_vals, m)
  }

  # Rescale so mean multiplier = 1.0 (preserves total read count)
  global_mean <- mean(all_vals)
  for (chr_name in names(mults)) {
    mults[[chr_name]] <- mults[[chr_name]] / global_mean
  }

  mults
}


# ---- Read-position generation (internal) ----------------------------------

#' Compute trisomy read-removal probability (Nguyen et al. 2023).
#'
#' For non-target chromosomes i != h, reads are removed with probability p
#' so that the relative coverage of the trisomy chromosome h increases to
#' C'h = Ch * (1 + k*f), where k = 0.5 for non-mosaic, k = 0.25 for mosaic.
#'
#' Derivation: keeping fraction (1 - p) on non-target chromosomes is equivalent
#' to multiplying their coverage by 1/(1 + k*f). Solving: p = k*f / (1 + k*f).
#' @noRd
.trisomy_removal_prob <- function(fetal_fraction, mosaic = FALSE) {
  k <- if (isTRUE(mosaic)) 0.25 else 0.5
  kf <- k * fetal_fraction
  kf / (1 + kf)
}

.sample_positions <- function(profile, seed,
                              reads_per_bin = READS_PER_BIN_DEFAULT,
                              gc_multipliers = NULL) {
  set.seed(seed)

  total_bins <- sum(BINS_PER_CHR)
  total_reads <- round(total_bins * reads_per_bin)
  genome_bp <- sum(as.numeric(COMPRESSED_CHR_LENGTHS))

  # Trisomy parameters
  has_trisomy <- !is.null(profile$trisomy_chr)
  if (has_trisomy) {
    removal_prob <- .trisomy_removal_prob(profile$fetal_fraction,
                                          profile$mosaic)
  }

  # GC bias: if multipliers are provided, oversample then thin.
  # max_mult is the maximum GC multiplier across the genome, used
  # to compute the oversampling factor.
  has_gc <- !is.null(gc_multipliers)
  if (has_gc) {
    max_mult <- max(vapply(gc_multipliers, max, numeric(1)))
  }

  positions_by_chr <- vector("list", length(CHR_NAMES))
  names(positions_by_chr) <- CHR_NAMES

  for (chr_i in seq_along(CHR_NAMES)) {
    chr_name <- CHR_NAMES[chr_i]
    chr_len <- COMPRESSED_CHR_LENGTHS[chr_name]
    mult <- profile$chr_mult[chr_name]
    if (is.na(mult) || mult == 0) {
      positions_by_chr[[chr_i]] <- integer(0)
      next
    }

    # Expected reads proportional to chr length (in compressed space)
    n_reads <- round(total_reads * (chr_len / genome_bp) * mult)
    if (n_reads == 0L) {
      positions_by_chr[[chr_i]] <- integer(0)
      next
    }

    # Max valid position: chr_len - READ_LEN + 1 (1-based)
    max_pos <- chr_len - READ_LEN + 1L
    if (max_pos < 1L) {
      positions_by_chr[[chr_i]] <- integer(0)
      next
    }

    if (has_gc) {
      # Oversample by max_mult, then rejection-sample per bin
      n_oversample <- min(as.integer(ceiling(n_reads * max_mult)), max_pos)
      positions <- sort(sample.int(max_pos, n_oversample, replace = FALSE))

      # Determine which bin each position falls in (0-based bin index)
      bin_idx <- (positions - 1L) %/% COMPRESSED_BINSIZE + 1L
      gc_m <- gc_multipliers[[chr_name]]
      # Acceptance probability = gc_mult[bin] / max_mult
      accept_prob <- gc_m[bin_idx] / max_mult
      keep_gc <- stats::runif(length(positions)) < accept_prob
      positions <- positions[keep_gc]
    } else {
      if (n_reads > max_pos) n_reads <- max_pos
      positions <- sort(sample.int(max_pos, n_reads, replace = FALSE))
    }

    n_reads <- length(positions)
    if (n_reads == 0L) {
      positions_by_chr[[chr_i]] <- integer(0)
      next
    }

    # Nguyen et al. trisomy simulation: remove reads from non-target
    # chromosomes with probability p = k*f/(1+k*f) so relative coverage
    # of the trisomy chromosome reflects an extra fetal copy.
    if (has_trisomy && chr_name != profile$trisomy_chr) {
      keep <- stats::runif(length(positions)) >= removal_prob
      positions <- positions[keep]
    }
    positions_by_chr[[chr_i]] <- positions
  }

  positions_by_chr
}


# ---- Core generator (exported) -------------------------------------------

#' Generate a synthetic BAM cohort for testing
#'
#' Creates coordinate-sorted, indexed BAM files in a specified directory
#' using compressed chromosome lengths: each 100kb GRCh37 bin is represented
#' as a 100bp region, yielding identical bin-count structure at a fraction
#' of the file size (~435KB per BAM).
#'
#' The cohort contains approximately 75% euploid females, 25% euploid males,
#' and 3 trisomy females (T21, T18, T13). Each sample uses approximately
#' 3 reads per bin (~91k reads total). Results are deterministic (sample *i*
#' uses `set.seed(42 + i)`).
#'
#' Trisomy simulation follows Nguyen et al. 2023
#' (doi:10.1101/2023.11.24.568620): instead of naively multiplying reads
#' on the trisomy chromosome, reads are stochastically removed from
#' non-target chromosomes so that the relative coverage of the trisomy
#' chromosome reflects an extra copy in the fetal genome. The removal
#' probability is `p = k*f / (1 + k*f)` where `k = 0.5` (non-mosaic) or
#' `k = 0.25` (mosaic).
#'
#' @param out_dir Directory to write BAM files into (created if needed).
#' @param n_samples Total number of samples to generate (default 50L).
#'   Must be at least 6 (3 trisomy + 1 male + 2 female minimum).
#' @param fetal_fraction Fetal DNA fraction for trisomy samples (default 0.10).
#'   Higher values produce stronger trisomy signal. Typical NIPT range: 0.05-0.20.
#' @param reads_per_bin Average reads per bin for euploid autosomes
#'   (default 3). Higher values produce larger BAMs but better signal-to-noise
#'   for trisomy detection. At `reads_per_bin = 3`, BAMs are ~435KB;
#'   at 10, ~1.3MB; at 30, ~4.4MB.
#' @param gc_bias_strength Numeric scaling factor for simulated GC bias
#'   amplitude (default 0.3). Each sample gets a random quadratic GC bias
#'   curve; higher values produce more between-sample variance. Set to 0
#'   for uniform coverage (no GC effect).
#' @param mosaic Logical; if `TRUE`, simulate 50% placental mosaicism (halved
#'   trisomy signal strength). Default `FALSE`.
#' @param verbose Logical; emit progress messages via [message()].
#'
#' @return A data frame (manifest) with columns: `sample_id`, `sex`,
#'   `trisomy`, `fetal_fraction`, `n_reads`, `bam_file`.
#'
#' @details
#' The generated BAMs are NIPTeR-compatible (no unmapped reads, unique
#' positions per chromosome) and can be used with both the NIPTeR binning
#' layer ([nipter_bin_bam()]) and the WisecondorX native pipeline
#' ([rwisecondorx_newref()], [rwisecondorx_predict()]) when binned with
#' `binsize = COMPRESSED_BINSIZE` (100).
#'
#' These BAMs are structurally valid and fast to generate, but the underlying
#' coverage model is still synthetic. Use real negative cohorts for
#' calibration, threshold selection, and any real conformance work.
#'
#' The manifest is also written to `manifest.tsv` in `out_dir`.
#'
#' @examples
#' \dontrun{
#' manifest <- generate_cohort(tempdir())
#' head(manifest)
#' }
#'
#' @export
generate_cohort <- function(out_dir, n_samples = 50L, fetal_fraction = 0.10,
                            reads_per_bin = READS_PER_BIN_DEFAULT,
                            gc_bias_strength = 0.3,
                            mosaic = FALSE, verbose = TRUE) {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  Rduckhts::rduckhts_load(con)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  profiles <- .build_cohort_profiles(n_samples, fetal_fraction, mosaic)

  # Generate GC landscape (same for all samples) and per-sample bias
  gc_landscape <- if (gc_bias_strength > 0) .gc_landscape() else NULL

  manifest <- data.frame(
    sample_id      = character(0),
    sex            = character(0),
    trisomy        = character(0),
    fetal_fraction = numeric(0),
    n_reads        = integer(0),
    bam_file       = character(0),
    stringsAsFactors = FALSE
  )

  for (idx in seq_along(profiles)) {
    prof <- profiles[[idx]]
    seed <- 42L + idx
    bam_path <- file.path(out_dir, paste0(prof$id, ".bam"))

    if (verbose) {
      message(sprintf("  [%d/%d] %s (sex=%s)",
                      idx, length(profiles), prof$id, prof$sex))
    }

    # Compute per-sample GC bias multipliers
    gc_mults <- if (!is.null(gc_landscape)) {
      .gc_bias_multipliers(gc_landscape, sample_seed = 10000L + idx,
                           gc_bias_strength = gc_bias_strength)
    } else {
      NULL
    }

    positions_by_chr <- .sample_positions(
      prof,
      seed,
      reads_per_bin = reads_per_bin,
      gc_multipliers = gc_mults
    )
    n_reads <- sum(lengths(positions_by_chr))

    write_synthetic_bam_cpp(
      path = bam_path,
      chr_names = CHR_NAMES,
      chr_lengths = as.integer(COMPRESSED_CHR_LENGTHS),
      positions_by_chr = positions_by_chr,
      read_len = READ_LEN,
      mapq = 60L,
      threads = 1L
    )
    Rduckhts::rduckhts_bam_index(con, bam_path, threads = 1L)

    trisomy <- "none"
    if (grepl("trisomy_21", prof$id)) trisomy <- "T21"
    if (grepl("trisomy_18", prof$id)) trisomy <- "T18"
    if (grepl("trisomy_13", prof$id)) trisomy <- "T13"

    ff <- if (!is.null(prof$fetal_fraction)) prof$fetal_fraction else NA_real_

    manifest <- rbind(manifest, data.frame(
      sample_id      = prof$id,
      sex            = prof$sex,
      trisomy        = trisomy,
      fetal_fraction = ff,
      n_reads        = n_reads,
      bam_file       = basename(bam_path),
      stringsAsFactors = FALSE
    ))
  }

  # Write manifest
  manifest_path <- file.path(out_dir, "manifest.tsv")
  utils::write.table(manifest, manifest_path, sep = "\t", row.names = FALSE,
                     quote = FALSE)

  if (verbose) {
    total_size <- sum(file.size(
      file.path(out_dir, list.files(out_dir, pattern = "\\.bam$"))
    ))
    message(sprintf("Done. %d BAMs (%.1f MB). Manifest: %s",
                    nrow(manifest), total_size / 1024 / 1024, manifest_path))
  }

  manifest
}
