#' Bin a BAM/CRAM file — NIPTeR style
#'
#' Replicates NIPTeR's `bin_bam_sample()`: counts reads in fixed-width bins
#' across chromosomes 1-22, X, and Y. The original NIPTeR counts all mapped
#' reads with no MAPQ filter; real-world NIPT pipelines typically pre-filter
#' with `samtools view --min-MQ 40 -F 1024` before binning. Both modes are
#' supported through the `mapq`, `exclude_flags`, and `require_flags`
#' parameters.
#'
#' The result is reshaped into a `NIPTeRSample` object whose structure
#' parallels NIPTeR's `NIPTSample`: autosomal reads as a chromosome-by-bin
#' matrix and sex chromosome reads as a separate two-row matrix.
#'
#' When `separate_strands = TRUE`, forward (`+`) and reverse (`-`) reads are
#' counted independently, producing two matrices per chromosome set (class
#' `"SeparatedStrands"`). This doubles the predictor pool for
#' [nipter_regression()] — see NIPTeR documentation for details.
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param binsize Bin size in base pairs. Default 50000 (NIPTeR's fixed bin
#'   size). Must match the binsize used when building any control group.
#' @param mapq Minimum mapping quality. Default `0L` (NIPTeR original: all
#'   mapped reads). Set to `40L` to match common NIPT pipeline pre-filtering
#'   (`samtools view --min-MQ 40`).
#' @param require_flags Integer bitmask; only reads with all bits set are kept
#'   (samtools `-f`). Default `0L` (no requirement).
#' @param exclude_flags Integer bitmask; reads with any bit set are dropped
#'   (samtools `-F`). Default `0L`. Set to `1024L` (`0x400`) to exclude
#'   reads marked as duplicates (`samtools view -F 1024`).
#' @param rmdup Duplicate removal strategy. `"none"` (default, NIPTeR standard)
#'   or `"flag"` (equivalent to `exclude_flags = 1024L`; for BAMs already
#'   processed by Picard / sambamba). NIPTeR does not perform streaming
#'   deduplication.
#' @param separate_strands Logical; when `TRUE`, produces a
#'   `SeparatedStrands` object with independent forward/reverse count matrices.
#'   Default `FALSE` (`CombinedStrands`).
#' @param con Optional open DBI connection with duckhts already loaded.
#' @param reference Optional FASTA reference path for CRAM inputs.
#'
#' @return An object of class `c("NIPTeRSample", <strand_type>)`:
#'
#'   **`CombinedStrands`** (default): `autosomal_chromosome_reads` is a
#'   list of one 22-row integer matrix (rows `"1"`–`"22"`);
#'   `sex_chromosome_reads` is a list of one 2-row matrix (rows `"X"`, `"Y"`).
#'
#'   **`SeparatedStrands`** (`separate_strands = TRUE`):
#'   `autosomal_chromosome_reads` is a list of two matrices — element 1 is
#'   forward (rows `"1F"`–`"22F"`), element 2 is reverse (rows `"1R"`–`"22R"`);
#'   `sex_chromosome_reads` similarly contains forward (`"XF"`, `"YF"`) and
#'   reverse (`"XR"`, `"YR"`) matrices.
#'
#' @seealso [bam_convert()], [nipter_bin_bam_bed()], [nipter_gc_correct()],
#'   [nipter_regression()]
#'
#' @examples
#' \dontrun{
#' # NIPTeR original defaults: all mapped reads, no dedup
#' sample <- nipter_bin_bam("sample.bam", binsize = 50000L)
#'
#' # Common NIPT pipeline: MAPQ >= 40, exclude duplicate-flagged reads
#' sample <- nipter_bin_bam("sample.dm.bam", binsize = 50000L,
#'                          mapq = 40L, exclude_flags = 1024L)
#'
#' # SeparatedStrands for regression with doubled predictor pool
#' sample_ss <- nipter_bin_bam("sample.bam", separate_strands = TRUE)
#'
#' sample$autosomal_chromosome_reads[[1]]["21", ]   # chr21 bin counts
#' }
#'
#' @export
nipter_bin_bam <- function(bam,
                           binsize          = 50000L,
                           mapq             = 0L,
                           require_flags    = 0L,
                           exclude_flags    = 0L,
                           rmdup            = c("none", "flag"),
                           separate_strands = FALSE,
                           con              = NULL,
                           reference        = NULL) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.logical(separate_strands), length(separate_strands) == 1L)
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  stopifnot(is.numeric(binsize),       length(binsize)       == 1L, binsize       >= 1L)
  stopifnot(is.numeric(mapq),          length(mapq)          == 1L, mapq          >= 0L)
  stopifnot(is.numeric(require_flags), length(require_flags) == 1L, require_flags >= 0L)
  stopifnot(is.numeric(exclude_flags), length(exclude_flags) == 1L, exclude_flags >= 0L)
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }
  binsize       <- as.integer(binsize)
  mapq          <- as.integer(mapq)
  require_flags <- as.integer(require_flags)
  exclude_flags <- as.integer(exclude_flags)

  bins <- bam_convert(bam,
                      binsize          = binsize,
                      mapq             = mapq,
                      require_flags    = require_flags,
                      exclude_flags    = exclude_flags,
                      rmdup            = rmdup,
                      separate_strands = separate_strands,
                      con              = con,
                      reference        = reference)

  name <- sub("\\.cram$|\\.bam$", "", basename(bam), ignore.case = TRUE)

  if (isTRUE(separate_strands)) {
    .stranded_bins_to_nipter_sample(bins$fwd, bins$rev, binsize, name)
  } else {
    .bins_to_nipter_sample(bins, binsize, name)
  }
}


#' Write NIPTeR-style bin counts to a bgzipped BED file
#'
#' Bins a BAM/CRAM file with NIPTeR defaults and writes the result to a
#' bgzipped, tabix-indexed BED file.
#'
#' When `separate_strands = FALSE` (default), the output has five columns:
#' `chrom`, `start`, `end`, `count`, `corrected_count`.
#'
#' When `separate_strands = TRUE`, the output has nine columns:
#' `chrom`, `start`, `end`, `count`, `count_fwd`, `count_rev`,
#' `corrected_count`, `corrected_fwd`, `corrected_rev`. `count` is the total
#' (forward + reverse); `corrected_count` is the total of the per-strand
#' corrected values.
#'
#' `corrected_count` (and `corrected_fwd`/`corrected_rev`) is `NA` until a
#' GC-corrected sample is supplied via the `corrected` parameter.
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param bed Path for the output `.bed.gz` file. The tabix index is written
#'   alongside as `paste0(bed, ".tbi")`.
#' @param binsize Bin size in base pairs (default 50000).
#' @param mapq Minimum mapping quality (default `0L`). Set to `40L` to match
#'   common NIPT pre-filtering (`samtools view --min-MQ 40`).
#' @param require_flags Integer bitmask; only reads with all bits set are kept
#'   (samtools `-f`). Default `0L`.
#' @param exclude_flags Integer bitmask; reads with any bit set are dropped
#'   (samtools `-F`). Default `0L`. Set to `1024L` to exclude duplicate-flagged
#'   reads (`samtools view -F 1024`).
#' @param rmdup Duplicate removal strategy: `"none"` (default) or `"flag"`.
#' @param separate_strands Logical; when `TRUE`, includes separate
#'   `count_fwd` and `count_rev` columns for forward- and reverse-strand
#'   counts. Default `FALSE`.
#' @param corrected Optional `NIPTeRSample` already processed by
#'   `nipter_gc_correct()`. When supplied, its corrected counts populate the
#'   `corrected_count` column; otherwise the column is `NA`.
#' @param con Optional open DBI connection with duckhts already loaded.
#' @param reference Optional FASTA reference path for CRAM inputs.
#' @param index Logical; create a tabix index (default `TRUE`).
#'
#' @return `bed` (invisibly).
#'
#' @seealso `nipter_bin_bam()`, `nipter_gc_correct()`, [bam_convert_bed()]
#'
#' @examples
#' \dontrun{
#' nipter_bin_bam_bed("sample.bam", "sample.nipter.bed.gz")
#'
#' # With pre-filtering matching a typical NIPT pipeline
#' nipter_bin_bam_bed("sample.dm.bam", "sample.nipter.bed.gz",
#'                    mapq = 40L, exclude_flags = 1024L)
#'
#' # Strand-separated output (9 columns)
#' nipter_bin_bam_bed("sample.bam", "sample.stranded.bed.gz",
#'                    separate_strands = TRUE)
#' }
#'
#' @export
nipter_bin_bam_bed <- function(bam,
                               bed,
                               binsize          = 50000L,
                               mapq             = 0L,
                               require_flags    = 0L,
                               exclude_flags    = 0L,
                               rmdup            = c("none", "flag"),
                               separate_strands = FALSE,
                               corrected        = NULL,
                               con              = NULL,
                               reference        = NULL,
                               index            = TRUE) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  stopifnot(is.character(bed), length(bed) == 1L, nzchar(bed))
  stopifnot(is.logical(separate_strands), length(separate_strands) == 1L)
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }
  stopifnot(is.numeric(binsize),       length(binsize)       == 1L, binsize       >= 1L)
  stopifnot(is.numeric(mapq),          length(mapq)          == 1L, mapq          >= 0L)
  stopifnot(is.numeric(require_flags), length(require_flags) == 1L, require_flags >= 0L)
  stopifnot(is.numeric(exclude_flags), length(exclude_flags) == 1L, exclude_flags >= 0L)
  stopifnot(is.logical(index), length(index) == 1L)
  binsize       <- as.integer(binsize)
  mapq          <- as.integer(mapq)
  require_flags <- as.integer(require_flags)
  exclude_flags <- as.integer(exclude_flags)

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  sample <- nipter_bin_bam(bam, binsize = binsize, mapq = mapq,
                           require_flags = require_flags,
                           exclude_flags = exclude_flags,
                           rmdup = rmdup, separate_strands = separate_strands,
                           con = con, reference = reference)

  if (isTRUE(separate_strands)) {
    df <- .nipter_bed_separated_strands(sample, binsize, corrected)
  } else {
    df <- .nipter_bed_combined_strands(sample, binsize, corrected)
  }

  tmp <- tempfile(fileext = ".bed")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.table(df, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE)

  Rduckhts::rduckhts_bgzip(con, tmp,
                           output_path = bed,
                           threads     = 1L,
                           keep        = TRUE,
                           overwrite   = TRUE)

  if (isTRUE(index)) {
    Rduckhts::rduckhts_tabix_index(con, bed, preset = "bed", threads = 1L)
  }

  invisible(bed)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Reshape bam_convert() list output into a NIPTeRSample.
# NIPTeRSample stores reads as matrices: rows = chromosomes, cols = bins.
# All chromosomes share the same n_bins (widest chromosome determines width;
# narrower chromosomes are zero-padded on the right).
.bins_to_nipter_sample <- function(bins, binsize, name) {
  auto_keys <- as.character(1:22)
  sex_keys  <- c("23", "24")   # internal: X=23, Y=24

  n_bins <- max(
    vapply(bins[c(auto_keys, sex_keys)], length, integer(1L)),
    na.rm = TRUE
  )
  if (is.infinite(n_bins) || n_bins == 0L) {
    stop("No reads found in chromosomes 1-22/X/Y.", call. = FALSE)
  }

  .pad <- function(v) {
    if (is.null(v)) return(integer(n_bins))
    length(v) <- n_bins   # zero-pads (NA → convert to 0)
    v[is.na(v)] <- 0L
    v
  }

  auto_mat <- do.call(rbind, lapply(bins[auto_keys], .pad))
  rownames(auto_mat) <- auto_keys
  colnames(auto_mat) <- as.character(seq_len(n_bins))

  sex_mat <- do.call(rbind, lapply(bins[sex_keys], .pad))
  rownames(sex_mat) <- c("X", "Y")
  colnames(sex_mat) <- as.character(seq_len(n_bins))

  structure(
    list(
      autosomal_chromosome_reads           = list(auto_mat),
      sex_chromosome_reads                 = list(sex_mat),
      correction_status_autosomal          = "Uncorrected",
      correction_status_sex                = "Uncorrected",
      sample_name                          = name
    ),
    class = c("NIPTeRSample", "CombinedStrands")
  )
}


# SeparatedStrands variant: fwd_bins and rev_bins are each structured like
# bam_convert()'s default return (one integer vector per chromosome key).
# Produces two matrices per chromosome set: [[1]] = forward, [[2]] = reverse.
# Rownames follow NIPTeR convention: "1F".."22F"/"1R".."22R" for autosomes,
# "XF"/"YF" and "XR"/"YR" for sex chromosomes.
.stranded_bins_to_nipter_sample <- function(fwd_bins, rev_bins, binsize, name) {
  auto_keys <- as.character(1:22)
  sex_keys  <- c("23", "24")
  all_keys  <- c(auto_keys, sex_keys)

  n_bins <- max(
    vapply(fwd_bins[all_keys], length, integer(1L)),
    vapply(rev_bins[all_keys], length, integer(1L)),
    na.rm = TRUE
  )
  if (is.infinite(n_bins) || n_bins == 0L) {
    stop("No reads found in chromosomes 1-22/X/Y.", call. = FALSE)
  }

  .pad <- function(v) {
    if (is.null(v)) return(integer(n_bins))
    length(v) <- n_bins
    v[is.na(v)] <- 0L
    v
  }

  col_names <- as.character(seq_len(n_bins))

  # Forward autosomal: rows "1F".."22F"
  fwd_auto <- do.call(rbind, lapply(fwd_bins[auto_keys], .pad))
  rownames(fwd_auto) <- paste0(auto_keys, "F")
  colnames(fwd_auto) <- col_names

  # Reverse autosomal: rows "1R".."22R"
  rev_auto <- do.call(rbind, lapply(rev_bins[auto_keys], .pad))
  rownames(rev_auto) <- paste0(auto_keys, "R")
  colnames(rev_auto) <- col_names

  # Forward sex: rows "XF", "YF"
  fwd_sex <- do.call(rbind, lapply(fwd_bins[sex_keys], .pad))
  rownames(fwd_sex) <- c("XF", "YF")
  colnames(fwd_sex) <- col_names

  # Reverse sex: rows "XR", "YR"
  rev_sex <- do.call(rbind, lapply(rev_bins[sex_keys], .pad))
  rownames(rev_sex) <- c("XR", "YR")
  colnames(rev_sex) <- col_names

  structure(
    list(
      autosomal_chromosome_reads  = list(fwd_auto, rev_auto),
      sex_chromosome_reads        = list(fwd_sex, rev_sex),
      correction_status_autosomal = "Uncorrected",
      correction_status_sex       = "Uncorrected",
      sample_name                 = name
    ),
    class = c("NIPTeRSample", "SeparatedStrands")
  )
}


# CombinedStrands BED data.frame: 5 columns
# chrom, start, end, count, corrected_count
.nipter_bed_combined_strands <- function(sample, binsize, corrected = NULL) {
  raw_auto <- sample$autosomal_chromosome_reads[[1L]]
  raw_sex  <- sample$sex_chromosome_reads[[1L]]
  combined <- rbind(raw_auto, raw_sex)   # 24 rows, n_bins columns

  chr_names <- c(as.character(1:22), "X", "Y")
  n_bins    <- ncol(combined)
  starts    <- as.integer(seq(0L, by = binsize, length.out = n_bins))

  df <- data.frame(
    chrom           = rep(chr_names, each = n_bins),
    start           = rep(starts, times = 24L),
    end             = rep(starts + binsize, times = 24L),
    count           = as.integer(t(combined)),
    corrected_count = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(corrected)) {
    corr_auto <- corrected$autosomal_chromosome_reads[[1L]]
    corr_sex  <- corrected$sex_chromosome_reads[[1L]]
    df$corrected_count <- as.numeric(t(rbind(corr_auto, corr_sex)))
  }

  df
}


# SeparatedStrands BED data.frame: 9 columns
# chrom, start, end, count, count_fwd, count_rev,
# corrected_count, corrected_fwd, corrected_rev
.nipter_bed_separated_strands <- function(sample, binsize, corrected = NULL) {
  fwd_auto <- sample$autosomal_chromosome_reads[[1L]]  # "1F".."22F"
  rev_auto <- sample$autosomal_chromosome_reads[[2L]]   # "1R".."22R"
  fwd_sex  <- sample$sex_chromosome_reads[[1L]]          # "XF","YF"
  rev_sex  <- sample$sex_chromosome_reads[[2L]]          # "XR","YR"

  fwd_combined <- rbind(fwd_auto, fwd_sex)   # 24 rows
  rev_combined <- rbind(rev_auto, rev_sex)    # 24 rows
  total        <- fwd_combined + rev_combined

  chr_names <- c(as.character(1:22), "X", "Y")
  n_bins    <- ncol(total)
  starts    <- as.integer(seq(0L, by = binsize, length.out = n_bins))

  df <- data.frame(
    chrom           = rep(chr_names, each = n_bins),
    start           = rep(starts, times = 24L),
    end             = rep(starts + binsize, times = 24L),
    count           = as.integer(t(total)),
    count_fwd       = as.integer(t(fwd_combined)),
    count_rev       = as.integer(t(rev_combined)),
    corrected_count = NA_real_,
    corrected_fwd   = NA_real_,
    corrected_rev   = NA_real_,
    stringsAsFactors = FALSE
  )

  # Corrected counts from a GC-corrected SeparatedStrands sample
  if (!is.null(corrected) && inherits(corrected, "SeparatedStrands")) {
    corr_fwd <- rbind(corrected$autosomal_chromosome_reads[[1L]],
                      corrected$sex_chromosome_reads[[1L]])
    corr_rev <- rbind(corrected$autosomal_chromosome_reads[[2L]],
                      corrected$sex_chromosome_reads[[2L]])
    df$corrected_count <- as.numeric(t(corr_fwd + corr_rev))
    df$corrected_fwd   <- as.numeric(t(corr_fwd))
    df$corrected_rev   <- as.numeric(t(corr_rev))
  } else if (!is.null(corrected)) {
    corr_auto <- corrected$autosomal_chromosome_reads[[1L]]
    corr_sex  <- corrected$sex_chromosome_reads[[1L]]
    df$corrected_count <- as.numeric(t(rbind(corr_auto, corr_sex)))
  }

  df
}
