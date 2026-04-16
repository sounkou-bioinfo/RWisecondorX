# bed_reader.R — Read bin counts from bgzipped BED files
#
# Provides bed_to_sample() and bed_to_nipter_sample() for loading pre-computed
# bin counts back into the in-memory formats expected by the WisecondorX and
# NIPTeR analysis pipelines.
#
# This closes the round-trip: bam_convert_bed() / nipter_bin_bam_bed() write
# bgzipped BED files; these functions read them back. All HTS I/O goes through
# Rduckhts's read_tabix() DuckDB table function, which returns all columns as
# VARCHAR — no BED-schema type coercion issues with doubles in columns 7+.


#' Read a WisecondorX-format BED file into a sample list
#'
#' Reads a 4-column bgzipped BED file (as written by [bam_convert_bed()]) and
#' returns a named list of integer vectors suitable for [rwisecondorx_newref()],
#' [rwisecondorx_predict()], [scale_sample()], or [bam_convert_npz()].
#'
#' The BED file must have 4 tab-delimited columns: `chrom`, `start`, `end`,
#' `count` (no header). Coordinates are 0-based half-open intervals. The bin
#' size is inferred from the first row (`end - start`) unless explicitly
#' provided.
#'
#' @param bed Path to a bgzipped (or plain) BED file with a `.tbi` index.
#' @param binsize Optional integer; bin size in base pairs. If `NULL` (default),
#'   inferred from the first row of the BED file.
#' @param con Optional open DBI connection with duckhts already loaded.
#'
#' @return A named list with one integer vector per chromosome key
#'   (`"1"`--`"22"`, `"23"` for X, `"24"` for Y). Each vector has length
#'   `max_bin + 1` for that chromosome. Chromosomes absent from the BED file
#'   are `NULL`. This is the same format returned by [bam_convert()].
#'
#' @seealso [bam_convert_bed()], [bam_convert()], [rwisecondorx_newref()],
#'   [rwisecondorx_predict()], [bed_to_nipter_sample()]
#'
#' @examples
#' \dontrun{
#' # Write bin counts to BED, then read them back
#' bam_convert_bed("sample.bam", "sample.bed.gz", binsize = 5000L)
#' bins <- bed_to_sample("sample.bed.gz")
#'
#' # Use directly with the native WisecondorX pipeline
#' samples <- lapply(bed_files, bed_to_sample)
#' ref <- rwisecondorx_newref(samples, binsize = 100000L, nipt = TRUE)
#' }
#'
#' @export
bed_to_sample <- function(bed, binsize = NULL, con = NULL) {
  stopifnot(is.character(bed), length(bed) == 1L, nzchar(bed))
  stopifnot(file.exists(bed))

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  bed_path <- normalizePath(bed, mustWork = TRUE)
  bed_sql <- gsub("'", "''", bed_path)

  # read_tabix() returns all columns as VARCHAR (column0, column1, ...).
  # 4-column BED: chrom, start, end, count
  sql <- sprintf("
    SELECT
      column0                    AS chrom,
      CAST(column1 AS INTEGER)   AS start_pos,
      CAST(column2 AS INTEGER)   AS end_pos,
      CAST(column3 AS INTEGER)   AS count
    FROM read_tabix('%s')
  ", bed_sql)

  rows <- DBI::dbGetQuery(con, sql)

  if (nrow(rows) == 0L) {
    return(stats::setNames(vector("list", 24L), as.character(1:24)))
  }

  # Infer binsize from first row if not provided
  if (is.null(binsize)) {
    binsize <- as.integer(rows$end_pos[1L] - rows$start_pos[1L])
  }
  binsize <- as.integer(binsize)
  stopifnot(binsize > 0L)

  .bed_rows_to_bins(rows, binsize)
}


#' Read a NIPTeR-format BED file into a NIPTeRSample
#'
#' Reads a 5-column or 9-column bgzipped BED file (as written by
#' [nipter_bin_bam_bed()]) and returns a `NIPTeRSample` object suitable for
#' all NIPTeR statistical functions: [nipter_gc_correct()],
#' [nipter_chi_correct()], [nipter_z_score()], [nipter_ncv_score()],
#' [nipter_regression()], and [nipter_predict_sex()].
#'
#' A 5-column BED (`chrom`, `start`, `end`, `count`, `corrected_count`)
#' produces a `CombinedStrands` sample. A 9-column BED (`chrom`, `start`,
#' `end`, `count`, `count_fwd`, `count_rev`, `corrected_count`,
#' `corrected_fwd`, `corrected_rev`) produces a `SeparatedStrands` sample
#' with independent forward/reverse count matrices. The number of columns is
#' detected automatically.
#'
#' When the corrected columns contain non-NA values (i.e. the BED was written
#' with a GC-corrected `corrected` argument), the returned sample's count
#' matrices are replaced with the corrected values and the correction status is
#' set to `"GC Corrected"`. For `SeparatedStrands`, the per-strand corrected
#' values (`corrected_fwd`, `corrected_rev`) are used to populate the
#' forward and reverse matrices independently.
#'
#' @param bed Path to a bgzipped (or plain) BED file with a `.tbi` index.
#' @param name Optional sample name. If `NULL` (default), derived from the BED
#'   file basename (e.g. `"sample"` from `"sample.bed.gz"`).
#' @param binsize Optional integer; bin size in base pairs. If `NULL` (default),
#'   inferred from the first row of the BED file.
#' @param con Optional open DBI connection with duckhts already loaded.
#'
#' @return An object of class `c("NIPTeRSample", <strand_type>)`:
#'
#'   **`CombinedStrands`** (from 5-column BED): same structure as
#'   [nipter_bin_bam()] with `separate_strands = FALSE`.
#'
#'   **`SeparatedStrands`** (from 9-column BED): same structure as
#'   [nipter_bin_bam()] with `separate_strands = TRUE`.
#'
#' @seealso [nipter_bin_bam_bed()], [nipter_bin_bam()], [bed_to_sample()],
#'   [nipter_gc_correct()], [nipter_z_score()]
#'
#' @examples
#' \dontrun{
#' # Write NIPTeR-style bin counts to BED, then read back
#' nipter_bin_bam_bed("sample.bam", "sample.bed.gz")
#' sample <- bed_to_nipter_sample("sample.bed.gz")
#'
#' # Use with NIPTeR scoring
#' samples <- lapply(bed_files, bed_to_nipter_sample)
#' cg <- nipter_as_control_group(samples)
#' z21 <- nipter_z_score(samples[[1]], cg, chromo_focus = 21)
#' }
#'
#' @export
bed_to_nipter_sample <- function(bed, name = NULL, binsize = NULL, con = NULL) {
  stopifnot(is.character(bed), length(bed) == 1L, nzchar(bed))
  stopifnot(file.exists(bed))

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  bed_path <- normalizePath(bed, mustWork = TRUE)
  bed_sql <- gsub("'", "''", bed_path)

  if (is.null(name)) {
    name <- sub("\\.bed(\\.gz)?$", "", basename(bed), ignore.case = TRUE)
  }

  # Detect column count by probing for column8 (0-indexed).
  # 9-col BED (SeparatedStrands): column0..column8 all present.
  # 5-col BED (CombinedStrands): column0..column4 only; column8 does not exist.
  # Checking column8 (not column5) avoids false positives from files with 6-8
  # columns that this package does not produce.
  probe_sql <- sprintf(
    "SELECT column8 FROM read_tabix('%s') LIMIT 1",
    bed_sql
  )
  is_separated <- tryCatch({
    probe <- DBI::dbGetQuery(con, probe_sql)
    nrow(probe) > 0L &&
      !is.na(probe$column8[1L]) &&
      nzchar(probe$column8[1L])
  }, error = function(e) FALSE)

  if (is_separated) {
    .bed_to_nipter_separated(con, bed_sql, name, binsize)
  } else {
    .bed_to_nipter_combined(con, bed_sql, name, binsize)
  }
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Convert BED rows to the named-list-of-integer-vectors format used by
# bam_convert() and expected by rwisecondorx_newref/predict.
.bed_rows_to_bins <- function(rows, binsize) {
  chr_keys <- as.character(1:24)
  result <- stats::setNames(vector("list", 24L), chr_keys)

  # Normalize chromosome names
  chrom <- .normalize_chr_name(rows$chrom)

  # Compute bin index from start position
  bin <- as.integer(rows$start_pos / binsize)

  for (chr in chr_keys) {
    idx <- which(chrom == chr)
    if (length(idx) == 0L) next

    chr_bin <- bin[idx]
    chr_count <- rows$count[idx]
    max_bin <- max(chr_bin)
    counts <- integer(max_bin + 1L)
    counts[chr_bin + 1L] <- as.integer(chr_count)
    result[[chr]] <- counts
  }

  result
}


# Core: build a CombinedStrandsSample from a pre-fetched data frame.
# rows must have columns: chrom, start_pos, end_pos, count, corrected_count.
# name and binsize are required; binsize=NULL infers from first row.
.rows_to_nipter_combined <- function(rows, name, binsize) {
  if (nrow(rows) == 0L)
    stop(sprintf("BED for sample '%s' has no rows.", name), call. = FALSE)
  if (is.null(binsize))
    binsize <- as.integer(rows$end_pos[1L] - rows$start_pos[1L])
  binsize <- as.integer(binsize)

  auto_keys  <- as.character(1:22)
  sex_labels <- c("X", "Y")
  sex_keys   <- c("23", "24")

  chrom   <- .normalize_chr_name(rows$chrom)
  bin_idx <- as.integer(rows$start_pos / binsize)

  n_bins <- 0L
  for (key in c(auto_keys, sex_keys)) {
    sel <- which(chrom == key)
    if (length(sel) > 0L)
      n_bins <- max(n_bins, max(bin_idx[sel]) + 1L)
  }

  col_names <- as.character(seq_len(n_bins))

  auto_mat <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(auto_keys, col_names))
  for (i in seq_along(auto_keys)) {
    sel <- which(chrom == auto_keys[i])
    if (length(sel) > 0L)
      auto_mat[i, bin_idx[sel] + 1L] <- as.integer(rows$count[sel])
  }

  sex_mat <- matrix(0L, nrow = 2L, ncol = n_bins,
                    dimnames = list(sex_labels, col_names))
  for (i in seq_along(sex_keys)) {
    sel <- which(chrom == sex_keys[i])
    if (length(sel) > 0L)
      sex_mat[i, bin_idx[sel] + 1L] <- as.integer(rows$count[sel])
  }

  # Determine correction status and matrices
  has_corr <- !is.null(rows$corrected_count) && !all(is.na(rows$corrected_count))

  if (has_corr) {
    corr_auto <- matrix(0, nrow = 22L, ncol = n_bins,
                        dimnames = list(auto_keys, col_names))
    corr_sex  <- matrix(0, nrow = 2L, ncol = n_bins,
                        dimnames = list(sex_labels, col_names))
    for (i in seq_along(auto_keys)) {
      sel <- which(chrom == auto_keys[i])
      if (length(sel) > 0L)
        corr_auto[i, bin_idx[sel] + 1L] <- rows$corrected_count[sel]
    }
    for (i in seq_along(sex_keys)) {
      sel <- which(chrom == sex_keys[i])
      if (length(sel) > 0L)
        corr_sex[i, bin_idx[sel] + 1L] <- rows$corrected_count[sel]
    }
    storage.mode(corr_auto) <- "double"
    storage.mode(corr_sex)  <- "double"
    CombinedStrandsSample(
      sample_name = name,
      binsize     = binsize,
      auto_matrix = corr_auto,
      sex_matrix_ = corr_sex,
      correction  = NIPTCorrectionRecord(
        autosomal = "GC Corrected",
        sex       = "GC Corrected"
      )
    )
  } else {
    CombinedStrandsSample(
      sample_name = name,
      binsize     = binsize,
      auto_matrix = auto_mat,
      sex_matrix_ = sex_mat
    )
  }
}

# Wrapper: fetch rows then delegate to .rows_to_nipter_combined.
.bed_to_nipter_combined <- function(con, bed_sql, name, binsize) {
  sql <- sprintf("
    SELECT
      column0                          AS chrom,
      CAST(column1 AS INTEGER)         AS start_pos,
      CAST(column2 AS INTEGER)         AS end_pos,
      CAST(column3 AS INTEGER)         AS count,
      TRY_CAST(column4 AS DOUBLE)      AS corrected_count
    FROM read_tabix('%s')
  ", bed_sql)
  rows <- DBI::dbGetQuery(con, sql)
  if (nrow(rows) == 0L) {
    stop(sprintf("BED file '%s' contains no rows; cannot construct NIPTeRSample.",
                 bed_sql), call. = FALSE)
  }
  .rows_to_nipter_combined(rows, name, binsize)
}


# Core: build a SeparatedStrandsSample from a pre-fetched data frame.
# rows must have columns: chrom, start_pos, end_pos, count_fwd, count_rev,
# corrected_fwd, corrected_rev (all typed).
.rows_to_nipter_sep <- function(rows, name, binsize) {
  if (nrow(rows) == 0L)
    stop(sprintf("BED for sample '%s' has no rows.", name), call. = FALSE)
  if (is.null(binsize))
    binsize <- as.integer(rows$end_pos[1L] - rows$start_pos[1L])
  binsize <- as.integer(binsize)

  auto_keys <- as.character(1:22)
  sex_keys  <- c("23", "24")

  chrom   <- .normalize_chr_name(rows$chrom)
  bin_idx <- as.integer(rows$start_pos / binsize)

  n_bins <- 0L
  for (key in c(auto_keys, sex_keys)) {
    sel <- which(chrom == key)
    if (length(sel) > 0L)
      n_bins <- max(n_bins, max(bin_idx[sel]) + 1L)
  }

  col_names <- as.character(seq_len(n_bins))

  fwd_auto <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(paste0(auto_keys, "F"), col_names))
  rev_auto <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(paste0(auto_keys, "R"), col_names))
  fwd_sex  <- matrix(0L, nrow = 2L, ncol = n_bins,
                     dimnames = list(c("XF", "YF"), col_names))
  rev_sex  <- matrix(0L, nrow = 2L, ncol = n_bins,
                     dimnames = list(c("XR", "YR"), col_names))

  for (i in seq_along(auto_keys)) {
    sel <- which(chrom == auto_keys[i])
    if (length(sel) > 0L) {
      fwd_auto[i, bin_idx[sel] + 1L] <- as.integer(rows$count_fwd[sel])
      rev_auto[i, bin_idx[sel] + 1L] <- as.integer(rows$count_rev[sel])
    }
  }
  for (i in seq_along(sex_keys)) {
    sel <- which(chrom == sex_keys[i])
    if (length(sel) > 0L) {
      fwd_sex[i, bin_idx[sel] + 1L] <- as.integer(rows$count_fwd[sel])
      rev_sex[i, bin_idx[sel] + 1L] <- as.integer(rows$count_rev[sel])
    }
  }

  has_corr <- !is.null(rows$corrected_fwd) && !all(is.na(rows$corrected_fwd)) &&
              !is.null(rows$corrected_rev) && !all(is.na(rows$corrected_rev))

  if (has_corr) {
    corr_fwd_auto <- matrix(0, nrow = 22L, ncol = n_bins,
                            dimnames = list(paste0(auto_keys, "F"), col_names))
    corr_rev_auto <- matrix(0, nrow = 22L, ncol = n_bins,
                            dimnames = list(paste0(auto_keys, "R"), col_names))
    corr_fwd_sex  <- matrix(0, nrow = 2L, ncol = n_bins,
                            dimnames = list(c("XF", "YF"), col_names))
    corr_rev_sex  <- matrix(0, nrow = 2L, ncol = n_bins,
                            dimnames = list(c("XR", "YR"), col_names))
    for (i in seq_along(auto_keys)) {
      sel <- which(chrom == auto_keys[i])
      if (length(sel) > 0L) {
        corr_fwd_auto[i, bin_idx[sel] + 1L] <- rows$corrected_fwd[sel]
        corr_rev_auto[i, bin_idx[sel] + 1L] <- rows$corrected_rev[sel]
      }
    }
    for (i in seq_along(sex_keys)) {
      sel <- which(chrom == sex_keys[i])
      if (length(sel) > 0L) {
        corr_fwd_sex[i, bin_idx[sel] + 1L] <- rows$corrected_fwd[sel]
        corr_rev_sex[i, bin_idx[sel] + 1L] <- rows$corrected_rev[sel]
      }
    }
    storage.mode(corr_fwd_auto) <- "double"; storage.mode(corr_rev_auto) <- "double"
    storage.mode(corr_fwd_sex)  <- "double"; storage.mode(corr_rev_sex)  <- "double"
    SeparatedStrandsSample(
      sample_name = name,
      binsize     = binsize,
      auto_fwd    = corr_fwd_auto,
      auto_rev    = corr_rev_auto,
      sex_fwd     = corr_fwd_sex,
      sex_rev     = corr_rev_sex,
      correction  = NIPTCorrectionRecord(
        autosomal = "GC Corrected",
        sex       = "GC Corrected"
      )
    )
  } else {
    SeparatedStrandsSample(
      sample_name = name,
      binsize     = binsize,
      auto_fwd    = fwd_auto,
      auto_rev    = rev_auto,
      sex_fwd     = fwd_sex,
      sex_rev     = rev_sex
    )
  }
}

# Wrapper: fetch rows then delegate to .rows_to_nipter_sep.
.bed_to_nipter_separated <- function(con, bed_sql, name, binsize) {
  sql <- sprintf("
    SELECT
      column0                          AS chrom,
      CAST(column1 AS INTEGER)         AS start_pos,
      CAST(column2 AS INTEGER)         AS end_pos,
      CAST(column3 AS INTEGER)         AS count,
      CAST(column4 AS INTEGER)         AS count_fwd,
      CAST(column5 AS INTEGER)         AS count_rev,
      TRY_CAST(column6 AS DOUBLE)      AS corrected_count,
      TRY_CAST(column7 AS DOUBLE)      AS corrected_fwd,
      TRY_CAST(column8 AS DOUBLE)      AS corrected_rev
    FROM read_tabix('%s')
  ", bed_sql)
  rows <- DBI::dbGetQuery(con, sql)
  if (nrow(rows) == 0L) {
    stop(sprintf("BED file '%s' contains no rows; cannot construct NIPTeRSample.",
                 bed_sql), call. = FALSE)
  }
  .rows_to_nipter_sep(rows, name, binsize)
}
