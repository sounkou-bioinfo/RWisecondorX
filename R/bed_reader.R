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

  # Detect column count by probing column5 (0-indexed).
  # 5-col BED: column0..column4 → column5 does not exist (DuckDB error).
  # 9-col BED: column0..column8 → column5 will have count_rev.
  probe_sql <- sprintf(
    "SELECT column5 FROM read_tabix('%s') LIMIT 1",
    bed_sql
  )
  is_separated <- tryCatch({
    probe <- DBI::dbGetQuery(con, probe_sql)
    nrow(probe) > 0L &&
      !is.na(probe$column5[1L]) &&
      nzchar(probe$column5[1L])
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
  chrom <- rows$chrom
  chrom <- sub("^[Cc][Hh][Rr]", "", chrom)
  chrom[chrom == "X" | chrom == "x"] <- "23"
  chrom[chrom == "Y" | chrom == "y"] <- "24"

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


# Read a 5-column CombinedStrands BED into a NIPTeRSample.
# Columns: chrom, start, end, count, corrected_count
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

  if (is.null(binsize)) {
    binsize <- as.integer(rows$end_pos[1L] - rows$start_pos[1L])
  }
  binsize <- as.integer(binsize)

  auto_keys <- as.character(1:22)
  sex_labels <- c("X", "Y")
  sex_keys <- c("23", "24")

  # Normalize chromosome names to keys
  chrom <- sub("^[Cc][Hh][Rr]", "", rows$chrom)
  chrom[chrom == "X" | chrom == "x"] <- "23"
  chrom[chrom == "Y" | chrom == "y"] <- "24"

  # Determine n_bins: widest chromosome
  bin_idx <- as.integer(rows$start_pos / binsize)
  n_bins <- 0L
  for (key in c(auto_keys, sex_keys)) {
    sel <- which(chrom == key)
    if (length(sel) > 0L) {
      n_bins <- max(n_bins, max(bin_idx[sel]) + 1L)
    }
  }

  col_names <- as.character(seq_len(n_bins))

  # Build autosomal matrix (22 x n_bins)
  auto_mat <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(auto_keys, col_names))
  for (i in seq_along(auto_keys)) {
    sel <- which(chrom == auto_keys[i])
    if (length(sel) > 0L) {
      auto_mat[i, bin_idx[sel] + 1L] <- as.integer(rows$count[sel])
    }
  }

  # Build sex matrix (2 x n_bins)
  sex_mat <- matrix(0L, nrow = 2L, ncol = n_bins,
                    dimnames = list(sex_labels, col_names))
  for (i in seq_along(sex_keys)) {
    sel <- which(chrom == sex_keys[i])
    if (length(sel) > 0L) {
      sex_mat[i, bin_idx[sel] + 1L] <- as.integer(rows$count[sel])
    }
  }

  sample_obj <- structure(
    list(
      autosomal_chromosome_reads  = list(auto_mat),
      sex_chromosome_reads        = list(sex_mat),
      correction_status_autosomal = "Uncorrected",
      correction_status_sex       = "Uncorrected",
      sample_name                 = name
    ),
    class = c("NIPTeRSample", "CombinedStrands")
  )

  # Apply corrected counts if present (not all NA)
  if (!all(is.na(rows$corrected_count))) {
    corr_auto <- matrix(0, nrow = 22L, ncol = n_bins,
                        dimnames = list(auto_keys, col_names))
    corr_sex  <- matrix(0, nrow = 2L, ncol = n_bins,
                        dimnames = list(sex_labels, col_names))
    for (i in seq_along(auto_keys)) {
      sel <- which(chrom == auto_keys[i])
      if (length(sel) > 0L) {
        corr_auto[i, bin_idx[sel] + 1L] <- rows$corrected_count[sel]
      }
    }
    for (i in seq_along(sex_keys)) {
      sel <- which(chrom == sex_keys[i])
      if (length(sel) > 0L) {
        corr_sex[i, bin_idx[sel] + 1L] <- rows$corrected_count[sel]
      }
    }
    sample_obj$autosomal_chromosome_reads <- list(corr_auto)
    sample_obj$sex_chromosome_reads       <- list(corr_sex)
    sample_obj$correction_status_autosomal <- "GC Corrected"
    sample_obj$correction_status_sex       <- "GC Corrected"
  }

  sample_obj
}


# Read a 9-column SeparatedStrands BED into a NIPTeRSample.
# Columns: chrom, start, end, count, count_fwd, count_rev,
#           corrected_count, corrected_fwd, corrected_rev
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

  if (is.null(binsize)) {
    binsize <- as.integer(rows$end_pos[1L] - rows$start_pos[1L])
  }
  binsize <- as.integer(binsize)

  auto_keys <- as.character(1:22)
  sex_keys <- c("23", "24")

  # Normalize chromosome names to keys
  chrom <- sub("^[Cc][Hh][Rr]", "", rows$chrom)
  chrom[chrom == "X" | chrom == "x"] <- "23"
  chrom[chrom == "Y" | chrom == "y"] <- "24"

  bin_idx <- as.integer(rows$start_pos / binsize)

  # Determine n_bins
  n_bins <- 0L
  for (key in c(auto_keys, sex_keys)) {
    sel <- which(chrom == key)
    if (length(sel) > 0L) {
      n_bins <- max(n_bins, max(bin_idx[sel]) + 1L)
    }
  }

  col_names <- as.character(seq_len(n_bins))

  # Build four raw-count matrices: fwd_auto, rev_auto, fwd_sex, rev_sex
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

  sample_obj <- structure(
    list(
      autosomal_chromosome_reads  = list(fwd_auto, rev_auto),
      sex_chromosome_reads        = list(fwd_sex, rev_sex),
      correction_status_autosomal = "Uncorrected",
      correction_status_sex       = "Uncorrected",
      sample_name                 = name
    ),
    class = c("NIPTeRSample", "SeparatedStrands")
  )

  # Apply per-strand corrected counts if present
  has_corr <- !all(is.na(rows$corrected_fwd)) &&
    !all(is.na(rows$corrected_rev))

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

    sample_obj$autosomal_chromosome_reads <- list(corr_fwd_auto, corr_rev_auto)
    sample_obj$sex_chromosome_reads       <- list(corr_fwd_sex, corr_rev_sex)
    sample_obj$correction_status_autosomal <- "GC Corrected"
    sample_obj$correction_status_sex       <- "GC Corrected"
  }

  sample_obj
}
