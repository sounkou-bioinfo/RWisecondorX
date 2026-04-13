#' Count reads per bin from a BAM or CRAM file
#'
#' Core read-counting engine shared by the WisecondorX and NIPTeR layers.
#' Reads aligned reads from a BAM or CRAM file and returns per-bin read counts
#' for chromosomes 1-22 and the sex chromosomes (X mapped to key "23", Y to "24").
#'
#' Filter parameters let callers replicate the exact behaviour of both tools.
#' WisecondorX defaults: `mapq = 1`, `rmdup = "streaming"`,
#' `filter_improper_pairs = TRUE`.  NIPTeR defaults: `mapq = 0`,
#' `rmdup = "none"`, `filter_improper_pairs = FALSE`.
#'
#' When `rmdup = "streaming"` the WisecondorX larp/larp2 state machine is
#' reproduced exactly: improper pairs are invisible to the dedup logic (they
#' do not update larp2), unpaired reads update larp but not larp2, and larp
#' is never reset between chromosomes.  Bin assignment uses truncating integer
#' division matching Python's `int(pos / binsize)`.
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param binsize Bin size in base pairs. Default 5000 (WisecondorX); use
#'   50000 for NIPTeR-style workflows.
#' @param mapq Minimum mapping quality to retain a read. Default `1L`
#'   (WisecondorX). Set to `0L` to disable MAPQ filtering (NIPTeR).
#' @param rmdup Duplicate removal strategy:
#'   `"streaming"` — WisecondorX larp/larp2 consecutive-position dedup
#'   (default; not meaningful when `filter_improper_pairs = FALSE`);
#'   `"flag"` — exclude reads with SAM flag `0x400` (pre-marked by Picard /
#'   sambamba); `"none"` — no duplicate removal.
#' @param filter_improper_pairs When `TRUE` (default, WisecondorX behaviour)
#'   paired reads that are not properly paired (`FLAG & 0x2 == 0`) are
#'   excluded. Set to `FALSE` to include all mapped reads regardless of pair
#'   status (NIPTeR behaviour).
#' @param con An optional open DBI connection with the duckhts extension already
#'   loaded. If `NULL` (default), a temporary in-memory DuckDB connection is
#'   created for this call.
#' @param reference Optional FASTA reference path for CRAM inputs.
#'
#' @return A named list with one integer vector per chromosome key (`"1"`–`"22"`,
#'   `"23"` for X, `"24"` for Y). Each vector contains per-bin read counts
#'   (bin 0 = positions 0 to binsize-1). Chromosomes absent from the BAM are
#'   `NULL`.
#'
#' @seealso [bam_convert_bed()], [bam_convert_npz()], `nipter_bin_bam()`
#'
#' @examples
#' \dontrun{
#' # WisecondorX defaults
#' bins <- bam_convert("sample.bam", binsize = 5000L, rmdup = "streaming")
#'
#' # NIPTeR defaults
#' bins <- bam_convert("sample.bam", binsize = 50000L, mapq = 0L,
#'                     rmdup = "none", filter_improper_pairs = FALSE)
#' }
#'
#' @export
bam_convert <- function(bam,
                        binsize               = 5000L,
                        mapq                  = 1L,
                        rmdup                 = c("streaming", "flag", "none"),
                        filter_improper_pairs = TRUE,
                        con                   = NULL,
                        reference             = NULL) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  stopifnot(is.numeric(mapq),    length(mapq)    == 1L, mapq    >= 0L)
  stopifnot(is.logical(filter_improper_pairs), length(filter_improper_pairs) == 1L)
  binsize <- as.integer(binsize)
  mapq    <- as.integer(mapq)

  own_con <- is.null(con)
  if (own_con) {
    if (!requireNamespace("Rduckhts", quietly = TRUE)) {
      stop("Rduckhts is required. Install it with: remotes::install_github('RGenomicsETL/duckhts/r/Rduckhts')",
           call. = FALSE)
    }
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  sql <- .convert_sql(bam, binsize, mapq, rmdup,
                      filter_improper_pairs = filter_improper_pairs,
                      reference = reference)
  rows <- DBI::dbGetQuery(con, sql)

  .rows_to_bins(rows, binsize)
}


#' Convert BAM/CRAM to a bgzipped BED bin-count file
#'
#' Runs [bam_convert()] and writes the per-bin read counts to a bgzipped,
#' tab-delimited BED file (four columns: `chrom`, `start`, `end`, `count`),
#' then creates a tabix index (`.tbi`) alongside it via Rduckhts.
#' This is the language-agnostic intermediate format for the RWisecondorX
#' native pipeline; the file can be queried directly with duckhts/DuckDB or
#' any tabix-aware tool.
#'
#' Coordinates are 0-based half-open intervals matching the BED convention
#' (`start = bin_index * binsize`, `end = start + binsize`). Chromosomes are
#' written as `1`–`22`, `X`, `Y` (no `chr` prefix) in numeric order. All
#' bins are written, including those with a count of zero.
#'
#' bgzip and tabix indexing are performed via [Rduckhts::rduckhts_bgzip()] and
#' [Rduckhts::rduckhts_tabix_index()], so no external tools are required.
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param bed Path for the output `.bed.gz` file (created or overwritten).
#'   The tabix index is written to `paste0(bed, ".tbi")`.
#' @param binsize Bin size in base pairs (default 5000).
#' @param rmdup Duplicate-removal strategy passed to [bam_convert()].
#' @param con Optional open DBI connection with duckhts already loaded.
#'   If `NULL` (default), a temporary in-memory DuckDB connection is created.
#' @param reference Optional FASTA reference path for CRAM inputs.
#' @param index Logical; when `TRUE` (default) a tabix index is created
#'   alongside the bgzipped output.
#'
#' @return `bed` (invisibly).
#'
#' @seealso [bam_convert()], [bam_convert_npz()], `wisecondorx_convert()`
#'
#' @examples
#' \dontrun{
#' bam_convert_bed("sample.bam", "sample.bed.gz", binsize = 5000, rmdup = "streaming")
#' }
#'
#' @export
bam_convert_bed <- function(bam,
                            bed,
                            binsize   = 5000L,
                            rmdup     = c("streaming", "none", "flag"),
                            con       = NULL,
                            reference = NULL,
                            index     = TRUE) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  stopifnot(is.character(bed), length(bed) == 1L, nzchar(bed))
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  stopifnot(is.logical(index), length(index) == 1L)
  binsize <- as.integer(binsize)

  # Manage connection lifecycle here so the same con is used for bam_convert,
  # rduckhts_bgzip, and rduckhts_tabix_index.
  own_con <- is.null(con)
  if (own_con) {
    if (!requireNamespace("Rduckhts", quietly = TRUE)) {
      stop(
        "Rduckhts is required. Install it with: ",
        "remotes::install_github('RGenomicsETL/duckhts/r/Rduckhts')",
        call. = FALSE
      )
    }
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  bins <- bam_convert(bam, binsize = binsize, rmdup = rmdup, con = con,
                      reference = reference)

  frames <- vector("list", 24L)
  for (i in seq_len(24L)) {
    key    <- as.character(i)
    counts <- bins[[key]]
    if (is.null(counts)) next
    n_bins <- length(counts)
    starts <- as.integer(seq(0L, by = binsize, length.out = n_bins))
    frames[[i]] <- data.frame(
      chrom = .bin_chr_name(key),
      start = starts,
      end   = starts + binsize,
      count = counts,
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, Filter(Negate(is.null), frames))
  if (is.null(df) || nrow(df) == 0L) {
    stop("bam_convert returned no data for ", bam, call. = FALSE)
  }

  # Write uncompressed BED to a temp file, then bgzip → tabix via Rduckhts.
  tmp <- tempfile(fileext = ".bed")
  on.exit(unlink(tmp), add = TRUE)
  write.table(df, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE)

  Rduckhts::rduckhts_bgzip(con, tmp,
                           output_path = bed,
                           threads     = 1L,
                           keep        = TRUE,
                           overwrite   = TRUE)

  if (isTRUE(index)) {
    Rduckhts::rduckhts_tabix_index(con, bed,
                                   preset  = "bed",
                                   threads = 1L)
  }

  invisible(bed)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.bin_chr_name <- function(key) {
  switch(key, "23" = "X", "24" = "Y", key)
}

.valid_chr_map <- c(
  as.character(1:22),
  "X" = "23", "x" = "23",
  "Y" = "24", "y" = "24"
)

.convert_sql <- function(bam_path, binsize, mapq, rmdup,
                         filter_improper_pairs = TRUE,
                         reference = NULL) {
  read_bam_call <- .read_bam_call(bam_path, reference = reference)

  # Optional proper-pair filter (WisecondorX behaviour).
  # When FALSE all mapped reads are included regardless of pair status (NIPTeR).
  pair_clause <- if (isTRUE(filter_improper_pairs)) {
    "AND NOT ((flag & 1) != 0 AND (flag & 2) = 0)"
  } else {
    ""
  }

  # rmdup strategy
  dedup_where <- switch(rmdup,
    streaming = .dedup_where_streaming(),
    none      = "TRUE",
    flag      = "(flag & 1024) = 0"
  )

  # Bin assignment: (pos - 1) converts duckhts 1-based POS to 0-based.
  # Integer division `//` matches Python's int(pos / binsize) (truncation).
  sprintf("
WITH raw AS (
    SELECT
        rname,
        pos - 1          AS pos0,
        pnext - 1        AS pnext0,
        mapq,
        flag,
        file_offset,
        (flag & 1)       AS is_paired
    FROM %s
    WHERE rname IS NOT NULL
      AND mapq >= %d
      %s
),
with_lag AS (
    SELECT *,
        LAG(pos0) OVER (ORDER BY file_offset) AS prev_pos,
        LAST_VALUE(CASE WHEN is_paired != 0 THEN pnext0 END IGNORE NULLS)
            OVER (ORDER BY file_offset ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING)
            AS prev_pnext
    FROM raw
),
deduped AS (
    SELECT * FROM with_lag
    WHERE %s
)
SELECT
    rname,
    (pos0 // %d)::INTEGER AS bin,
    COUNT(*)::INTEGER     AS n
FROM deduped
GROUP BY rname, bin
ORDER BY
    TRY_CAST(regexp_replace(
        CASE upper(regexp_replace(rname, '^[Cc][Hh][Rr]', ''))
            WHEN 'X' THEN '23'
            WHEN 'Y' THEN '24'
            WHEN 'M' THEN '25'
            ELSE regexp_replace(rname, '^[Cc][Hh][Rr]', '')
        END,
    '[^0-9]', '', 'g') AS INTEGER),
    bin
", read_bam_call, mapq, pair_clause, dedup_where, binsize)
}

.read_bam_call <- function(bam_path, reference = NULL) {
  bam_path <- gsub("'", "''", bam_path)

  if (is.null(reference)) {
    return(sprintf("read_bam('%s')", bam_path))
  }

  reference <- gsub("'", "''", reference)
  sprintf("read_bam('%s', reference := '%s')", bam_path, reference)
}

.dedup_where_streaming <- function() {
  # Exact replication of WisecondorX's larp/larp2 dedup.
  # Drop a read if its (pos, pnext) matches the previous proper read's values,
  # where pnext comparison uses the previous PAIRED read's pnext (larp2).
  paste(
    "NOT (",
    "  prev_pos IS NOT NULL",
    "  AND",
    "  pos0 = prev_pos",
    "  AND (is_paired = 0 OR pnext0 = prev_pnext)",
    ")"
  )
}

.rows_to_bins <- function(rows, binsize) {
  chr_keys <- as.character(1:24)
  result <- stats::setNames(vector("list", 24L), chr_keys)

  if (nrow(rows) == 0L) return(result)

  # Normalize chromosome names: strip "chr" prefix, map X->23, Y->24
  # read_bam() returns RNAME (uppercase); tolower the column lookup to be safe.
  rname <- rows[[grep("^rname$", names(rows), ignore.case = TRUE, value = TRUE)]]
  rname <- sub("^[Cc][Hh][Rr]", "", rname)
  rname[rname == "X" | rname == "x"] <- "23"
  rname[rname == "Y" | rname == "y"] <- "24"

  for (chr in chr_keys) {
    idx <- which(rname == chr)
    if (length(idx) == 0L) next

    chr_rows <- rows[idx, , drop = FALSE]
    max_bin <- max(chr_rows$bin)
    counts <- integer(max_bin + 1L)
    counts[chr_rows$bin + 1L] <- chr_rows$n   # bin is 0-based -> +1 for R index
    result[[chr]] <- counts
  }

  result
}
