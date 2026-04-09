#' Convert BAM/CRAM to WisecondorX bin counts
#'
#' Replicates the WisecondorX `convert` step: reads aligned reads from a
#' BAM or CRAM file, applies the same filtering and duplicate-removal strategy
#' as the upstream Python implementation, and returns per-bin read counts for
#' chromosomes 1-22, X (as 23), and Y (as 24).
#'
#' The default `rmdup = "streaming"` exactly reproduces the pysam larp/larp2
#' deduplication used by WisecondorX. Key subtleties preserved:
#' - Improper pairs are invisible to the dedup state machine (they never
#'   update larp/larp2 in pysam's `continue` branch).
#' - Unpaired reads update larp but NOT larp2.
#' - larp is never reset between chromosomes.
#' - Bin assignment matches Python's `int(pos / binsize)` (truncating division).
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param reference Optional FASTA reference path for CRAM inputs. Passed to
#'   `read_bam(..., reference := ...)`. Leave `NULL` for BAM inputs.
#' @param binsize Bin size in base pairs (default 5000, matching WisecondorX
#'   convert default). The reference bin size should be a multiple of this value.
#' @param rmdup Duplicate removal strategy. One of:
#'   \describe{
#'     \item{`"streaming"`}{Consecutive-position dedup matching WisecondorX's
#'       larp/larp2 streaming state machine (default). Recommended when the BAM
#'       has not been pre-processed with a dedup tool.}
#'     \item{`"none"`}{No duplicate removal. Use for NIPT where read depth is low
#'       and duplicate removal is not recommended (corresponds to WisecondorX
#'       `--normdup` flag).}
#'     \item{`"flag"`}{Use the SAM 0x400 duplicate flag. Use when the BAM has
#'       been processed with Picard, sambamba, or similar tools.}
#'   }
#' @param con An optional open DBI connection with the duckhts extension already
#'   loaded. If `NULL` (default), a temporary in-memory DuckDB connection is
#'   created for this call.
#'
#' @return A named list with one integer vector per chromosome key ("1"-"22",
#'   "23" for X, "24" for Y). Each vector has length
#'   `floor(chr_length / binsize) + 1` and contains per-bin read counts.
#'   Chromosomes with no reads present in the BAM are `NULL`.
#'
#' @seealso
#' WisecondorX paper: Huijsdens-van Amsterdam et al. (2018).
#' Conformance reference: `wisecondorx_convert_conformance.py` in the duckhts
#' repository, which validates exact bin-for-bin agreement on real NIPT data.
#'
#' @examples
#' \dontrun{
#' bins <- bam_convert("sample.bam", binsize = 5000, rmdup = "streaming")
#' bins[["11"]]  # bin counts for chromosome 11
#' }
#'
#' @export
bam_convert <- function(bam,
                        binsize = 5000L,
                        rmdup   = c("streaming", "none", "flag"),
                        con     = NULL,
                        reference = NULL) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }
  stopifnot(is.numeric(binsize), length(binsize) == 1L, binsize >= 1L)
  binsize <- as.integer(binsize)

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

  sql <- .convert_sql(bam, binsize, rmdup, reference = reference)
  rows <- DBI::dbGetQuery(con, sql)

  .rows_to_bins(rows, binsize)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.valid_chr_map <- c(
  as.character(1:22),
  "X" = "23", "x" = "23",
  "Y" = "24", "y" = "24"
)

.convert_sql <- function(bam_path, binsize, rmdup, reference = NULL) {
  read_bam_call <- .read_bam_call(bam_path, reference = reference)

  # Base filtering: keep proper pairs and unpaired reads; remove improper pairs.
  # This matches pysam's `continue` on `not read.is_proper_pair`.
  proper_filter <- "NOT ((flag & 1) != 0 AND (flag & 2) = 0)"

  # Apply rmdup strategy
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
      AND %s
),
with_lag AS (
    SELECT *,
        -- prev_pos: last pos of ANY proper read (paired or unpaired)
        -- mirrors pysam larp which is always updated
        LAG(pos0) OVER (ORDER BY file_offset) AS prev_pos,

        -- prev_pnext: last pnext of a PAIRED proper read only.
        -- In pysam, unpaired reads update larp but NOT larp2, so larp2 retains
        -- the last paired-read pnext even when unpaired reads appear between.
        LAST_VALUE(CASE WHEN is_paired != 0 THEN pnext0 END IGNORE NULLS)
            OVER (ORDER BY file_offset ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING)
            AS prev_pnext
    FROM raw
),
deduped AS (
    SELECT * FROM with_lag
    WHERE %s
),
filtered AS (
    SELECT * FROM deduped WHERE mapq >= 1
)
SELECT
    rname,
    (pos0 // %d)::INTEGER AS bin,
    COUNT(*)::INTEGER     AS n
FROM filtered
GROUP BY rname, bin
ORDER BY
    -- Numeric chromosome order: strip 'chr' prefix, map X->23 Y->24, sort as integer.
    -- Lexicographic ORDER BY rname would give 1,10,11,...,2,20,... which is wrong.
    TRY_CAST(regexp_replace(
        CASE upper(regexp_replace(rname, '^[Cc][Hh][Rr]', ''))
            WHEN 'X' THEN '23'
            WHEN 'Y' THEN '24'
            WHEN 'M' THEN '25'
            ELSE regexp_replace(rname, '^[Cc][Hh][Rr]', '')
        END,
    '[^0-9]', '', 'g') AS INTEGER),
    bin
", read_bam_call, proper_filter, dedup_where, binsize)
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
