#' Count reads per bin from a BAM or CRAM file
#'
#' Core read-counting engine shared by the WisecondorX and NIPTeR layers.
#' Reads aligned reads from a BAM or CRAM file and returns per-bin read counts
#' for chromosomes 1-22 and the sex chromosomes (X mapped to key `"23"`, Y to
#' `"24"`).
#'
#' Read filtering mirrors the `samtools view` convention: `mapq` sets the
#' minimum mapping quality; `require_flags` is a bitmask of flags that must
#' **all** be set (equivalent to `samtools view -f`); `exclude_flags` is a
#' bitmask of flags that must **all** be clear (equivalent to
#' `samtools view -F`).  Use the duckhts UDF `sam_flag_bits(flag)` to inspect
#' named flag fields, or `sam_flag_has(flag, bit)` to test individual bits.
#'
#' `rmdup` controls duplicate removal independently of the flag filters:
#' `"streaming"` applies the WisecondorX larp/larp2 consecutive-position state
#' machine and also enforces the WisecondorX improper-pair rule (paired reads
#' that are not properly paired are excluded from both counting and the dedup
#' state — this is intrinsic to the algorithm, not a flag option); `"flag"`
#' additionally excludes reads with SAM flag `0x400` set (Picard / sambamba
#' pre-marked duplicates); `"none"` applies no deduplication.
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param binsize Bin size in base pairs. Default `5000L` (WisecondorX); use
#'   `50000L` for NIPTeR-style workflows.
#' @param mapq Minimum mapping quality. Default `1L` (WisecondorX / samtools
#'   default). Set to `0L` to retain all reads regardless of MAPQ (NIPTeR).
#' @param require_flags Integer bitmask. Only reads for which
#'   `(FLAG & require_flags) == require_flags` are retained. `0L` (default)
#'   imposes no requirement. Example: `require_flags = 0x2L` keeps only
#'   properly paired reads.
#' @param exclude_flags Integer bitmask. Reads for which
#'   `(FLAG & exclude_flags) != 0` are dropped. `0L` (default) drops nothing.
#'   Example: `exclude_flags = 0xF04L` drops unmapped, secondary, QC-fail and
#'   supplementary reads (common samtools pre-filter).
#' @param rmdup Duplicate removal strategy. `"streaming"` (default) applies the
#'   WisecondorX larp/larp2 algorithm (also excludes improper pairs).
#'   `"flag"` drops reads with SAM flag `0x400`. `"none"` keeps all reads that
#'   pass the other filters.
#' @param separate_strands Logical; when `TRUE`, returns per-strand counts
#'   (forward `+` and reverse `-`). The return value changes to a list of two
#'   named lists (`fwd` and `rev`), each structured like the default return.
#'   Used by `nipter_bin_bam(separate_strands = TRUE)` for the NIPTeR
#'   SeparatedStrands object. Default `FALSE`.
#' @param con Optional open DBI connection with duckhts already loaded. If
#'   `NULL` (default) a temporary in-memory DuckDB connection is created.
#' @param reference Optional FASTA reference path for CRAM inputs.
#'
#' @return When `separate_strands = FALSE` (default): a named list with one
#'   integer vector per chromosome key (`"1"`–`"22"`, `"23"` for X, `"24"` for
#'   Y). Each vector contains per-bin read counts (bin 0 = positions 0 to
#'   `binsize - 1`). Chromosomes absent from the BAM are `NULL`.
#'
#'   When `separate_strands = TRUE`: a list with two elements, `fwd` and `rev`,
#'   each structured like the default return.
#'
#' @seealso [bam_convert_bed()], [bam_convert_npz()], `nipter_bin_bam()`
#'
#' @examples
#' \dontrun{
#' # WisecondorX defaults
#' bins <- bam_convert("sample.bam")
#'
#' # NIPTeR defaults — all mapped reads, no dedup, 50 kb bins
#' bins <- bam_convert("sample.bam", binsize = 50000L, mapq = 0L,
#'                     rmdup = "none")
#'
#' # Pre-filtered BAM: skip unmapped + secondary + supplementary, flag dedup
#' bins <- bam_convert("sample.bam",
#'                     exclude_flags = bitwOr(0x4L, bitwOr(0x100L, 0x800L)),
#'                     rmdup = "flag")
#' }
#'
#' @export
bam_convert <- function(bam,
                        binsize          = 5000L,
                        mapq             = 1L,
                        require_flags    = 0L,
                        exclude_flags    = 0L,
                        rmdup            = c("streaming", "flag", "none"),
                        separate_strands = FALSE,
                        con              = NULL,
                        reference        = NULL) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }
  stopifnot(is.numeric(binsize),       length(binsize)       == 1L, binsize       >= 1L)
  stopifnot(is.numeric(mapq),          length(mapq)          == 1L, mapq          >= 0L)
  stopifnot(is.numeric(require_flags), length(require_flags) == 1L, require_flags >= 0L)
  stopifnot(is.numeric(exclude_flags), length(exclude_flags) == 1L, exclude_flags >= 0L)
  binsize       <- as.integer(binsize)
  mapq          <- as.integer(mapq)
  require_flags <- as.integer(require_flags)
  exclude_flags <- as.integer(exclude_flags)
  stopifnot(is.logical(separate_strands), length(separate_strands) == 1L)

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  sql <- .convert_sql(bam, binsize, mapq, require_flags, exclude_flags,
                      rmdup, reference = reference,
                      separate_strands = separate_strands)
  rows <- DBI::dbGetQuery(con, sql)

  if (isTRUE(separate_strands)) {
    fwd_rows <- rows[rows$strand == "+", , drop = FALSE]
    rev_rows <- rows[rows$strand == "-", , drop = FALSE]
    list(fwd = .rows_to_bins(fwd_rows, binsize),
         rev = .rows_to_bins(rev_rows, binsize))
  } else {
    .rows_to_bins(rows, binsize)
  }
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
#' @param mapq Minimum mapping quality (default `1L`).
#' @param require_flags Integer bitmask; only reads with all bits set are kept
#'   (samtools `-f`). Default `0L` (no requirement).
#' @param exclude_flags Integer bitmask; reads with any bit set are dropped
#'   (samtools `-F`). Default `0L` (nothing dropped).
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
                            binsize       = 5000L,
                            mapq          = 1L,
                            require_flags = 0L,
                            exclude_flags = 0L,
                            rmdup         = c("streaming", "flag", "none"),
                            con           = NULL,
                            reference     = NULL,
                            index         = TRUE) {
  rmdup <- match.arg(rmdup)
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  stopifnot(is.character(bed), length(bed) == 1L, nzchar(bed))
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

  # Manage connection lifecycle here so the same con is used for bam_convert,
  # rduckhts_bgzip, and rduckhts_tabix_index.
  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  bins <- bam_convert(bam, binsize = binsize, mapq = mapq,
                      require_flags = require_flags, exclude_flags = exclude_flags,
                      rmdup = rmdup, con = con, reference = reference)

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
  utils::write.table(df, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
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

.convert_sql <- function(bam_path, binsize, mapq, require_flags, exclude_flags,
                         rmdup, reference = NULL, separate_strands = FALSE) {
  read_bam_call <- .read_bam_call(bam_path, reference = reference)

  # User-supplied flag filters (samtools -f / -F style).
  req_clause <- if (require_flags > 0L)
    sprintf("AND (flag & %d) = %d", require_flags, require_flags) else ""
  exc_clause <- if (exclude_flags > 0L)
    sprintf("AND (flag & %d) = 0", exclude_flags) else ""

  # Improper-pair filter: intrinsic to WisecondorX streaming dedup algorithm.
  # Applied only when rmdup = "streaming"; NOT a user flag option.
  improper_clause <- if (rmdup == "streaming")
    "AND NOT ((flag & 1) != 0 AND (flag & 2) = 0)" else ""

  # rmdup strategy
  dedup_where <- switch(rmdup,
    streaming = .dedup_where_streaming(),
    none      = "TRUE",
    flag      = "(flag & 1024) = 0"
  )

  # Strand separation support
  strand_select <- if (isTRUE(separate_strands))
    ",\n        CASE WHEN (flag & 16) = 0 THEN '+' ELSE '-' END AS strand" else ""
  strand_group <- if (isTRUE(separate_strands)) ", strand" else ""

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
      %s
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
    (pos0 // %d)::INTEGER AS bin%s,
    COUNT(*)::INTEGER     AS n
FROM deduped
GROUP BY rname, bin%s
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
", read_bam_call, mapq, req_clause, exc_clause, improper_clause, dedup_where,
   binsize, strand_select, strand_group)
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
