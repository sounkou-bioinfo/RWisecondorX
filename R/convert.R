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
#' pre-marked duplicates); `"none"` applies no deduplication. The heavy lifting
#' is delegated to the native `bam_bin_counts(...)` kernel bundled in
#' `Rduckhts`; `RWisecondorX` reshapes that fixed-bin output into the
#' chromosome-keyed objects expected by the WisecondorX and NIPTeR layers.
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
#' @return When `separate_strands = FALSE` (default): a list-like
#'   `WisecondorXSample` S7 object with one integer vector per chromosome key
#'   (`"1"`–`"22"`, `"23"` for X, `"24"` for Y). Each vector contains per-bin
#'   read counts (bin 0 = positions 0 to `binsize - 1`). Chromosomes present in
#'   the BAM header are returned as dense vectors padded with trailing zeros up
#'   to the chromosome span implied by the header and `binsize`; chromosomes
#'   absent from the header are `NULL`.
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

  rows <- .bam_bin_count_rows(
    con = con,
    bam = bam,
    binsize = binsize,
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    rmdup = rmdup,
    reference = reference
  )
  chr_lengths <- .bam_chr_lengths(con, bam)

  if (isTRUE(separate_strands)) {
    list(
      fwd = .count_rows_to_bins(rows, "count_fwd",
                                chr_lengths = chr_lengths,
                                binsize = binsize),
      rev = .count_rows_to_bins(rows, "count_rev",
                                chr_lengths = chr_lengths,
                                binsize = binsize)
    )
  } else {
    .as_wcx_sample(
      .count_rows_to_bins(rows, "count_total",
                          chr_lengths = chr_lengths,
                          binsize = binsize)
    )
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
#' header-defined bins are written, including those with a count of zero.
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
#' @param metadata Optional named list or named atomic vector of provenance
#'   metadata to write as leading `##RWX_<key>=<value>` lines before the BED
#'   body. The data rows remain headerless. Default `NULL` writes no metadata.
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
                            index         = TRUE,
                            metadata      = NULL) {
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

  rows <- .bam_bin_count_rows(
    con = con,
    bam = bam,
    binsize = binsize,
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    rmdup = rmdup,
    reference = reference,
    stats = "gc,mq",
    include_unmapped = TRUE
  )
  chr_lengths <- .bam_chr_lengths(con, bam)
  bins <- .as_wcx_sample(
    .count_rows_to_bins(rows, "count_total",
                        chr_lengths = chr_lengths,
                        binsize = binsize)
  )

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
  .write_tabix_body(
    df,
    tmp,
    metadata = .merge_tabix_metadata(
      metadata,
      .merge_tabix_metadata(
        .merge_tabix_metadata(
          .bam_alignment_count_metadata(con, bam),
          .bam_read_counts_metadata(rows)
        ),
        list(read_counts_source = "samtools_idxstats+bam_bin_counts(gc,mq)")
      )
    )
  )

  Rduckhts::rduckhts_bgzip(con, tmp,
                           output_path = bed,
                           threads     = 1L,
                           keep        = TRUE,
                           overwrite   = TRUE)

  if (isTRUE(index)) {
    Rduckhts::rduckhts_tabix_index(con, bed,
                                   preset  = "bed",
                                   threads = 1L,
                                   comment_char = "#")
  }

  invisible(bed)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.bin_chr_name <- function(key) {
  switch(key, "23" = "X", "24" = "Y", key)
}

#' Query BAM header for chromosome lengths, keyed "1"-"24"
#'
#' Returns a named list: keys are "1"-"24" (matching bam_convert() convention),
#' values are integer chromosome lengths. Only chromosomes 1-22/X/Y present in
#' the BAM header are included; contigs, decoys, etc. are ignored.
#'
#' @param con Open DBI connection with duckhts loaded.
#' @param bam Path to BAM/CRAM.
#' @return Named list of integer lengths.
#' @keywords internal
.bam_chr_lengths <- function(con, bam) {
  hdr <- Rduckhts::rduckhts_hts_header(con, bam)
  sq  <- hdr[hdr$record_type == "SQ", c("id", "length"), drop = FALSE]

  result <- list()
  for (i in seq_len(nrow(sq))) {
    nm  <- sq$id[i]
    len <- sq$length[i]
    key <- .normalize_chr_name(nm)
    if (key %in% as.character(1:24)) {
      result[[key]] <- as.integer(len)
    }
  }
  result
}

.chr_n_bins_dense <- function(chr_length, binsize) {
  if (is.null(chr_length) || is.na(chr_length) || chr_length <= 0L) return(0L)
  as.integer((as.double(chr_length) + as.double(binsize) - 1) %/% as.double(binsize))
}

.count_rows_to_bins <- function(rows, value_col, chr_lengths = NULL, binsize = NULL) {
  chr_keys <- as.character(1:24)
  result <- stats::setNames(vector("list", 24L), chr_keys)
  if (!is.null(chr_lengths)) {
    stopifnot(!is.null(binsize), length(binsize) == 1L, binsize >= 1L)
    for (chr in chr_keys) {
      chr_length <- chr_lengths[[chr]]
      if (!is.null(chr_length)) {
        n_bins <- .chr_n_bins_dense(chr_length, binsize)
        result[[chr]] <- integer(n_bins)
      }
    }
  }

  if (nrow(rows) == 0L) return(result)

  chrom <- rows[[grep("^chrom$", names(rows), ignore.case = TRUE, value = TRUE)]]
  chrom <- .normalize_chr_name(chrom)
  values <- as.integer(rows[[value_col]])
  bins <- as.integer(rows$bin_id)

  for (chr in chr_keys) {
    idx <- which(chrom == chr)
    if (length(idx) == 0L) next

    chr_bins <- bins[idx]
    max_bin <- max(chr_bins)
    counts <- result[[chr]]
    if (is.null(counts) || length(counts) < (max_bin + 1L)) {
      counts <- integer(max_bin + 1L)
    }
    counts[chr_bins + 1L] <- values[idx]
    result[[chr]] <- counts
  }

  result
}

.bam_bin_count_rows <- function(con,
                                bam,
                                binsize,
                                mapq,
                                require_flags,
                                exclude_flags,
                                rmdup,
                                reference = NULL,
                                stats = NULL,
                                include_unmapped = FALSE) {
  Rduckhts::rduckhts_bam_bin_counts(
    con,
    path = bam,
    bin_width = binsize,
    include_unmapped = include_unmapped,
    reference = reference,
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    rmdup = rmdup,
    stats = stats
  )
}

.weighted_mean_or_na <- function(values, weights) {
  ok <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  stats::weighted.mean(values[ok], weights[ok])
}

.bam_alignment_count_metadata <- function(con, bam) {
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam), file.exists(bam))

  out_tsv <- tempfile(fileext = ".idxstats.tsv")
  on.exit(unlink(out_tsv), add = TRUE)

  result <- Rduckhts::rduckhts_samtools_idxstats(
    con = con,
    path = bam,
    output = out_tsv,
    overwrite = TRUE
  )
  if (!nrow(result) || !isTRUE(result$success[[1L]]) || !file.exists(out_tsv)) {
    stop(
      "Failed to compute idxstats-style alignment counts for ", bam,
      if (nrow(result) && "error_message" %in% names(result) &&
          !is.na(result$error_message[[1L]]) &&
          nzchar(result$error_message[[1L]])) {
        paste0(": ", result$error_message[[1L]])
      } else {
        "."
      },
      call. = FALSE
    )
  }

  idxstats <- utils::read.delim(
    out_tsv,
    sep = "\t",
    header = FALSE,
    check.names = FALSE,
    col.names = c("chromosome", "chromosome_length", "mapped_reads", "unmapped_reads")
  )

  mapped_total <- sum(as.numeric(idxstats$mapped_reads), na.rm = TRUE)
  unmapped_total <- sum(as.numeric(idxstats$unmapped_reads), na.rm = TRUE)

  list(
    read_counts_input_total = mapped_total + unmapped_total,
    read_counts_mapped_total = mapped_total,
    read_counts_unmapped_total = unmapped_total
  )
}

.bam_read_counts_metadata <- function(rows) {
  if (!nrow(rows)) {
    return(NULL)
  }

  out <- list()
  if ("count_total" %in% names(rows)) {
    out$read_counts_binned_post_sum <- sum(as.numeric(rows$count_total), na.rm = TRUE)
    out$read_counts_binned_nonzero_bins_post <- sum(as.numeric(rows$count_total) > 0, na.rm = TRUE)
  }
  if ("count_fwd" %in% names(rows)) {
    out$read_counts_binned_fwd_sum <- sum(as.numeric(rows$count_fwd), na.rm = TRUE)
  }
  if ("count_rev" %in% names(rows)) {
    out$read_counts_binned_rev_sum <- sum(as.numeric(rows$count_rev), na.rm = TRUE)
  }
  if ("count_pre" %in% names(rows)) {
    out$read_counts_binned_pre_sum <- sum(as.numeric(rows$count_pre), na.rm = TRUE)
  }
  if (all(c("gc_perc_pre", "count_pre") %in% names(rows))) {
    out$gc_read_perc_pre <- .weighted_mean_or_na(
      as.numeric(rows$gc_perc_pre),
      as.numeric(rows$count_pre)
    )
  }
  if (all(c("gc_perc_post", "count_total") %in% names(rows))) {
    out$gc_read_perc_post <- .weighted_mean_or_na(
      as.numeric(rows$gc_perc_post),
      as.numeric(rows$count_total)
    )
  }
  if (all(c("mean_mapq_post", "count_total") %in% names(rows))) {
    out$mean_mapq_post <- .weighted_mean_or_na(
      as.numeric(rows$mean_mapq_post),
      as.numeric(rows$count_total)
    )
  }
  if (!is.null(out$read_counts_binned_pre_sum) &&
      !is.null(out$read_counts_binned_post_sum) &&
      is.finite(out$read_counts_binned_pre_sum) &&
      out$read_counts_binned_pre_sum > 0 &&
      is.finite(out$read_counts_binned_post_sum)) {
    out$read_counts_binned_retained_fraction <-
      as.numeric(out$read_counts_binned_post_sum) /
      as.numeric(out$read_counts_binned_pre_sum)
  }

  out[!vapply(out, is.null, logical(1L))]
}
