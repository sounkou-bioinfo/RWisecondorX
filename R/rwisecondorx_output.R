# WisecondorX native R implementation — output generation
#
# Ports predict_output.py from upstream WisecondorX.
# Writes BED files and statistics files in the same format as the upstream
# Python implementation.


#' Write prediction output files
#'
#' Writes bins BED, segments BED, aberrations BED, and statistics text file.
#'
#' @param prediction A `WisecondorXPrediction` object.
#' @param outprefix Path prefix for output files.
#'
#' @return `outprefix` (invisibly).
#'
#' @export
write_wisecondorx_output <- function(prediction, outprefix) {
  stopifnot(inherits(prediction, "WisecondorXPrediction"))
  stopifnot(is.character(outprefix), length(outprefix) == 1L, nzchar(outprefix))
  .write_prediction_output(prediction, outprefix)
  invisible(outprefix)
}


#' Internal output writer
#' @keywords internal
.write_prediction_output <- function(prediction, outprefix) {
  .write_bins_bed(prediction, outprefix)
  .write_segments_bed(prediction, outprefix)
  .write_aberrations_bed(prediction, outprefix)
  .write_statistics(prediction, outprefix)
}


#' Write bins BED file
#' @keywords internal
.write_bins_bed <- function(prediction, outprefix) {
  path <- paste0(outprefix, "_bins.bed")
  binsize <- prediction$binsize
  lines <- character(0)
  lines <- c(lines, "chr\tstart\tend\tid\tratio\tzscore")

  for (chr_idx in seq_along(prediction$results_r)) {
    r <- prediction$results_r[[chr_idx]]
    z <- prediction$results_z[[chr_idx]]
    chr_name <- .idx_to_chr_name(chr_idx)
    n <- length(r)
    if (n == 0L) next

    feat <- 1L
    for (i in seq_len(n)) {
      ri <- if (r[i] == 0) "nan" else as.character(round(r[i], 6))
      zi <- if (z[i] == 0) "nan" else as.character(round(z[i], 6))
      feat_str <- sprintf("%s:%d-%d", chr_name, feat, feat + binsize - 1L)
      lines <- c(lines, paste(chr_name, feat, feat + binsize - 1L,
                               feat_str, ri, zi, sep = "\t"))
      feat <- feat + binsize
    }
  }

  writeLines(lines, path)
  invisible(path)
}


#' Write segments BED file
#' @keywords internal
.write_segments_bed <- function(prediction, outprefix) {
  path <- paste0(outprefix, "_segments.bed")
  binsize <- prediction$binsize
  rc <- prediction$results_c
  lines <- c("chr\tstart\tend\tratio\tzscore")

  if (nrow(rc) > 0L) {
    for (i in seq_len(nrow(rc))) {
      chr_name <- .idx_to_chr_name(rc$chr[i])
      start_bp <- as.integer(rc$start[i] * binsize + 1L)
      end_bp   <- as.integer(rc$end[i] * binsize)
      lines <- c(lines, paste(chr_name, start_bp, end_bp,
                               round(rc$ratio[i], 6),
                               round(rc$zscore[i], 6), sep = "\t"))
    }
  }

  writeLines(lines, path)
  invisible(path)
}


#' Write aberrations BED file
#' @keywords internal
.write_aberrations_bed <- function(prediction, outprefix) {
  path <- paste0(outprefix, "_aberrations.bed")
  binsize <- prediction$binsize
  ab <- prediction$aberrations
  lines <- c("chr\tstart\tend\tratio\tzscore\ttype")

  if (nrow(ab) > 0L) {
    for (i in seq_len(nrow(ab))) {
      chr_name <- .idx_to_chr_name(ab$chr[i])
      start_bp <- as.integer(ab$start[i] * binsize + 1L)
      end_bp   <- as.integer(ab$end[i] * binsize)
      lines <- c(lines, paste(chr_name, start_bp, end_bp,
                               round(ab$ratio[i], 6),
                               round(ab$zscore[i], 6),
                               ab$type[i], sep = "\t"))
    }
  }

  writeLines(lines, path)
  invisible(path)
}


#' Write statistics file
#' @keywords internal
.write_statistics <- function(prediction, outprefix) {
  path <- paste0(outprefix, "_statistics.txt")
  stats_df <- prediction$statistics
  lines <- c("chr\tratio.mean\tratio.median\tzscore")

  for (i in seq_len(nrow(stats_df))) {
    lines <- c(lines, paste(stats_df$chr[i], stats_df$ratio_mean[i],
                             stats_df$ratio_median[i], stats_df$zscore[i],
                             sep = "\t"))
  }

  # Append summary lines
  lines <- c(lines,
    sprintf("Gender based on --yfrac (or manually overridden by --gender): %s",
            prediction$gender),
    sprintf("Number of reads: %d", prediction$n_reads),
    sprintf("Standard deviation of the ratios per chromosome: %s",
            round(stats::sd(stats_df$ratio_mean, na.rm = TRUE), 5)),
    sprintf("Median segment variance per bin (doi: 10.1093/nar/gky1263): %s",
            .get_msv(prediction$results_c, prediction$results_r)),
    sprintf("Copy number profile abnormality (CPA) score (doi: 10.1186/s13073-020-00735-4): %s",
            .get_cpa(prediction$results_c, prediction$binsize))
  )

  writeLines(lines, path)
  invisible(path)
}


#' Convert 1-based chromosome index to display name
#' @keywords internal
.idx_to_chr_name <- function(idx) {
  if (idx == 23L) "X"
  else if (idx == 24L) "Y"
  else as.character(idx)
}
