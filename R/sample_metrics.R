# Helpers for per-sample SeqFF / Y-unique metrics used by the preprocessing
# and BED-writing CLIs. These stay internal; the public surface remains the
# existing BED writers plus generic `metadata =`.

.null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

.too_few_reads_for_metrics_flag <- function(count_pre,
                                            count_post,
                                            min_reads = 1e5) {
  vals <- suppressWarnings(as.numeric(c(count_pre, count_post)))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(NA)
  }
  any(vals < as.numeric(min_reads))
}

.read_metrics_table <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))
  stopifnot(file.exists(path))

  first <- readLines(path, n = 1L, warn = FALSE)
  delim <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    ","
  } else if (length(first) && gregexpr(",", first, fixed = TRUE)[[1L]][1L] > 0L &&
             (!grepl("\t", first, fixed = TRUE) ||
                lengths(regmatches(first, gregexpr(",", first, fixed = TRUE))) >=
                  lengths(regmatches(first, gregexpr("\t", first, fixed = TRUE))))) {
    ","
  } else {
    "\t"
  }

  utils::read.delim(
    path,
    sep = delim,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.match_sample_col <- function(df,
                              candidates = c(
                                "sample_name", "sample", "Sample_name", "Sample",
                                "sample_id", "SampleID"
                              )) {
  hit <- candidates[candidates %in% names(df)][1L]
  if (is.na(hit) || !nzchar(hit)) NULL else hit
}

.match_metric_row <- function(path, sample_name) {
  df <- .read_metrics_table(path)
  sample_col <- .match_sample_col(df)
  if (is.null(sample_col)) {
    if (nrow(df) == 1L) {
      return(df[1L, , drop = FALSE])
    }
    stop(
      "Cannot infer sample column in metrics table: ", path,
      call. = FALSE
    )
  }

  hit <- which(as.character(df[[sample_col]]) == sample_name)
  if (!length(hit)) {
    stop(
      "Sample '", sample_name, "' not found in metrics table: ", path,
      call. = FALSE
    )
  }
  if (length(hit) > 1L) {
    stop(
      "Sample '", sample_name, "' appears multiple times in metrics table: ", path,
      call. = FALSE
    )
  }
  df[hit, , drop = FALSE]
}

.record_to_list <- function(record) {
  if (is.null(record)) {
    return(NULL)
  }
  if (is.list(record) && !length(record)) {
    return(list())
  }
  if (is.data.frame(record)) {
    if (nrow(record) != 1L) {
      stop("Metric records supplied as data frames must have exactly one row.", call. = FALSE)
    }
    out <- as.list(record[1L, , drop = FALSE])
    names(out) <- names(record)
    return(out)
  }
  if (is.atomic(record) && !is.null(names(record))) {
    return(as.list(record))
  }
  if (is.list(record) && !is.null(names(record))) {
    return(record)
  }
  stop("Metric record must be a named vector, named list, one-row data frame, or NULL.",
       call. = FALSE)
}

.metric_scalar <- function(record, candidates, type = c("numeric", "character")) {
  type <- match.arg(type)
  rec <- .record_to_list(record)
  if (is.null(rec) || !length(rec)) {
    return(NULL)
  }
  hit <- candidates[candidates %in% names(rec)][1L]
  if (is.na(hit) || !nzchar(hit)) {
    return(NULL)
  }
  value <- rec[[hit]]
  if (!length(value)) {
    return(NULL)
  }
  value <- value[[1L]]
  if (type == "numeric") {
    out <- suppressWarnings(as.numeric(value))
    if (!length(out)) return(NULL)
    out
  } else {
    if (isTRUE(is.na(value))) {
      NA_character_
    } else {
      as.character(value)
    }
  }
}

.read_y_unique_regions_file <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))
  stopifnot(file.exists(path))

  df <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df) && all(c("Chromosome", "Start", "End") %in% names(df))) {
    if (!"GeneName" %in% names(df)) {
      df$GeneName <- NA_character_
    }
    return(df[, c("Chromosome", "Start", "End", "GeneName"), drop = FALSE])
  }

  df <- utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(df) < 3L) {
    stop("Y-unique regions file must have at least 3 columns: ", path, call. = FALSE)
  }
  df <- df[, seq_len(min(4L, ncol(df))), drop = FALSE]
  names(df)[1:3] <- c("Chromosome", "Start", "End")
  if (ncol(df) < 4L) {
    df$GeneName <- NA_character_
  } else {
    names(df)[4L] <- "GeneName"
  }
  df
}

.format_y_unique_genes <- function(regions) {
  if (is.null(regions) || !nrow(regions) || !"GeneName" %in% names(regions)) {
    return(NULL)
  }
  vals <- unique(as.character(regions$GeneName))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) NULL else paste(vals, collapse = ";")
}

.format_y_unique_intervals <- function(regions) {
  if (is.null(regions) || !nrow(regions)) {
    return(NULL)
  }
  gene <- if ("GeneName" %in% names(regions)) as.character(regions$GeneName) else rep(NA_character_, nrow(regions))
  paste(
    sprintf(
      "%s:%s-%s%s",
      as.character(regions$Chromosome),
      as.integer(regions$Start),
      as.integer(regions$End),
      ifelse(!is.na(gene) & nzchar(gene), paste0("(", gene, ")"), "")
    ),
    collapse = ";"
  )
}

.seqff_metric_record <- function(record = NULL,
                                 pre = NULL,
                                 post = NULL,
                                 source = NULL) {
  rec <- .record_to_list(record)
  pre_rec <- .record_to_list(pre)
  post_rec <- .record_to_list(post)

  out <- list(
    seqff_pre = .metric_scalar(rec, c("SeqFF_pre", "seqff_pre"), "numeric"),
    enet_pre = .metric_scalar(rec, c("Enet_pre", "enet_pre"), "numeric"),
    wrsc_pre = .metric_scalar(rec, c("WRSC_pre", "wrsc_pre"), "numeric"),
    seqff_post = .metric_scalar(rec, c("SeqFF_post", "seqff_post"), "numeric"),
    enet_post = .metric_scalar(rec, c("Enet_post", "enet_post"), "numeric"),
    wrsc_post = .metric_scalar(rec, c("WRSC_post", "wrsc_post"), "numeric"),
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("seqff_source", "metrics_seqff_source"), "character")
    )
  )

  if (!is.null(pre_rec)) {
    out$seqff_pre <- .metric_scalar(pre_rec, c("SeqFF", "seqff", "SeqFF_pre", "seqff_pre"), "numeric")
    out$enet_pre <- .metric_scalar(pre_rec, c("Enet", "enet", "Enet_pre", "enet_pre"), "numeric")
    out$wrsc_pre <- .metric_scalar(pre_rec, c("WRSC", "wrsc", "WRSC_pre", "wrsc_pre"), "numeric")
  }
  if (!is.null(post_rec)) {
    out$seqff_post <- .metric_scalar(post_rec, c("SeqFF", "seqff", "SeqFF_post", "seqff_post"), "numeric")
    out$enet_post <- .metric_scalar(post_rec, c("Enet", "enet", "Enet_post", "enet_post"), "numeric")
    out$wrsc_post <- .metric_scalar(post_rec, c("WRSC", "wrsc", "WRSC_post", "wrsc_post"), "numeric")
  }

  out
}

.y_unique_metric_record <- function(record = NULL,
                                    pre = NULL,
                                    post = NULL,
                                    regions = NULL,
                                    regions_file = NULL,
                                    source = NULL) {
  rec <- .record_to_list(record)
  pre_rec <- .record_to_list(pre)
  post_rec <- .record_to_list(post)

  out <- list(
    ratio_pre = .metric_scalar(rec, c("y_unique_ratio_pre"), "numeric"),
    reads_pre = .metric_scalar(rec, c("y_unique_reads_pre"), "numeric"),
    total_pre = .metric_scalar(rec, c("total_nuclear_reads_pre"), "numeric"),
    ratio_post = .metric_scalar(rec, c("y_unique_ratio_post"), "numeric"),
    reads_post = .metric_scalar(rec, c("y_unique_reads_post"), "numeric"),
    total_post = .metric_scalar(rec, c("total_nuclear_reads_post"), "numeric"),
    regions_file = .null_coalesce(
      regions_file,
      .metric_scalar(rec, c("y_unique_regions_file", "YUniqueRegionsFile"), "character")
    ),
    genes = .metric_scalar(rec, c("y_unique_genes", "YUniqueGenes"), "character"),
    intervals = .metric_scalar(rec, c("y_unique_regions", "YUniqueRegions"), "character"),
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("y_unique_source", "metrics_y_unique_source"), "character")
    )
  )

  if (!is.null(pre_rec)) {
    out$ratio_pre <- .metric_scalar(pre_rec, c("ratio", "y_unique_ratio", "YUniqueRatio", "Ratio"), "numeric")
    out$reads_pre <- .metric_scalar(pre_rec, c("y_unique_reads", "YUniqueReads"), "numeric")
    out$total_pre <- .metric_scalar(pre_rec, c("total_nuclear_reads", "TotalNuclearReads"), "numeric")
  }
  if (!is.null(post_rec)) {
    out$ratio_post <- .metric_scalar(post_rec, c("ratio", "y_unique_ratio", "YUniqueRatio", "Ratio"), "numeric")
    out$reads_post <- .metric_scalar(post_rec, c("y_unique_reads", "YUniqueReads"), "numeric")
    out$total_post <- .metric_scalar(post_rec, c("total_nuclear_reads", "TotalNuclearReads"), "numeric")
  }

  regions_df <- NULL
  if (!is.null(post_rec) && is.data.frame(post_rec$regions)) {
    regions_df <- post_rec$regions
  } else if (!is.null(pre_rec) && is.data.frame(pre_rec$regions)) {
    regions_df <- pre_rec$regions
  } else if (is.data.frame(regions)) {
    regions_df <- regions
  } else if (!is.null(out$regions_file) && nzchar(out$regions_file) && file.exists(out$regions_file)) {
    regions_df <- .read_y_unique_regions_file(out$regions_file)
  }

  if (!is.null(regions_df)) {
    out$genes <- .null_coalesce(out$genes, .format_y_unique_genes(regions_df))
    out$intervals <- .null_coalesce(out$intervals, .format_y_unique_intervals(regions_df))
    if (is.null(out$regions_file) && !is.null(regions_file) && file.exists(regions_file)) {
      out$regions_file <- normalizePath(regions_file, winslash = "/", mustWork = TRUE)
    }
  } else if (!is.null(regions_file) && file.exists(regions_file)) {
    out$regions_file <- normalizePath(regions_file, winslash = "/", mustWork = TRUE)
  }

  out
}

.native_bin_metric_record <- function(record = NULL,
                                      rows = NULL,
                                      source = NULL) {
  if (!is.null(rows)) {
    stopifnot(is.data.frame(rows))
    record <- .bam_bin_stats_metadata(rows)
  }
  rec <- .record_to_list(record)

  out <- list(
    count_pre = .metric_scalar(
      rec,
      c("native_count_pre_sum", "count_pre_sum", "TotalMappedReads"),
      "numeric"
    ),
    count_post = .metric_scalar(
      rec,
      c("native_count_post_sum", "count_total_sum", "TotalUniqueReads"),
      "numeric"
    ),
    count_fwd = .metric_scalar(rec, c("native_count_fwd_sum", "count_fwd_sum"), "numeric"),
    count_rev = .metric_scalar(rec, c("native_count_rev_sum", "count_rev_sum"), "numeric"),
    n_nonzero_bins_post = .metric_scalar(
      rec,
      c("native_n_nonzero_bins_post", "n_nonzero_bins_post"),
      "numeric"
    ),
    gc_pre = .metric_scalar(
      rec,
      c("gc_read_perc_pre", "GCPCTBeforeFiltering", "gc_perc_pre_weighted_mean"),
      "numeric"
    ),
    gc_post = .metric_scalar(
      rec,
      c("gc_read_perc_post", "GCPCTAfterFiltering", "gc_perc_post_weighted_mean"),
      "numeric"
    ),
    mean_mapq_post = .metric_scalar(
      rec,
      c("mean_mapq_post", "mean_mapq_post_weighted_mean"),
      "numeric"
    ),
    source = .null_coalesce(
      source,
      .metric_scalar(rec, c("native_stats_source", "native_source"), "character")
    )
  )

  if (is.null(out$source) && !is.null(rows)) {
    out$source <- "bam_bin_counts(gc,mq)"
  }

  out
}

.sample_metrics_row <- function(sample_name,
                                bam = NULL,
                                seqff = NULL,
                                y_unique = NULL,
                                native = NULL,
                                filters_pre = NULL,
                                filters_post = NULL) {
  stopifnot(is.character(sample_name), length(sample_name) == 1L, nzchar(sample_name))
  seqff <- .null_coalesce(seqff, list())
  y_unique <- .null_coalesce(y_unique, list())
  native <- .null_coalesce(native, list())
  filters_pre <- .null_coalesce(.record_to_list(filters_pre), list())
  filters_post <- .null_coalesce(.record_to_list(filters_post), list())
  unique_pct <- if (!is.null(native$count_pre) &&
                    is.finite(native$count_pre) &&
                    native$count_pre > 0 &&
                    !is.null(native$count_post) &&
                    is.finite(native$count_post)) {
    as.numeric(native$count_post) / as.numeric(native$count_pre)
  } else {
    NA_real_
  }
  too_few_reads_for_metrics <- .too_few_reads_for_metrics_flag(
    count_pre = native$count_pre,
    count_post = native$count_post
  )

  data.frame(
    sample_name = sample_name,
    sample = sample_name,
    bam = if (is.null(bam)) NA_character_ else bam,
    native_count_pre_sum = unname(.null_coalesce(native$count_pre, NA_real_)),
    native_count_post_sum = unname(.null_coalesce(native$count_post, NA_real_)),
    native_count_fwd_sum = unname(.null_coalesce(native$count_fwd, NA_real_)),
    native_count_rev_sum = unname(.null_coalesce(native$count_rev, NA_real_)),
    native_n_nonzero_bins_post = unname(.null_coalesce(native$n_nonzero_bins_post, NA_real_)),
    gc_read_perc_pre = unname(.null_coalesce(native$gc_pre, NA_real_)),
    gc_read_perc_post = unname(.null_coalesce(native$gc_post, NA_real_)),
    mean_mapq_post = unname(.null_coalesce(native$mean_mapq_post, NA_real_)),
    GCPCTBeforeFiltering = unname(.null_coalesce(native$gc_pre, NA_real_)),
    GCPCTAfterFiltering = unname(.null_coalesce(native$gc_post, NA_real_)),
    TotalMappedReads = unname(.null_coalesce(native$count_pre, NA_real_)),
    TotalUniqueReads = unname(.null_coalesce(native$count_post, NA_real_)),
    TotalUniqueReadsMapped = unname(.null_coalesce(native$count_post, NA_real_)),
    TotalUniqueReadsPercent = unname(unique_pct),
    TotalUniqueReadsPercentMapped = unname(unique_pct),
    too_few_reads_for_metrics = unname(too_few_reads_for_metrics),
    SeqFF_pre = unname(.null_coalesce(seqff$seqff_pre, NA_real_)),
    Enet_pre = unname(.null_coalesce(seqff$enet_pre, NA_real_)),
    WRSC_pre = unname(.null_coalesce(seqff$wrsc_pre, NA_real_)),
    SeqFF_post = unname(.null_coalesce(seqff$seqff_post, NA_real_)),
    Enet_post = unname(.null_coalesce(seqff$enet_post, NA_real_)),
    WRSC_post = unname(.null_coalesce(seqff$wrsc_post, NA_real_)),
    NonFiltered_SeqFF = unname(.null_coalesce(seqff$seqff_pre, NA_real_)),
    NonFiltered_Enet = unname(.null_coalesce(seqff$enet_pre, NA_real_)),
    NonFiltered_WRSC = unname(.null_coalesce(seqff$wrsc_pre, NA_real_)),
    NonFiltered_ConsensusFF = unname(.null_coalesce(seqff$seqff_pre, NA_real_)),
    Filtered_SeqFF = unname(.null_coalesce(seqff$seqff_post, NA_real_)),
    Filtered_Enet = unname(.null_coalesce(seqff$enet_post, NA_real_)),
    Filtered_WRSC = unname(.null_coalesce(seqff$wrsc_post, NA_real_)),
    Filtered_ConsensusFF = unname(.null_coalesce(seqff$seqff_post, NA_real_)),
    y_unique_ratio_pre = unname(.null_coalesce(y_unique$ratio_pre, NA_real_)),
    y_unique_reads_pre = unname(.null_coalesce(y_unique$reads_pre, NA_real_)),
    total_nuclear_reads_pre = unname(.null_coalesce(y_unique$total_pre, NA_real_)),
    y_unique_ratio_post = unname(.null_coalesce(y_unique$ratio_post, NA_real_)),
    y_unique_reads_post = unname(.null_coalesce(y_unique$reads_post, NA_real_)),
    total_nuclear_reads_post = unname(.null_coalesce(y_unique$total_post, NA_real_)),
    YUniqueRatio = unname(.null_coalesce(y_unique$ratio_pre, NA_real_)),
    YUniqueRatioFiltered = unname(.null_coalesce(y_unique$ratio_post, NA_real_)),
    y_unique_regions_file = .null_coalesce(y_unique$regions_file, NA_character_),
    metrics_pre_mapq = unname(.null_coalesce(.metric_scalar(filters_pre, c("mapq"), "numeric"), NA_real_)),
    metrics_pre_require_flags = unname(.null_coalesce(.metric_scalar(filters_pre, c("require_flags"), "numeric"), NA_real_)),
    metrics_pre_exclude_flags = unname(.null_coalesce(.metric_scalar(filters_pre, c("exclude_flags"), "numeric"), NA_real_)),
    metrics_post_mapq = unname(.null_coalesce(.metric_scalar(filters_post, c("mapq"), "numeric"), NA_real_)),
    metrics_post_require_flags = unname(.null_coalesce(.metric_scalar(filters_post, c("require_flags"), "numeric"), NA_real_)),
    metrics_post_exclude_flags = unname(.null_coalesce(.metric_scalar(filters_post, c("exclude_flags"), "numeric"), NA_real_)),
    seqff_source = .null_coalesce(seqff$source, NA_character_),
    y_unique_source = .null_coalesce(y_unique$source, NA_character_),
    native_stats_source = .null_coalesce(native$source, NA_character_),
    stringsAsFactors = FALSE
  )
}

.sample_metrics_with_nipter_status <- function(row,
                                               written = NULL,
                                               error = NULL) {
  stopifnot(is.data.frame(row), nrow(row) == 1L)

  if (!"nipter_bed_status" %in% names(row)) {
    row$nipter_bed_status <- NA_character_
  }
  if (!"nipter_bed_written" %in% names(row)) {
    row$nipter_bed_written <- NA
  }
  if (!"nipter_bed_error" %in% names(row)) {
    row$nipter_bed_error <- NA_character_
  }

  if (!is.null(written)) {
    row$nipter_bed_written <- as.logical(written)[[1L]]
  }
  if (!is.null(error)) {
    err <- as.character(error)[[1L]]
    row$nipter_bed_error <- if (is.na(err) || !nzchar(err)) NA_character_ else err
  }

  if (isTRUE(row$nipter_bed_written[[1L]])) {
    row$nipter_bed_status <- "ok"
    row$nipter_bed_error <- NA_character_
  } else if (identical(row$nipter_bed_written[[1L]], FALSE) ||
             (!is.na(row$nipter_bed_error[[1L]]) &&
                nzchar(row$nipter_bed_error[[1L]]))) {
    row$nipter_bed_status <- "failed"
  }

  row
}

.sample_metrics_tabix_metadata <- function(seqff = NULL,
                                           y_unique = NULL,
                                           native = NULL,
                                           filters_pre = NULL,
                                           filters_post = NULL) {
  seqff <- .null_coalesce(seqff, list())
  y_unique <- .null_coalesce(y_unique, list())
  native <- .null_coalesce(native, list())
  filters_pre <- .null_coalesce(.record_to_list(filters_pre), list())
  filters_post <- .null_coalesce(.record_to_list(filters_post), list())

  meta <- list(
    metrics_pre_mapq = .metric_scalar(filters_pre, c("mapq"), "numeric"),
    metrics_pre_require_flags = .metric_scalar(filters_pre, c("require_flags"), "numeric"),
    metrics_pre_exclude_flags = .metric_scalar(filters_pre, c("exclude_flags"), "numeric"),
    metrics_post_mapq = .metric_scalar(filters_post, c("mapq"), "numeric"),
    metrics_post_require_flags = .metric_scalar(filters_post, c("require_flags"), "numeric"),
    metrics_post_exclude_flags = .metric_scalar(filters_post, c("exclude_flags"), "numeric"),
    ff_seqff_pre = seqff$seqff_pre,
    ff_enet_pre = seqff$enet_pre,
    ff_wrsc_pre = seqff$wrsc_pre,
    ff_consensus_pre = seqff$seqff_pre,
    ff_seqff_post = seqff$seqff_post,
    ff_enet_post = seqff$enet_post,
    ff_wrsc_post = seqff$wrsc_post,
    ff_consensus_post = seqff$seqff_post,
    seqff_source = seqff$source,
    y_unique_ratio_pre = y_unique$ratio_pre,
    y_unique_reads_pre = y_unique$reads_pre,
    y_unique_total_nuclear_reads_pre = y_unique$total_pre,
    y_unique_ratio_post = y_unique$ratio_post,
    y_unique_reads_post = y_unique$reads_post,
    y_unique_total_nuclear_reads_post = y_unique$total_post,
    y_unique_source = y_unique$source,
    y_unique_regions_file = y_unique$regions_file,
    y_unique_genes = y_unique$genes,
    y_unique_regions = y_unique$intervals,
    gc_read_perc_pre = native$gc_pre,
    gc_read_perc_post = native$gc_post,
    mean_mapq_post = native$mean_mapq_post,
    GCPCTBeforeFiltering = native$gc_pre,
    GCPCTAfterFiltering = native$gc_post,
    total_read_starts_pre = native$count_pre,
    total_read_starts_post = native$count_post,
    total_read_starts_retained_fraction = if (!is.null(native$count_pre) &&
      is.finite(native$count_pre) &&
      native$count_pre > 0 &&
      !is.null(native$count_post) &&
      is.finite(native$count_post)) {
      as.numeric(native$count_post) / as.numeric(native$count_pre)
    } else {
      NULL
    },
    native_count_fwd_sum = native$count_fwd,
    native_count_rev_sum = native$count_rev,
    native_n_nonzero_bins_post = native$n_nonzero_bins_post,
    too_few_reads_for_metrics = .too_few_reads_for_metrics_flag(
      count_pre = native$count_pre,
      count_post = native$count_post
    ),
    too_few_reads_for_metrics_min_reads = 1e5,
    native_stats_source = native$source
  )

  meta <- meta[!vapply(meta, is.null, logical(1))]
  .normalize_tabix_metadata(meta)
}
