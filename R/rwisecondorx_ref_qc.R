# Native reference QC for WisecondorX references
#
# Mirrors the upstream Python ref_qc.py heuristics for in-memory
# WisecondorXReference objects built by rwisecondorx_newref().

.ref_qc_outlier_n_sigma <- 3
.ref_qc_chry_min_usable_pct_warn <- 50
.ref_qc_chry_min_usable_pct_fail <- 20


#' QC a native WisecondorX reference
#'
#' Native R implementation of the upstream WisecondorX `ref_qc.py` heuristics
#' for `WisecondorXReference` objects. Inspects within-sample reference-bin
#' distances, compatibility between autosomal and gonosomal mask prefixes, and
#' (for male references) chrY usability metrics.
#'
#' The function returns a structured QC report and can optionally write the
#' report as JSON for later inspection.
#'
#' @param reference A `WisecondorXReference` object from
#'   [rwisecondorx_newref()].
#' @param min_ref_bins Integer; minimum number of usable normalizing reference
#'   bins retained for a target bin before the report flags a warning. This is
#'   not the number of cohort samples. Default `150L`, matching the upstream
#'   Python QC heuristic.
#' @param output_json Optional path for writing the QC report as JSON. When
#'   supplied, requires the `jsonlite` package.
#'
#' @return A named list with class `"WisecondorXReferenceQC"` containing:
#' \describe{
#'   \item{overall_verdict}{`"PASS"`, `"WARN"`, or `"FAIL"`.}
#'   \item{worst_severity}{Integer severity code: `0L` pass, `1L` warn, `2L`
#'     fail.}
#'   \item{compat_issues}{Character vector of autosomal-prefix compatibility
#'     issues between autosomal and gonosomal sub-references.}
#'   \item{metrics}{Named list of per-branch metrics (`A`, `F`, `M`) including
#'     verdict, message, mean/std of per-bin mean distances, outlier counts, and
#'     counts of target bins that fall below the normalizing-reference-bin
#'     threshold. Male reports also include chrY metrics.}
#' }
#'
#' @seealso [rwisecondorx_newref()], [rwisecondorx_predict()]
#'
#' @export
rwisecondorx_ref_qc <- function(reference,
                                min_ref_bins = 150L,
                                output_json = NULL) {
  reference <- .as_wcx_reference(reference)
  stopifnot(.is_wcx_reference(reference))
  stopifnot(is.numeric(min_ref_bins), length(min_ref_bins) == 1L,
            min_ref_bins >= 1L)
  if (!is.null(output_json)) {
    stopifnot(is.character(output_json), length(output_json) == 1L,
              nzchar(output_json))
  }

  suffixes <- .reference_qc_suffixes(reference)
  if (!length(suffixes)) {
    stop("Reference object does not contain any bins_per_chr components.")
  }

  min_ref_bins <- as.integer(min_ref_bins)
  compat_issues <- .reference_qc_compat_issues(reference)
  worst <- if (length(compat_issues)) 2L else 0L

  qc <- list(
    overall_verdict = "PASS",
    worst_severity = 0L,
    compat_issues = compat_issues,
    metrics = list(),
    branch_summary = data.frame(stringsAsFactors = FALSE),
    readiness = list(
      reference_mask_compatible = !length(compat_issues),
      any_branch_fail = NA,
      ready_for_prediction = NA,
      evaluated_branches = character(0)
    )
  )

  for (label in names(suffixes)) {
    suffix <- suffixes[[label]]
    metrics <- .reference_qc_metrics(reference, suffix, min_ref_bins)
    if (is.null(metrics)) {
      next
    }

    verdict <- if (identical(label, "M")) {
      .reference_qc_verdict_m(metrics, min_ref_bins)
    } else {
      .reference_qc_verdict_f(metrics, min_ref_bins)
    }

    worst <- max(worst, verdict$severity)

    branch_report <- list(
      n_bins = metrics$n_bins,
      n_valid = metrics$n_valid,
      n_missing_mean_distance = as.integer(metrics$n_bins - metrics$n_valid),
      verdict = verdict$verdict,
      severity = as.integer(verdict$severity),
      message = verdict$message,
      min_reference_bins_threshold = as.integer(min_ref_bins),
      mean_of_means = metrics$mean_of_means,
      std_of_means = metrics$std_of_means,
      n_mean_outlier = metrics$n_mean_outlier,
      outlier_pct = metrics$outlier_pct,
      n_low_refs = metrics$n_low_refs,
      n_target_bins_below_reference_bin_threshold = metrics$n_low_refs,
      has_low_reference_bins = isTRUE(metrics$n_low_refs > 0L),
      has_outlier_bins = isTRUE(metrics$n_mean_outlier > 0L),
      branch_ready_for_prediction = isTRUE(verdict$severity < 2L)
    )
    if (!is.null(metrics$chrY)) {
      branch_report$chrY <- metrics$chrY
      if (isTRUE(metrics$chrY$n_valid > 0L) && is.finite(metrics$chrY$usable_pct)) {
        branch_report$has_chrY_low_usable_pct_warn <-
          metrics$chrY$usable_pct < .ref_qc_chry_min_usable_pct_warn
        branch_report$has_chrY_low_usable_pct_fail <-
          metrics$chrY$usable_pct < .ref_qc_chry_min_usable_pct_fail
      } else {
        branch_report$has_chrY_low_usable_pct_warn <- NA
        branch_report$has_chrY_low_usable_pct_fail <- NA
      }
    }
    qc$metrics[[label]] <- branch_report
  }

  qc$worst_severity <- as.integer(worst)
  qc$overall_verdict <- c("PASS", "WARN", "FAIL")[qc$worst_severity + 1L]
  if (length(qc$metrics)) {
    qc$branch_summary <- do.call(
      rbind,
      lapply(names(qc$metrics), function(branch) {
        m <- qc$metrics[[branch]]
        data.frame(
          branch = branch,
          verdict = m$verdict,
          severity = m$severity,
          n_bins = m$n_bins,
          n_valid = m$n_valid,
          n_missing_mean_distance = m$n_missing_mean_distance,
          n_low_refs = m$n_low_refs,
          has_low_reference_bins = m$has_low_reference_bins,
          has_outlier_bins = m$has_outlier_bins,
          branch_ready_for_prediction = m$branch_ready_for_prediction,
          stringsAsFactors = FALSE
        )
      })
    )
  }
  qc$readiness <- list(
    reference_mask_compatible = !length(compat_issues),
    any_branch_fail = any(vapply(qc$metrics, function(x) identical(x$verdict, "FAIL"), logical(1L))),
    ready_for_prediction = (qc$worst_severity < 2L) && !length(compat_issues),
    evaluated_branches = names(qc$metrics)
  )
  qc <- .as_wcx_reference_qc(qc)

  if (!is.null(output_json)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("output_json requires the 'jsonlite' package to be installed.")
    }
    out_dir <- dirname(output_json)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    jsonlite::write_json(unclass(qc), output_json,
                         auto_unbox = TRUE, pretty = TRUE,
                         null = "null")
  }

  qc
}


.reference_qc_suffixes <- function(reference) {
  keys <- names(reference)
  suffixes <- character(0)
  if ("bins_per_chr.F" %in% keys) {
    suffixes <- c(suffixes, F = ".F")
  }
  if ("bins_per_chr.M" %in% keys) {
    suffixes <- c(suffixes, M = ".M")
  }
  if ("bins_per_chr" %in% keys && !length(suffixes)) {
    suffixes <- c(suffixes, A = "")
  }
  suffixes
}


.reference_qc_compat_issues <- function(reference) {
  required_keys <- c("mask", "bins_per_chr", "masked_bins_per_chr_cum")
  if (!all(required_keys %in% names(reference))) {
    return(character(0))
  }

  a_bins <- reference$bins_per_chr
  a_cum <- reference$masked_bins_per_chr_cum
  if (length(a_bins) < 22L || length(a_cum) < 22L) {
    return(character(0))
  }

  a_autosomal_span <- as.integer(sum(a_bins[seq_len(22L)]))
  a_autosomal_masked <- as.integer(a_cum[22L])
  a_prefix <- reference$mask[seq_len(a_autosomal_span)]
  issues <- character(0)

  for (suffix in c(".F", ".M")) {
    bins_key <- paste0("bins_per_chr", suffix)
    cum_key <- paste0("masked_bins_per_chr_cum", suffix)
    mask_key <- paste0("mask", suffix)
    if (!all(c(bins_key, cum_key, mask_key) %in% names(reference))) {
      next
    }

    s_bins <- reference[[bins_key]]
    s_cum <- reference[[cum_key]]
    s_mask <- reference[[mask_key]]
    if (length(s_bins) < 22L || length(s_cum) < 22L) {
      next
    }

    s_autosomal_span <- as.integer(sum(s_bins[seq_len(22L)]))
    s_autosomal_masked <- as.integer(s_cum[22L])
    if (s_autosomal_masked != a_autosomal_masked) {
      issues <- c(
        issues,
        sprintf(
          "%s autosomal masked-bin count (%d) differs from autosomal reference A (%d)",
          suffix,
          s_autosomal_masked,
          a_autosomal_masked
        )
      )
    }

    s_prefix <- s_mask[seq_len(s_autosomal_span)]
    if (s_autosomal_span != a_autosomal_span || !identical(s_prefix, a_prefix)) {
      issues <- c(
        issues,
        sprintf(
          "%s autosomal mask prefix is not aligned to the finalized autosomal reference A",
          suffix
        )
      )
    }
  }

  issues
}


.reference_qc_per_bin_stats <- function(indexes, distances) {
  indexes <- as.matrix(indexes)
  distances <- as.matrix(distances)

  valid <- indexes > 0L & is.finite(distances)
  n_refs <- rowSums(valid)
  sum_dist <- rowSums(ifelse(valid, distances, 0))
  mean_dist <- rep(NA_real_, nrow(distances))
  has_refs <- n_refs > 0L
  mean_dist[has_refs] <- sum_dist[has_refs] / n_refs[has_refs]

  list(mean_dist = mean_dist, n_refs = as.integer(n_refs))
}


.reference_qc_chrY_metrics <- function(reference,
                                       suffix,
                                       mean_dist,
                                       n_refs,
                                       cutoff_outlier,
                                       min_ref_bins) {
  if (!identical(suffix, ".M")) {
    return(NULL)
  }

  cum_key <- paste0("masked_bins_per_chr_cum", suffix)
  if (!(cum_key %in% names(reference))) {
    return(NULL)
  }

  masked_bins_cum <- reference[[cum_key]]
  if (length(masked_bins_cum) < 24L) {
    return(NULL)
  }

  start_idx <- as.integer(masked_bins_cum[23L]) + 1L
  end_idx <- as.integer(masked_bins_cum[24L])
  if (start_idx > end_idx) {
    return(list(n_bins = 0L))
  }

  mean_y <- mean_dist[start_idx:end_idx]
  n_refs_y <- n_refs[start_idx:end_idx]
  valid <- is.finite(mean_y)
  n_valid <- sum(valid)
  if (!n_valid) {
    return(list(n_bins = length(mean_y), n_valid = 0L, mean_of_means = NA_real_))
  }

  usable <- valid & (n_refs_y >= min_ref_bins)
  n_usable <- sum(usable)
  std_of_means <- if (n_valid > 1L) stats::sd(mean_y[valid]) else 0

  list(
    n_bins = length(mean_y),
    n_valid = as.integer(n_valid),
    mean_of_means = mean(mean_y[valid]),
    std_of_means = std_of_means,
    n_mean_outlier = as.integer(sum(mean_y[valid] >= cutoff_outlier)),
    n_low_refs = as.integer(sum(n_refs_y < min_ref_bins)),
    n_usable = as.integer(n_usable),
    usable_pct = 100 * n_usable / max(length(mean_y), 1L)
  )
}


.reference_qc_metrics <- function(reference, suffix, min_ref_bins) {
  idx_key <- paste0("indexes", suffix)
  dist_key <- paste0("distances", suffix)
  if (!(idx_key %in% names(reference)) || !(dist_key %in% names(reference))) {
    return(NULL)
  }

  indexes <- as.matrix(reference[[idx_key]])
  distances <- as.matrix(reference[[dist_key]])
  n_bins <- nrow(indexes)
  if (!n_bins) {
    return(list(n_bins = 0L))
  }

  stats <- .reference_qc_per_bin_stats(indexes, distances)
  valid <- is.finite(stats$mean_dist)
  n_valid <- sum(valid)
  if (!n_valid) {
    return(list(n_bins = n_bins, n_valid = 0L))
  }

  mean_of_means <- mean(stats$mean_dist[valid])
  std_of_means <- if (n_valid > 1L) stats::sd(stats$mean_dist[valid]) else 0
  cutoff_outlier <- mean_of_means + .ref_qc_outlier_n_sigma * std_of_means

  list(
    n_bins = n_bins,
    n_valid = as.integer(n_valid),
    mean_of_means = mean_of_means,
    std_of_means = std_of_means,
    n_mean_outlier = as.integer(sum(stats$mean_dist[valid] >= cutoff_outlier)),
    outlier_pct = 100 * sum(stats$mean_dist[valid] >= cutoff_outlier) / n_valid,
    n_low_refs = as.integer(sum(stats$n_refs < min_ref_bins)),
    chrY = .reference_qc_chrY_metrics(reference, suffix, stats$mean_dist,
                                      stats$n_refs, cutoff_outlier,
                                      min_ref_bins)
  )
}


.reference_qc_verdict_f <- function(metrics, min_ref_bins) {
  if (is.null(metrics) || isTRUE(metrics$n_valid == 0L)) {
    return(list(verdict = "FAIL", severity = 2L, message = "no data"))
  }
  if (metrics$n_low_refs > 0L) {
    return(list(
      verdict = "WARN",
      severity = 1L,
      message = sprintf(
        "usable normalizing reference bins per target < %d for %d target bins",
        min_ref_bins,
        metrics$n_low_refs
      )
    ))
  }
  if (metrics$std_of_means > 10) {
    return(list(
      verdict = "FAIL",
      severity = 2L,
      message = sprintf("std(per-bin mean dist) = %.2f (high)",
                        metrics$std_of_means)
    ))
  }
  if (metrics$std_of_means > 2) {
    return(list(
      verdict = "WARN",
      severity = 1L,
      message = sprintf("std(per-bin mean dist) = %.2f", metrics$std_of_means)
    ))
  }
  if (metrics$outlier_pct > 1) {
    return(list(
      verdict = "WARN",
      severity = 1L,
      message = sprintf("outlier bins = %.2f%%", metrics$outlier_pct)
    ))
  }

  list(verdict = "PASS", severity = 0L, message = "")
}


.reference_qc_verdict_m <- function(metrics, min_ref_bins) {
  if (is.null(metrics) || isTRUE(metrics$n_valid == 0L)) {
    return(list(verdict = "FAIL", severity = 2L, message = "no data"))
  }

  verdict <- "PASS"
  severity <- 0L
  message <- ""
  update <- function(new_verdict, new_severity, new_message) {
    if (new_severity > severity) {
      verdict <<- new_verdict
      severity <<- new_severity
      message <<- new_message
    }
  }

  if (metrics$n_low_refs > 0L) {
    update("WARN", 1L,
           sprintf(
             "usable normalizing reference bins per target < %d for %d target bins",
             min_ref_bins,
             metrics$n_low_refs
           ))
  }
  if (metrics$mean_of_means > 10) {
    update("FAIL", 2L,
           sprintf("mean(per-bin mean dist) = %.2f (heavy tail)",
                   metrics$mean_of_means))
  } else if (metrics$mean_of_means > 2) {
    update("WARN", 1L,
           sprintf("mean(per-bin mean dist) = %.2f", metrics$mean_of_means))
  }

  chrY <- metrics$chrY
  if (!is.null(chrY) && isTRUE(chrY$n_valid > 0L) &&
      is.finite(chrY$usable_pct)) {
    if (chrY$usable_pct < .ref_qc_chry_min_usable_pct_fail) {
      update("FAIL", 2L,
             sprintf("usable chrY bins = %.1f%% (<%d%%)",
                     chrY$usable_pct, .ref_qc_chry_min_usable_pct_fail))
    } else if (chrY$usable_pct < .ref_qc_chry_min_usable_pct_warn) {
      update("WARN", 1L,
             sprintf("usable chrY bins = %.1f%% (<%d%%)",
                     chrY$usable_pct, .ref_qc_chry_min_usable_pct_warn))
    }
  }
  if (!is.null(chrY) && isTRUE(chrY$n_valid > 0L) &&
      is.finite(chrY$mean_of_means)) {
    if (chrY$mean_of_means > 100) {
      update("FAIL", 2L,
             sprintf("chrY mean distance = %.1f (very poor chrY)",
                     chrY$mean_of_means))
    } else if (chrY$mean_of_means > 5) {
      update("WARN", 1L,
             sprintf("chrY mean distance = %.1f", chrY$mean_of_means))
    }
  }
  if (metrics$outlier_pct > 1) {
    update("WARN", 1L,
           sprintf("outlier bins = %.2f%%", metrics$outlier_pct))
  }

  list(verdict = verdict, severity = severity, message = message)
}
