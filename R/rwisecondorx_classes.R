# WisecondorX S7 classes
#
# These classes intentionally subclass base `list` so the existing
# WisecondorX algorithms can keep using `[[`, `$`, `names()`, `vapply()`, and
# other list-oriented code paths while still gaining a typed object boundary
# with validators.

.wcx_list_class <- S7::new_S3_class(
  "list",
  constructor = function(.data) .data
)

.wcx_chr_keys <- as.character(1:24)

.is_wcx_sample <- function(x) {
  inherits(x, "WisecondorXSample") || S7::S7_inherits(x, WisecondorXSample)
}

.is_wcx_reference <- function(x) {
  inherits(x, "WisecondorXReference") || S7::S7_inherits(x, WisecondorXReference)
}

.is_wcx_prediction <- function(x) {
  inherits(x, "WisecondorXPrediction") || S7::S7_inherits(x, WisecondorXPrediction)
}

.is_wcx_reference_qc <- function(x) {
  inherits(x, "WisecondorXReferenceQC") || S7::S7_inherits(x, WisecondorXReferenceQC)
}

.validate_named_numeric_vectors <- function(x, keys, label) {
  if (!is.list(x)) {
    return(sprintf("%s must be a list.", label))
  }
  if (is.null(names(x)) || !identical(names(x), keys)) {
    return(sprintf(
      "%s must be a named list with keys %s.",
      label, paste(keys, collapse = ", ")
    ))
  }
  ok <- vapply(x, function(el) {
    is.null(el) || (is.atomic(el) && is.numeric(el) && is.null(dim(el)))
  }, logical(1L))
  if (!all(ok)) {
    bad <- paste(names(x)[!ok], collapse = ", ")
    return(sprintf("%s elements must be NULL or numeric vectors. Bad keys: %s",
                   label, bad))
  }
  NULL
}

.validate_reference_branch <- function(x, suffix, expected_chr_count) {
  req <- c(
    paste0("mask", suffix),
    paste0("bins_per_chr", suffix),
    paste0("masked_bins_per_chr", suffix),
    paste0("masked_bins_per_chr_cum", suffix),
    paste0("pca_components", suffix),
    paste0("pca_mean", suffix),
    paste0("indexes", suffix),
    paste0("distances", suffix),
    paste0("null_ratios", suffix)
  )
  if (!all(req %in% names(x))) {
    return(sprintf("Reference branch '%s' is missing required keys.", suffix %||% "A"))
  }

  bins_per_chr <- x[[paste0("bins_per_chr", suffix)]]
  masked_bins_per_chr <- x[[paste0("masked_bins_per_chr", suffix)]]
  masked_bins_per_chr_cum <- x[[paste0("masked_bins_per_chr_cum", suffix)]]
  mask <- x[[paste0("mask", suffix)]]
  pca_components <- x[[paste0("pca_components", suffix)]]
  pca_mean <- x[[paste0("pca_mean", suffix)]]
  indexes <- x[[paste0("indexes", suffix)]]
  distances <- x[[paste0("distances", suffix)]]
  null_ratios <- x[[paste0("null_ratios", suffix)]]

  if (!is.numeric(bins_per_chr) || length(bins_per_chr) != expected_chr_count) {
    return(sprintf("bins_per_chr%s must be a numeric vector of length %d.",
                   suffix, expected_chr_count))
  }
  if (!is.numeric(masked_bins_per_chr) || length(masked_bins_per_chr) != expected_chr_count) {
    return(sprintf("masked_bins_per_chr%s must align with bins_per_chr%s.",
                   suffix, suffix))
  }
  if (!is.numeric(masked_bins_per_chr_cum) || length(masked_bins_per_chr_cum) != expected_chr_count) {
    return(sprintf("masked_bins_per_chr_cum%s must align with bins_per_chr%s.",
                   suffix, suffix))
  }
  if (!is.logical(mask) || length(mask) != sum(as.integer(bins_per_chr))) {
    return(sprintf("mask%s must have length sum(bins_per_chr%s).", suffix, suffix))
  }
  if (!is.matrix(indexes) || !is.matrix(distances) || !identical(dim(indexes), dim(distances))) {
    return(sprintf("indexes%s and distances%s must be matrices with identical dimensions.",
                   suffix, suffix))
  }
  if (nrow(indexes) != sum(as.integer(masked_bins_per_chr))) {
    return(sprintf("indexes%s row count must equal sum(masked_bins_per_chr%s).",
                   suffix, suffix))
  }
  if (!is.matrix(null_ratios) || nrow(null_ratios) != nrow(indexes)) {
    return(sprintf("null_ratios%s must be a matrix aligned to indexes%s rows.",
                   suffix, suffix))
  }
  if (!is.matrix(pca_components) || !is.numeric(pca_mean) ||
      ncol(pca_components) != length(pca_mean)) {
    return(sprintf("pca_components%s and pca_mean%s must describe the same feature space.",
                   suffix, suffix))
  }
  NULL
}

#' Abstract base class for WisecondorX list-like objects
#'
#' @export
WisecondorXObject <- S7::new_class(
  "WisecondorXObject",
  parent = .wcx_list_class
)

#' WisecondorX binned sample
#'
#' A list-like S7 object keyed by chromosome `"1"`–`"24"`, where each element
#' is a numeric vector of bin counts.
#'
#' @export
WisecondorXSample <- S7::new_class(
  "WisecondorXSample",
  parent = WisecondorXObject,
  validator = function(self) {
    .validate_named_numeric_vectors(self, .wcx_chr_keys, "WisecondorXSample")
  }
)

#' WisecondorX reference panel
#'
#' A list-like S7 object produced by [rwisecondorx_newref()]. It contains the
#' autosomal branch plus optional female (`.F`) and male (`.M`) gonosomal
#' branches.
#'
#' @export
WisecondorXReference <- S7::new_class(
  "WisecondorXReference",
  parent = WisecondorXObject,
  validator = function(self) {
    required <- c("binsize", "is_nipt", "trained_cutoff", "has_female", "has_male")
    if (!all(required %in% names(self))) {
      return("WisecondorXReference is missing core metadata keys.")
    }
    if (!is.numeric(self$binsize) || length(self$binsize) != 1L || self$binsize < 1L) {
      return("Reference 'binsize' must be a positive scalar.")
    }
    if (!is.logical(self$is_nipt) || length(self$is_nipt) != 1L) {
      return("Reference 'is_nipt' must be TRUE or FALSE.")
    }
    if (!is.logical(self$has_female) || length(self$has_female) != 1L) {
      return("Reference 'has_female' must be TRUE or FALSE.")
    }
    if (!is.logical(self$has_male) || length(self$has_male) != 1L) {
      return("Reference 'has_male' must be TRUE or FALSE.")
    }
    if (!is.numeric(self$trained_cutoff) || length(self$trained_cutoff) != 1L ||
        !is.finite(self$trained_cutoff) || self$trained_cutoff <= 0 ||
        self$trained_cutoff >= 1) {
      return("Reference 'trained_cutoff' must be a finite scalar in (0, 1).")
    }

    base_issue <- .validate_reference_branch(self, "", 22L)
    if (!is.null(base_issue)) return(base_issue)

    if (isTRUE(self$has_female)) {
      female_issue <- .validate_reference_branch(self, ".F", 23L)
      if (!is.null(female_issue)) return(female_issue)
    }
    if (isTRUE(self$has_male)) {
      male_issue <- .validate_reference_branch(self, ".M", 24L)
      if (!is.null(male_issue)) return(male_issue)
    }

    NULL
  }
)

#' WisecondorX prediction result
#'
#' A list-like S7 object produced by [rwisecondorx_predict()].
#'
#' @export
WisecondorXPrediction <- S7::new_class(
  "WisecondorXPrediction",
  parent = WisecondorXObject,
  validator = function(self) {
    required <- c(
      "results_r", "results_z", "results_w", "results_c", "aberrations",
      "statistics", "gender", "ref_gender", "n_reads", "binsize", "bins_per_chr"
    )
    if (!all(required %in% names(self))) {
      return("WisecondorXPrediction is missing required result fields.")
    }
    if (!is.list(self$results_r) || !is.list(self$results_z) || !is.list(self$results_w)) {
      return("Prediction ratio, z-score, and weight outputs must be lists.")
    }
    if (!identical(length(self$results_r), length(self$results_z)) ||
        !identical(length(self$results_r), length(self$results_w))) {
      return("Prediction ratio, z-score, and weight lists must align.")
    }
    if (!is.data.frame(self$results_c) || !is.data.frame(self$aberrations) ||
        !is.data.frame(self$statistics)) {
      return("Prediction segment, aberration, and statistics outputs must be data frames.")
    }
    if (!is.character(self$gender) || length(self$gender) != 1L || !(self$gender %in% c("F", "M"))) {
      return("Prediction 'gender' must be 'F' or 'M'.")
    }
    if (!is.character(self$ref_gender) || length(self$ref_gender) != 1L || !(self$ref_gender %in% c("F", "M"))) {
      return("Prediction 'ref_gender' must be 'F' or 'M'.")
    }
    if (!is.numeric(self$n_reads) || length(self$n_reads) != 1L || self$n_reads < 0) {
      return("Prediction 'n_reads' must be a non-negative scalar.")
    }
    if (!is.numeric(self$binsize) || length(self$binsize) != 1L || self$binsize < 1L) {
      return("Prediction 'binsize' must be a positive scalar.")
    }
    if (!is.numeric(self$bins_per_chr)) {
      return("Prediction 'bins_per_chr' must be numeric.")
    }
    NULL
  }
)

#' WisecondorX reference QC report
#'
#' A structured QC report returned by [rwisecondorx_ref_qc()].
#'
#' @export
WisecondorXReferenceQC <- S7::new_class(
  "WisecondorXReferenceQC",
  parent = WisecondorXObject,
  validator = function(self) {
    required <- c("overall_verdict", "worst_severity", "compat_issues", "metrics")
    if (!all(required %in% names(self))) {
      return("WisecondorXReferenceQC is missing required report fields.")
    }
    if (!is.character(self$overall_verdict) || length(self$overall_verdict) != 1L ||
        !(self$overall_verdict %in% c("PASS", "WARN", "FAIL"))) {
      return("QC 'overall_verdict' must be PASS, WARN, or FAIL.")
    }
    if (!is.numeric(self$worst_severity) || length(self$worst_severity) != 1L) {
      return("QC 'worst_severity' must be a scalar.")
    }
    if (!is.character(self$compat_issues)) {
      return("QC 'compat_issues' must be a character vector.")
    }
    if (!is.list(self$metrics)) {
      return("QC 'metrics' must be a named list.")
    }
    NULL
  }
)

.as_wcx_sample <- function(x) {
  if (.is_wcx_sample(x)) x else WisecondorXSample(x)
}

.as_wcx_reference <- function(x) {
  if (.is_wcx_reference(x)) x else WisecondorXReference(x)
}

.as_wcx_prediction <- function(x) {
  if (.is_wcx_prediction(x)) x else WisecondorXPrediction(x)
}

.as_wcx_reference_qc <- function(x) {
  if (.is_wcx_reference_qc(x)) x else WisecondorXReferenceQC(x)
}
