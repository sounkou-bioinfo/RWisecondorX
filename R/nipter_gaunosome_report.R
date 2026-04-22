# Internal helpers for batch gaunosome reporting.

.normalize_named_list_by_sample <- function(values, samples, arg) {
  if (is.null(values)) {
    return(NULL)
  }
  if (!is.list(values) || !length(values)) {
    stop(sprintf("'%s' must be NULL or a non-empty list.", arg),
         call. = FALSE)
  }

  sample_names <- vapply(samples, .sample_name, character(1L))
  if (is.null(names(values))) {
    if (length(values) != length(sample_names)) {
      stop(
        sprintf(
          "Unnamed '%s' lists must have length equal to the number of samples.",
          arg
        ),
        call. = FALSE
      )
    }
    names(values) <- sample_names
  } else {
    if (!all(sample_names %in% names(values))) {
      stop(
        sprintf("'%s' names must cover every sample name.", arg),
        call. = FALSE
      )
    }
    values <- values[sample_names]
  }

  names(values) <- sample_names
  values
}

.as_gaunosome_sample_list <- function(samples) {
  if (.is_nipt_control_group_object(samples)) {
    samples <- samples$samples
  }
  if (!is.list(samples) || !length(samples)) {
    stop("'samples' must be a non-empty list of NIPTeRSample/NIPTSample objects.",
         call. = FALSE)
  }
  ok <- vapply(samples, .is_nipt_sample_object, logical(1L))
  if (!all(ok)) {
    stop("'samples' must contain only NIPTeRSample/NIPTSample objects.",
         call. = FALSE)
  }
  samples
}

.gaunosome_report_from_scores <- function(scores, focus_chromosomes) {
  sample_names <- vapply(scores, `[[`, character(1L), "sample_name")
  names(scores) <- sample_names
  summary_rows <- lapply(scores, function(score) {
    data.frame(
      sample_name = score$sample_name,
      score$summary,
      stringsAsFactors = FALSE
    )
  })

  .as_nipt_gaunosome_report(list(
    scores = scores,
    summary = do.call(rbind, summary_rows),
    sample_names = sample_names,
    focus_chromosomes = unique(match.arg(
      focus_chromosomes,
      c("X", "Y"),
      several.ok = TRUE
    ))
  ))
}

#' Score multiple samples against a typed gaunosome reference
#'
#' Batch wrapper around \code{\link{nipter_gaunosome_score}}. Scores each
#' sample against one \code{\link{NIPTReferenceModel}} and returns a typed
#' cohort-level report with both the per-sample objects and one flattened
#' summary table.
#'
#' @param samples A non-empty list of \code{NIPTeRSample}/\code{NIPTSample}
#'   objects, or a \code{NIPTControlGroup} whose \code{$samples} should be
#'   scored.
#' @param reference A \code{NIPTReferenceModel} prepared with
#'   \code{\link{nipter_build_gaunosome_models}}.
#' @param y_unique_ratios Optional named numeric vector keyed by sample name,
#'   used when the reference includes a \code{"y_unique"} sex model.
#' @param min_controls Minimum number of non-outlier same-sex controls required
#'   for the z-score component. Default \code{2L}.
#' @param sample_predictors Optional named list keyed by sample name. Each
#'   value should be a named list of extra regression predictors for that
#'   sample, used when the fitted models include extra columns such as
#'   \code{gc_read_perc_post}.
#' @param focus_chromosomes Character vector; any subset of \code{c("X", "Y")}.
#'
#' @return A typed \code{NIPTGaunosomeReport}.
#'
#' @export
nipter_gaunosome_report <- function(samples,
                                    reference,
                                    y_unique_ratios = NULL,
                                    min_controls = 2L,
                                    sample_predictors = NULL,
                                    focus_chromosomes = c("X", "Y")) {
  sample_list <- .as_gaunosome_sample_list(samples)
  y_unique_ratios <- .normalize_named_numeric_by_sample(
    y_unique_ratios,
    sample_list,
    arg = "y_unique_ratios"
  )
  sample_predictors <- .normalize_named_list_by_sample(
    sample_predictors,
    sample_list,
    arg = "sample_predictors"
  )
  focus_chromosomes <- unique(match.arg(
    focus_chromosomes,
    c("X", "Y"),
    several.ok = TRUE
  ))

  scores <- lapply(sample_list, function(sample) {
    sample_name <- .sample_name(sample)
    nipter_gaunosome_score(
      sample = sample,
      reference = reference,
      y_unique_ratio = if (is.null(y_unique_ratios)) NULL else y_unique_ratios[[sample_name]],
      min_controls = min_controls,
      sample_predictors = if (is.null(sample_predictors)) NULL else sample_predictors[[sample_name]],
      focus_chromosomes = focus_chromosomes
    )
  })

  .gaunosome_report_from_scores(scores, focus_chromosomes = focus_chromosomes)
}

#' Write gaunosome summary output
#'
#' Writes the flattened gaunosome summary table to a tab-separated text file.
#' Accepts either a single \code{NIPTGaunosomeScore} or a batch
#' \code{NIPTGaunosomeReport}; single scores are written as one-sample reports.
#'
#' @param x A \code{NIPTGaunosomeScore} or \code{NIPTGaunosomeReport}.
#' @param outprefix Path prefix for the output file.
#'
#' @return The written summary-table path, invisibly.
#'
#' @export
write_nipter_gaunosome_output <- function(x, outprefix) {
  stopifnot(is.character(outprefix), length(outprefix) == 1L, nzchar(outprefix))
  report <- if (.is_nipt_gaunosome_report(x)) {
    x
  } else if (.is_nipt_gaunosome_score(x)) {
    .gaunosome_report_from_scores(
      list(x),
      focus_chromosomes = unique(x$summary$chromosome)
    )
  } else {
    stop("'x' must be a NIPTGaunosomeScore or NIPTGaunosomeReport.",
         call. = FALSE)
  }

  path <- paste0(outprefix, "_gaunosomes.tsv")
  utils::write.table(
    report$summary,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  invisible(path)
}
