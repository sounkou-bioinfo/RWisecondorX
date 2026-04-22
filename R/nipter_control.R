#' Build a NIPTeR control group from a list of binned samples
#'
#' Collects multiple \code{NIPTeRSample} objects into a
#' \code{NIPTeRControlGroup}, the reference cohort used by
#' \code{\link{nipter_z_score}}, \code{\link{nipter_ncv_score}},
#' \code{\link{nipter_chi_correct}}, and \code{\link{nipter_regression}}.
#'
#' All samples must share the same strand type (\code{"CombinedStrands"} or
#' \code{"SeparatedStrands"}). Duplicate sample names are silently removed.
#'
#' @param samples A list of \code{NIPTeRSample} objects.
#' @param description Label for the group (default
#'   \code{"General control group"}).
#' @param sample_sex Optional character vector of known sex labels for the
#'   control samples. Accepted values are \code{"female"}, \code{"male"},
#'   \code{"ambiguous"}, and \code{"unknown"}. When unnamed, values are matched
#'   in sample order; when named, names must match \code{sample_name}.
#' @param sex_source Optional string describing where \code{sample_sex} came
#'   from, e.g. \code{"explicit"}, \code{"consensus_gmm"}, or
#'   \code{"laboratory_lims"}.
#'
#' @return An object of class \code{c("NIPTeRControlGroup", <strand_type>)}.
#'
#' @seealso [nipter_diagnose_control_group()], [nipter_match_control_group()],
#'   [nipter_z_score()], [nipter_chi_correct()]
#'
#' @examples
#' \dontrun{
#' samples <- lapply(bam_files, nipter_bin_bam)
#' cg <- nipter_as_control_group(samples)
#' }
#'
#' @export
nipter_as_control_group <- function(samples,
                                    description = "General control group",
                                    sample_sex = NULL,
                                    sex_source = NULL) {
  stopifnot(is.list(samples), length(samples) >= 2L)
  stopifnot(is.character(description), length(description) == 1L,
            nzchar(description))

  # Validate that all entries are NIPTeRSample (S3) or NIPTSample (S7)
  ok <- vapply(samples, .is_nipt_sample_object, logical(1L))
  if (!all(ok)) {
    stop("All elements of 'samples' must be NIPTeRSample or NIPTSample objects.",
         call. = FALSE)
  }

  # Validate same strand type using .strand_type_of() which handles both
  strand_types <- vapply(samples, .strand_type_of, character(1L))
  if (length(unique(strand_types)) > 1L) {
    stop(
      "All samples must have the same strand type. Found: ",
      paste(unique(strand_types), collapse = ", "),
      call. = FALSE
    )
  }
  strand_type_val <- strand_types[1L]   # "combined" or "separated"

  # Remove duplicate sample names (keep first occurrence)
  names_vec <- vapply(samples, .sample_name, character(1L))
  dups <- duplicated(names_vec)
  if (any(dups)) {
    message("Removing ", sum(dups), " duplicate sample(s) by name.")
    samples <- samples[!dups]
    if (!is.null(sample_sex)) {
      if (is.null(names(sample_sex))) {
        sample_sex <- sample_sex[!dups]
      } else {
        sample_sex <- sample_sex[names_vec[!dups]]
      }
    }
  }

  sample_sex <- .normalize_sample_sex(sample_sex, samples)

  # Return S7 control group
  if (identical(strand_type_val, "separated")) {
    SeparatedControlGroup(
      samples = samples,
      description = description,
      sample_sex = sample_sex,
      sex_source = sex_source
    )
  } else {
    CombinedControlGroup(
      samples = samples,
      description = description,
      sample_sex = sample_sex,
      sex_source = sex_source
    )
  }
}


#' Diagnose a NIPTeR control group
#'
#' Computes per-chromosome Z-scores across all samples in a control group and
#' flags outliers (|Z| > 3). The Shapiro-Wilk test is applied to each
#' chromosome's fraction distribution.
#'
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param collapse_strands Logical; for \code{SeparatedStrands} control groups,
#'   collapse forward and reverse strand fractions into 22 autosomal fractions
#'   (\code{TRUE}, default) or keep the upstream-style 44-row strand-resolved
#'   diagnostics (\code{FALSE}). Ignored for \code{CombinedStrands}.
#' @param z_cutoff Absolute Z-score cutoff used to flag aberrant rows.
#'   Default \code{3}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{z_scores}{Chromosome-by-sample matrix of Z-scores (rows =
#'     chromosomes 1-22, or strand-resolved rows 1F-22F/1R-22R when
#'     \code{collapse_strands = FALSE} on a \code{SeparatedStrands}
#'     control group; columns = samples).}
#'   \item{aberrant_scores}{A \code{data.frame} with columns
#'     \code{chromosome}, \code{sample_name}, \code{z_score} for all
#'     \code{|Z| > z_cutoff}, or \code{NULL} if none.}
#'   \item{statistics}{A matrix with rows per chromosome and columns
#'     \code{mean}, \code{SD}, \code{shapiro_p_value}.}
#' }
#'
#' @seealso [nipter_as_control_group()]
#'
#' @examples
#' \dontrun{
#' diag <- nipter_diagnose_control_group(cg)
#' diag$aberrant_scores
#' }
#'
#' @export
nipter_diagnose_control_group <- function(control_group,
                                         collapse_strands = TRUE,
                                         z_cutoff = 3) {
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.logical(collapse_strands), length(collapse_strands) == 1L)
  stopifnot(is.numeric(z_cutoff), length(z_cutoff) == 1L,
            is.finite(z_cutoff), z_cutoff > 0)

  fracs <- .diagnose_control_group_fractions(
    control_group,
    collapse_strands = collapse_strands
  )

  # Z-score each chromosome across samples
  z_mat <- t(apply(fracs, 1L, scale))
  sample_names <- control_names(control_group)
  colnames(z_mat) <- sample_names
  rownames(z_mat) <- rownames(fracs)

  # Find aberrant |Z| > cutoff
  aberrant <- which(abs(z_mat) > z_cutoff, arr.ind = TRUE)
  if (nrow(aberrant) > 0L) {
    ab_df <- data.frame(
      chromosome  = rownames(z_mat)[aberrant[, 1L]],
      sample_name = colnames(z_mat)[aberrant[, 2L]],
      z_score     = z_mat[aberrant],
      stringsAsFactors = FALSE
    )
  } else {
    ab_df <- NULL
  }

  # Per-chromosome stats
  chr_means <- rowMeans(fracs)
  chr_sds   <- apply(fracs, 1L, stats::sd)
  chr_shap  <- apply(fracs, 1L, function(x) {
    if (length(unique(x)) < 3L) return(NA_real_)
    stats::shapiro.test(x)$p.value
  })

  stats_mat <- cbind(mean = chr_means, SD = chr_sds,
                     shapiro_p_value = chr_shap)
  rownames(stats_mat) <- rownames(fracs)


  list(
    z_scores        = z_mat,
    aberrant_scores = ab_df,
    statistics      = stats_mat
  )
}

.diagnose_control_group_fractions <- function(control_group,
                                              collapse_strands = TRUE) {
  separated <- identical(.strand_type_of(control_group), "separated")
  if (separated && !isTRUE(collapse_strands)) {
    fracs <- fractions_for_regression(control_group)
  } else {
    fracs <- fractions_auto(control_group)
  }
  stopifnot(is.matrix(fracs), ncol(fracs) == n_controls(control_group))
  fracs
}


#' Drop samples from a NIPTeR control group by sample name
#'
#' Convenience helper for quickly pruning a control cohort before rebuilding
#' downstream reference artifacts.
#'
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param sample_names Character vector of sample names to drop.
#'
#' @return A \code{NIPTeRControlGroup} of the same strand type with the
#'   requested samples removed.
#'
#' @seealso [nipter_as_control_group()], [nipter_diagnose_control_group()]
#'
#' @export
nipter_drop_control_group_samples <- function(control_group, sample_names) {
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.character(sample_names))

  sample_names <- unique(sample_names[nzchar(sample_names)])
  if (!length(sample_names)) {
    return(control_group)
  }

  keep <- !control_names(control_group) %in% sample_names
  if (sum(keep) < 2L) {
    stop("Dropping the requested samples would leave fewer than 2 controls.",
         call. = FALSE)
  }

  samples <- control_group$samples[keep]
  sample_sex <- .control_sample_sex(control_group)
  if (!is.null(sample_sex)) {
    sample_sex <- sample_sex[control_names(control_group)[keep]]
  }

  nipter_as_control_group(
    samples,
    description = control_group$description,
    sample_sex = sample_sex,
    sex_source = .control_sex_source(control_group)
  )
}

.match_sample_qc_col <- function(df,
                                 provided,
                                 candidates,
                                 label) {
  nms <- names(df)
  if (!length(nms)) {
    stop("sample_qc must have named columns.", call. = FALSE)
  }

  if (!is.null(provided) && nzchar(provided)) {
    hit <- match(tolower(provided), tolower(nms), nomatch = 0L)
    if (hit < 1L) {
      stop(
        "Could not find sample_qc column '", provided,
        "' for ", label, ".",
        call. = FALSE
      )
    }
    return(nms[[hit]])
  }

  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) {
    stop(
      "Could not infer a sample_qc ", label, " column. Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  nms[[hit[[1L]]]]
}

.append_filter_reason <- function(current, new_reason) {
  ifelse(
    is.na(current) | !nzchar(current),
    new_reason,
    paste(current, new_reason, sep = ";")
  )
}

.flag_mad_outliers <- function(values,
                               cutoff,
                               label) {
  stopifnot(is.numeric(values))
  stopifnot(is.numeric(cutoff), length(cutoff) == 1L, is.finite(cutoff), cutoff > 0)

  finite <- is.finite(values)
  if (sum(finite) < 3L) {
    stop(
      "Need at least 3 finite ", label, " values to apply MAD-based filtering.",
      call. = FALSE
    )
  }

  center <- stats::median(values[finite], na.rm = TRUE)
  spread <- stats::mad(values[finite], center = center, constant = 1.4826,
                       na.rm = TRUE)
  if (!is.finite(spread) || spread <= .Machine$double.eps) {
    return(rep(FALSE, length(values)))
  }

  out <- rep(FALSE, length(values))
  out[!finite] <- TRUE
  out[finite] <- abs(values[finite] - center) / spread > cutoff
  out
}

#' Filter control samples by read-depth and GC QC
#'
#' Applies hard sample-level QC gates before reference/model building. This is
#' intended for cohort curation steps such as excluding controls with too few
#' unique reads or aberrant GC content.
#'
#' @param control_group A \code{NIPTeRControlGroup}.
#' @param sample_qc Data frame containing at least a sample-name column and, if
#'   the corresponding filters are enabled, QC columns for total unique reads
#'   and/or GC content.
#' @param sample_col Optional sample-name column in \code{sample_qc}. When
#'   \code{NULL}, common names such as \code{sample_name} and \code{Sample}
#'   are inferred.
#' @param total_unique_reads_col Optional total-unique-reads column in
#'   \code{sample_qc}. Required when either read-depth threshold is supplied.
#'   If \code{NULL}, \code{read_counts_binned_post_sum} is inferred.
#' @param gc_col Optional GC column in \code{sample_qc}. Required when
#'   \code{gc_mad_cutoff} is supplied. If \code{NULL},
#'   \code{gc_read_perc_post} is inferred.
#' @param min_total_unique_reads Optional minimum allowed total unique reads.
#'   Samples below this threshold are dropped.
#' @param max_total_unique_reads Optional maximum allowed total unique reads.
#'   Samples above this threshold are dropped.
#' @param gc_mad_cutoff Optional robust MAD cutoff for the GC column. Samples
#'   more than this many MADs from the cohort median are dropped.
#'
#' @return A list with:
#' \describe{
#'   \item{control_group}{The filtered control group.}
#'   \item{retained_samples}{Character vector of retained sample names.}
#'   \item{excluded_samples}{Character vector of excluded sample names.}
#'   \item{exclusion_table}{Data frame describing dropped samples and reasons.}
#'   \item{matched_qc}{The QC rows matched to the control-group sample order.}
#'   \item{settings}{Named list of applied QC filter settings.}
#' }
#'
#' @seealso [nipter_drop_control_group_samples()], [nipter_build_reference()]
#'
#' @export
nipter_filter_control_group_qc <- function(control_group,
                                           sample_qc,
                                           sample_col = NULL,
                                           total_unique_reads_col = NULL,
                                           gc_col = NULL,
                                           min_total_unique_reads = NULL,
                                           max_total_unique_reads = NULL,
                                           gc_mad_cutoff = NULL) {
  stopifnot(.is_nipt_control_group_object(control_group))

  thresholds_active <- !is.null(min_total_unique_reads) ||
    !is.null(max_total_unique_reads) ||
    !is.null(gc_mad_cutoff)
  if (!thresholds_active) {
    return(list(
      control_group = control_group,
      retained_samples = control_names(control_group),
      excluded_samples = character(),
      exclusion_table = data.frame(
        sample_name = character(),
        exclusion_reason = character(),
        total_unique_reads = numeric(),
        gc = numeric(),
        stringsAsFactors = FALSE
      ),
      matched_qc = NULL,
      settings = list(
        min_total_unique_reads = min_total_unique_reads,
        max_total_unique_reads = max_total_unique_reads,
        gc_mad_cutoff = gc_mad_cutoff
      )
    ))
  }

  if (is.null(sample_qc)) {
    stop(
      "sample_qc must be supplied when QC thresholds are enabled.",
      call. = FALSE
    )
  }

  stopifnot(is.data.frame(sample_qc))
  qc <- as.data.frame(sample_qc, stringsAsFactors = FALSE)
  sample_col <- .match_sample_qc_col(
    qc,
    sample_col,
    c("sample_name", "sample", "Sample_name", "Sample", "sample_id", "SampleID"),
    "sample-name"
  )
  qc_names <- as.character(qc[[sample_col]])
  if (anyNA(qc_names) || any(!nzchar(qc_names))) {
    stop("sample_qc sample names must be non-missing and non-empty.",
         call. = FALSE)
  }
  if (anyDuplicated(qc_names)) {
    dup <- unique(qc_names[duplicated(qc_names)])
    stop(
      "Duplicate sample names in sample_qc: ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }

  sample_names <- control_names(control_group)
  missing_qc <- setdiff(sample_names, qc_names)
  if (length(missing_qc)) {
    stop(
      "sample_qc is missing rows for control samples: ",
      paste(missing_qc, collapse = ", "),
      call. = FALSE
    )
  }

  qc <- qc[match(sample_names, qc_names), , drop = FALSE]
  rownames(qc) <- sample_names
  reason <- rep(NA_character_, length(sample_names))

  reads_vals <- rep(NA_real_, length(sample_names))
  if (!is.null(min_total_unique_reads) || !is.null(max_total_unique_reads)) {
    total_unique_reads_col <- .match_sample_qc_col(
      qc,
      total_unique_reads_col,
      c("read_counts_binned_post_sum"),
      "total-unique-reads"
    )
    reads_vals <- suppressWarnings(as.numeric(qc[[total_unique_reads_col]]))
    if (any(!is.finite(reads_vals))) {
      bad <- sample_names[!is.finite(reads_vals)]
      stop(
        "Non-numeric or missing total unique reads for samples: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }

    if (!is.null(min_total_unique_reads)) {
      below <- reads_vals < as.numeric(min_total_unique_reads)
      reason[below] <- .append_filter_reason(reason[below], "below_min_unique_reads")
    }
    if (!is.null(max_total_unique_reads)) {
      above <- reads_vals > as.numeric(max_total_unique_reads)
      reason[above] <- .append_filter_reason(reason[above], "above_max_unique_reads")
    }
  }

  gc_vals <- rep(NA_real_, length(sample_names))
  if (!is.null(gc_mad_cutoff)) {
    gc_col <- .match_sample_qc_col(
      qc,
      gc_col,
      c("gc_read_perc_post"),
      "GC"
    )
    gc_vals <- suppressWarnings(as.numeric(qc[[gc_col]]))
    if (any(!is.finite(gc_vals))) {
      bad <- sample_names[!is.finite(gc_vals)]
      stop(
        "Non-numeric or missing GC values for samples: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
    gc_outliers <- .flag_mad_outliers(
      gc_vals,
      cutoff = as.numeric(gc_mad_cutoff),
      label = "GC"
    )
    reason[gc_outliers] <- .append_filter_reason(reason[gc_outliers], "gc_mad_outlier")
  }

  excluded <- sample_names[!is.na(reason) & nzchar(reason)]
  filtered_group <- if (length(excluded)) {
    nipter_drop_control_group_samples(control_group, excluded)
  } else {
    control_group
  }

  list(
    control_group = filtered_group,
    retained_samples = control_names(filtered_group),
    excluded_samples = excluded,
    exclusion_table = data.frame(
      sample_name = excluded,
      exclusion_reason = reason[match(excluded, sample_names)],
      total_unique_reads = reads_vals[match(excluded, sample_names)],
      gc = gc_vals[match(excluded, sample_names)],
      stringsAsFactors = FALSE
    ),
    matched_qc = qc,
    settings = list(
      sample_col = sample_col,
      total_unique_reads_col = if (exists("total_unique_reads_col")) total_unique_reads_col else NULL,
      gc_col = if (exists("gc_col")) gc_col else NULL,
      min_total_unique_reads = min_total_unique_reads,
      max_total_unique_reads = max_total_unique_reads,
      gc_mad_cutoff = gc_mad_cutoff
    )
  )
}

.flag_control_group_outliers <- function(diagnostics,
                                         collapse_strands = TRUE,
                                         max_aberrant_chromosomes = 2L,
                                         outlier_rule = c("any_aberrant_score", "bidirectional_or_multichromosome")) {
  stopifnot(is.list(diagnostics))
  stopifnot(is.logical(collapse_strands), length(collapse_strands) == 1L)
  stopifnot(is.numeric(max_aberrant_chromosomes),
            length(max_aberrant_chromosomes) == 1L,
            max_aberrant_chromosomes >= 0)
  outlier_rule <- match.arg(outlier_rule)

  ab <- diagnostics$aberrant_scores
  if (is.null(ab) || !nrow(ab)) {
    return(character())
  }

  if (identical(outlier_rule, "any_aberrant_score")) {
    return(unique(as.character(ab$sample_name)))
  }

  if (isTRUE(collapse_strands)) {
    chrom_summary <- stats::aggregate(
      chromosome ~ sample_name,
      data = ab,
      FUN = function(x) length(unique(x))
    )
    return(chrom_summary$sample_name[
      chrom_summary$chromosome > as.integer(max_aberrant_chromosomes)
    ])
  }

  ab$chromosome_core <- sub("[FR]$", "", ab$chromosome)
  sample_split <- split(ab, ab$sample_name)

  flagged <- vapply(sample_split, function(df) {
    per_chr <- table(df$chromosome_core)
    any(per_chr > 1L) ||
      length(unique(df$chromosome_core)) > as.integer(max_aberrant_chromosomes)
  }, logical(1L))

  names(flagged)[flagged]
}


#' Iteratively prune aberrant controls using chi-correction diagnostics
#'
#' Mirrors the legacy production cleaning loop for NIPTeR control cohorts:
#' chi-correct the current control group, diagnose chromosomal-fraction
#' outliers, drop flagged samples, and repeat until no further outliers remain
#' or the minimum control size would be violated.
#'
#' For \code{SeparatedStrands} control groups with
#' \code{collapse_strands = FALSE} (default), the pruning rule matches the
#' legacy script: drop a sample if both strands of any chromosome are
#' aberrant, or if more than \code{max_aberrant_chromosomes} distinct
#' chromosomes are aberrant.
#'
#' @param control_group A \code{NIPTeRControlGroup}.
#' @param sample Optional \code{NIPTeRSample} used as the dummy sample passed
#'   into [nipter_chi_correct()]. Defaults to the first control sample.
#' @param chi_cutoff Numeric chi-correction cutoff passed to
#'   [nipter_chi_correct()]. Default \code{3.5}.
#' @param collapse_strands Logical; keep 44 strand-resolved diagnostics for
#'   \code{SeparatedStrands} (\code{FALSE}, default) or collapse to 22
#'   autosomal fractions.
#' @param z_cutoff Absolute Z-score cutoff passed to
#'   [nipter_diagnose_control_group()].
#' @param max_aberrant_chromosomes Maximum number of distinct aberrant
#'   chromosomes allowed before a sample is dropped. Default \code{2L}.
#' @param outlier_rule Character scalar controlling how aberrant samples are
#'   dropped. \code{"any_aberrant_score"} (default) removes any sample
#'   appearing in \code{abberant_scores}, matching the upstream NIPTeR
#'   vignette. \code{"bidirectional_or_multichromosome"} drops a sample only
#'   when both strands of one chromosome are aberrant or when more than
#'   \code{max_aberrant_chromosomes} distinct chromosomes are aberrant.
#' @param min_controls Minimum allowed retained control count. The iteration
#'   stops before dropping samples if doing so would leave fewer controls than
#'   this threshold.
#' @param max_iterations Maximum number of pruning iterations. Default
#'   \code{100L}.
#' @param verbose Logical; emit per-iteration messages.
#'
#' @return A named list with:
#' \describe{
#'   \item{control_group}{The final uncorrected control group after dropping
#'     flagged samples.}
#'   \item{chi_corrected_control_group}{The final chi-corrected control group
#'     used for the terminal diagnostic pass.}
#'   \item{diagnostics}{The terminal diagnostic list from
#'     [nipter_diagnose_control_group()].}
#'   \item{dropped_samples}{Character vector of unique dropped sample names in
#'     drop order.}
#'   \item{iteration_log}{Data frame recording flagged and retained control
#'     counts per iteration.}
#'   \item{converged}{Logical; \code{TRUE} when the loop ended without flagged
#'     samples.}
#'   \item{stop_reason}{Text reason for the stopping condition.}
#' }
#'
#' @seealso [nipter_chi_correct()], [nipter_diagnose_control_group()],
#'   [nipter_drop_control_group_samples()]
#'
#' @export
nipter_prune_control_group_outliers <- function(control_group,
                                                sample = NULL,
                                                chi_cutoff = 3.5,
                                                collapse_strands = FALSE,
                                                z_cutoff = 3,
                                                max_aberrant_chromosomes = 2L,
                                                outlier_rule = c("any_aberrant_score", "bidirectional_or_multichromosome"),
                                                min_controls = 10L,
                                                max_iterations = 100L,
                                                verbose = FALSE) {
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.null(sample) || .is_nipt_sample_object(sample))
  stopifnot(is.numeric(chi_cutoff), length(chi_cutoff) == 1L)
  stopifnot(is.logical(collapse_strands), length(collapse_strands) == 1L)
  stopifnot(is.numeric(z_cutoff), length(z_cutoff) == 1L, z_cutoff > 0)
  stopifnot(is.numeric(max_aberrant_chromosomes),
            length(max_aberrant_chromosomes) == 1L,
            max_aberrant_chromosomes >= 0)
  outlier_rule <- match.arg(outlier_rule)
  stopifnot(is.numeric(min_controls), length(min_controls) == 1L, min_controls >= 2L)
  stopifnot(is.numeric(max_iterations), length(max_iterations) == 1L,
            max_iterations >= 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  sample <- if (is.null(sample)) control_group$samples[[1L]] else sample
  if (!identical(.strand_type_of(sample), .strand_type_of(control_group))) {
    stop("sample and control_group must have the same strand type.",
         call. = FALSE)
  }

  current_raw <- control_group
  dropped <- character()
  iter_log <- vector("list", as.integer(max_iterations))
  stop_reason <- "max_iterations"
  converged <- FALSE

  for (iter in seq_len(as.integer(max_iterations))) {
    current_chi <- nipter_chi_correct(
      sample = sample,
      control_group = current_raw,
      chi_cutoff = chi_cutoff
    )[["control_group"]]
    diag <- nipter_diagnose_control_group(
      current_chi,
      collapse_strands = collapse_strands,
      z_cutoff = z_cutoff
    )
    flagged_chi <- .flag_control_group_outliers(
      diag,
      collapse_strands = collapse_strands,
      max_aberrant_chromosomes = max_aberrant_chromosomes,
      outlier_rule = outlier_rule
    )
    flagged <- flagged_chi

    iter_log[[iter]] <- data.frame(
      iteration = iter,
      n_controls = n_controls(current_raw),
      n_flagged_post_chi = length(flagged_chi),
      flagged_samples_post_chi = paste(flagged_chi, collapse = ";"),
      n_flagged = length(flagged),
      flagged_samples = paste(flagged, collapse = ";"),
      stringsAsFactors = FALSE
    )

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "nipter_prune_control_group_outliers iter %d: %d controls, %d flagged post-chi",
          iter, n_controls(current_raw), length(flagged_chi)
        )
      )
    }

    if (!length(flagged)) {
      stop_reason <- "no_flagged_samples"
      converged <- TRUE
      iter_log <- iter_log[seq_len(iter)]
      return(list(
        control_group = current_raw,
        chi_corrected_control_group = current_chi,
        diagnostics = diag,
        dropped_samples = dropped,
        iteration_log = do.call(rbind, iter_log),
        converged = converged,
        stop_reason = stop_reason
      ))
    }

    if ((n_controls(current_raw) - length(flagged)) < min_controls) {
      stop_reason <- "min_controls"
      iter_log <- iter_log[seq_len(iter)]
      return(list(
        control_group = current_raw,
        chi_corrected_control_group = current_chi,
        diagnostics = diag,
        dropped_samples = dropped,
        iteration_log = do.call(rbind, iter_log),
        converged = converged,
        stop_reason = stop_reason
      ))
    }

    current_raw <- nipter_drop_control_group_samples(current_raw, flagged)
    sample <- current_raw$samples[[1L]]
    dropped <- c(dropped, flagged)
  }

  current_chi <- nipter_chi_correct(
    sample = sample,
    control_group = current_raw,
    chi_cutoff = chi_cutoff
  )[["control_group"]]
  diag <- nipter_diagnose_control_group(
    current_chi,
    collapse_strands = collapse_strands,
    z_cutoff = z_cutoff
  )

  list(
    control_group = current_raw,
    chi_corrected_control_group = current_chi,
    diagnostics = diag,
    dropped_samples = dropped,
    iteration_log = do.call(rbind, iter_log),
    converged = converged,
    stop_reason = stop_reason
  )
}


#' Select best-matching controls for a sample
#'
#' Ranks control samples by similarity of chromosomal fractions to a test
#' sample. Uses sum-of-squared differences of collapsed chromosomal fractions
#' (chromosomes 1-12, 14-17, 19-20, 22 by default — excluding trisomy
#' chromosomes 13, 18, 21). The distance computation is accelerated by an
#' Rcpp+OpenMP kernel.
#'
#' @param sample A \code{NIPTeRSample} object to match against.
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param n Integer; number of best-matching controls to return.
#' @param mode \code{"subset"} (default) returns a new
#'   \code{NIPTeRControlGroup}; \code{"report"} returns a named numeric vector
#'   of sum-of-squares scores.
#' @param exclude_chromosomes Integer vector of chromosomes to exclude from the
#'   distance calculation (default \code{c(13, 18, 21)}).
#' @param include_chromosomes Integer vector of chromosomes to include. If
#'   \code{NULL} (default), uses all autosomal chromosomes minus
#'   \code{exclude_chromosomes}.
#' @param cpus Integer; number of OpenMP threads for the SSD computation.
#'   Default \code{1L}.
#'
#' @return A \code{NIPTeRControlGroup} (if \code{mode = "subset"}) or a named
#'   numeric vector of sum-of-squares distances (if \code{mode = "report"}).
#'
#' @seealso [nipter_as_control_group()], [nipter_match_matrix()]
#'
#' @export
nipter_match_control_group <- function(sample,
                                       control_group,
                                       n,
                                       mode = c("subset", "report"),
                                       exclude_chromosomes = c(13L, 18L, 21L),
                                       include_chromosomes = NULL,
                                       cpus = 1L) {
  mode <- match.arg(mode)
  stopifnot(.is_nipt_sample_object(sample))
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1L)

  # Strand-type compatibility guard
  sample_st <- .strand_type_of(sample)
  cg_st     <- .strand_type_of(control_group)
  if (!identical(sample_st, cg_st)) {
    stop(sprintf(
      "Strand type mismatch: sample is '%s' but control_group is '%s'.",
      sample_st, cg_st
    ), call. = FALSE)
  }

  cpus <- as.integer(cpus)

  # Determine comparison chromosomes (0-based row indices for Rcpp)
  if (is.null(include_chromosomes)) {
    compare_chroms <- setdiff(1:22, exclude_chromosomes)
  } else {
    compare_chroms <- as.integer(include_chromosomes)
  }
  compare_idx <- as.integer(compare_chroms) - 1L  # 0-based for Rcpp

  # Pre-extract the full 22×N fractions matrix once
  fracs_mat  <- .control_group_fractions_collapsed(control_group)   # 22 × N
  names_vec  <- colnames(fracs_mat)
  query_frac <- .sample_chr_fractions_collapsed(sample)              # 22-element

  # Rcpp kernel: one query vs N controls (OpenMP-parallelized)
  scores <- nipter_ssd_scores_cpp(fracs_mat, query_frac, compare_idx, cpus)
  names(scores) <- names_vec

  # Sort ascending (most similar first)
  scores <- sort(scores)

  if (mode == "report") return(scores)

  # Subset mode: return top n as a new control group
  n <- min(n, length(scores))
  keep_names <- names(scores)[seq_len(n)]
  all_samples <- control_group$samples
  keep_idx   <- match(keep_names, vapply(all_samples, .sample_name, character(1L)))
  sname <- .sample_name(sample)
  keep_sex <- NULL
  if (.is_s7_nipt_control_group(control_group) &&
      !is.null(control_group$sample_sex)) {
    keep_sex <- control_group$sample_sex[keep_names]
  }
  nipter_as_control_group(
    all_samples[keep_idx],
    description = sprintf("Fitted to %s", sname),
    sample_sex = keep_sex,
    sex_source = if (.is_s7_nipt_control_group(control_group))
      control_group$sex_source else NULL
  )
}


#' Compute the full pairwise SSD matrix for a control group
#'
#' Returns the symmetric N×N matrix of sum-of-squared-differences between all
#' pairs of control samples' chromosomal fractions. The diagonal is zero.
#' Row means of this matrix are the per-sample "matching score" used for QC
#' in the production matching loop (see \emph{Details}).
#'
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param exclude_chromosomes Integer vector of chromosomes to exclude from the
#'   distance calculation (default \code{c(13, 18, 21)}).
#' @param include_chromosomes Integer vector of chromosomes to include. If
#'   \code{NULL} (default), uses all autosomal chromosomes minus
#'   \code{exclude_chromosomes}.
#' @param cpus Integer; OpenMP threads. Default \code{1L}.
#'
#' @return A numeric N×N matrix with sample names as row and column names.
#'
#' @details
#' The production NIPT pipeline uses this matrix to identify outlier controls
#' before scoring: each sample's mean SSD against all others is computed, and
#' samples with mean SSD more than 3 SD above the group mean are iteratively
#' removed. This function replaces the \code{lapply} over \code{match_control_group}
#' calls in \code{CoverageProjectionSCA_Reports.R} with a single vectorized
#' Rcpp kernel.
#'
#' @seealso [nipter_match_control_group()]
#'
#' @export
nipter_match_matrix <- function(control_group,
                                exclude_chromosomes = c(13L, 18L, 21L),
                                include_chromosomes = NULL,
                                cpus = 1L) {
  stopifnot(.is_nipt_control_group_object(control_group))
  cpus <- as.integer(cpus)

  if (is.null(include_chromosomes)) {
    compare_chroms <- setdiff(1:22, exclude_chromosomes)
  } else {
    compare_chroms <- as.integer(include_chromosomes)
  }
  compare_idx <- as.integer(compare_chroms) - 1L  # 0-based for Rcpp

  fracs_mat <- .control_group_fractions_collapsed(control_group)  # 22 × N
  ssd_mat   <- nipter_ssd_matrix_cpp(fracs_mat, compare_idx, cpus)

  nm <- control_names(control_group)
  rownames(ssd_mat) <- nm
  colnames(ssd_mat) <- nm
  ssd_mat
}


#' Load a NIPTeR control group from a directory of TSV.bgz files
#'
#' Reads all \code{.bed.gz} (or \code{.tsv.bgz}) files in \code{bed_dir} using
#' \code{rduckhts_tabix_multi()} — a single multi-file DuckDB scan — and
#' constructs a \code{NIPTeRControlGroup} from the results. This is much faster
#' than \code{lapply(files, bed_to_nipter_sample)} for large cohorts because all
#' files are read in one pass.
#'
#' @param bed_dir Character; directory containing one \code{.bed.gz} or
#'   \code{.tsv.bgz} file per control sample, each produced by
#'   \code{\link{nipter_bin_bam_bed}}.
#' @param pattern Glob pattern for file discovery (default
#'   \code{"*.bed.gz"}).
#' @param binsize Optional integer; bin size in base pairs. If \code{NULL}
#'   (default), inferred from the first row of the first file.
#' @param autosomal_source Which autosomal counts to realize in each imported
#'   sample. `"auto"` (default) uses corrected autosomal columns when present,
#'   otherwise raw counts. `"raw"` always uses raw count columns.
#'   `"corrected"` requires corrected columns.
#' @param sex_counts Which sex-chromosome counts to realize in each imported
#'   sample. `"match"` (default) follows `autosomal_source`. `"raw"` always
#'   uses raw count columns. `"corrected"` requires corrected sex columns.
#' @param description Label for the resulting control group (default
#'   \code{"General control group"}).
#' @param sample_sex Optional character vector of known sex labels for the
#'   samples in \code{bed_dir}. Names must match the inferred sample names when
#'   supplied as a named vector.
#' @param sex_source Optional string describing the provenance of
#'   \code{sample_sex}.
#' @param con Optional open DBI connection with duckhts loaded.
#'
#' @details
#' The column count (5 or 9) is detected automatically from the first file:
#' 5-column BEDs produce a \code{CombinedStrands} control group;
#' 9-column BEDs (written by \code{nipter_bin_bam_bed(separate_strands = TRUE)})
#' produce a \code{SeparatedStrands} control group. All files in the directory
#' must have the same column count. When corrected BED columns are present,
#' `autosomal_source` and `sex_counts` control whether the imported samples use
#' raw or corrected values for each compartment before constructing the control
#' group.
#'
#' @return A \code{NIPTeRControlGroup}.
#'
#' @seealso [nipter_bin_bam_bed()], [nipter_as_control_group()],
#'   [bed_to_nipter_sample()]
#'
#' @examples
#' \dontrun{
#' # Bin all controls to BED once
#' for (bam in bam_files) {
#'   nipter_bin_bam_bed(bam, file.path("controls/", sub(".bam$", ".bed.gz", basename(bam))))
#' }
#' # Load them all at once
#' cg <- nipter_control_group_from_beds("controls/")
#' }
#'
#' @export
nipter_control_group_from_beds <- function(bed_dir,
                                           pattern     = "*.bed.gz",
                                           binsize     = NULL,
                                           autosomal_source = c("auto", "raw", "corrected"),
                                           sex_counts  = c("match", "raw", "corrected"),
                                           description = "General control group",
                                           sample_sex  = NULL,
                                           sex_source  = NULL,
                                           con         = NULL) {
  stopifnot(is.character(bed_dir), length(bed_dir) == 1L,
            nzchar(bed_dir), dir.exists(bed_dir))
  autosomal_source <- match.arg(autosomal_source)
  sex_counts <- match.arg(sex_counts)

  files <- Sys.glob(file.path(bed_dir, pattern))
  if (length(files) == 0L) {
    stop("No files matching '", pattern, "' found in '", bed_dir, "'.",
         call. = FALSE)
  }

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  # Register all files as a single multi-reader DuckDB table.
  # rduckhts_tabix_multi() builds UNION ALL BY NAME of per-file read_tabix()
  # queries and adds a `filename` column identifying the source.  The tabix
  # reader returns all data columns as VARCHAR, so we cast in the SELECT below
  # via TRY_CAST (handles literal "NA" written by write.table()).
  tbl <- paste0("nipter_cg_beds_",
                as.hexmode(sample.int(.Machine$integer.max, 1L)))
  Rduckhts::rduckhts_tabix_multi(con, tbl, files, overwrite = TRUE)
  on.exit(
    tryCatch(DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS \"%s\"", tbl)),
             error = function(e) NULL),
    add = TRUE
  )

  # Auto-detect strand type: 9-column BEDs (SeparatedStrands) have column8;
  # 5-column BEDs (CombinedStrands) do not.
  probe_sql <- sprintf('SELECT column8 FROM "%s" LIMIT 1', tbl)
  is_separated <- tryCatch({
    probe <- DBI::dbGetQuery(con, probe_sql)
    nrow(probe) > 0L && !is.na(probe[[1L]][1L]) && nzchar(probe[[1L]][1L])
  }, error = function(e) FALSE)

  # Row order across files is non-deterministic (UNION ALL does not guarantee
  # order). Matrix assignment uses bin_idx derived from start_pos, so row order
  # within a sample does not affect correctness.
  if (is_separated) {
    rows <- DBI::dbGetQuery(con, sprintf(
      'SELECT
         filename,
         column0                       AS chrom,
         CAST(column1 AS INTEGER)      AS start_pos,
         CAST(column2 AS INTEGER)      AS end_pos,
         CAST(column3 AS INTEGER)      AS count,
         CAST(column4 AS INTEGER)      AS count_fwd,
         CAST(column5 AS INTEGER)      AS count_rev,
         TRY_CAST(column6 AS DOUBLE)   AS corrected_count,
         TRY_CAST(column7 AS DOUBLE)   AS corrected_fwd,
         TRY_CAST(column8 AS DOUBLE)   AS corrected_rev
       FROM "%s"', tbl
    ))
  } else {
    rows <- DBI::dbGetQuery(con, sprintf(
      'SELECT
         filename,
         column0                       AS chrom,
         CAST(column1 AS INTEGER)      AS start_pos,
         CAST(column2 AS INTEGER)      AS end_pos,
         CAST(column3 AS INTEGER)      AS count,
         TRY_CAST(column4 AS DOUBLE)   AS corrected_count
       FROM "%s"', tbl
    ))
  }

  if (nrow(rows) == 0L) {
    stop("All BED files in '", bed_dir, "' appear to be empty.", call. = FALSE)
  }

  binsize <- .infer_bed_binsize(rows, binsize, bed_dir)

  # Derive sample name from file path (strip directory + extension)
  rows$sample_name <- sub("\\.bed(\\.gz)?$|\\.tsv(\\.bgz)?$", "",
                          basename(rows$filename), ignore.case = TRUE)

  # Compute the global maximum bin index so all sample matrices share the same
  # column width.  Narrower chromosomes are zero-padded on the right.
  bin_idx_all <- as.integer(rows$start_pos / binsize)
  n_bins_global <- max(bin_idx_all) + 1L

  sample_names <- unique(rows$sample_name)

  if (is_separated) {
    samples <- lapply(sample_names, function(nm) {
      sub_rows <- rows[rows$sample_name == nm, , drop = FALSE]
      .rows_to_nipter_sep(
        sub_rows,
        nm,
        binsize,
        n_bins = n_bins_global,
        autosomal_source = autosomal_source,
        sex_source = sex_counts
      )
    })
  } else {
    samples <- lapply(sample_names, function(nm) {
      sub_rows <- rows[rows$sample_name == nm, , drop = FALSE]
      .rows_to_nipter_combined(
        sub_rows,
        nm,
        binsize,
        n_bins = n_bins_global,
        autosomal_source = autosomal_source,
        sex_source = sex_counts
      )
    })
  }

  nipter_as_control_group(
    samples,
    description = description,
    sample_sex = sample_sex,
    sex_source = sex_source
  )
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Compute chromosomal fractions for one NIPTeRSample.
#
# CombinedStrands: returns a named 22-element numeric vector (names "1"–"22").
#   Each element = rowSums(auto)[chr] / sum(auto).
#
# SeparatedStrands: returns a named 44-element numeric vector
#   (names "1F".."22F","1R".."22R").  Each element =
#   rowSums(strand_mat)[row] / sum(strand_mat) / 2.
#   This matches NIPTeR's chrfractions.SeparatedStrands, which computes
#   sapply(auto_list, function(x) (rowSums(x) / sum(x)) / 2)  and then
#   flattens the 22×2 result into a 44-vector via sapply over samples.
.sample_chr_fractions <- function(sample) {
  if (.strand_type_of(sample) == "separated") {
    # SeparatedStrands: per-strand fractions, each strand normalised to its own
    # total, then divided by 2 (matches NIPTeR's chrfractions.SeparatedStrands).
    auto_list <- .sample_autosomal_reads(sample)
    fracs <- unlist(lapply(auto_list, function(mat) {
      s <- sum(mat)
      if (s == 0) return(stats::setNames(rep(0, nrow(mat)), rownames(mat)))
      (rowSums(mat) / s) / 2
    }))
    return(fracs)
  }
  # CombinedStrands: use autosomal_matrix() which works for both S3 and S7
  auto <- autosomal_matrix(sample)
  chr_sums <- rowSums(auto)
  total <- sum(chr_sums)
  if (total == 0) return(stats::setNames(rep(0, 22L), as.character(1:22)))
  chr_sums / total
}

# Compute collapsed (22-element) chromosomal fractions for any NIPTeRSample.
# For SeparatedStrands, sums forward + reverse per chromosome (matching
# NIPTeR's retrieve_fractions_of_interest.SeparatedStrands which sums
# frac[paste0(chr,"F"),] + frac[paste0(chr,"R"),]).
# Used by nipter_z_score(), nipter_ncv_score(), nipter_match_control_group().
.sample_chr_fractions_collapsed <- function(sample) {
  if (.strand_type_of(sample) == "separated") {
    frac44 <- .sample_chr_fractions(sample)
    keys_f <- paste0(1:22, "F")
    keys_r <- paste0(1:22, "R")
    collapsed <- frac44[keys_f] + frac44[keys_r]
    names(collapsed) <- as.character(1:22)
    return(collapsed)
  }
  .sample_chr_fractions(sample)
}

# Compute a fractions matrix for a control group.
# CombinedStrands: 22 x n_samples. SeparatedStrands: 44 x n_samples.
.control_group_fractions <- function(control_group) {
  if (.is_s7_nipt_control_group(control_group)) {
    return(fractions_for_regression(control_group))
  }
  frac_list <- lapply(control_group$samples, .sample_chr_fractions)
  mat <- do.call(cbind, frac_list)
  colnames(mat) <- vapply(control_group$samples,
                          function(s) s$sample_name, character(1L))
  mat
}

# Compute collapsed 22 x n_samples fractions matrix for any control group.
.control_group_fractions_collapsed <- function(control_group) {
  if (.is_s7_nipt_control_group(control_group)) {
    return(fractions_auto(control_group))
  }
  frac_list <- lapply(control_group$samples, .sample_chr_fractions_collapsed)
  mat <- do.call(cbind, frac_list)
  colnames(mat) <- vapply(control_group$samples,
                          function(s) s$sample_name, character(1L))
  mat
}

# Compute total reads per chromosome (22-element vector, always collapsed)
# for a NIPTeRSample. Uses autosomal reads only.
# For SeparatedStrands, sums forward + reverse via autosomal_matrix().
.sample_chr_reads <- function(sample) {
  mat <- autosomal_matrix(sample)
  stats::setNames(rowSums(mat), rownames(mat))
}

.sample_reference_frame_row <- function(sample) {
  chroms <- c(as.character(1:22), "X", "Y")
  auto_counts <- .sample_chr_reads(sample)
  sex_counts <- stats::setNames(rowSums(sex_matrix(sample)), c("X", "Y"))
  chr_counts <- c(auto_counts[as.character(1:22)], sex_counts[c("X", "Y")])
  auto_total <- sum(auto_counts)
  chr_fracs <- chr_counts / max(auto_total, .Machine$double.eps)

  row <- c(
    Sample_name = .sample_name(sample),
    stats::setNames(as.list(as.numeric(chr_counts)),
                    paste0("NChrReads_", chroms)),
    stats::setNames(as.list(as.numeric(chr_fracs)),
                    paste0("FrChrReads_", chroms))
  )
  as.data.frame(row, stringsAsFactors = FALSE)
}


#' Build a chromosome-level reference frame from a NIPTeR control group
#'
#' Produces the per-sample count and fraction table that downstream sex-aware
#' NIPT model building actually needs. This keeps application-level gaunosome
#' modelling data out of the raw \code{NIPTeRSample} class while providing a
#' stable training frame for future sex-chromosome Z-score, NCV, and
#' regression models.
#'
#' @param control_group A \code{NIPTeRControlGroup} object.
#' @param sample_sex Optional character vector overriding the sex labels stored
#'   on \code{control_group}. Accepted values are \code{"female"},
#'   \code{"male"}, \code{"ambiguous"}, and \code{"unknown"}.
#'
#' @return A typed \code{NIPTReferenceFrame} data frame with one row per sample
#'   and columns: \code{Sample_name}, optional \code{SampleSex},
#'   \code{NChrReads_*}, and \code{FrChrReads_*} for chromosomes \code{1:22},
#'   \code{X}, and \code{Y}.
#'
#' @export
nipter_reference_frame <- function(control_group, sample_sex = NULL) {
  stopifnot(.is_nipt_control_group_object(control_group))

  samples <- control_group$samples
  sample_names <- vapply(samples, .sample_name, character(1L))
  sample_sex <- .normalize_sample_sex(
    if (is.null(sample_sex)) control_group$sample_sex else sample_sex,
    samples
  )

  chroms <- c(as.character(1:22), "X", "Y")
  rows <- lapply(samples, .sample_reference_frame_row)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  count_cols <- paste0("NChrReads_", chroms)
  frac_cols  <- paste0("FrChrReads_", chroms)
  for (col in count_cols) out[[col]] <- as.numeric(out[[col]])
  for (col in frac_cols) out[[col]] <- as.numeric(out[[col]])

  if (!is.null(sample_sex)) {
    out$SampleSex <- unname(sample_sex[out$Sample_name])
    out <- out[, c("Sample_name", "SampleSex", count_cols, frac_cols),
               drop = FALSE]
  } else {
    out <- out[, c("Sample_name", count_cols, frac_cols), drop = FALSE]
  }

  .as_nipt_reference_frame(out)
}
