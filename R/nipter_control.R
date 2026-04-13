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
                                    description = "General control group") {
  stopifnot(is.list(samples), length(samples) >= 2L)
  stopifnot(is.character(description), length(description) == 1L,
            nzchar(description))

  # Validate that all entries are NIPTeRSample

  ok <- vapply(samples, inherits, logical(1L), what = "NIPTeRSample")
  if (!all(ok)) {
    stop("All elements of 'samples' must be NIPTeRSample objects.", call. = FALSE)
  }

  # Validate same strand type
  strand_types <- vapply(samples, function(s) class(s)[2L], character(1L))
  if (length(unique(strand_types)) > 1L) {
    stop(
      "All samples must have the same strand type. Found: ",
      paste(unique(strand_types), collapse = ", "),
      call. = FALSE
    )
  }
  strand_type <- strand_types[1L]

  # Remove duplicate sample names (keep first occurrence)
  names_vec <- vapply(samples, function(s) s$sample_name, character(1L))
  dups <- duplicated(names_vec)
  if (any(dups)) {
    message("Removing ", sum(dups), " duplicate sample(s) by name.")
    samples <- samples[!dups]
  }

  # Collect correction statuses
  auto_status <- unique(unlist(lapply(samples, `[[`,
                                      "correction_status_autosomal")))
  sex_status  <- unique(unlist(lapply(samples, `[[`,
                                      "correction_status_sex")))

  structure(
    list(
      samples                       = samples,
      correction_status_autosomal   = auto_status,
      correction_status_sex         = sex_status,
      description                   = description
    ),
    class = c("NIPTeRControlGroup", strand_type)
  )
}


#' Diagnose a NIPTeR control group
#'
#' Computes per-chromosome Z-scores across all samples in a control group and
#' flags outliers (|Z| > 3). The Shapiro-Wilk test is applied to each
#' chromosome's fraction distribution.
#'
#' @param control_group A \code{NIPTeRControlGroup} object.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{z_scores}{Chromosome-by-sample matrix of Z-scores (rows =
#'     chromosomes 1-22, columns = samples).}
#'   \item{aberrant_scores}{A \code{data.frame} with columns
#'     \code{chromosome}, \code{sample_name}, \code{z_score} for all
#'     |Z| > 3, or \code{NULL} if none.}
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
nipter_diagnose_control_group <- function(control_group) {
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))

  fracs <- .control_group_fractions_collapsed(control_group)
  # fracs: 22 x n_samples matrix

  # Z-score each chromosome across samples
  z_mat <- t(apply(fracs, 1L, scale))
  sample_names <- vapply(control_group$samples, function(s) s$sample_name,
                         character(1L))
  colnames(z_mat) <- sample_names
  rownames(z_mat) <- as.character(1:22)

  # Find aberrant |Z| > 3
  aberrant <- which(abs(z_mat) > 3, arr.ind = TRUE)
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
  rownames(stats_mat) <- as.character(1:22)


  list(
    z_scores        = z_mat,
    aberrant_scores = ab_df,
    statistics      = stats_mat
  )
}


#' Select best-matching controls for a sample
#'
#' Ranks control samples by similarity of chromosomal fractions to a test
#' sample. Uses sum-of-squared differences of control-chromosome fractions
#' (chromosomes 1-12, 14-17, 19-20, 22 by default — excluding trisomy
#' chromosomes 13, 18, 21).
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
#'
#' @return A \code{NIPTeRControlGroup} (if \code{mode = "subset"}) or a named
#'   numeric vector of sum-of-squares distances (if \code{mode = "report"}).
#'
#' @seealso [nipter_as_control_group()]
#'
#' @export
nipter_match_control_group <- function(sample,
                                       control_group,
                                       n,
                                       mode = c("subset", "report"),
                                       exclude_chromosomes = c(13L, 18L, 21L),
                                       include_chromosomes = NULL) {
  mode <- match.arg(mode)
  stopifnot(inherits(sample, "NIPTeRSample"))
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1L)

  # Determine comparison chromosomes
  if (is.null(include_chromosomes)) {
    compare_chroms <- setdiff(1:22, exclude_chromosomes)
  } else {
    compare_chroms <- as.integer(include_chromosomes)
  }
  compare_keys <- as.character(compare_chroms)

  # Get sample fractions (collapsed to 22-row for SeparatedStrands)
  sample_frac <- .sample_chr_fractions_collapsed(sample)[compare_keys]

  # Get control fractions and compute SSD
  n_controls <- length(control_group$samples)
  names_vec  <- character(n_controls)
  scores     <- numeric(n_controls)

  for (i in seq_len(n_controls)) {
    ctrl <- control_group$samples[[i]]
    names_vec[i] <- ctrl$sample_name
    ctrl_frac    <- .sample_chr_fractions_collapsed(ctrl)[compare_keys]
    scores[i]    <- sum((sample_frac - ctrl_frac)^2)
  }
  names(scores) <- names_vec

  # Sort ascending (most similar first)
  scores <- sort(scores)

  if (mode == "report") return(scores)

  # Subset mode: return top n as a new control group
  n <- min(n, length(scores))
  keep_names <- names(scores)[seq_len(n)]
  keep_idx   <- match(keep_names, vapply(control_group$samples,
                                         function(s) s$sample_name,
                                         character(1L)))
  nipter_as_control_group(
    control_group$samples[keep_idx],
    description = sprintf("Fitted to %s", sample$sample_name)
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
  if (inherits(sample, "SeparatedStrands")) {
    auto <- sample$autosomal_chromosome_reads  # list of 2 matrices (fwd, rev)
    # Each matrix: rows = "1F".."22F" (or "1R".."22R"), cols = bins
    fracs <- unlist(lapply(auto, function(mat) {
      s <- sum(mat)
      if (s == 0) return(stats::setNames(rep(0, nrow(mat)), rownames(mat)))
      (rowSums(mat) / s) / 2
    }))
    return(fracs)
  }
  # CombinedStrands
  auto <- sample$autosomal_chromosome_reads[[1L]]
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
  if (inherits(sample, "SeparatedStrands")) {
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
  n <- length(control_group$samples)
  frac_list <- lapply(control_group$samples, .sample_chr_fractions)
  mat <- do.call(cbind, frac_list)
  colnames(mat) <- vapply(control_group$samples,
                          function(s) s$sample_name, character(1L))
  mat
}

# Compute collapsed 22 x n_samples fractions matrix for any control group.
.control_group_fractions_collapsed <- function(control_group) {
  n <- length(control_group$samples)
  frac_list <- lapply(control_group$samples, .sample_chr_fractions_collapsed)
  mat <- do.call(cbind, frac_list)
  colnames(mat) <- vapply(control_group$samples,
                          function(s) s$sample_name, character(1L))
  mat
}

# Compute total reads per chromosome (22-element vector, always collapsed)
# for a NIPTeRSample. Uses autosomal reads only.
# For SeparatedStrands, sums forward + reverse via Reduce("+", auto).
.sample_chr_reads <- function(sample) {
  auto <- sample$autosomal_chromosome_reads
  if (inherits(sample, "SeparatedStrands")) {
    summed <- Reduce("+", auto)
    rownames(summed) <- as.character(1:22)
    return(stats::setNames(rowSums(summed), rownames(summed)))
  }
  stats::setNames(rowSums(auto[[1L]]), rownames(auto[[1L]]))
}
