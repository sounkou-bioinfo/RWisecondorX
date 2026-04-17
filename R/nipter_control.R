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
  stopifnot(.is_nipt_control_group_object(control_group))

  fracs <- .control_group_fractions_collapsed(control_group)
  # fracs: 22 x n_samples matrix

  # Z-score each chromosome across samples
  z_mat <- t(apply(fracs, 1L, scale))
  sample_names <- control_names(control_group)
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
