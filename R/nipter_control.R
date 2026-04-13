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
  stopifnot(inherits(sample, "NIPTeRSample"))
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1L)
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
  keep_idx   <- match(keep_names, vapply(control_group$samples,
                                         function(s) s$sample_name,
                                         character(1L)))
  nipter_as_control_group(
    control_group$samples[keep_idx],
    description = sprintf("Fitted to %s", sample$sample_name)
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
  stopifnot(inherits(control_group, "NIPTeRControlGroup"))
  cpus <- as.integer(cpus)

  if (is.null(include_chromosomes)) {
    compare_chroms <- setdiff(1:22, exclude_chromosomes)
  } else {
    compare_chroms <- as.integer(include_chromosomes)
  }
  compare_idx <- as.integer(compare_chroms) - 1L  # 0-based for Rcpp

  fracs_mat <- .control_group_fractions_collapsed(control_group)  # 22 × N
  ssd_mat   <- nipter_ssd_matrix_cpp(fracs_mat, compare_idx, cpus)

  nm <- vapply(control_group$samples, function(s) s$sample_name, character(1L))
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
#' @param description Label for the resulting control group (default
#'   \code{"General control group"}).
#' @param con Optional open DBI connection with duckhts loaded.
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
                                           description = "General control group",
                                           con         = NULL) {
  stopifnot(is.character(bed_dir), length(bed_dir) == 1L,
            nzchar(bed_dir), dir.exists(bed_dir))

  files <- Sys.glob(file.path(bed_dir, pattern))
  if (length(files) == 0L) {
    stop("No files matching '", pattern, "' found in '", bed_dir, "'.",
         call. = FALSE)
  }

  own_con <- is.null(con)
  if (own_con) {
    if (!requireNamespace("Rduckhts", quietly = TRUE)) {
      stop("Rduckhts is required.", call. = FALSE)
    }
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  # Register all files as a single multi-reader DuckDB table.
  # rduckhts_tabix_multi() builds UNION ALL BY NAME of per-file read_tabix()
  # queries and adds a `filename` column identifying the source.  The tabix
  # reader returns all data columns as VARCHAR (column0..column4 for a
  # 5-column BED), so we must NOT pass column_types here — the "NA" strings
  # written by write.table() would fail a DOUBLE cast.  Cast in the SELECT
  # below with TRY_CAST, mirroring the pattern in bed_reader.R.
  tbl <- paste0("nipter_cg_beds_", as.integer(proc.time()[[3L]] * 1e6))
  Rduckhts::rduckhts_tabix_multi(con, tbl, files, overwrite = TRUE)
  on.exit(
    tryCatch(DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS \"%s\"", tbl)),
             error = function(e) NULL),
    add = TRUE
  )

  # Columns from read_tabix() without header_names: column0..column4 (VARCHAR).
  # column0 = chrom, column1 = start, column2 = end,
  # column3 = count, column4 = corrected_count ("NA" for uncorrected bins).
  #
  # Row order across files is non-deterministic (UNION ALL does not guarantee
  # order). This is safe because matrix assignment uses position-indexed access:
  #   auto_mat[chr_row, bin_idx[sel] + 1L] <- count
  # Each bin_idx is derived from start_pos, not from the position of the row
  # in the data frame. Inter-file row interleaving does not affect correctness.
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

  if (nrow(rows) == 0L) {
    stop("All BED files in '", bed_dir, "' appear to be empty.", call. = FALSE)
  }

  # Infer binsize from first row if not supplied
  if (is.null(binsize)) {
    binsize <- as.integer(rows$end_pos[1L] - rows$start_pos[1L])
  }
  binsize <- as.integer(binsize)

  # Derive sample name from file path (strip directory + extension)
  rows$sample_name <- sub("\\.bed(\\.gz)?$|\\.tsv(\\.bgz)?$", "",
                          basename(rows$filename), ignore.case = TRUE)

  # Compute the global maximum bin index so all sample matrices share the same
  # column width.  Narrower chromosomes are zero-padded on the right.
  bin_idx_all <- as.integer(rows$start_pos / binsize)
  n_bins_global <- max(bin_idx_all) + 1L

  sample_names <- unique(rows$sample_name)
  samples <- lapply(sample_names, function(nm) {
    sub_rows <- rows[rows$sample_name == nm, , drop = FALSE]
    .bed_rows_to_nipter_combined(sub_rows, nm, binsize, n_bins_global)
  })

  nipter_as_control_group(samples, description = description)
}


# Internal helper: build a NIPTeRSample from a data frame of rows already
# fetched from the multi-reader table.
# Columns: chrom, start_pos, end_pos, count, corrected_count (all typed).
# n_bins_global: shared column width so all samples in a control group have
#   identical matrix dimensions (narrower chromosomes are zero-padded).
.bed_rows_to_nipter_combined <- function(rows, name, binsize, n_bins_global) {
  auto_keys  <- as.character(1:22)
  sex_labels <- c("X", "Y")
  sex_keys   <- c("23", "24")

  chrom <- sub("^[Cc][Hh][Rr]", "", rows$chrom)
  chrom[chrom == "X" | chrom == "x"] <- "23"
  chrom[chrom == "Y" | chrom == "y"] <- "24"

  bin_idx <- as.integer(rows$start_pos / binsize)
  n_bins  <- n_bins_global

  col_names <- as.character(seq_len(n_bins))

  auto_mat <- matrix(0L, nrow = 22L, ncol = n_bins,
                     dimnames = list(auto_keys, col_names))
  sex_mat  <- matrix(0L, nrow = 2L, ncol = n_bins,
                     dimnames = list(sex_labels, col_names))

  for (i in seq_along(auto_keys)) {
    sel <- which(chrom == auto_keys[i])
    if (length(sel) > 0L)
      auto_mat[i, bin_idx[sel] + 1L] <- as.integer(rows$count[sel])
  }
  for (i in seq_along(sex_keys)) {
    sel <- which(chrom == sex_keys[i])
    if (length(sel) > 0L)
      sex_mat[i, bin_idx[sel] + 1L] <- as.integer(rows$count[sel])
  }

  obj <- structure(
    list(autosomal_chromosome_reads  = list(auto_mat),
         sex_chromosome_reads        = list(sex_mat),
         correction_status_autosomal = "Uncorrected",
         correction_status_sex       = "Uncorrected",
         sample_name                 = name),
    class = c("NIPTeRSample", "CombinedStrands")
  )

  if (!all(is.na(rows$corrected_count))) {
    corr_auto <- matrix(0, nrow = 22L, ncol = n_bins,
                        dimnames = list(auto_keys, col_names))
    corr_sex  <- matrix(0, nrow = 2L, ncol = n_bins,
                        dimnames = list(sex_labels, col_names))
    for (i in seq_along(auto_keys)) {
      sel <- which(chrom == auto_keys[i])
      if (length(sel) > 0L)
        corr_auto[i, bin_idx[sel] + 1L] <- rows$corrected_count[sel]
    }
    for (i in seq_along(sex_keys)) {
      sel <- which(chrom == sex_keys[i])
      if (length(sel) > 0L)
        corr_sex[i, bin_idx[sel] + 1L] <- rows$corrected_count[sel]
    }
    obj$autosomal_chromosome_reads  <- list(corr_auto)
    obj$sex_chromosome_reads        <- list(corr_sex)
    obj$correction_status_autosomal <- "GC Corrected"
    obj$correction_status_sex       <- "GC Corrected"
  }

  obj
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
