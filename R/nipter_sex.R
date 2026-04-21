#' Build a sex prediction model from a NIPTeR control group
#'
#' Fits a two-component Gaussian mixture model (GMM) on sex chromosome
#' fractions derived from a \code{NIPTeRControlGroup}. The model
#' distinguishes male and female samples based on Y-chromosome read fraction
#' (univariate) or X+Y chromosome fractions (bivariate).
#'
#' This mirrors the approach used in clinical NIPT pipelines (see
#' \emph{Details}), ported to R using \code{mclust::Mclust()} in place of
#' Python's \code{sklearn.GaussianMixture}.
#'
#' @param control_group A \code{NIPTeRControlGroup} object. Ideally
#'   chi-squared corrected via \code{\link{nipter_chi_correct}}.
#' @param method Character; the feature space for the GMM:
#'   \describe{
#'     \item{\code{"y_fraction"}}{Univariate: Y-chromosome read fraction
#'       relative to total autosomal reads.}
#'     \item{\code{"xy_fraction"}}{Bivariate: X and Y chromosome fractions.}
#'   }
#'
#' @return A list-like \code{NIPTeRSexModel} S7 object with elements:
#' \describe{
#'   \item{model}{The \code{mclust::Mclust} fitted object.}
#'   \item{method}{Character; the method used.}
#'   \item{male_cluster}{Integer (1 or 2); which cluster is male
#'     (higher Y fraction).}
#'   \item{classifications}{Named character vector of \code{"male"}/
#'     \code{"female"} labels for each control sample.}
#'   \item{fractions}{Matrix of the input features used for fitting
#'     (samples as rows). Columns are \code{"y_fraction"} or
#'     \code{c("x_fraction", "y_fraction")}.}
#' }
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Compute sex chromosome read fractions for every sample in the
#'     control group. Y fraction = \code{sum(Y bins) / sum(autosomal bins)};
#'     X fraction = \code{sum(X bins) / sum(autosomal bins)}.
#'   \item Fit a two-component Gaussian mixture with equal mixing proportions
#'     (\code{mclust::Mclust(data, G = 2,
#'     control = mclust::emControl(equalPro = TRUE))}).
#'   \item Identify the male cluster as the component with the higher median
#'     Y fraction.
#' }
#'
#' This follows the clinical NIPT pipeline pattern of building sex models
#' from control cohorts. The user's pipeline builds three models (Y-unique
#' ratio from samtools, XY fractions, Y fraction) and takes a majority vote;
#' \code{nipter_predict_sex()} implements the consensus when given multiple
#' models.
#'
#' @seealso [nipter_predict_sex()], [nipter_chi_correct()]
#'
#' @importFrom mclust Mclust emControl
#'
#' @examples
#' \dontrun{
#' # Build sex model from chi-corrected control group
#' sex_model <- nipter_sex_model(chi_cg, method = "y_fraction")
#'
#' # Bivariate X+Y model
#' sex_model_xy <- nipter_sex_model(chi_cg, method = "xy_fraction")
#' }
#'
#' @export
nipter_sex_model <- function(control_group,
                             method = c("y_fraction", "xy_fraction")) {
  method <- match.arg(method)
  stopifnot(.is_nipt_control_group_object(control_group))

  n_samples <- length(control_group$samples)
  if (n_samples < 4L) {
    stop("At least 4 samples are needed to fit a sex model.", call. = FALSE)
  }

  # Compute sex chromosome fractions
  fracs <- .sex_chr_fractions(control_group)
  # fracs: n_samples x 2 matrix with columns "x_fraction", "y_fraction"

  # Select features for the GMM
  if (method == "y_fraction") {
    fit_data <- fracs[, "y_fraction"]
  } else {
    fit_data <- fracs[, c("x_fraction", "y_fraction"), drop = FALSE]
  }

  # Fit 2-component GMM with equal mixing proportions.
  # NOTE: mclust::Mclust() internally calls mclustBIC() without namespace
  # qualification, so a bare mclust::Mclust() fails when the package is loaded
  # via requireNamespace() but not attached. We evaluate the call inside the
  # mclust namespace to work around this.
  gmm <- .mclust_fit(fit_data, G = 2L, equalPro = TRUE)

  if (is.null(gmm)) {
    stop("mclust::Mclust() failed to converge. ",
         "The control group may lack sufficient sex chromosome signal.",
         call. = FALSE)
  }

  # Identify the male cluster (higher median Y fraction)
  y_fracs <- fracs[, "y_fraction"]
  median_1 <- stats::median(y_fracs[gmm$classification == 1L], na.rm = TRUE)
  median_2 <- stats::median(y_fracs[gmm$classification == 2L], na.rm = TRUE)
  male_cluster <- if (median_1 > median_2) 1L else 2L

  # Build labels
  sample_names <- vapply(control_group$samples, .sample_name, character(1L))
  labels <- ifelse(gmm$classification == male_cluster, "male", "female")
  names(labels) <- sample_names

  .as_nipter_sex_model(list(
    model           = gmm,
    method          = method,
    male_cluster    = male_cluster,
    classifications = labels,
    fractions       = fracs
  ))
}


#' Predict fetal sex for a NIPTeR sample
#'
#' Uses one or more \code{NIPTeRSexModel} objects to classify a sample as
#' male or female. When multiple models are supplied, a majority vote
#' determines the consensus call.
#'
#' @param sample A \code{NIPTeRSample} object.
#' @param ... One or more \code{NIPTeRSexModel} objects, a single list of such
#'   models, or one \code{NIPTReferenceModel} containing pre-built
#'   \code{sex_models}.
#' @param y_unique_ratio Optional numeric scalar; a pre-computed Y-unique
#'   ratio (from \code{\link{nipter_y_unique_ratio}}) for the sample. This
#'   is only used when one of the models has \code{method = "y_unique"}.
#'   If a \code{"y_unique"} model is present, this argument is required.
#'
#' @return A list-like \code{NIPTeRSexPrediction} S7 object with elements:
#' \describe{
#'   \item{prediction}{Character; \code{"male"} or \code{"female"} (consensus
#'     when multiple models).}
#'   \item{model_predictions}{Named character vector of per-model predictions
#'     (names are the method names).}
#'   \item{fractions}{Named numeric vector with \code{x_fraction} and
#'     \code{y_fraction} for the sample.}
#'   \item{sample_name}{Character; the sample name.}
#' }
#'
#' @details
#' For each model, the sample's sex chromosome fractions are computed and
#' classified using \code{predict(model$model, newdata = ...)}. The
#' classification is mapped to \code{"male"}/\code{"female"} using the
#' model's \code{male_cluster}.
#'
#' For \code{"y_unique"} models, the \code{y_unique_ratio} argument is
#' used as the input feature instead of chromosome fractions derived from
#' binned read counts.
#'
#' With multiple models, the consensus is the label that appears most
#' frequently (majority vote). In case of a tie, \code{"female"} is
#' returned (conservative default for NIPT, as false-male calls could mask
#' sex chromosome aneuploidies).
#'
#' @seealso [nipter_sex_model()], [nipter_sex_model_y_unique()],
#'   [nipter_y_unique_ratio()]
#'
#' @examples
#' \dontrun{
#' # Single model
#' pred <- nipter_predict_sex(test_sample, sex_model)
#' pred$prediction
#'
#' # Consensus from three models (majority vote)
#' yr <- nipter_y_unique_ratio("sample.bam")
#' pred <- nipter_predict_sex(test_sample, model_y, model_xy, model_yu,
#'                            y_unique_ratio = yr$ratio)
#' pred$prediction
#' }
#'
#' @export
nipter_predict_sex <- function(sample, ..., y_unique_ratio = NULL) {
  stopifnot(.is_nipt_sample_object(sample))

  models <- list(...)
  # Allow a single list of models

  if (length(models) == 1L && .is_nipt_reference_model(models[[1L]])) {
    models <- models[[1L]]$sex_models
  } else if (length(models) == 1L && is.list(models[[1L]]) &&
             !.is_nipter_sex_model(models[[1L]])) {
    models <- models[[1L]]
  }

  if (length(models) == 0L) {
    stop("At least one NIPTeRSexModel must be provided.", call. = FALSE)
  }

  for (i in seq_along(models)) {
    if (!.is_nipter_sex_model(models[[i]])) {
      stop("All models must be NIPTeRSexModel objects.", call. = FALSE)
    }
  }

  # Validate y_unique_ratio if provided
  if (!is.null(y_unique_ratio)) {
    stopifnot(is.numeric(y_unique_ratio), length(y_unique_ratio) == 1L)
  }

  # Compute sample's sex chromosome fractions
  sample_fracs <- .sample_sex_fractions(sample)

  model_preds <- character(length(models))
  names(model_preds) <- make.unique(vapply(models, function(m) m$method, character(1L)))

  for (i in seq_along(models)) {
    mdl <- models[[i]]

    if (mdl$method == "y_unique") {
      if (is.null(y_unique_ratio)) {
        stop(
          "A y_unique sex model was supplied but 'y_unique_ratio' is NULL. ",
          "Provide the sample's post-filter Y-unique ratio explicitly.",
          call. = FALSE
        )
      }
      newdata <- y_unique_ratio
    } else if (mdl$method == "y_fraction") {
      newdata <- sample_fracs["y_fraction"]
    } else {
      newdata <- matrix(sample_fracs[c("x_fraction", "y_fraction")],
                        nrow = 1L)
      colnames(newdata) <- c("x_fraction", "y_fraction")
    }

    pred <- stats::predict(mdl$model, newdata = newdata)
    cluster <- pred$classification[1L]
    model_preds[i] <- if (cluster == mdl$male_cluster) "male" else "female"
  }

  # Consensus: majority vote (tie goes to "female" — conservative for NIPT)
  n_male   <- sum(model_preds == "male")
  n_female <- sum(model_preds == "female")
  if (n_male == n_female) {
    warning(
      sprintf(
        "Sex-model vote tied for sample '%s' (%d male vs %d female); returning 'female' conservatively.",
        .sample_name(sample), n_male, n_female
      ),
      call. = FALSE
    )
  }
  consensus <- if (n_male > n_female) "male" else "female"

  .as_nipter_sex_prediction(list(
    prediction = consensus,
    model_predictions = model_preds,
    fractions = sample_fracs,
    sample_name = .sample_name(sample)
  ))
}


#' Build a sex prediction model from Y-unique ratios
#'
#' Fits a two-component Gaussian mixture model (GMM) on Y-unique region
#' read ratios, which are computed from BAM files by
#' \code{\link{nipter_y_unique_ratio}}.
#'
#' This is a companion to \code{\link{nipter_sex_model}}, which operates on
#' binned \code{NIPTeRSample} fractions. The Y-unique model operates at the
#' BAM level and does not require prior binning. The resulting model object
#' is compatible with \code{\link{nipter_predict_sex}} for majority-vote
#' consensus.
#'
#' @param ratios Named numeric vector of Y-unique ratios (one per sample).
#'   Typically obtained by calling \code{\link{nipter_y_unique_ratio}} on
#'   each BAM in the control cohort.
#'
#' @return A list-like \code{NIPTeRSexModel} S7 object with elements:
#' \describe{
#'   \item{model}{The \code{mclust::Mclust} fitted object.}
#'   \item{method}{\code{"y_unique"}.}
#'   \item{male_cluster}{Integer (1 or 2); which cluster is male
#'     (higher Y-unique ratio).}
#'   \item{classifications}{Named character vector of \code{"male"}/
#'     \code{"female"} labels for each input sample.}
#'   \item{fractions}{The input \code{ratios} vector (named).}
#' }
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Fit a two-component Gaussian mixture with equal mixing proportions
#'     (\code{mclust::Mclust(ratios, G = 2,
#'     control = mclust::emControl(equalPro = TRUE))}).
#'   \item Identify the male cluster as the component with the higher
#'     median Y-unique ratio.
#' }
#'
#' @seealso [nipter_y_unique_ratio()], [nipter_predict_sex()],
#'   [nipter_sex_model()]
#'
#' @examples
#' \dontrun{
#' # Compute Y-unique ratios for a cohort
#' bams <- list.files("bams/", pattern = "\\.bam$", full.names = TRUE)
#' ratios <- vapply(bams, function(b) nipter_y_unique_ratio(b)$ratio,
#'                  numeric(1L))
#' names(ratios) <- basename(bams)
#'
#' # Build sex model
#' model_yu <- nipter_sex_model_y_unique(ratios)
#' model_yu$classifications
#' }
#'
#' @export
nipter_sex_model_y_unique <- function(ratios) {
  stopifnot(is.numeric(ratios), length(ratios) >= 4L)
  if (is.null(names(ratios))) {
    names(ratios) <- paste0("sample_", seq_along(ratios))
  }

  gmm <- .mclust_fit(ratios, G = 2L, equalPro = TRUE)

  if (is.null(gmm)) {
    stop("mclust::Mclust() failed to converge. ",
         "The cohort may lack sufficient Y-unique signal.",
         call. = FALSE)
  }

  # Male cluster = higher median Y-unique ratio
  median_1 <- stats::median(ratios[gmm$classification == 1L], na.rm = TRUE)
  median_2 <- stats::median(ratios[gmm$classification == 2L], na.rm = TRUE)
  male_cluster <- if (median_1 > median_2) 1L else 2L

  labels <- ifelse(gmm$classification == male_cluster, "male", "female")
  names(labels) <- names(ratios)

  .as_nipter_sex_model(list(
    model           = gmm,
    method          = "y_unique",
    male_cluster    = male_cluster,
    classifications = labels,
    fractions       = ratios
  ))
}


#' Compute Y-unique region read ratio from a BAM file
#'
#' Counts reads overlapping a Y-chromosome interval set and divides by total
#' nuclear genome reads (chromosomes 1--22, X, Y). This ratio is a strong
#' univariate sex predictor used in the clinical NIPT pipeline mirrored here.
#'
#' For GRCh37, the default interval set is the bundled legacy
#' \code{extdata/grch37_Y_chrom_blacklist.bed}, which matches the historical
#' \code{samtools stats -t human_g1k_v37.fasta_Y_chrom_blacklist.bed} path used
#' by the production shell pipeline. A smaller 7-gene interval file can still
#' be supplied explicitly through \code{regions_file}, but it is no longer the
#' default because it does not reproduce the legacy \code{YUniqueRatioFiltered}
#' semantics.
#'
#' The BAM must be indexed (\code{.bai} or \code{.csi}).
#'
#' Read counting uses DuckDB/duckhts BED coverage summaries via
#' [Rduckhts::rduckhts_bam_bed_coverage()]. Interval chromosome names are
#' matched to the BAM header so both \code{Y} and \code{chrY}-style contigs
#' work. Both the Y-target BED and the whole-nuclear BED apply the same MAPQ
#' and flag filters.
#'
#' Y-target reads are summed across the supplied intervals exactly as listed in
#' \code{regions_file}. This intentionally preserves the historical
#' interval-wise summation semantics: if two target intervals overlap, a read
#' contributing to both intervals is counted twice.
#'
#' @param bam Path to an indexed BAM or CRAM file.
#' @param mapq Minimum mapping quality. Default \code{1L}.
#' @param exclude_flags Integer bitmask; reads with any of these flags set
#'   are dropped (samtools \code{-F}). Default \code{0L}.
#' @param require_flags Integer bitmask; only reads with all bits set are
#'   kept (samtools \code{-f}). Default \code{0L}.
#' @param regions_file Path to a Y-unique regions file. Headered TSV files
#'   with columns \code{Chromosome}, \code{Start}, \code{End},
#'   \code{GeneName} are accepted, as are headerless BED-like files whose
#'   first 3 columns are chromosome/start/end and optional 4th column is a
#'   gene or interval label. Defaults to the bundled GRCh37 file. Supply a
#'   custom file for GRCh38 or other assemblies.
#' @param con Optional open DBI connection with duckhts loaded. If
#'   \code{NULL} (default) a temporary in-memory DuckDB connection is
#'   created.
#' @param reference Optional FASTA reference path for CRAM inputs.
#' @param decompression_threads Integer; number of htslib decompression worker
#'   threads to use for BAM/CRAM decoding in the BED-coverage path. Default
#'   \code{0L}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{ratio}{Numeric; Y-unique reads / total nuclear reads. \code{0}
#'     if total nuclear reads is zero.}
#'   \item{y_unique_reads}{Integer; reads overlapping the supplied Y interval
#'     set.}
#'   \item{total_nuclear_reads}{Integer; total reads on chromosomes 1--22,
#'     X, Y.}
#'   \item{regions}{Data frame of the intervals used (columns:
#'     \code{Chromosome}, \code{Start}, \code{End}, \code{GeneName}).}
#'   \item{regions_file}{Normalized path to the regions file used.}
#' }
#'
#' @seealso [nipter_sex_model()], [nipter_predict_sex()]
#'
#' @examples
#' \dontrun{
#' yr <- nipter_y_unique_ratio("sample.bam", mapq = 1L)
#' yr$ratio
#' yr$y_unique_reads
#' yr$total_nuclear_reads
#' }
#'
#' @export
nipter_y_unique_ratio <- function(bam,
                                  mapq          = 1L,
                                  require_flags = 0L,
                                  exclude_flags = 0L,
                                  regions_file  = NULL,
                                  con           = NULL,
                                  reference     = NULL,
                                  decompression_threads = 0L) {
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  stopifnot(is.numeric(mapq),          length(mapq)          == 1L, mapq          >= 0L)
  stopifnot(is.numeric(require_flags), length(require_flags) == 1L, require_flags >= 0L)
  stopifnot(is.numeric(exclude_flags), length(exclude_flags) == 1L, exclude_flags >= 0L)
  stopifnot(is.numeric(decompression_threads),
            length(decompression_threads) == 1L,
            decompression_threads >= 0L)
  mapq          <- as.integer(mapq)
  require_flags <- as.integer(require_flags)
  exclude_flags <- as.integer(exclude_flags)
  decompression_threads <- as.integer(decompression_threads)

  # Load regions
  if (is.null(regions_file)) {
    regions_file <- system.file(
      "extdata",
      "grch37_Y_chrom_blacklist.bed",
      package = "RWisecondorX",
      mustWork = TRUE
    )
  }
  stopifnot(file.exists(regions_file))
  regions <- .read_y_unique_regions_file(regions_file)
  regions_file <- normalizePath(regions_file, winslash = "/", mustWork = TRUE)

  # Manage connection
  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  hdr <- Rduckhts::rduckhts_hts_header(con, bam)
  hdr <- hdr[
    hdr$record_type == "SQ" &
      !is.na(hdr$id) &
      nzchar(as.character(hdr$id)) &
      !is.na(hdr$length),
    ,
    drop = FALSE
  ]
  header_chr <- as.character(hdr$id)
  header_chr_norm <- .normalize_chr_name(header_chr, xy_to_numeric = FALSE)

  region_chr_query <- vapply(regions$Chromosome, function(chr) {
    chr_norm <- .normalize_chr_name(chr, xy_to_numeric = FALSE)
    hit <- match(chr_norm, header_chr_norm)
    if (is.na(hit)) NA_character_ else header_chr[[hit]]
  }, character(1L))
  keep_regions <- !is.na(region_chr_query) & nzchar(region_chr_query)
  regions_used <- regions[keep_regions, c("Chromosome", "Start", "End", "GeneName"), drop = FALSE]
  if (nrow(regions_used)) {
    regions_bed <- data.frame(
      chrom = region_chr_query[keep_regions],
      start = as.integer(regions_used$Start),
      end = as.integer(regions_used$End),
      name = ifelse(
        is.na(regions_used$GeneName) | !nzchar(as.character(regions_used$GeneName)),
        sprintf(
          "%s:%s-%s",
          regions_used$Chromosome,
          as.integer(regions_used$Start),
          as.integer(regions_used$End)
        ),
        as.character(regions_used$GeneName)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    regions_bed <- data.frame(
      chrom = character(),
      start = integer(),
      end = integer(),
      name = character()
    )
  }

  nuclear_keep <- header_chr_norm %in% c(as.character(seq_len(22L)), "X", "Y")
  nuclear_bed <- data.frame(
    chrom = header_chr[nuclear_keep],
    start = 0L,
    end = as.integer(hdr$length[nuclear_keep]),
    name = paste0("nuclear_", header_chr[nuclear_keep]),
    stringsAsFactors = FALSE
  )

  y_bed_path <- tempfile("y_unique_regions_", fileext = ".bed")
  nuclear_bed_path <- tempfile("nuclear_regions_", fileext = ".bed")
  on.exit(unlink(c(y_bed_path, nuclear_bed_path), force = TRUE), add = TRUE)

  utils::write.table(
    regions_bed,
    file = y_bed_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  utils::write.table(
    nuclear_bed,
    file = nuclear_bed_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  coverage_args <- list(
    con = con,
    path = bam,
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    reference = reference,
    decompression_threads = decompression_threads,
    strand_outputs = FALSE
  )

  if (nrow(regions_bed)) {
    y_cov <- do.call(
      Rduckhts::rduckhts_bam_bed_coverage,
      c(coverage_args, list(bed_path = y_bed_path))
    )
    y_unique_reads <- as.integer(sum(y_cov$numreads_post, na.rm = TRUE))
  } else {
    y_unique_reads <- 0L
  }

  if (nrow(nuclear_bed)) {
    total_cov <- do.call(
      Rduckhts::rduckhts_bam_bed_coverage,
      c(coverage_args, list(bed_path = nuclear_bed_path))
    )
    total_nuclear_reads <- as.integer(sum(total_cov$numreads_post, na.rm = TRUE))
  } else {
    total_nuclear_reads <- 0L
  }

  ratio <- if (total_nuclear_reads > 0L) {
    as.numeric(y_unique_reads) / as.numeric(total_nuclear_reads)
  } else 0

  list(
    ratio               = ratio,
    y_unique_reads      = y_unique_reads,
    total_nuclear_reads = total_nuclear_reads,
    regions             = regions_used,
    regions_file        = regions_file
  )
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Compute X and Y chromosome fractions for all samples in a control group.
# Returns an n_samples x 2 matrix with columns "x_fraction", "y_fraction".
# Each fraction = sum(sex_chr bins) / sum(autosomal bins).
.sex_chr_fractions <- function(control_group) {
  n <- length(control_group$samples)
  mat <- matrix(NA_real_, nrow = n, ncol = 2L)
  colnames(mat) <- c("x_fraction", "y_fraction")

  for (i in seq_len(n)) {
    mat[i, ] <- .sample_sex_fractions(control_group$samples[[i]])
  }

  rownames(mat) <- vapply(control_group$samples, .sample_name, character(1L))
  mat
}


# Compute X and Y chromosome fractions for a single NIPTeRSample.
# Returns named numeric vector c(x_fraction = ..., y_fraction = ...).
# Fractions are relative to total autosomal reads (sum of all autosomal bins).
.sample_sex_fractions <- function(sample) {
  # Autosomal total
  auto <- .sample_autosomal_reads(sample)
  if (.strand_type_of(sample) == "separated") {
    auto_total <- sum(Reduce("+", auto))
  } else {
    auto_total <- sum(auto[[1L]])
  }

  if (auto_total == 0) {
    return(c(x_fraction = 0, y_fraction = 0))
  }

  # Sex chromosome totals
  sex <- .sample_sex_reads(sample)
  if (.strand_type_of(sample) == "separated") {
    sex_totals <- c(X = 0, Y = 0)
    for (mat in sex) {
      rs <- rowSums(mat)
      key <- sub("[FR]$", "", names(rs))
      agg <- tapply(rs, key, sum)
      sex_totals[names(agg)] <- sex_totals[names(agg)] + agg
    }
    x_total <- sex_totals[["X"]]
    y_total <- sex_totals[["Y"]]
  } else {
    x_total <- sum(sex[[1L]]["X", ])
    y_total <- sum(sex[[1L]]["Y", ])
  }

  c(x_fraction = x_total / auto_total,
    y_fraction = y_total / auto_total)
}


# Wrapper around mclust::Mclust() that evaluates the call inside the mclust
# namespace.  mclust::Mclust() internally calls mclustBIC() without its own
# namespace prefix, so the bare mclust::Mclust() path fails when the package
# is loaded via requireNamespace() but not attached to the search path.
# Evaluating in the namespace makes internal symbol lookup succeed, which is
# the CRAN-friendly alternative to library()/require() inside a function.
.mclust_fit <- function(data, G = 2L, equalPro = TRUE) {
  .mclust_fit_in_namespace(data = data, G = G, equalPro = equalPro)
}
