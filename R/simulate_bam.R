#' Simulate a trisomy BAM by thinning a donor alignment
#'
#' Creates a BAM by copying all reads on the target chromosome unchanged and
#' randomly removing reads from all other mapped chromosomes according to the
#' Nguyen et al. fetal-fraction model. This preserves the donor alignment's
#' read names, flags, CIGAR strings, mapping qualities, and coordinate order
#' while enriching the target chromosome relative to the rest of the genome.
#'
#' This is a BAM-level simulation primitive intended for realistic fixture and
#' benchmark generation from real negative donor alignments. It is more faithful
#' than the compressed synthetic cohort generator in [generate_cohort()], but it
#' is still a simulation and does not replace real-data validation.
#'
#' @param bam Input BAM or CRAM path.
#' @param out_bam Output BAM path.
#' @param trisomy_chr Target chromosome to enrich, e.g. `"21"` or `"chr21"`.
#' @param fetal_fraction Numeric fetal fraction in `(0, 1]`. Default `0.10`.
#' @param mosaic Logical; if `TRUE`, use the 50% mosaic enrichment model
#'   (`k = 0.25`) instead of the non-mosaic model (`k = 0.5`).
#' @param seed Integer RNG seed for deterministic thinning.
#' @param threads Integer thread count for htslib I/O and BAM indexing.
#' @param reference Optional FASTA path required when the input is CRAM.
#' @param index Logical; if `TRUE` (default), create a `.bai` index for
#'   `out_bam`.
#' @param overwrite Logical; overwrite `out_bam` and its `.bai` if present.
#'
#' @return A named list with input/output record counts and thinning metadata.
#'
#' @seealso [generate_cohort()], [bam_convert()]
#'
#' @examples
#' \dontrun{
#' simulate_trisomy_bam(
#'   bam = "negative_sample.bam",
#'   out_bam = "simulated_t21.bam",
#'   trisomy_chr = "21",
#'   fetal_fraction = 0.12,
#'   seed = 42L
#' )
#' }
#'
#' @export
simulate_trisomy_bam <- function(bam,
                                 out_bam,
                                 trisomy_chr,
                                 fetal_fraction = 0.10,
                                 mosaic = FALSE,
                                 seed = 1L,
                                 threads = 1L,
                                 reference = NULL,
                                 index = TRUE,
                                 overwrite = FALSE) {
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam), file.exists(bam))
  stopifnot(is.character(out_bam), length(out_bam) == 1L, nzchar(out_bam))
  stopifnot(is.character(trisomy_chr), length(trisomy_chr) == 1L, nzchar(trisomy_chr))
  stopifnot(is.numeric(fetal_fraction), length(fetal_fraction) == 1L,
            is.finite(fetal_fraction), fetal_fraction > 0, fetal_fraction <= 1)
  seed <- as.integer(seed)
  threads <- as.integer(threads)
  stopifnot(!is.na(seed), !is.na(threads), threads >= 1L)
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L,
              nzchar(reference), file.exists(reference))
  }

  bam_norm <- normalizePath(bam, winslash = "/", mustWork = TRUE)
  out_norm <- normalizePath(out_bam, winslash = "/", mustWork = FALSE)
  if (identical(bam_norm, out_norm)) {
    stop("Input and output paths must differ.", call. = FALSE)
  }
  out_dir <- dirname(out_bam)
  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir, call. = FALSE)
  }

  index_path <- paste0(out_bam, ".bai")
  if (!overwrite && (file.exists(out_bam) || file.exists(index_path))) {
    stop("Output BAM or BAI already exists. Use overwrite = TRUE to replace it.",
         call. = FALSE)
  }
  if (overwrite) {
    unlink(c(out_bam, index_path), force = TRUE)
  }

  target_chr <- .normalize_chr_name(trisomy_chr, xy_to_numeric = FALSE)
  removal_prob <- .trisomy_removal_prob(fetal_fraction = fetal_fraction,
                                        mosaic = mosaic)
  other_keep_prob <- 1 - removal_prob
  reference_path <- if (is.null(reference)) "" else reference

  res <- simulate_trisomy_bam_cpp(
    input_path = bam,
    output_path = out_bam,
    trisomy_chr = target_chr,
    other_keep_prob = other_keep_prob,
    seed = seed,
    threads = threads,
    reference = reference_path
  )

  if (isTRUE(index)) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    Rduckhts::rduckhts_bam_index(con, out_bam, threads = threads)
  }

  res$output_bam <- out_bam
  res$index_bam <- isTRUE(index)
  res$index_path <- if (isTRUE(index)) index_path else NA_character_
  res$fetal_fraction <- fetal_fraction
  res$mosaic <- isTRUE(mosaic)
  res$seed <- seed
  res$threads <- threads
  res$removal_prob <- removal_prob
  res
}
