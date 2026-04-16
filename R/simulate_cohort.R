.alignment_stem <- function(path) {
  sub("\\.(bam|cram)$", "", basename(path), ignore.case = TRUE)
}

.alignment_index_path <- function(path) {
  candidates <- c(
    paste0(path, ".bai"),
    paste0(path, ".crai"),
    sub("\\.(bam|cram)$", ".bai", path, ignore.case = TRUE),
    sub("\\.(bam|cram)$", ".crai", path, ignore.case = TRUE)
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) NA_character_ else hits[[1]]
}

.resolve_alignment_inputs <- function(bams = NULL, bam_dir = NULL,
                                      pattern = "\\.(bam|cram)$") {
  if (is.null(bams) == is.null(bam_dir)) {
    stop("Supply exactly one of `bams` or `bam_dir`.", call. = FALSE)
  }

  if (!is.null(bams)) {
    stopifnot(is.character(bams), length(bams) >= 1L)
    if (any(!nzchar(bams)) || any(!file.exists(bams))) {
      stop("All `bams` entries must be existing BAM/CRAM paths.", call. = FALSE)
    }
    return(normalizePath(bams, winslash = "/", mustWork = TRUE))
  }

  stopifnot(is.character(bam_dir), length(bam_dir) == 1L, nzchar(bam_dir),
            dir.exists(bam_dir))
  paths <- list.files(bam_dir, pattern = pattern, full.names = TRUE)
  paths <- sort(paths)
  if (length(paths) == 0L) {
    stop("No BAM/CRAM files matched `pattern` in `bam_dir`.", call. = FALSE)
  }
  normalizePath(paths, winslash = "/", mustWork = TRUE)
}

.build_simulation_grid <- function(trisomy_chromosomes, fetal_fraction, mosaic) {
  stopifnot(is.character(trisomy_chromosomes), length(trisomy_chromosomes) >= 1L,
            all(nzchar(trisomy_chromosomes)))
  stopifnot(is.numeric(fetal_fraction), length(fetal_fraction) >= 1L,
            all(is.finite(fetal_fraction)), all(fetal_fraction > 0),
            all(fetal_fraction <= 1))
  stopifnot(is.logical(mosaic), length(mosaic) >= 1L, all(!is.na(mosaic)))

  expand.grid(
    trisomy_chr = unique(.normalize_chr_name(trisomy_chromosomes,
                                             xy_to_numeric = FALSE)),
    fetal_fraction = unique(as.numeric(fetal_fraction)),
    mosaic = unique(as.logical(mosaic)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

.format_sim_tag <- function(x) {
  txt <- formatC(x, format = "f", digits = 3)
  txt <- sub("0+$", "", txt)
  txt <- sub("\\.$", "", txt)
  gsub("\\.", "p", txt)
}

.simulation_output_name <- function(donor_id, trisomy_chr, fetal_fraction, mosaic) {
  suffix <- if (isTRUE(mosaic)) "_mosaic" else ""
  sprintf("%s_sim_T%s_ff%s%s.bam",
          donor_id,
          .normalize_chr_name(trisomy_chr, xy_to_numeric = FALSE),
          .format_sim_tag(fetal_fraction),
          suffix)
}

#' Simulate a trisomy cohort from donor BAMs/CRAMs
#'
#' Generates a set of positive BAMs from one or more negative donor alignments
#' by applying [simulate_trisomy_bam()] across a simulation grid. The result is
#' a manifest that ties every simulated BAM back to its donor, chromosome,
#' fetal fraction, and mosaic setting.
#'
#' The donor BAMs themselves serve as the negative cohort; by default they are
#' recorded in the manifest but not copied into `out_dir`.
#'
#' @param bams Optional character vector of donor BAM/CRAM paths.
#' @param bam_dir Optional directory containing donor BAM/CRAM files.
#'   Supply exactly one of `bams` or `bam_dir`.
#' @param out_dir Output directory for simulated positive BAMs and the manifest.
#' @param donor_ids Optional donor IDs, one per input BAM/CRAM. Defaults to
#'   filename stems with duplicates made unique.
#' @param trisomy_chromosomes Character vector of target chromosomes.
#'   Default `c("21", "18", "13")`.
#' @param fetal_fraction Numeric vector of fetal fractions. A separate
#'   simulation is produced for each distinct value.
#' @param mosaic Logical vector of mosaic settings. A separate simulation is
#'   produced for each distinct value. Default `FALSE`.
#' @param seed Integer base seed. Individual simulations use deterministic
#'   seeds derived from this base.
#' @param threads Integer thread count for htslib I/O and indexing.
#' @param reference Optional FASTA path applied to all CRAM inputs.
#' @param index Logical; create `.bai` indexes for simulated BAMs. Default
#'   `TRUE`.
#' @param include_donors Logical; include the donor negatives as rows in the
#'   manifest. Default `TRUE`.
#' @param pattern File-matching regex used only when `bam_dir` is supplied.
#' @param overwrite Logical; overwrite existing simulated BAMs and manifest.
#'
#' @return A data frame manifest. The same table is written to
#'   `<out_dir>/manifest.tsv`.
#'
#' @seealso [simulate_trisomy_bam()], [generate_cohort()]
#'
#' @examples
#' \dontrun{
#' manifest <- simulate_trisomy_cohort(
#'   bam_dir = "negatives/",
#'   out_dir = "simulated_positives/",
#'   trisomy_chromosomes = c("21", "18"),
#'   fetal_fraction = c(0.08, 0.12)
#' )
#' }
#'
#' @export
simulate_trisomy_cohort <- function(bams = NULL,
                                    bam_dir = NULL,
                                    out_dir,
                                    donor_ids = NULL,
                                    trisomy_chromosomes = c("21", "18", "13"),
                                    fetal_fraction = 0.10,
                                    mosaic = FALSE,
                                    seed = 1L,
                                    threads = 1L,
                                    reference = NULL,
                                    index = TRUE,
                                    include_donors = TRUE,
                                    pattern = "\\.(bam|cram)$",
                                    overwrite = FALSE) {
  stopifnot(is.character(out_dir), length(out_dir) == 1L, nzchar(out_dir))
  seed <- as.integer(seed)
  threads <- as.integer(threads)
  stopifnot(!is.na(seed), !is.na(threads), threads >= 1L)
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L,
              nzchar(reference), file.exists(reference))
  }

  input_bams <- .resolve_alignment_inputs(bams = bams, bam_dir = bam_dir,
                                          pattern = pattern)
  if (is.null(donor_ids)) {
    donor_ids <- make.unique(vapply(input_bams, .alignment_stem, character(1)))
  } else {
    stopifnot(is.character(donor_ids), length(donor_ids) == length(input_bams),
              all(nzchar(donor_ids)))
    if (anyDuplicated(donor_ids)) {
      stop("`donor_ids` must be unique.", call. = FALSE)
    }
  }

  sim_grid <- .build_simulation_grid(trisomy_chromosomes = trisomy_chromosomes,
                                     fetal_fraction = fetal_fraction,
                                     mosaic = mosaic)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  manifest_path <- file.path(out_dir, "manifest.tsv")
  if (!overwrite && file.exists(manifest_path)) {
    stop("Manifest already exists: ", manifest_path,
         ". Use overwrite = TRUE to replace it.", call. = FALSE)
  }

  manifest_rows <- vector("list", length(input_bams) * nrow(sim_grid) +
                            if (isTRUE(include_donors)) length(input_bams) else 0L)
  row_i <- 1L

  if (isTRUE(index)) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  sim_i <- 0L
  for (i in seq_along(input_bams)) {
    donor_bam <- input_bams[i]
    donor_id <- donor_ids[i]

    if (isTRUE(include_donors)) {
      manifest_rows[[row_i]] <- data.frame(
        cohort_role = "donor_negative",
        donor_id = donor_id,
        source_bam = donor_bam,
        output_bam = donor_bam,
        index_path = .alignment_index_path(donor_bam),
        simulated = FALSE,
        trisomy_chr = NA_character_,
        fetal_fraction = NA_real_,
        mosaic = NA,
        seed = NA_integer_,
        input_records = NA_real_,
        output_records = NA_real_,
        target_input = NA_real_,
        target_output = NA_real_,
        other_input = NA_real_,
        other_output = NA_real_,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }

    for (j in seq_len(nrow(sim_grid))) {
      sim_i <- sim_i + 1L
      sim_seed <- seed + sim_i - 1L
      sim_chr <- sim_grid$trisomy_chr[j]
      sim_ff <- sim_grid$fetal_fraction[j]
      sim_mosaic <- sim_grid$mosaic[j]
      out_bam <- file.path(
        out_dir,
        .simulation_output_name(
          donor_id = donor_id,
          trisomy_chr = sim_chr,
          fetal_fraction = sim_ff,
          mosaic = sim_mosaic
        )
      )

      sim_res <- simulate_trisomy_bam(
        bam = donor_bam,
        out_bam = out_bam,
        trisomy_chr = sim_chr,
        fetal_fraction = sim_ff,
        mosaic = sim_mosaic,
        seed = sim_seed,
        threads = threads,
        reference = reference,
        index = FALSE,
        overwrite = overwrite
      )

      if (isTRUE(index)) {
        Rduckhts::rduckhts_bam_index(con, out_bam, threads = threads)
      }

      manifest_rows[[row_i]] <- data.frame(
        cohort_role = "simulated_positive",
        donor_id = donor_id,
        source_bam = donor_bam,
        output_bam = sim_res$output_bam,
        index_path = if (isTRUE(index)) paste0(sim_res$output_bam, ".bai") else NA_character_,
        simulated = TRUE,
        trisomy_chr = sim_res$target_chr,
        fetal_fraction = sim_res$fetal_fraction,
        mosaic = sim_res$mosaic,
        seed = sim_res$seed,
        input_records = sim_res$input_records,
        output_records = sim_res$output_records,
        target_input = sim_res$target_input,
        target_output = sim_res$target_output,
        other_input = sim_res$other_input,
        other_output = sim_res$other_output,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  manifest <- do.call(rbind, manifest_rows[seq_len(row_i - 1L)])
  utils::write.table(manifest, manifest_path, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  attr(manifest, "manifest_path") <- manifest_path
  manifest
}
