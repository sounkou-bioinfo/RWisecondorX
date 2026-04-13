# cohort.R — Synthetic BAM cohort generator for WisecondorX / NIPTeR testing
#
# Generates coordinate-sorted, indexed BAM files using "compressed" chromosome
# lengths: each 100kb genomic bin maps to 100bp, so bin COUNT structure is
# identical to real GRCh37 at 100kb resolution, but BAMs are tiny (~435KB each).
#
# Cohort composition (50 samples):
#   35 euploid female (XX, zero Y reads)
#   12 euploid male   (XY, ~50% X/Y relative to autosomes)
#    1 trisomy 21 female (1.5x chr21)
#    1 trisomy 18 female (1.5x chr18)
#    1 trisomy 13 female (1.5x chr13)
#
# Deterministic: sample i uses set.seed(42 + i).
# NIPTeR-compatible: no unmapped reads, unique positions per chromosome.
# Requires: samtools on PATH.


# ---- GRCh37 bin structure (package-internal constants) --------------------

# Number of 100kb bins per chromosome (ceiling(chr_length / 100000))
BINS_PER_CHR <- c(
  "1"  = 2493L, "2"  = 2432L, "3"  = 1981L, "4"  = 1912L,
  "5"  = 1810L, "6"  = 1712L, "7"  = 1592L, "8"  = 1464L,
  "9"  = 1413L, "10" = 1356L, "11" = 1351L, "12" = 1339L,
  "13" = 1152L, "14" = 1074L, "15" = 1026L, "16" =  904L,
  "17" =  812L, "18" =  781L, "19" =  592L, "20" =  631L,
  "21" =  482L, "22" =  514L, "X"  = 1553L, "Y"  =  594L
)

#' Compressed bin size for synthetic cohort BAMs
#'
#' Each 100kb genomic bin is compressed to 100bp, yielding identical bin
#' counts to GRCh37 at 100kb resolution.
#'
#' @export
COMPRESSED_BINSIZE <- 100L

# Compressed chromosome lengths: bins * 100bp per bin
COMPRESSED_CHR_LENGTHS <- BINS_PER_CHR * COMPRESSED_BINSIZE

# SAM chromosome names (no "chr" prefix, matching GRCh37/b37 convention)
CHR_NAMES <- names(BINS_PER_CHR)

# Read properties
READ_LEN <- 36L
READ_SEQ <- paste(rep("A", READ_LEN), collapse = "")
READ_QUAL <- paste(rep("F", READ_LEN), collapse = "")
READ_CIGAR <- paste0(READ_LEN, "M")

# Average reads per bin for euploid autosomes
READS_PER_BIN <- 3.0


# ---- Sample profile constructors (internal) ------------------------------

.make_profile <- function(id, sex, chr_mult) {
  list(id = id, sex = sex, chr_mult = chr_mult)
}

.euploid_female <- function(id) {
  mult <- stats::setNames(rep(1.0, 24), CHR_NAMES)
  mult["Y"] <- 0.0
  .make_profile(id, "F", mult)
}

.euploid_male <- function(id) {
  mult <- stats::setNames(rep(1.0, 24), CHR_NAMES)
  mult["X"] <- 0.5
  mult["Y"] <- 0.5
  .make_profile(id, "M", mult)
}

.trisomy_female <- function(id, trisomy_chr) {
  mult <- stats::setNames(rep(1.0, 24), CHR_NAMES)
  mult["Y"] <- 0.0
  mult[as.character(trisomy_chr)] <- 1.5
  .make_profile(id, "F", mult)
}


# ---- Cohort definition ---------------------------------------------------

.build_cohort_profiles <- function() {
  profiles <- vector("list", 50L)
  idx <- 1L

  for (i in 1:35) {
    profiles[[idx]] <- .euploid_female(sprintf("euploid_f_%02d", i))
    idx <- idx + 1L
  }
  for (i in 1:12) {
    profiles[[idx]] <- .euploid_male(sprintf("euploid_m_%02d", i))
    idx <- idx + 1L
  }

  profiles[[idx]] <- .trisomy_female("trisomy_21", "21"); idx <- idx + 1L
  profiles[[idx]] <- .trisomy_female("trisomy_18", "18"); idx <- idx + 1L
  profiles[[idx]] <- .trisomy_female("trisomy_13", "13")

  profiles
}


# ---- SAM generation (internal) -------------------------------------------

.sam_header <- function() {
  hd <- "@HD\tVN:1.6\tSO:coordinate"
  sq <- sprintf("@SQ\tSN:%s\tLN:%d", CHR_NAMES, COMPRESSED_CHR_LENGTHS)
  paste(c(hd, sq), collapse = "\n")
}

.sam_records <- function(profile, seed) {
  set.seed(seed)

  total_bins <- sum(BINS_PER_CHR)
  total_reads <- round(total_bins * READS_PER_BIN)
  genome_bp <- sum(as.numeric(COMPRESSED_CHR_LENGTHS))

  all_records <- character(0)
  read_counter <- 0L

  for (chr_name in CHR_NAMES) {
    chr_len <- COMPRESSED_CHR_LENGTHS[chr_name]
    mult <- profile$chr_mult[chr_name]
    if (is.na(mult) || mult == 0) next

    # Expected reads proportional to chr length (in compressed space)
    n_reads <- round(total_reads * (chr_len / genome_bp) * mult)
    if (n_reads == 0L) next

    # Max valid position: chr_len - READ_LEN + 1 (1-based)
    max_pos <- chr_len - READ_LEN + 1L
    if (max_pos < 1L) next
    if (n_reads > max_pos) n_reads <- max_pos

    positions <- sort(sample.int(max_pos, n_reads, replace = FALSE))

    recs <- sprintf("r%d\t0\t%s\t%d\t60\t%s\t*\t0\t0\t%s\t%s",
                    read_counter + seq_len(n_reads),
                    chr_name, positions,
                    READ_CIGAR, READ_SEQ, READ_QUAL)
    all_records <- c(all_records, recs)
    read_counter <- read_counter + n_reads
  }

  all_records
}


# ---- Core generator (exported) -------------------------------------------

#' Generate a synthetic BAM cohort for testing
#'
#' Creates 50 coordinate-sorted, indexed BAM files in a specified directory
#' using compressed chromosome lengths: each 100kb GRCh37 bin is represented
#' as a 100bp region, yielding identical bin-count structure at a fraction
#' of the file size (~435KB per BAM).
#'
#' The cohort contains 35 euploid females, 12 euploid males, and 3 trisomy
#' females (T21, T18, T13). Each sample uses approximately 3 reads per bin
#' (~91k reads total). Results are deterministic (sample *i* uses
#' `set.seed(42 + i)`).
#'
#' @param out_dir Directory to write BAM files into (created if needed).
#' @param verbose Logical; emit progress messages via [message()].
#'
#' @return A data frame (manifest) with columns: `sample_id`, `sex`,
#'   `trisomy`, `n_reads`, `bam_file`.
#'
#' @details
#' The generated BAMs are NIPTeR-compatible (no unmapped reads, unique
#' positions per chromosome) and can be used with both the NIPTeR binning
#' layer ([nipter_bin_bam()]) and the WisecondorX native pipeline
#' ([rwisecondorx_newref()], [rwisecondorx_predict()]) when binned with
#' `binsize = COMPRESSED_BINSIZE` (100).
#'
#' The manifest is also written to `manifest.tsv` in `out_dir`.
#'
#' Requires `samtools` on PATH.
#'
#' @examples
#' \dontrun{
#' manifest <- generate_cohort(tempdir())
#' head(manifest)
#' }
#'
#' @export
generate_cohort <- function(out_dir, verbose = TRUE) {
  samtools <- Sys.which("samtools")
  if (nchar(samtools) == 0L) {
    stop("samtools not found on PATH", call. = FALSE)
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  profiles <- .build_cohort_profiles()
  header_text <- .sam_header()

  manifest <- data.frame(
    sample_id = character(0),
    sex       = character(0),
    trisomy   = character(0),
    n_reads   = integer(0),
    bam_file  = character(0),
    stringsAsFactors = FALSE
  )

  for (idx in seq_along(profiles)) {
    prof <- profiles[[idx]]
    seed <- 42L + idx
    bam_path <- file.path(out_dir, paste0(prof$id, ".bam"))

    if (verbose) {
      message(sprintf("  [%d/%d] %s (sex=%s)",
                      idx, length(profiles), prof$id, prof$sex))
    }

    records <- .sam_records(prof, seed)
    sam_text <- paste(c(header_text, records), collapse = "\n")

    # Pipe SAM through samtools to produce sorted BAM
    tmp_bam <- tempfile(fileext = ".bam")

    con <- pipe(sprintf("%s view -b -o %s -",
                        shQuote(samtools), shQuote(tmp_bam)), "wb")
    writeLines(sam_text, con)
    close(con)

    # Sort (should be sorted by construction, but be safe)
    system2(samtools, c("sort", "-o", shQuote(bam_path), shQuote(tmp_bam)),
            stdout = FALSE, stderr = FALSE)
    unlink(tmp_bam)

    # Index
    system2(samtools, c("index", shQuote(bam_path)),
            stdout = FALSE, stderr = FALSE)

    trisomy <- "none"
    if (grepl("trisomy_21", prof$id)) trisomy <- "T21"
    if (grepl("trisomy_18", prof$id)) trisomy <- "T18"
    if (grepl("trisomy_13", prof$id)) trisomy <- "T13"

    manifest <- rbind(manifest, data.frame(
      sample_id = prof$id,
      sex       = prof$sex,
      trisomy   = trisomy,
      n_reads   = length(records),
      bam_file  = basename(bam_path),
      stringsAsFactors = FALSE
    ))
  }

  # Write manifest
  manifest_path <- file.path(out_dir, "manifest.tsv")
  utils::write.table(manifest, manifest_path, sep = "\t", row.names = FALSE,
                     quote = FALSE)

  if (verbose) {
    total_size <- sum(file.size(
      file.path(out_dir, list.files(out_dir, pattern = "\\.bam$"))
    ))
    message(sprintf("Done. %d BAMs (%.1f MB). Manifest: %s",
                    nrow(manifest), total_size / 1024 / 1024, manifest_path))
  }

  manifest
}
