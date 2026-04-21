library(tinytest)
library(RWisecondorX)

manifest_enable <- tolower(Sys.getenv("NIPTER_REAL_BAM_ENABLE", unset = "false"))
if (!manifest_enable %in% c("1", "true", "yes")) {
  exit_file("NIPTER_REAL_BAM_ENABLE not set; skipping internal real-BAM validation")
}

manifest_path <- Sys.getenv(
  "NIPTER_REAL_BAM_LIST",
  unset = "/mnt/data/BixCTF/NiptSeqNeo/all_bam_list_sample_500.txt"
)
if (!nzchar(manifest_path) || !file.exists(manifest_path)) {
  exit_file("real BAM manifest not available")
}

.read_manifest_paths <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) {
    stop("Manifest contains no BAM paths: ", path, call. = FALSE)
  }
  lines
}

.index_path <- function(path) {
  candidates <- c(
    paste0(path, ".bai"),
    sub("\\.bam$", ".bai", path, ignore.case = TRUE)
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits)) hits[[1L]] else NA_character_
}

paths <- .read_manifest_paths(manifest_path)
limit <- suppressWarnings(as.integer(Sys.getenv("NIPTER_REAL_BAM_LIMIT",
                                                unset = "50")))
if (is.na(limit) || limit < 10L) {
  stop("NIPTER_REAL_BAM_LIMIT must be >= 10 when set.", call. = FALSE)
}
paths <- paths[seq_len(min(length(paths), limit))]

missing <- paths[!file.exists(paths)]
if (length(missing)) {
  stop("Manifest contains missing BAMs:\n  ",
       paste(head(missing, 10L), collapse = "\n  "), call. = FALSE)
}

indexes <- vapply(paths, .index_path, character(1L))
expect_true(all(!is.na(indexes)),
            info = "all manifest BAMs have BAI indexes")

samples <- lapply(seq_along(paths), function(i) {
  path <- paths[[i]]
  nipter_bin_bam(
    bam = path,
    binsize = 50000L,
    mapq = 40L,
    exclude_flags = 1024L,
    rmdup = "none"
  )
})

sample_ids <- sprintf(
  "%03d_%s",
  seq_along(paths),
  sub("\\.(bam|cram)$", "", basename(paths), ignore.case = TRUE)
)
samples <- lapply(seq_along(samples), function(i) {
  RWisecondorX:::.nipt_sample_dollar_assign(samples[[i]], "sample_name", sample_ids[[i]])
})
names(samples) <- sample_ids
cg <- nipter_as_control_group(
  samples,
  description = sprintf("Real BAM manifest cohort (%d samples)", length(samples))
)
qc <- nipter_control_group_qc(cg, include_bins = FALSE)

expect_identical(n_controls(cg), length(paths),
                 info = "control group includes every manifest BAM")
expect_identical(length(qc$sample_names), length(paths),
                 info = "QC report includes every manifest BAM")
expect_true(all(as.character(1:22) %in% qc$chromosome_summary$chromosome),
            info = "QC chromosome summary covers all autosomes")
expect_true(nrow(qc$sample_summary) == length(paths),
            info = "QC sample summary has one row per manifest BAM")
expect_true(all(qc$sample_summary$sample_name %in% qc$sample_names),
            info = "QC sample summary names align with sample_names")
expect_true(all(vapply(cg$samples, function(s) {
  m <- autosomal_matrix(s)
  is.matrix(m) && is.numeric(m) && nrow(m) == 22L
}, logical(1L))),
info = "all real BAMs bin into numeric 22-row autosomal matrices")
