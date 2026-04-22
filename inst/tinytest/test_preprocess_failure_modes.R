library(tinytest)
library(RWisecondorX)

exit_if_not(nzchar(Sys.which("samtools")))

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.run_rscript <- function(args, env = character()) {
  out <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = args,
    stdout = TRUE,
    stderr = TRUE,
    env = env
  )
  list(
    status = as.integer(attr(out, "status") %||% 0L),
    output = out
  )
}

.read_one_row <- function(path) {
  stopifnot(file.exists(path))
  utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.write_toy_fasta <- function(path, chr_len = 100000L) {
  chroms <- c(as.character(seq_len(22L)), "X", "Y")
  seq <- paste(rep("ACGT", length.out = chr_len), collapse = "")
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  for (chr in chroms) {
    writeLines(paste0(">", chr), con)
    writeLines(seq, con)
  }
  invisible(path)
}

.write_toy_sam <- function(path, chr_len = 100000L, empty = FALSE) {
  chroms <- c(as.character(seq_len(22L)), "X", "Y")
  header <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    sprintf("@SQ\tSN:%s\tLN:%d", chroms, chr_len)
  )
  reads <- if (isTRUE(empty)) {
    character()
  } else {
    c(
      paste(
        "r1", "0", "1", "100", "60", "50M", "*", "0", "0",
        paste(rep("A", 50L), collapse = ""),
        paste(rep("I", 50L), collapse = ""),
        sep = "\t"
      ),
      paste(
        "r2", "0", "1", "60000", "60", "50M", "*", "0", "0",
        paste(rep("C", 50L), collapse = ""),
        paste(rep("I", 50L), collapse = ""),
        sep = "\t"
      ),
      paste(
        "r3", "0", "2", "100", "60", "50M", "*", "0", "0",
        paste(rep("G", 50L), collapse = ""),
        paste(rep("I", 50L), collapse = ""),
        sep = "\t"
      )
    )
  }
  writeLines(c(header, reads), path)
  invisible(path)
}

.make_toy_bam <- function(dir, stem = "tiny", empty = FALSE) {
  fasta <- file.path(dir, "toy_grch37.fa")
  sam <- file.path(dir, paste0(stem, ".sam"))
  bam_unsorted <- file.path(dir, paste0(stem, ".unsorted.bam"))
  bam <- file.path(dir, paste0(stem, ".bam"))

  .write_toy_fasta(fasta)
  .write_toy_sam(sam, empty = empty)

  system2("samtools", c("faidx", fasta))
  system2("samtools", c("view", "-bS", "-o", bam_unsorted, sam))
  system2("samtools", c("sort", "-o", bam, bam_unsorted))
  system2("samtools", c("index", bam))

  list(
    fasta = fasta,
    bam = bam
  )
}

.make_malformed_bam <- function(path) {
  writeBin(as.raw(c(0x42, 0x41, 0x44, 0x21, 0x00, 0xff)), path)
  invisible(path)
}

.script_path <- function(name) {
  root <- tinytest::get_call_wd()
  candidates <- c(
    file.path(root, "inst", "scripts", name),
    file.path(root, "scripts", name),
    file.path(root, "..", "scripts", name)
  )
  hit <- candidates[file.exists(candidates)][1L]
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

.canonical_failure_fields <- c(
  "sample_processing_status",
  "sample_failure_stage",
  "sample_failure_step",
  "sample_failure_message",
  "nipter_bed_status"
)

tmpdir <- tempfile("preprocess_failure_modes_")
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

convert_script <- .script_path("convert_sample.R")
preprocess_script <- .script_path("preprocess_cohort.R")

# ---------------------------------------------------------------------------
# convert_sample.R: malformed BAM still writes a failure QC row
# ---------------------------------------------------------------------------

malformed_dir <- file.path(tmpdir, "malformed")
dir.create(malformed_dir, recursive = TRUE, showWarnings = FALSE)
toy_ref_for_malformed <- .make_toy_bam(malformed_dir, stem = "refonly", empty = TRUE)$fasta
malformed_bam <- file.path(malformed_dir, "broken.bam")
.make_malformed_bam(malformed_bam)
malformed_qc <- file.path(malformed_dir, "broken.sample_qc.tsv")

malformed_run <- .run_rscript(c(
  convert_script,
  "--mode", "nipter",
  "--bam", malformed_bam,
  "--out", file.path(malformed_dir, "broken.bed.gz"),
  "--fasta", toy_ref_for_malformed,
  "--sample-qc-out", malformed_qc
))

expect_true(
  malformed_run$status != 0L,
  info = "convert_sample.R fails with a malformed BAM input"
)
expect_true(
  file.exists(malformed_qc),
  info = "convert_sample.R still writes sample_qc.tsv on malformed BAM failure"
)

malformed_row <- .read_one_row(malformed_qc)
expect_identical(malformed_row$sample_processing_status, "failed")
expect_identical(malformed_row$sample_failure_stage, "samtools_stats")
expect_identical(malformed_row$sample_failure_step, "summary")
expect_identical(malformed_row$nipter_bed_status, "failed")
expect_true(
  nzchar(malformed_row$sample_failure_message),
  info = "malformed BAM failure row preserves the failure message"
)

# ---------------------------------------------------------------------------
# convert_sample.R: empty BAMs fail at the same semantic stage as tiny BAMs,
# but with zero valid bins recorded in the message.
# ---------------------------------------------------------------------------

empty_dir <- file.path(tmpdir, "empty")
dir.create(empty_dir, recursive = TRUE, showWarnings = FALSE)
empty <- .make_toy_bam(empty_dir, stem = "empty", empty = TRUE)

empty_convert_qc <- file.path(empty_dir, "empty.sample_qc.tsv")
empty_convert_run <- .run_rscript(c(
  convert_script,
  "--mode", "nipter",
  "--bam", empty$bam,
  "--out", file.path(empty_dir, "empty.bed.gz"),
  "--fasta", empty$fasta,
  "--sample-qc-out", empty_convert_qc
))

expect_true(
  empty_convert_run$status != 0L,
  info = "convert_sample.R fails on an empty but structurally valid BAM"
)
expect_true(
  file.exists(empty_convert_qc),
  info = "convert_sample.R writes sample_qc.tsv for an empty valid BAM failure"
)

empty_convert_row <- .read_one_row(empty_convert_qc)
expect_identical(empty_convert_row$sample_processing_status, "failed")
expect_identical(empty_convert_row$sample_failure_stage, "nipter_bed")
expect_identical(empty_convert_row$sample_failure_step, "gc_correction")
expect_identical(empty_convert_row$nipter_bed_status, "failed")
expect_true(
  grepl("valid bins = 0", empty_convert_row$sample_failure_message),
  info = "empty BAM failure row records the zero-valid-bin GC-correction failure"
)

# ---------------------------------------------------------------------------
# tiny valid BAM: convert_sample.R and preprocess_cohort.R agree on the
# canonical failure fields when NIPTeR GC correction cannot proceed.
# ---------------------------------------------------------------------------

tiny_dir <- file.path(tmpdir, "tiny")
dir.create(tiny_dir, recursive = TRUE, showWarnings = FALSE)
tiny <- .make_toy_bam(tiny_dir, stem = "tiny", empty = FALSE)

tiny_convert_qc <- file.path(tiny_dir, "tiny.sample_qc.tsv")
tiny_convert_run <- .run_rscript(c(
  convert_script,
  "--mode", "nipter",
  "--bam", tiny$bam,
  "--out", file.path(tiny_dir, "tiny.bed.gz"),
  "--fasta", tiny$fasta,
  "--sample-qc-out", tiny_convert_qc
))

expect_true(
  tiny_convert_run$status != 0L,
  info = "convert_sample.R fails on a tiny but structurally valid BAM when NIPTeR GC correction lacks enough support"
)
expect_true(
  file.exists(tiny_convert_qc),
  info = "convert_sample.R writes sample_qc.tsv for a tiny valid BAM failure"
)

tiny_convert_row <- .read_one_row(tiny_convert_qc)
expect_identical(tiny_convert_row$sample_processing_status, "failed")
expect_identical(tiny_convert_row$sample_failure_stage, "nipter_bed")
expect_identical(tiny_convert_row$sample_failure_step, "gc_correction")
expect_identical(tiny_convert_row$nipter_bed_status, "failed")
expect_true(
  grepl("Too few valid autosomal bins", tiny_convert_row$sample_failure_message),
  info = "convert_sample.R exposes the NIPTeR GC-correction failure reason"
)

manifest <- file.path(tiny_dir, "tiny_manifest.txt")
writeLines(c("bam", tiny$bam), manifest)
cohort_out_root <- file.path(tiny_dir, "cohort_out")

tiny_cohort_run <- .run_rscript(c(
  preprocess_script,
  "--bam-list", manifest,
  "--out-root", cohort_out_root,
  "--fasta", tiny$fasta,
  "--jobs", "1",
  "--threads", "2",
  "--seqff=FALSE",
  "--overwrite=TRUE"
))

expect_identical(
  tiny_cohort_run$status,
  0L,
  info = "preprocess_cohort.R completes even when one sample fails downstream"
)

cohort_qc <- file.path(cohort_out_root, "sample_qc", "sample_qc.tsv")
cohort_failures <- file.path(cohort_out_root, "sample_qc", "sample_failures.tsv")
cohort_sample_qc <- file.path(cohort_out_root, "sample_qc", "tiny.sample_qc.tsv")

expect_true(file.exists(cohort_qc))
expect_true(file.exists(cohort_failures))
expect_true(file.exists(cohort_sample_qc))

tiny_cohort_row <- .read_one_row(cohort_sample_qc)
tiny_cohort_summary <- .read_one_row(cohort_qc)
tiny_cohort_failures <- .read_one_row(cohort_failures)

for (field in .canonical_failure_fields) {
  expect_identical(
    tiny_cohort_row[[field]],
    tiny_convert_row[[field]],
    info = sprintf(
      "convert_sample.R and preprocess_cohort.R agree on '%s' for the same tiny BAM failure",
      field
    )
  )
  expect_identical(
    tiny_cohort_summary[[field]],
    tiny_cohort_row[[field]],
    info = sprintf(
      "cohort sample_qc.tsv and per-sample QC row agree on '%s'",
      field
    )
  )
  expect_identical(
    tiny_cohort_failures[[field]],
    tiny_cohort_row[[field]],
    info = sprintf(
      "sample_failures.tsv and per-sample QC row agree on '%s'",
      field
    )
  )
}

expect_true(
  grepl("Too few valid autosomal bins", tiny_cohort_row$sample_failure_message),
  info = "preprocess_cohort.R preserves the concrete NIPTeR failure message for a tiny BAM"
)
