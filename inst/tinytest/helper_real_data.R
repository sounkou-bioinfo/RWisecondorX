.real_manifest_path <- function() {
  candidates <- c(
    Sys.getenv("RWISECONDORX_REAL_BAM_LIST", unset = NA_character_),
    Sys.getenv("NIPTER_REAL_BAM_LIST", unset = NA_character_)
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
  if (length(candidates)) candidates[[1L]] else NULL
}

.read_real_manifest <- function(path = .real_manifest_path()) {
  if (is.null(path) || !file.exists(path)) {
    return(character())
  }
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) >= 1L && identical(lines[[1L]], "bam")) {
    lines <- lines[-1L]
  }
  lines[file.exists(lines)]
}

.first_real_bam <- function() {
  direct <- Sys.getenv("RWISECONDORX_TEST_BAM", unset = NA_character_)
  if (!is.na(direct) && nzchar(direct) && file.exists(direct)) {
    return(direct)
  }
  paths <- .read_real_manifest()
  if (length(paths)) paths[[1L]] else NULL
}

.real_reference_fasta <- function() {
  candidates <- c(
    Sys.getenv("RWISECONDORX_REAL_FASTA", unset = NA_character_),
    Sys.getenv("NIPTER_REAL_FASTA", unset = NA_character_)
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
  if (length(candidates)) candidates[[1L]] else NULL
}
