# Optional metadata helpers for our generic tabix-backed TSV/BED artifacts.
#
# Metadata is stored as leading comment lines:
#   ##RWX_<key>=<value>
#
# The data body remains headerless, so existing read_tabix()-based consumers
# keep working unchanged. The metadata is intended for provenance only.

.rwx_tabix_metadata_prefix <- "##RWX_"

.normalize_tabix_metadata_prefix <- function(prefix) {
  stopifnot(is.character(prefix), length(prefix) == 1L, nzchar(prefix))
  if (!startsWith(prefix, "#")) {
    prefix <- paste0("##", prefix)
  }
  prefix
}

.stringify_tabix_metadata_value <- function(value, key) {
  if (is.null(value)) {
    stop("Metadata value for key '", key, "' cannot be NULL.", call. = FALSE)
  }
  if (length(value) != 1L) {
    stop("Metadata value for key '", key, "' must have length 1.", call. = FALSE)
  }
  if (is.list(value)) {
    stop("Metadata value for key '", key, "' must be a scalar atomic value.", call. = FALSE)
  }
  out <- if (isTRUE(is.na(value))) {
    "NA"
  } else if (is.logical(value)) {
    if (isTRUE(value)) "true" else "false"
  } else {
    as.character(value)
  }
  if (grepl("[\r\n]", out)) {
    stop("Metadata value for key '", key, "' cannot contain newlines.", call. = FALSE)
  }
  out
}

.normalize_tabix_metadata <- function(metadata) {
  if (is.null(metadata)) {
    return(stats::setNames(character(), character()))
  }
  if (!(is.list(metadata) || is.atomic(metadata))) {
    stop("metadata must be NULL, a named list, or a named atomic vector.", call. = FALSE)
  }

  nms <- names(metadata)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) {
    stop("metadata must be named.", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("metadata keys must be unique.", call. = FALSE)
  }
  if (any(grepl("[=\t\r\n]", nms))) {
    stop("metadata keys may not contain '=', tabs, or newlines.", call. = FALSE)
  }

  vals <- vapply(
    seq_along(metadata),
    function(i) .stringify_tabix_metadata_value(metadata[[i]], nms[[i]]),
    character(1L)
  )
  stats::setNames(vals, nms)
}

.merge_tabix_metadata <- function(base = NULL, extra = NULL) {
  base_md <- .normalize_tabix_metadata(base)
  extra_md <- .normalize_tabix_metadata(extra)

  if (!length(base_md)) {
    return(extra_md)
  }
  if (!length(extra_md)) {
    return(base_md)
  }

  overlap <- intersect(names(base_md), names(extra_md))
  if (length(overlap)) {
    extra_md <- extra_md[setdiff(names(extra_md), overlap)]
  }
  c(base_md, extra_md)
}

.tabix_metadata_lines <- function(metadata, prefix = .rwx_tabix_metadata_prefix) {
  md <- .normalize_tabix_metadata(metadata)
  if (!length(md)) {
    return(character())
  }
  prefix <- .normalize_tabix_metadata_prefix(prefix)
  sprintf("%s%s=%s", prefix, names(md), unname(md))
}

.write_tabix_body <- function(df, path, metadata = NULL) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))

  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  meta_lines <- .tabix_metadata_lines(metadata)
  if (length(meta_lines)) {
    writeLines(meta_lines, con = con, sep = "\n", useBytes = TRUE)
  }

  utils::write.table(
    df,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  invisible(path)
}

#' Read RWisecondorX tabix metadata headers
#'
#' Reads optional provenance metadata stored in leading comment lines of a
#' generic tabix-indexed TSV/BED file. Metadata lines use the form
#' `##RWX_<key>=<value>`.
#'
#' This is intended for our bgzipped, tabix-indexed TSV artifacts such as
#' RWisecondorX and NIPTeR BED outputs. The data rows remain headerless; this
#' function only inspects the optional metadata comment block.
#'
#' @param path Path to a tabix-indexed TSV/BED file.
#' @param prefix Metadata line prefix to match. Defaults to `##RWX_`.
#' @param con Optional open DBI connection with duckhts already loaded.
#'
#' @return A named character vector of parsed metadata values. Returns an empty
#'   named character vector when no matching metadata lines are present.
#'
#' @export
tabix_metadata <- function(path,
                           prefix = .rwx_tabix_metadata_prefix,
                           con = NULL) {
  stopifnot(is.character(path), length(path) == 1L, nzchar(path))
  stopifnot(file.exists(path))
  prefix <- .normalize_tabix_metadata_prefix(prefix)

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  hdr <- Rduckhts::rduckhts_hts_header(
    con,
    normalizePath(path, winslash = "/", mustWork = TRUE),
    format = "tabix",
    mode = "raw"
  )
  if (!nrow(hdr) || !"raw" %in% names(hdr)) {
    return(stats::setNames(character(), character()))
  }

  lines <- hdr$raw[startsWith(hdr$raw, prefix)]
  if (!length(lines)) {
    return(stats::setNames(character(), character()))
  }

  payload <- substring(lines, nchar(prefix) + 1L)
  eq_pos <- regexpr("=", payload, fixed = TRUE)
  if (any(eq_pos < 1L)) {
    stop("Malformed metadata line(s) in tabix header.", call. = FALSE)
  }

  keys <- substr(payload, 1L, eq_pos - 1L)
  vals <- substring(payload, eq_pos + 1L)
  if (anyDuplicated(keys)) {
    stop("Duplicate metadata keys in tabix header: ",
         paste(unique(keys[duplicated(keys)]), collapse = ", "),
         call. = FALSE)
  }

  stats::setNames(vals, keys)
}
