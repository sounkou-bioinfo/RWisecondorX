# NIPTeR S7 class hierarchy
#
# CombinedStrandsSample  -+
#                          +- NIPTSample (abstract)
# SeparatedStrandsSample -+
#
# CombinedControlGroup  -+
#                         +- NIPTControlGroup (abstract)
# SeparatedControlGroup -+
#
# Key generics (dispatch on NIPTSample or NIPTControlGroup):
#   autosomal_matrix(x)         -> 22 x n_bins double matrix
#   sex_matrix(x)               -> 2  x n_bins double matrix (rows: X, Y)
#   strand_type(x)              -> "combined" or "separated"
#   n_bins(x)                   -> integer
#   fractions_auto(cg)          -> 22 x N fraction matrix (cached)
#   fractions_for_regression(cg)-> 22 or 44 x N matrix
#
# Correction state is tracked as typed correction-step objects inside
# NIPTCorrectionRecord. Human-readable labels are only produced at the
# compatibility/reporting boundary.

.numeric_matrix_property <- function(n_rows, label) {
  S7::new_property(
    S7::class_any,
    validator = function(value) {
      if (!is.matrix(value) || !is.numeric(value) || nrow(value) != n_rows) {
        sprintf("%s must be a %d-row numeric matrix", label, n_rows)
      }
    }
  )
}

.nipt_model_list_class <- S7::new_S3_class(
  "list",
  constructor = function(.data) .data
)

.nipt_model_df_class <- S7::new_S3_class(
  "data.frame",
  constructor = function(.data) .data
)

.nipt_correction_codes <- c("raw", "gc", "chi")
.nipt_correction_labels <- c(
  raw = "Uncorrected",
  gc = "GC Corrected",
  chi = "Chi square corrected"
)

.validate_correction_code <- function(value) {
  if (length(value) != 1L || !value %in% .nipt_correction_codes) {
    sprintf("code must be one of: %s",
            paste(.nipt_correction_codes, collapse = ", "))
  }
}


# ---- Correction steps --------------------------------------------------------

#' Typed correction step in a NIPTeR correction path
#'
#' @param code Correction step code. One of \code{"raw"}, \code{"gc"}, or
#'   \code{"chi"}.
#'
#' @export
NIPTCorrectionStep <- S7::new_class(
  "NIPTCorrectionStep",
  properties = list(
    code = S7::new_property(S7::class_character, validator = .validate_correction_code)
  )
)

.is_nipt_correction_step <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTCorrectionStep), error = function(e) FALSE)
}

.as_nipt_correction_step <- function(x) {
  if (.is_nipt_correction_step(x)) return(x)
  if (is.character(x) && length(x) == 1L) {
    if (x %in% names(.nipt_correction_labels)) {
      return(NIPTCorrectionStep(code = x))
    }
    label_match <- names(.nipt_correction_labels)[match(x, .nipt_correction_labels)]
    if (!is.na(label_match)) {
      return(NIPTCorrectionStep(code = unname(label_match)))
    }
  }
  stop("Cannot coerce value to NIPTCorrectionStep.", call. = FALSE)
}

.nipt_raw_correction_step <- function() NIPTCorrectionStep(code = "raw")
.nipt_gc_correction_step <- function() NIPTCorrectionStep(code = "gc")
.nipt_chi_correction_step <- function() NIPTCorrectionStep(code = "chi")

.nipt_raw_correction_path <- function() {
  NIPTCorrectionPath(steps = list(.nipt_raw_correction_step()))
}

.nipt_gc_correction_path <- function() {
  NIPTCorrectionPath(steps = list(.nipt_gc_correction_step()))
}

.nipt_gc_both_correction_record <- function() {
  NIPTCorrectionRecord(
    autosomal = .nipt_gc_correction_path(),
    sex = .nipt_gc_correction_path()
  )
}


#' Ordered correction path for autosomal or sex-chromosome processing
#'
#' @param steps List of \code{NIPTCorrectionStep} objects. Defaults to a single
#'   \code{"raw"} step.
#'
#' @export
NIPTCorrectionPath <- S7::new_class(
  "NIPTCorrectionPath",
  properties = list(
    steps = S7::class_list
  ),
  validator = function(self) {
    if (!length(self@steps)) {
      return("steps must contain at least one correction step")
    }
    if (!all(vapply(self@steps, .is_nipt_correction_step, logical(1L)))) {
      return("all steps must be NIPTCorrectionStep objects")
    }
    NULL
  }
)

.is_nipt_correction_path <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTCorrectionPath), error = function(e) FALSE)
}

.as_nipt_correction_path <- function(x) {
  if (.is_nipt_correction_path(x)) return(x)
  if (is.null(x)) return(.nipt_raw_correction_path())
  if (.is_nipt_correction_step(x)) return(NIPTCorrectionPath(steps = list(x)))
  if (is.character(x)) {
    steps <- lapply(as.list(x), .as_nipt_correction_step)
    return(NIPTCorrectionPath(steps = steps))
  }
  if (is.list(x)) {
    steps <- lapply(x, .as_nipt_correction_step)
    return(NIPTCorrectionPath(steps = steps))
  }
  stop("Cannot coerce value to NIPTCorrectionPath.", call. = FALSE)
}

.correction_path_codes <- function(path) {
  vapply(path@steps, function(step) step@code, character(1L))
}

.correction_path_labels <- function(path) {
  codes <- .correction_path_codes(path)
  if (length(codes) > 1L) {
    codes <- codes[codes != "raw"]
  }
  unname(.nipt_correction_labels[codes])
}

.append_correction_step <- function(path, step) {
  path <- .as_nipt_correction_path(path)
  step <- .as_nipt_correction_step(step)
  existing <- .correction_path_codes(path)
  if (step@code %in% existing) {
    return(path)
  }
  steps <- path@steps
  if (step@code != "raw") {
    keep <- vapply(steps, function(x) x@code != "raw", logical(1L))
    steps <- steps[keep]
  }
  NIPTCorrectionPath(steps = c(steps, list(step)))
}


# ---- Correction record -------------------------------------------------------

#' Correction record for a NIPTeR sample
#'
#' Tracks the sequence of corrections applied to the autosomal and sex
#' chromosome matrices of a \code{NIPTSample} object.
#'
#' @param autosomal A \code{NIPTCorrectionPath} describing the autosomal
#'   correction history. Defaults to a single \code{"raw"} step.
#' @param sex A \code{NIPTCorrectionPath} describing the sex-chromosome
#'   correction history. Defaults to a single \code{"raw"} step.
#'
#' @export
NIPTCorrectionRecord <- S7::new_class(
  "NIPTCorrectionRecord",
  properties = list(
    autosomal = S7::new_property(
      NIPTCorrectionPath,
      default = quote(.nipt_raw_correction_path())
    ),
    sex = S7::new_property(
      NIPTCorrectionPath,
      default = quote(.nipt_raw_correction_path())
    )
  )
)


# ---- Abstract NIPTSample -----------------------------------------------------

#' Abstract base class for NIPTeR bin-count samples
#'
#' Never instantiated directly. Subclasses are \code{CombinedStrandsSample}
#' (reads are summed across strands) and \code{SeparatedStrandsSample} (forward
#' and reverse counts stored independently).
#'
#' @param sample_name Sample identifier.
#' @param binsize Positive integer bin width in base pairs.
#' @param correction A \code{NIPTCorrectionRecord} describing the applied
#'   corrections.
#'
#' @export
NIPTSample <- S7::new_class(
  "NIPTSample",
  abstract = TRUE,
  properties = list(
    sample_name = S7::class_character,
    binsize     = S7::new_property(
      S7::class_integer,
      validator = function(value) {
        if (length(value) != 1L || value < 1L) "binsize must be a positive integer"
      }
    ),
    correction = S7::new_property(
      NIPTCorrectionRecord,
      default = quote(NIPTCorrectionRecord())
    )
  )
)


# ---- Generics ----------------------------------------------------------------

#' Extract the autosomal read-count matrix from a NIPTSample
#'
#' Returns a 22 x n_bins numeric matrix (rows = chromosomes "1"-"22",
#' cols = bin indices).  For \code{SeparatedStrandsSample}, forward and
#' reverse counts are summed.
#'
#' @param x A \code{NIPTSample} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return A 22-row numeric matrix.
#' @export
autosomal_matrix <- S7::new_generic("autosomal_matrix", "x")

#' Extract the sex-chromosome read-count matrix from a NIPTSample
#'
#' Returns a 2 x n_bins numeric matrix (row 1 = X, row 2 = Y).
#' For \code{SeparatedStrandsSample}, forward and reverse counts are summed.
#'
#' @param x A \code{NIPTSample} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return A 2-row numeric matrix.
#' @export
sex_matrix <- S7::new_generic("sex_matrix", "x")

#' Return the strand type of a NIPTSample
#'
#' @param x A \code{NIPTSample} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return \code{"combined"} or \code{"separated"}.
#' @export
strand_type <- S7::new_generic("strand_type", "x")

#' Number of bins in a NIPTSample
#'
#' @param x A \code{NIPTSample} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return A positive integer.
#' @export
n_bins <- S7::new_generic("n_bins", "x")

S7::method(n_bins, NIPTSample) <- function(x) ncol(autosomal_matrix(x))


# ---- CombinedStrandsSample ---------------------------------------------------

#' NIPTeR sample with combined forward+reverse read counts
#'
#' @inheritParams NIPTSample
#' @param auto_matrix 22 x n_bins numeric matrix for autosomes.
#' @param sex_matrix_ 2 x n_bins numeric matrix for X/Y.
#'
#' @slot auto_matrix  22 x n_bins numeric matrix (rows = chr 1-22).
#' @slot sex_matrix_  2 x n_bins numeric matrix (row 1 = X, row 2 = Y).
#'
#' @export
CombinedStrandsSample <- S7::new_class(
  "CombinedStrandsSample",
  parent = NIPTSample,
  properties = list(
    auto_matrix = .numeric_matrix_property(22L, "auto_matrix"),
    sex_matrix_ = .numeric_matrix_property(2L, "sex_matrix_")
  )
)

S7::method(autosomal_matrix, CombinedStrandsSample) <- function(x) x@auto_matrix
S7::method(sex_matrix,       CombinedStrandsSample) <- function(x) x@sex_matrix_
S7::method(strand_type,      CombinedStrandsSample) <- function(x) "combined"


# ---- SeparatedStrandsSample --------------------------------------------------

#' NIPTeR sample with separate forward and reverse read counts
#'
#' @inheritParams NIPTSample
#' @param auto_fwd 22 x n_bins numeric forward-strand autosomal matrix.
#' @param auto_rev 22 x n_bins numeric reverse-strand autosomal matrix.
#' @param sex_fwd 2 x n_bins numeric forward-strand X/Y matrix.
#' @param sex_rev 2 x n_bins numeric reverse-strand X/Y matrix.
#'
#' @slot auto_fwd  22 x n_bins numeric matrix (forward strand).
#' @slot auto_rev  22 x n_bins numeric matrix (reverse strand).
#' @slot sex_fwd   2 x n_bins numeric matrix (forward strand; row1=X, row2=Y).
#' @slot sex_rev   2 x n_bins numeric matrix (reverse strand).
#'
#' @export
SeparatedStrandsSample <- S7::new_class(
  "SeparatedStrandsSample",
  parent = NIPTSample,
  properties = list(
    auto_fwd = .numeric_matrix_property(22L, "auto_fwd"),
    auto_rev = .numeric_matrix_property(22L, "auto_rev"),
    sex_fwd  = .numeric_matrix_property(2L, "sex_fwd"),
    sex_rev  = .numeric_matrix_property(2L, "sex_rev")
  )
)

S7::method(autosomal_matrix, SeparatedStrandsSample) <- function(x) {
  m <- x@auto_fwd + x@auto_rev
  rownames(m) <- as.character(1:22)
  m
}

S7::method(sex_matrix, SeparatedStrandsSample) <- function(x) {
  m <- x@sex_fwd + x@sex_rev
  rownames(m) <- c("X", "Y")
  m
}

S7::method(strand_type, SeparatedStrandsSample) <- function(x) "separated"


# ---- Compatibility helpers ---------------------------------------------------

.legacy_field <- function(x, name) {
  .subset2(unclass(x), name)
}

.legacy_set_field <- function(x, name, value) {
  attrs <- attributes(x)
  y <- unclass(x)
  y[[name]] <- value
  attributes(y) <- attrs
  y
}

.is_s7_nipt_sample <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTSample), error = function(e) FALSE)
}

.is_s7_combined_sample <- function(x) {
  tryCatch(S7::S7_inherits(x, CombinedStrandsSample), error = function(e) FALSE)
}

.is_s7_separated_sample <- function(x) {
  tryCatch(S7::S7_inherits(x, SeparatedStrandsSample), error = function(e) FALSE)
}

.is_s7_nipt_control_group <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTControlGroup), error = function(e) FALSE)
}

.is_s7_separated_control_group <- function(x) {
  tryCatch(S7::S7_inherits(x, SeparatedControlGroup), error = function(e) FALSE)
}

.is_nipt_sample_object <- function(x) {
  inherits(x, "NIPTSample") || inherits(x, "NIPTeRSample") || .is_s7_nipt_sample(x)
}

.is_nipt_control_group_object <- function(x) {
  inherits(x, "NIPTControlGroup") ||
    inherits(x, "NIPTeRControlGroup") ||
    .is_s7_nipt_control_group(x)
}

.sample_autosomal_reads <- function(x) {
  if (.is_s7_combined_sample(x)) {
    return(list(x@auto_matrix))
  }
  if (.is_s7_separated_sample(x)) {
    return(list(x@auto_fwd, x@auto_rev))
  }
  .legacy_field(x, "autosomal_chromosome_reads")
}

.sample_sex_reads <- function(x) {
  if (.is_s7_combined_sample(x)) {
    return(list(x@sex_matrix_))
  }
  if (.is_s7_separated_sample(x)) {
    return(list(x@sex_fwd, x@sex_rev))
  }
  .legacy_field(x, "sex_chromosome_reads")
}

.sample_name <- function(x) {
  if (.is_s7_nipt_sample(x)) x@sample_name
  else .legacy_field(x, "sample_name")
}

.sample_binsize <- function(x) {
  if (.is_s7_nipt_sample(x)) x@binsize else NA_integer_
}

.sample_correction_status <- function(x,
                                      scope = c("autosomal", "sex")) {
  scope <- match.arg(scope)
  if (.is_s7_nipt_sample(x)) {
    path <- if (scope == "autosomal") x@correction@autosomal else x@correction@sex
    .correction_path_labels(path)
  } else {
    if (scope == "autosomal") .legacy_field(x, "correction_status_autosomal")
    else .legacy_field(x, "correction_status_sex")
  }
}

.sample_correction_path <- function(x,
                                    scope = c("autosomal", "sex")) {
  scope <- match.arg(scope)
  if (.is_s7_nipt_sample(x)) {
    if (scope == "autosomal") x@correction@autosomal else x@correction@sex
  } else {
    .as_nipt_correction_path(.sample_correction_status(x, scope))
  }
}

.sample_with_reads <- function(x, autosomal = NULL, sex = NULL) {
  if (!.is_s7_nipt_sample(x)) {
    if (!is.null(autosomal)) {
      x <- .legacy_set_field(x, "autosomal_chromosome_reads", autosomal)
    }
    if (!is.null(sex)) {
      x <- .legacy_set_field(x, "sex_chromosome_reads", sex)
    }
    return(x)
  }

  if (!is.null(autosomal)) {
    if (.is_s7_combined_sample(x)) {
      x@auto_matrix <- autosomal[[1L]]
    } else {
      x@auto_fwd <- autosomal[[1L]]
      x@auto_rev <- autosomal[[2L]]
    }
  }
  if (!is.null(sex)) {
    if (.is_s7_combined_sample(x)) {
      x@sex_matrix_ <- sex[[1L]]
    } else {
      x@sex_fwd <- sex[[1L]]
      x@sex_rev <- sex[[2L]]
    }
  }
  x
}

.sample_with_correction_status <- function(x,
                                           autosomal = NULL,
                                           sex = NULL) {
  if (!.is_s7_nipt_sample(x)) {
    if (!is.null(autosomal)) {
      x <- .legacy_set_field(
        x, "correction_status_autosomal",
        .correction_path_labels(.as_nipt_correction_path(autosomal))
      )
    }
    if (!is.null(sex)) {
      x <- .legacy_set_field(
        x, "correction_status_sex",
        .correction_path_labels(.as_nipt_correction_path(sex))
      )
    }
    return(x)
  }

  correction <- x@correction
  if (!is.null(autosomal)) correction@autosomal <- .as_nipt_correction_path(autosomal)
  if (!is.null(sex)) correction@sex <- .as_nipt_correction_path(sex)
  x@correction <- correction
  x
}

.sample_append_correction_step <- function(x,
                                           scope = c("autosomal", "sex"),
                                           step) {
  scope <- match.arg(scope)
  current <- .sample_correction_path(x, scope)
  updated <- .append_correction_step(current, step)
  if (scope == "autosomal") {
    .sample_with_correction_status(x, autosomal = updated)
  } else {
    .sample_with_correction_status(x, sex = updated)
  }
}

.control_group_correction_status <- function(cg,
                                             scope = c("autosomal", "sex")) {
  scope <- match.arg(scope)
  labels <- unlist(lapply(
    cg@samples,
    .sample_correction_status,
    scope = scope
  ), use.names = FALSE)
  unique(labels)
}

.nipt_sample_dollar <- function(x, name) {
  switch(
    name,
    autosomal_chromosome_reads  = .sample_autosomal_reads(x),
    sex_chromosome_reads        = .sample_sex_reads(x),
    correction_status_autosomal = .sample_correction_status(x, "autosomal"),
    correction_status_sex       = .sample_correction_status(x, "sex"),
    sample_name                 = .sample_name(x),
    binsize                     = .sample_binsize(x),
    NULL
  )
}

.nipt_sample_subset2 <- function(x, i, ...) {
  if (is.character(i) && length(i) == 1L) {
    return(.nipt_sample_dollar(x, i))
  }
  stop("NIPTSample only supports named [[ access.", call. = FALSE)
}

.nipt_sample_dollar_assign <- function(x, name, value) {
  switch(
    name,
    autosomal_chromosome_reads = {
      x <- .sample_with_reads(x, autosomal = value)
      x
    },
    sex_chromosome_reads = {
      x <- .sample_with_reads(x, sex = value)
      x
    },
    correction_status_autosomal = {
      x <- .sample_with_correction_status(x, autosomal = .as_nipt_correction_path(value))
      x
    },
    correction_status_sex = {
      x <- .sample_with_correction_status(x, sex = .as_nipt_correction_path(value))
      x
    },
    sample_name = {
      if (.is_s7_nipt_sample(x)) {
        x@sample_name <- value
      } else {
        x <- .legacy_set_field(x, "sample_name", value)
      }
      x
    },
    binsize = {
      if (.is_s7_nipt_sample(x)) {
        x@binsize <- value
      } else {
        x <- .legacy_set_field(x, "binsize", value)
      }
      x
    },
    stop("Unknown field '", name, "' for NIPTSample.", call. = FALSE)
  )
}


# ---- Backwards-compat generics for old S3 objects ---------------------------
# These allow all internal helpers to call autosomal_matrix() on both old S3
# NIPTeRSample lists and new S7 objects during the migration period.

.s3_NIPTeRSample <- S7::new_S3_class("NIPTeRSample")

S7::method(autosomal_matrix, .s3_NIPTeRSample) <- function(x) {
  auto <- x$autosomal_chromosome_reads
  if (inherits(x, "SeparatedStrands")) {
    m <- auto[[1L]] + auto[[2L]]
    rownames(m) <- as.character(1:22)
    storage.mode(m) <- "double"
    m
  } else {
    auto[[1L]]
  }
}

S7::method(sex_matrix, .s3_NIPTeRSample) <- function(x) {
  sx <- x$sex_chromosome_reads
  if (inherits(x, "SeparatedStrands")) {
    m <- sx[[1L]] + sx[[2L]]
    rownames(m) <- c("X", "Y")
    storage.mode(m) <- "double"
    m
  } else {
    sx[[1L]]
  }
}

S7::method(strand_type, .s3_NIPTeRSample) <- function(x) {
  if (inherits(x, "SeparatedStrands")) "separated" else "combined"
}

S7::method(n_bins, .s3_NIPTeRSample) <- function(x) {
  ncol(x$autosomal_chromosome_reads[[1L]])
}


# ---- NIPTControlGroup --------------------------------------------------------

#' Extract collapsed 22 x N autosomal chromosome fractions from a control group
#'
#' @param cg A \code{NIPTControlGroup} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return A 22 x N numeric matrix of chromosome fractions.
#' @export
fractions_auto <- S7::new_generic("fractions_auto", "cg")

#' Extract fractions suitable for regression
#'
#' For \code{CombinedControlGroup}: 22 x N matrix.
#' For \code{SeparatedControlGroup}: 44 x N matrix (rows "1F".."22F","1R".."22R").
#'
#' @param cg A \code{NIPTControlGroup} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return A numeric matrix.
#' @export
fractions_for_regression <- S7::new_generic("fractions_for_regression", "cg")

#' Number of samples in a control group
#' @param cg A \code{NIPTControlGroup} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return An integer.
#' @export
n_controls <- S7::new_generic("n_controls", "cg")

#' Sample names of a control group
#' @param cg A \code{NIPTControlGroup} object.
#' @param ... Reserved for S7 method dispatch; currently unused.
#' @return A character vector.
#' @export
control_names <- S7::new_generic("control_names", "cg")


#' Abstract base class for NIPTeR control groups
#'
#' @param samples List of \code{NIPTeRSample} or \code{NIPTSample} objects.
#' @param description Human-readable label for the control group.
#' @param sample_sex Optional named character vector keyed by sample name with
#'   values in \code{female}, \code{male}, \code{ambiguous}, or
#'   \code{unknown}.
#' @param sex_source Optional scalar string describing the origin of
#'   \code{sample_sex}.
#' @param .cache Internal environment used for lazy cached summaries.
#'
#' @export
NIPTControlGroup <- S7::new_class(
  "NIPTControlGroup",
  abstract = TRUE,
  properties = list(
    samples     = S7::class_list,
    description = S7::new_property(S7::class_character,
                                   default = "General control group"),
    sample_sex  = S7::new_property(
      S7::class_any,
      default = quote(NULL),
      validator = function(value) {
        if (!is.null(value) && !is.character(value)) {
          "'sample_sex' must be NULL or a character vector"
        }
      }
    ),
    sex_source = S7::new_property(
      S7::class_any,
      default = quote(NULL),
      validator = function(value) {
        if (!is.null(value) &&
            (!is.character(value) || length(value) != 1L || !nzchar(value))) {
          "'sex_source' must be NULL or a non-empty string"
        }
      }
    ),
    .cache = S7::new_property(
      S7::class_environment,
      default = quote(new.env(parent = emptyenv()))
    )
  ),
  validator = function(self) {
    if (length(self@samples) < 2L) {
      return("need at least 2 samples")
    }

    ok <- vapply(self@samples, .is_nipt_sample_object, logical(1L))
    if (!all(ok)) {
      return("all elements of 'samples' must be NIPTeRSample or NIPTSample objects")
    }

    strand_types <- vapply(self@samples, .strand_type_of, character(1L))
    if (length(unique(strand_types)) > 1L) {
      return("all samples in a control group must have the same strand type")
    }

    binsizes <- vapply(self@samples, .sample_binsize, integer(1L))
    binsizes <- binsizes[!is.na(binsizes)]
    if (length(unique(binsizes)) > 1L) {
      return("all S7 samples in a control group must share the same binsize")
    }

    n_bins_vec <- vapply(self@samples, n_bins, integer(1L))
    if (length(unique(n_bins_vec)) > 1L) {
      return("all samples in a control group must share the same number of bins")
    }

    if (!is.null(self@sample_sex)) {
      sample_names <- vapply(self@samples, .sample_name, character(1L))
      if (!length(self@sample_sex)) {
        return("'sample_sex' must be NULL or a non-empty character vector")
      }
      if (is.null(names(self@sample_sex))) {
        return("'sample_sex' must be a named character vector keyed by sample_name")
      }
      if (!setequal(names(self@sample_sex), sample_names)) {
        return("'sample_sex' names must exactly match the control-group sample names")
      }
      allowed <- c("female", "male", "ambiguous", "unknown")
      if (!all(self@sample_sex %in% allowed)) {
        return("'sample_sex' values must be one of: female, male, ambiguous, unknown")
      }
    }

    if (!is.null(self@sex_source) &&
        (length(self@sex_source) != 1L || !nzchar(self@sex_source))) {
      return("'sex_source' must be NULL or a non-empty string")
    }

    NULL
  }
)

S7::method(n_controls, NIPTControlGroup) <- function(cg) length(cg@samples)
S7::method(control_names, NIPTControlGroup) <- function(cg)
  vapply(cg@samples, function(s) {
    if (.is_s7_nipt_sample(s)) s@sample_name else .legacy_field(s, "sample_name")
  }, character(1L))


#' NIPTeR control group for CombinedStrands samples
#'
#' @inheritParams NIPTControlGroup
#' @export
CombinedControlGroup <- S7::new_class(
  "CombinedControlGroup",
  parent = NIPTControlGroup
)

#' NIPTeR control group for SeparatedStrands samples
#'
#' @inheritParams NIPTControlGroup
#' @export
SeparatedControlGroup <- S7::new_class(
  "SeparatedControlGroup",
  parent = NIPTControlGroup
)

# fractions_auto: 22 x N collapsed fractions
.compute_fractions_auto <- function(cg) {
  mats <- lapply(cg@samples, autosomal_matrix)
  totals <- vapply(mats, sum, numeric(1L))
  m <- mapply(function(mat, tot) rowSums(mat) / tot, mats, totals)
  rownames(m) <- as.character(1:22)
  colnames(m) <- control_names(cg)
  m  # 22 x N
}

S7::method(fractions_auto, CombinedControlGroup) <- function(cg) {
  if (exists("fractions_auto", envir = cg@.cache, inherits = FALSE)) {
    return(cg@.cache$fractions_auto)
  }
  fracs <- .compute_fractions_auto(cg)
  cg@.cache$fractions_auto <- fracs
  fracs
}

S7::method(fractions_auto, SeparatedControlGroup) <- function(cg) {
  if (exists("fractions_auto", envir = cg@.cache, inherits = FALSE)) {
    return(cg@.cache$fractions_auto)
  }
  fracs <- .compute_fractions_auto(cg)
  cg@.cache$fractions_auto <- fracs
  fracs
}

S7::method(fractions_for_regression, CombinedControlGroup) <- function(cg) {
  fractions_auto(cg)
}

S7::method(fractions_for_regression, SeparatedControlGroup) <- function(cg) {
  if (exists("fractions_regression", envir = cg@.cache, inherits = FALSE)) {
    return(cg@.cache$fractions_regression)
  }
  mats_fwd <- lapply(cg@samples, function(s) {
    if (.is_s7_separated_sample(s)) s@auto_fwd
    else s$autosomal_chromosome_reads[[1L]]
  })
  mats_rev <- lapply(cg@samples, function(s) {
    if (.is_s7_separated_sample(s)) s@auto_rev
    else s$autosomal_chromosome_reads[[2L]]
  })
  totals <- vapply(cg@samples, function(s) sum(autosomal_matrix(s)), numeric(1L))
  fwd <- mapply(function(m, t) rowSums(m) / t, mats_fwd, totals)
  rev <- mapply(function(m, t) rowSums(m) / t, mats_rev, totals)
  rownames(fwd) <- paste0(1:22, "F")
  rownames(rev) <- paste0(1:22, "R")
  out <- rbind(fwd, rev)
  colnames(out) <- control_names(cg)
  cg@.cache$fractions_regression <- out
  out
}


# ---- Backwards-compat methods for old S3 NIPTeRControlGroup -----------------

.s3_NIPTeRControlGroup <- S7::new_S3_class("NIPTeRControlGroup")

S7::method(n_controls, .s3_NIPTeRControlGroup) <- function(cg)
  length(cg$samples)

S7::method(control_names, .s3_NIPTeRControlGroup) <- function(cg)
  vapply(cg$samples, function(s) s$sample_name, character(1L))

S7::method(fractions_auto, .s3_NIPTeRControlGroup) <- function(cg) {
  # Delegate to the S3 helper in nipter_control.R
  .control_group_fractions_collapsed(cg)
}

S7::method(fractions_for_regression, .s3_NIPTeRControlGroup) <- function(cg) {
  # Delegate to the S3 helper in nipter_control.R
  .control_group_fractions(cg)
}


.normalize_sample_sex <- function(sample_sex, samples) {
  if (is.null(sample_sex)) return(NULL)

  sample_names <- vapply(samples, .sample_name, character(1L))
  if (!length(sample_sex)) {
    stop("'sample_sex' must be NULL or a non-empty character vector.",
         call. = FALSE)
  }
  if (is.null(names(sample_sex))) {
    if (length(sample_sex) != length(sample_names)) {
      stop(
        "Unnamed 'sample_sex' vectors must have length equal to the number of samples.",
        call. = FALSE
      )
    }
    names(sample_sex) <- sample_names
  } else {
    if (!all(sample_names %in% names(sample_sex))) {
      stop(
        "'sample_sex' names must cover every control-group sample name.",
        call. = FALSE
      )
    }
    sample_sex <- sample_sex[sample_names]
  }

  names(sample_sex) <- sample_names
  sample_sex
}

.control_sample_sex <- function(cg) {
  if (.is_s7_nipt_control_group(cg)) cg@sample_sex else NULL
}

.control_sex_source <- function(cg) {
  if (.is_s7_nipt_control_group(cg)) cg@sex_source else NULL
}

.control_with_sample_sex <- function(cg, sample_sex = NULL, sex_source = NULL) {
  if (!.is_s7_nipt_control_group(cg)) {
    return(cg)
  }
  cg@sample_sex <- .normalize_sample_sex(sample_sex, cg@samples)
  if (!is.null(sex_source)) cg@sex_source <- sex_source
  cg
}

.control_group_with_samples <- function(cg, samples) {
  if (!.is_s7_nipt_control_group(cg)) {
    cg$samples <- samples
    return(cg)
  }
  .nipt_control_group_dollar_assign(cg, "samples", samples)
}

.nipt_control_group_dollar <- function(x, name) {
  switch(
    name,
    samples = x@samples,
    Description = x@description,
    description = x@description,
    sample_sex = x@sample_sex,
    sex_source = x@sex_source,
    correction_status_autosomal = .control_group_correction_status(
      x, "autosomal"
    ),
    correction_status_sex = .control_group_correction_status(
      x, "sex"
    ),
    NULL
  )
}

.nipt_control_group_subset2 <- function(x, i, ...) {
  if (is.character(i) && length(i) == 1L) {
    return(.nipt_control_group_dollar(x, i))
  }
  stop("NIPTControlGroup only supports named [[ access.", call. = FALSE)
}

.nipt_control_group_dollar_assign <- function(x, name, value) {
  switch(
    name,
    samples = {
      x@samples <- value
      x@.cache <- new.env(parent = emptyenv())
      if (!is.null(x@sample_sex)) {
        x@sample_sex <- .normalize_sample_sex(x@sample_sex, x@samples)
      }
      x
    },
    Description = {
      x@description <- value
      x
    },
    description = {
      x@description <- value
      x
    },
    sample_sex = {
      x@sample_sex <- .normalize_sample_sex(value, x@samples)
      x
    },
    sex_source = {
      x@sex_source <- value
      x
    },
    correction_status_autosomal = x,
    correction_status_sex = x,
    stop("Unknown field '", name, "' for NIPTControlGroup.", call. = FALSE)
  )
}


# ---- NCVTemplate -------------------------------------------------------------

#' Pre-computed NCV denominator template
#'
#' Returned by \code{nipter_prepare_ncv()} and consumed by
#' \code{nipter_ncv_score()}.
#'
#' @param focus_chromosome Integer chromosome identifier of the numerator.
#' @param denominators Character vector of denominator chromosome labels.
#' @param ctrl_mean Mean control ratio for the chosen denominator set.
#' @param ctrl_sd Standard deviation of the control ratios.
#' @param ctrl_cv Coefficient of variation of the control ratios.
#' @param shapiro_p Shapiro-Wilk p-value for the control-ratio distribution.
#' @param test_z_scores Numeric vector of held-out control z-scores.
#' @param test_sample_names Names corresponding to \code{test_z_scores}.
#' @param train_sample_names Control-sample names used to fit the template.
#'
#' @export
NCVTemplate <- S7::new_class(
  "NCVTemplate",
  properties = list(
    focus_chromosome   = S7::class_integer,
    denominators       = S7::class_character,
    ctrl_mean          = S7::class_double,
    ctrl_sd            = S7::class_double,
    ctrl_cv            = S7::class_double,
    shapiro_p          = S7::class_double,
    test_z_scores      = S7::class_double,
    test_sample_names  = S7::class_character,
    train_sample_names = S7::class_character
  )
)


# ---- Reference/model objects -------------------------------------------------

.nipt_reference_chroms <- c(as.character(1:22), "X", "Y")
.nipt_reference_count_cols <- paste0("NChrReads_", .nipt_reference_chroms)
.nipt_reference_frac_cols <- paste0("FrChrReads_", .nipt_reference_chroms)

.is_nipter_sex_model <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTeRSexModel), error = function(e) FALSE)
}

.is_nipter_sex_prediction <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTeRSexPrediction), error = function(e) FALSE)
}

.is_nipt_reference_frame <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTReferenceFrame), error = function(e) FALSE)
}

.is_nipt_reference_model <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTReferenceModel), error = function(e) FALSE)
}

.is_nipt_sex_score <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTSexScore), error = function(e) FALSE)
}

.is_nipt_sex_ncv_model <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTSexNCVModel), error = function(e) FALSE)
}

.is_nipt_sex_regression_model <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTSexRegressionModel), error = function(e) FALSE)
}

.is_nipt_sex_ncv_score <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTSexNCVScore), error = function(e) FALSE)
}

.is_nipt_sex_regression_score <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTSexRegressionScore), error = function(e) FALSE)
}

.is_nipt_gaunosome_score <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTGaunosomeScore), error = function(e) FALSE)
}

.is_nipt_gaunosome_report <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTGaunosomeReport), error = function(e) FALSE)
}

.is_nipt_control_group_qc <- function(x) {
  tryCatch(S7::S7_inherits(x, NIPTControlGroupQC), error = function(e) FALSE)
}

.validate_reference_model_map <- function(x, predicate, label,
                                          nested_models = FALSE) {
  if (!is.list(x)) {
    return(sprintf("%s must be a list.", label))
  }
  sexes <- c("female", "male")
  focuses <- c("X", "Y")
  if (!all(sexes %in% names(x))) {
    return(sprintf("%s must contain female and male entries.", label))
  }
  for (sex in sexes) {
    if (!is.list(x[[sex]])) {
      return(sprintf("%s[['%s']] must be a list.", label, sex))
    }
    present_focuses <- intersect(focuses, names(x[[sex]]))
    if (!length(present_focuses)) {
      return(sprintf("%s[['%s']] must contain at least one of X or Y.", label, sex))
    }
    for (focus in present_focuses) {
      entry <- x[[sex]][[focus]]
      if (nested_models) {
        if (!is.list(entry) || !length(entry) ||
            !all(vapply(entry, predicate, logical(1L)))) {
          return(sprintf(
            "%s[['%s']][['%s']] must be a non-empty list of typed model objects.",
            label, sex, focus
          ))
        }
      } else {
        if (!predicate(entry)) {
          return(sprintf(
            "%s[['%s']][['%s']] must be a typed model object.",
            label, sex, focus
          ))
        }
      }
    }
  }
  NULL
}

.validate_reference_frame_df <- function(x) {
  if (!is.data.frame(x)) {
    return("NIPTReferenceFrame must be a data.frame.")
  }
  required <- c("Sample_name", .nipt_reference_count_cols, .nipt_reference_frac_cols)
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    return(sprintf(
      "NIPTReferenceFrame is missing required columns: %s",
      paste(missing, collapse = ", ")
    ))
  }
  if (!is.character(x$Sample_name) || anyNA(x$Sample_name) || anyDuplicated(x$Sample_name)) {
    return("NIPTReferenceFrame$Sample_name must be a unique character vector.")
  }
  if ("SampleSex" %in% names(x)) {
    allowed <- c("female", "male", "ambiguous", "unknown")
    if (!is.character(x$SampleSex) || anyNA(x$SampleSex) ||
        !all(x$SampleSex %in% allowed)) {
      return("NIPTReferenceFrame$SampleSex must contain only female, male, ambiguous, or unknown.")
    }
  }
  if ("ConsensusGender" %in% names(x)) {
    if (!is.character(x$ConsensusGender) || anyNA(x$ConsensusGender) ||
        !all(x$ConsensusGender %in% c("female", "male", "ambiguous", "unknown"))) {
      return("NIPTReferenceFrame$ConsensusGender must contain only female, male, ambiguous, or unknown.")
    }
  }
  for (col in .nipt_reference_count_cols) {
    if (!is.numeric(x[[col]])) {
      return(sprintf("NIPTReferenceFrame column '%s' must be numeric.", col))
    }
  }
  for (col in .nipt_reference_frac_cols) {
    if (!is.numeric(x[[col]])) {
      return(sprintf("NIPTReferenceFrame column '%s' must be numeric.", col))
    }
  }
  optional_numeric <- c("RR_X", "RR_Y", "RR_X_SexClassMAD", "RR_Y_SexClassMAD",
                        "YUniqueRatio")
  for (col in optional_numeric) {
    if (col %in% names(x) && !is.numeric(x[[col]])) {
      return(sprintf("NIPTReferenceFrame column '%s' must be numeric.", col))
    }
  }
  if ("IsRefSexOutlier" %in% names(x) &&
      (!is.logical(x$IsRefSexOutlier) || anyNA(x$IsRefSexOutlier))) {
    return("NIPTReferenceFrame$IsRefSexOutlier must be a logical vector without NA values.")
  }
  NULL
}

#' Chromosome-level NIPT reference frame
#'
#' A typed tabular training frame derived from a \code{NIPTControlGroup}. Each
#' row corresponds to one reference sample and includes chromosome-level read
#' counts and fractions for chromosomes \code{1:22}, \code{X}, and \code{Y}.
#'
#' @param .data Data frame payload for the typed reference frame constructor.
#'
#' @export
NIPTReferenceFrame <- S7::new_class(
  "NIPTReferenceFrame",
  parent = .nipt_model_df_class,
  validator = function(self) {
    .validate_reference_frame_df(self)
  }
)

#' NIPTeR sex-prediction model
#'
#' A list-like S7 object wrapping one trained sex-classification model.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTeRSexModel <- S7::new_class(
  "NIPTeRSexModel",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("model", "method", "male_cluster", "classifications", "fractions")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTeRSexModel is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$method) || length(self$method) != 1L ||
        !(self$method %in% c("y_fraction", "xy_fraction", "y_unique"))) {
      return("NIPTeRSexModel$method must be one of y_fraction, xy_fraction, or y_unique.")
    }
    if (!is.numeric(self$male_cluster) || length(self$male_cluster) != 1L ||
        !(as.integer(self$male_cluster) %in% c(1L, 2L))) {
      return("NIPTeRSexModel$male_cluster must be 1 or 2.")
    }
    if (!is.character(self$classifications) || !length(self$classifications) ||
        is.null(names(self$classifications)) ||
        !all(self$classifications %in% c("male", "female"))) {
      return("NIPTeRSexModel$classifications must be a named male/female character vector.")
    }
    if (identical(self$method, "y_unique")) {
      if (!is.numeric(self$fractions) || !length(self$fractions)) {
        return("Y-unique NIPTeRSexModel$fractions must be a numeric vector.")
      }
      if (is.null(names(self$fractions))) {
        return("Y-unique NIPTeRSexModel$fractions must be named by sample.")
      }
    } else {
      if (!is.matrix(self$fractions) || !is.numeric(self$fractions) ||
          !identical(colnames(self$fractions), c("x_fraction", "y_fraction"))) {
        return("Fraction-based NIPTeRSexModel$fractions must be a numeric matrix with x_fraction and y_fraction columns.")
      }
      if (is.null(rownames(self$fractions))) {
        return("Fraction-based NIPTeRSexModel$fractions must have sample row names.")
      }
    }
    NULL
  }
)

#' NIPTeR sex prediction result
#'
#' A list-like S7 object containing the consensus sex call for one sample.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTeRSexPrediction <- S7::new_class(
  "NIPTeRSexPrediction",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("prediction", "model_predictions", "fractions", "sample_name")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTeRSexPrediction is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$prediction) || length(self$prediction) != 1L ||
        !(self$prediction %in% c("male", "female"))) {
      return("NIPTeRSexPrediction$prediction must be male or female.")
    }
    if (!is.character(self$model_predictions) || !length(self$model_predictions) ||
        is.null(names(self$model_predictions))) {
      return("NIPTeRSexPrediction$model_predictions must be a named character vector.")
    }
    valid_preds <- stats::na.omit(self$model_predictions)
    if (length(valid_preds) && !all(valid_preds %in% c("male", "female"))) {
      return("NIPTeRSexPrediction$model_predictions must contain male, female, or NA values.")
    }
    if (!is.numeric(self$fractions) ||
        !identical(names(self$fractions), c("x_fraction", "y_fraction"))) {
      return("NIPTeRSexPrediction$fractions must be a named numeric vector with x_fraction and y_fraction.")
    }
    if (!is.character(self$sample_name) || length(self$sample_name) != 1L || !nzchar(self$sample_name)) {
      return("NIPTeRSexPrediction$sample_name must be a non-empty string.")
    }
    NULL
  }
)

#' Built NIPT reference model
#'
#' A validated, serialisable package-level reference object that bundles the
#' control group, chromosome-level reference frame, and sex-prediction models.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTReferenceModel <- S7::new_class(
  "NIPTReferenceModel",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("control_group", "reference_frame", "sex_models",
                  "sample_sex_source", "build_date", "build_params")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTReferenceModel is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!.is_nipt_control_group_object(self$control_group)) {
      return("NIPTReferenceModel$control_group must be a NIPTControlGroup.")
    }
    if (!.is_nipt_reference_frame(self$reference_frame)) {
      return("NIPTReferenceModel$reference_frame must be a NIPTReferenceFrame.")
    }
    if (!is.list(self$sex_models)) {
      return("NIPTReferenceModel$sex_models must be a list.")
    }
    if (length(self$sex_models) &&
        !all(vapply(self$sex_models, .is_nipter_sex_model, logical(1L)))) {
      return("NIPTReferenceModel$sex_models must contain only NIPTeRSexModel objects.")
    }
    if (!is.null(self$sample_sex_source) &&
        (!is.character(self$sample_sex_source) || length(self$sample_sex_source) != 1L)) {
      return("NIPTReferenceModel$sample_sex_source must be NULL or a scalar character string.")
    }
    if (!is.character(self$build_date) || length(self$build_date) != 1L || !nzchar(self$build_date)) {
      return("NIPTReferenceModel$build_date must be a non-empty string.")
    }
    if (!is.list(self$build_params)) {
      return("NIPTReferenceModel$build_params must be a list.")
    }
    if ("sex_ncv_models" %in% names(self)) {
      msg <- .validate_reference_model_map(
        self$sex_ncv_models,
        predicate = .is_nipt_sex_ncv_model,
        label = "NIPTReferenceModel$sex_ncv_models"
      )
      if (!is.null(msg)) {
        return(msg)
      }
    }
    if ("sex_regression_models" %in% names(self)) {
      msg <- .validate_reference_model_map(
        self$sex_regression_models,
        predicate = .is_nipt_sex_regression_model,
        label = "NIPTReferenceModel$sex_regression_models",
        nested_models = TRUE
      )
      if (!is.null(msg)) {
        return(msg)
      }
    }
    NULL
  }
)

#' NIPTeR sex-chromosome score
#'
#' A list-like S7 object containing sex-matched X/Y z-scores and related
#' reference statistics for one sample.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTSexScore <- S7::new_class(
  "NIPTSexScore",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("sample_name", "predicted_sex", "sex_prediction",
                  "sample_metrics", "z_scores", "cv",
                  "reference_sizes", "reference_sample_names")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTSexScore is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$sample_name) || length(self$sample_name) != 1L ||
        !nzchar(self$sample_name)) {
      return("NIPTSexScore$sample_name must be a non-empty string.")
    }
    if (!is.character(self$predicted_sex) || length(self$predicted_sex) != 1L ||
        !(self$predicted_sex %in% c("female", "male"))) {
      return("NIPTSexScore$predicted_sex must be female or male.")
    }
    if (!.is_nipter_sex_prediction(self$sex_prediction)) {
      return("NIPTSexScore$sex_prediction must be a NIPTeRSexPrediction.")
    }
    if (!is.numeric(self$sample_metrics) ||
        !identical(names(self$sample_metrics),
                   c("FrChrReads_X", "FrChrReads_Y", "RR_X", "RR_Y"))) {
      return("NIPTSexScore$sample_metrics must be a named numeric vector with FrChrReads_X, FrChrReads_Y, RR_X, and RR_Y.")
    }
    expected_z <- c("Z_FrChrReads_X", "Z_FrChrReads_Y",
                    "Z_FrChrReads_X_XX", "Z_FrChrReads_Y_XX",
                    "Z_FrChrReads_X_XY", "Z_FrChrReads_Y_XY")
    if (!is.numeric(self$z_scores) || !identical(names(self$z_scores), expected_z)) {
      return("NIPTSexScore$z_scores must be a named numeric vector of the selected, XX, and XY sex-chromosome z-scores.")
    }
    expected_cv <- c("Z_FrChrReads_CV_X", "Z_FrChrReads_CV_Y",
                     "Z_FrChrReads_CV_X_XX", "Z_FrChrReads_CV_Y_XX",
                     "Z_FrChrReads_CV_X_XY", "Z_FrChrReads_CV_Y_XY")
    if (!is.numeric(self$cv) || !identical(names(self$cv), expected_cv)) {
      return("NIPTSexScore$cv must be a named numeric vector of the selected, XX, and XY coefficient-of-variation values.")
    }
    if (!is.numeric(self$reference_sizes) ||
        !identical(names(self$reference_sizes),
                   c("female", "male", "same_sex"))) {
      return("NIPTSexScore$reference_sizes must be a named numeric vector with female, male, and same_sex counts.")
    }
    if (!is.character(self$reference_sample_names)) {
      return("NIPTSexScore$reference_sample_names must be a character vector.")
    }
    NULL
  }
)

#' Sex-chromosome NCV model
#'
#' Typed denominator template for X/Y NCV scoring against a sex-matched
#' reference subset.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTSexNCVModel <- S7::new_class(
  "NIPTSexNCVModel",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("focus_chromosome", "reference_sex", "denominators",
                  "control_statistics", "reference_sample_names")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTSexNCVModel is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$focus_chromosome) || length(self$focus_chromosome) != 1L ||
        !(self$focus_chromosome %in% c("X", "Y"))) {
      return("NIPTSexNCVModel$focus_chromosome must be X or Y.")
    }
    if (!is.character(self$reference_sex) || length(self$reference_sex) != 1L ||
        !(self$reference_sex %in% c("female", "male"))) {
      return("NIPTSexNCVModel$reference_sex must be female or male.")
    }
    if (!is.character(self$denominators) || !length(self$denominators)) {
      return("NIPTSexNCVModel$denominators must be a non-empty character vector.")
    }
    expected_stats <- c("mean", "sd", "cv", "shapiro_p_value")
    if (!is.numeric(self$control_statistics) ||
        !identical(names(self$control_statistics), expected_stats)) {
      return("NIPTSexNCVModel$control_statistics must be a named numeric vector with mean, sd, cv, and shapiro_p_value.")
    }
    if (!is.character(self$reference_sample_names) || !length(self$reference_sample_names)) {
      return("NIPTSexNCVModel$reference_sample_names must be a non-empty character vector.")
    }
    NULL
  }
)

#' Sex-chromosome regression model
#'
#' Typed sex-matched ratio model for X/Y regression-based scoring.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTSexRegressionModel <- S7::new_class(
  "NIPTSexRegressionModel",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("fit", "focus_chromosome", "reference_sex", "response_column",
                  "predictors", "control_statistics", "reference_sample_names")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTSexRegressionModel is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!inherits(self$fit, "lm")) {
      return("NIPTSexRegressionModel$fit must be an lm object.")
    }
    if (!is.character(self$focus_chromosome) || length(self$focus_chromosome) != 1L ||
        !(self$focus_chromosome %in% c("X", "Y"))) {
      return("NIPTSexRegressionModel$focus_chromosome must be X or Y.")
    }
    if (!is.character(self$reference_sex) || length(self$reference_sex) != 1L ||
        !(self$reference_sex %in% c("female", "male"))) {
      return("NIPTSexRegressionModel$reference_sex must be female or male.")
    }
    if (!is.character(self$response_column) || length(self$response_column) != 1L ||
        !nzchar(self$response_column)) {
      return("NIPTSexRegressionModel$response_column must be a non-empty string.")
    }
    if (!is.character(self$predictors) || !length(self$predictors)) {
      return("NIPTSexRegressionModel$predictors must be a non-empty character vector.")
    }
    expected_stats <- c("mean_ratio", "sd_ratio", "cv",
                        "shapiro_p_value", "adj_r_squared")
    if (!is.numeric(self$control_statistics) ||
        !identical(names(self$control_statistics), expected_stats)) {
      return("NIPTSexRegressionModel$control_statistics must be a named numeric vector with mean_ratio, sd_ratio, cv, shapiro_p_value, and adj_r_squared.")
    }
    if (!is.character(self$reference_sample_names) || !length(self$reference_sample_names)) {
      return("NIPTSexRegressionModel$reference_sample_names must be a non-empty character vector.")
    }
    NULL
  }
)

#' Sex-chromosome NCV score
#'
#' Typed result for X/Y NCV scoring across male and female reference models.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTSexNCVScore <- S7::new_class(
  "NIPTSexNCVScore",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("sample_name", "focus_chromosome", "predicted_sex",
                  "sex_prediction", "sample_scores", "denominators",
                  "control_statistics", "reference_sizes")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTSexNCVScore is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$sample_name) || length(self$sample_name) != 1L || !nzchar(self$sample_name)) {
      return("NIPTSexNCVScore$sample_name must be a non-empty string.")
    }
    if (!is.character(self$focus_chromosome) || length(self$focus_chromosome) != 1L ||
        !(self$focus_chromosome %in% c("X", "Y"))) {
      return("NIPTSexNCVScore$focus_chromosome must be X or Y.")
    }
    if (!is.character(self$predicted_sex) || length(self$predicted_sex) != 1L ||
        !(self$predicted_sex %in% c("female", "male"))) {
      return("NIPTSexNCVScore$predicted_sex must be female or male.")
    }
    if (!.is_nipter_sex_prediction(self$sex_prediction)) {
      return("NIPTSexNCVScore$sex_prediction must be a NIPTeRSexPrediction.")
    }
    if (!is.numeric(self$sample_scores) ||
        !identical(names(self$sample_scores), c("female", "male", "selected"))) {
      return("NIPTSexNCVScore$sample_scores must be a named numeric vector with female, male, and selected scores.")
    }
    if (!is.list(self$denominators) ||
        !identical(sort(names(self$denominators)), c("female", "male"))) {
      return("NIPTSexNCVScore$denominators must be a list with female and male entries.")
    }
    if (!is.list(self$control_statistics) ||
        !identical(sort(names(self$control_statistics)), c("female", "male"))) {
      return("NIPTSexNCVScore$control_statistics must be a list with female and male entries.")
    }
    if (!is.numeric(self$reference_sizes) ||
        !identical(names(self$reference_sizes), c("female", "male", "same_sex"))) {
      return("NIPTSexNCVScore$reference_sizes must be a named numeric vector with female, male, and same_sex counts.")
    }
    NULL
  }
)

#' Sex-chromosome regression score
#'
#' Typed result for X/Y regression-based scoring across male and female
#' reference models.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTSexRegressionScore <- S7::new_class(
  "NIPTSexRegressionScore",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("sample_name", "focus_chromosome", "predicted_sex",
                  "sex_prediction", "scores", "ratios", "predictors",
                  "aggregate_scores", "reference_sizes")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTSexRegressionScore is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$sample_name) || length(self$sample_name) != 1L || !nzchar(self$sample_name)) {
      return("NIPTSexRegressionScore$sample_name must be a non-empty string.")
    }
    if (!is.character(self$focus_chromosome) || length(self$focus_chromosome) != 1L ||
        !(self$focus_chromosome %in% c("X", "Y"))) {
      return("NIPTSexRegressionScore$focus_chromosome must be X or Y.")
    }
    if (!is.character(self$predicted_sex) || length(self$predicted_sex) != 1L ||
        !(self$predicted_sex %in% c("female", "male"))) {
      return("NIPTSexRegressionScore$predicted_sex must be female or male.")
    }
    if (!.is_nipter_sex_prediction(self$sex_prediction)) {
      return("NIPTSexRegressionScore$sex_prediction must be a NIPTeRSexPrediction.")
    }
    expected_scores <- c("female", "male")
    if (!is.list(self$scores) || !identical(sort(names(self$scores)), expected_scores)) {
      return("NIPTSexRegressionScore$scores must be a list with female and male entries.")
    }
    if (!is.list(self$ratios) || !identical(sort(names(self$ratios)), expected_scores)) {
      return("NIPTSexRegressionScore$ratios must be a list with female and male entries.")
    }
    if (!is.list(self$predictors) || !identical(sort(names(self$predictors)), expected_scores)) {
      return("NIPTSexRegressionScore$predictors must be a list with female and male entries.")
    }
    if (!is.numeric(self$aggregate_scores) ||
        !identical(names(self$aggregate_scores), c("female_mean", "male_mean", "selected_mean"))) {
      return("NIPTSexRegressionScore$aggregate_scores must be a named numeric vector with female_mean, male_mean, and selected_mean.")
    }
    if (!is.numeric(self$reference_sizes) ||
        !identical(names(self$reference_sizes), c("female", "male", "same_sex"))) {
      return("NIPTSexRegressionScore$reference_sizes must be a named numeric vector with female, male, and same_sex counts.")
    }
    NULL
  }
)

#' Aggregate gaunosome score report
#'
#' Typed package-level report bundling sex, NCV, and regression-based X/Y
#' scoring against a built \code{NIPTReferenceModel}.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTGaunosomeScore <- S7::new_class(
  "NIPTGaunosomeScore",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("sample_name", "predicted_sex", "sex_prediction",
                  "sex_score", "ncv_scores", "regression_scores",
                  "summary")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTGaunosomeScore is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$sample_name) || length(self$sample_name) != 1L ||
        !nzchar(self$sample_name)) {
      return("NIPTGaunosomeScore$sample_name must be a non-empty string.")
    }
    if (!is.character(self$predicted_sex) || length(self$predicted_sex) != 1L ||
        !(self$predicted_sex %in% c("female", "male"))) {
      return("NIPTGaunosomeScore$predicted_sex must be female or male.")
    }
    if (!.is_nipter_sex_prediction(self$sex_prediction)) {
      return("NIPTGaunosomeScore$sex_prediction must be a NIPTeRSexPrediction.")
    }
    if (!.is_nipt_sex_score(self$sex_score)) {
      return("NIPTGaunosomeScore$sex_score must be a NIPTSexScore.")
    }
    if (!is.list(self$ncv_scores)) {
      return("NIPTGaunosomeScore$ncv_scores must be a list.")
    }
    if (length(self$ncv_scores) &&
        !all(vapply(self$ncv_scores, .is_nipt_sex_ncv_score, logical(1L)))) {
      return("NIPTGaunosomeScore$ncv_scores must contain only NIPTSexNCVScore objects.")
    }
    if (!is.list(self$regression_scores)) {
      return("NIPTGaunosomeScore$regression_scores must be a list.")
    }
    if (length(self$regression_scores) &&
        !all(vapply(self$regression_scores, .is_nipt_sex_regression_score, logical(1L)))) {
      return("NIPTGaunosomeScore$regression_scores must contain only NIPTSexRegressionScore objects.")
    }
    if (!is.data.frame(self$summary)) {
      return("NIPTGaunosomeScore$summary must be a data.frame.")
    }
    required_cols <- c(
      "chromosome", "predicted_sex", "z_score", "cv",
      "ncv_score_female", "ncv_score_male", "ncv_score_selected",
      "regression_score_female", "regression_score_male",
      "regression_score_selected"
    )
    missing_cols <- setdiff(required_cols, names(self$summary))
    if (length(missing_cols)) {
      return(sprintf(
        "NIPTGaunosomeScore$summary is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ))
    }
    NULL
  }
)

#' Aggregate gaunosome cohort report
#'
#' Typed package-level report bundling multiple \code{NIPTGaunosomeScore}
#' objects and their flattened per-sample summary table.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTGaunosomeReport <- S7::new_class(
  "NIPTGaunosomeReport",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("scores", "summary", "sample_names", "focus_chromosomes")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTGaunosomeReport is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.list(self$scores) || !length(self$scores)) {
      return("NIPTGaunosomeReport$scores must be a non-empty list.")
    }
    if (!all(vapply(self$scores, .is_nipt_gaunosome_score, logical(1L)))) {
      return("NIPTGaunosomeReport$scores must contain only NIPTGaunosomeScore objects.")
    }
    if (is.null(names(self$scores)) || any(!nzchar(names(self$scores)))) {
      return("NIPTGaunosomeReport$scores must be named by sample.")
    }
    if (!is.character(self$sample_names) || !length(self$sample_names)) {
      return("NIPTGaunosomeReport$sample_names must be a non-empty character vector.")
    }
    if (!identical(unname(self$sample_names), unname(names(self$scores)))) {
      return("NIPTGaunosomeReport$sample_names must align with names(NIPTGaunosomeReport$scores).")
    }
    if (!is.character(self$focus_chromosomes) || !length(self$focus_chromosomes) ||
        !all(self$focus_chromosomes %in% c("X", "Y"))) {
      return("NIPTGaunosomeReport$focus_chromosomes must contain X and/or Y.")
    }
    if (!is.data.frame(self$summary)) {
      return("NIPTGaunosomeReport$summary must be a data.frame.")
    }
    required_cols <- c(
      "sample_name", "chromosome", "predicted_sex", "z_score", "cv",
      "ncv_score_female", "ncv_score_male", "ncv_score_selected",
      "regression_score_female", "regression_score_male",
      "regression_score_selected"
    )
    missing_cols <- setdiff(required_cols, names(self$summary))
    if (length(missing_cols)) {
      return(sprintf(
        "NIPTGaunosomeReport$summary is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ))
    }
    NULL
  }
)

.validate_qc_summary_df <- function(x, required_cols, label) {
  if (!is.data.frame(x)) {
    return(sprintf("%s must be a data.frame.", label))
  }
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols)) {
    return(sprintf(
      "%s is missing required columns: %s",
      label,
      paste(missing_cols, collapse = ", ")
    ))
  }
  NULL
}

#' Control-group QC report
#'
#' Typed QC object bundling chromosome-level, sample-level, sex-stratified,
#' and optional chi/bin-level summaries for a \code{NIPTControlGroup}.
#'
#' @param .data Named list payload for the typed constructor.
#'
#' @export
NIPTControlGroupQC <- S7::new_class(
  "NIPTControlGroupQC",
  parent = .nipt_model_list_class,
  validator = function(self) {
    required <- c("sample_names", "chromosome_summary", "sample_summary",
                  "sex_summary", "bin_summary", "settings")
    missing <- setdiff(required, names(self))
    if (length(missing)) {
      return(sprintf("NIPTControlGroupQC is missing required fields: %s",
                     paste(missing, collapse = ", ")))
    }
    if (!is.character(self$sample_names) || !length(self$sample_names)) {
      return("NIPTControlGroupQC$sample_names must be a non-empty character vector.")
    }

    msg <- .validate_qc_summary_df(
      self$chromosome_summary,
      c("chromosome", "mean_fraction", "sd_fraction",
        "cv_fraction", "shapiro_p_value"),
      "NIPTControlGroupQC$chromosome_summary"
    )
    if (!is.null(msg)) return(msg)

    msg <- .validate_qc_summary_df(
      self$sample_summary,
      c("sample_name", "mean_ssd", "max_abs_z", "is_matching_outlier"),
      "NIPTControlGroupQC$sample_summary"
    )
    if (!is.null(msg)) return(msg)

    if (!identical(as.character(self$sample_summary$sample_name), self$sample_names)) {
      return("NIPTControlGroupQC$sample_summary$sample_name must align with sample_names.")
    }
    if (!is.logical(self$sample_summary$is_matching_outlier) ||
        anyNA(self$sample_summary$is_matching_outlier)) {
      return("NIPTControlGroupQC$sample_summary$is_matching_outlier must be logical without NA values.")
    }

    if (!is.null(self$sex_summary)) {
      msg <- .validate_qc_summary_df(
        self$sex_summary,
        c("sex_group", "n_samples",
          "rr_x_mean", "rr_x_sd", "rr_x_cv",
          "rr_y_mean", "rr_y_sd", "rr_y_cv",
          "rr_x_outliers", "rr_y_outliers",
          "can_build_ncv", "can_build_regression"),
        "NIPTControlGroupQC$sex_summary"
      )
      if (!is.null(msg)) return(msg)
      if (!is.logical(self$sex_summary$can_build_ncv) ||
          !is.logical(self$sex_summary$can_build_regression)) {
        return("NIPTControlGroupQC$sex_summary build flags must be logical.")
      }
    }

    if (!is.null(self$bin_summary)) {
      msg <- .validate_qc_summary_df(
        self$bin_summary,
        c("chromosome", "bin", "mean_scaled", "sd_scaled", "cv_scaled",
          "expected_count", "chi_square", "chi_z",
          "overdispersed", "correction_factor"),
        "NIPTControlGroupQC$bin_summary"
      )
      if (!is.null(msg)) return(msg)
      if (!is.logical(self$bin_summary$overdispersed) ||
          anyNA(self$bin_summary$overdispersed)) {
        return("NIPTControlGroupQC$bin_summary$overdispersed must be logical without NA values.")
      }
    }

    if (!is.list(self$settings)) {
      return("NIPTControlGroupQC$settings must be a list.")
    }
    NULL
  }
)

.as_nipt_reference_frame <- function(x) {
  if (.is_nipt_reference_frame(x)) x else NIPTReferenceFrame(x)
}

.as_nipter_sex_model <- function(x) {
  if (.is_nipter_sex_model(x)) x else NIPTeRSexModel(x)
}

.as_nipter_sex_prediction <- function(x) {
  if (.is_nipter_sex_prediction(x)) x else NIPTeRSexPrediction(x)
}

.as_nipt_reference_model <- function(x) {
  if (.is_nipt_reference_model(x)) x else NIPTReferenceModel(x)
}

.as_nipt_sex_score <- function(x) {
  if (.is_nipt_sex_score(x)) x else NIPTSexScore(x)
}

.as_nipt_sex_ncv_model <- function(x) {
  if (.is_nipt_sex_ncv_model(x)) x else NIPTSexNCVModel(x)
}

.as_nipt_sex_regression_model <- function(x) {
  if (.is_nipt_sex_regression_model(x)) x else NIPTSexRegressionModel(x)
}

.as_nipt_sex_ncv_score <- function(x) {
  if (.is_nipt_sex_ncv_score(x)) x else NIPTSexNCVScore(x)
}

.as_nipt_sex_regression_score <- function(x) {
  if (.is_nipt_sex_regression_score(x)) x else NIPTSexRegressionScore(x)
}

.as_nipt_gaunosome_score <- function(x) {
  if (.is_nipt_gaunosome_score(x)) x else NIPTGaunosomeScore(x)
}

.as_nipt_gaunosome_report <- function(x) {
  if (.is_nipt_gaunosome_report(x)) x else NIPTGaunosomeReport(x)
}

.as_nipt_control_group_qc <- function(x) {
  if (.is_nipt_control_group_qc(x)) x else NIPTControlGroupQC(x)
}

# ---- Strand-type helper used by guards in scoring functions -----------------

# Returns "separated" or "combined" for S7 NIPTSample/NIPTControlGroup or
# S3 NIPTeRSample/NIPTeRControlGroup objects.
.strand_type_of <- function(x) {
  if (.is_s7_nipt_sample(x)) return(strand_type(x))
  if (.is_s7_separated_control_group(x)) return("separated")
  if (.is_s7_nipt_control_group(x)) return("combined")
  if (inherits(x, "SeparatedStrands")) "separated" else "combined"
}
