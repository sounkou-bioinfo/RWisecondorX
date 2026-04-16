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
# Correction state is tracked as an ordered character vector in
# NIPTCorrectionRecord.  Each correction function appends to the record
# and scoring functions can assert that corrections have been applied.

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


# ---- Correction record -------------------------------------------------------

#' Correction record for a NIPTeR sample
#'
#' Tracks the sequence of corrections applied to the autosomal and sex
#' chromosome matrices of a \code{NIPTSample} object.
#'
#' @export
NIPTCorrectionRecord <- S7::new_class(
  "NIPTCorrectionRecord",
  properties = list(
    autosomal = S7::new_property(S7::class_character,
                                 default = "Uncorrected"),
    sex       = S7::new_property(S7::class_character,
                                 default = "Uncorrected")
  )
)


# ---- Abstract NIPTSample -----------------------------------------------------

#' Abstract base class for NIPTeR bin-count samples
#'
#' Never instantiated directly. Subclasses are \code{CombinedStrandsSample}
#' (reads are summed across strands) and \code{SeparatedStrandsSample} (forward
#' and reverse counts stored independently).
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
#' @return A 22-row numeric matrix.
#' @export
autosomal_matrix <- S7::new_generic("autosomal_matrix", "x")

#' Extract the sex-chromosome read-count matrix from a NIPTSample
#'
#' Returns a 2 x n_bins numeric matrix (row 1 = X, row 2 = Y).
#' For \code{SeparatedStrandsSample}, forward and reverse counts are summed.
#'
#' @param x A \code{NIPTSample} object.
#' @return A 2-row numeric matrix.
#' @export
sex_matrix <- S7::new_generic("sex_matrix", "x")

#' Return the strand type of a NIPTSample
#'
#' @param x A \code{NIPTSample} object.
#' @return \code{"combined"} or \code{"separated"}.
#' @export
strand_type <- S7::new_generic("strand_type", "x")

#' Number of bins in a NIPTSample
#'
#' @param x A \code{NIPTSample} object.
#' @return A positive integer.
#' @export
n_bins <- S7::new_generic("n_bins", "x")

S7::method(n_bins, NIPTSample) <- function(x) ncol(autosomal_matrix(x))


# ---- CombinedStrandsSample ---------------------------------------------------

#' NIPTeR sample with combined forward+reverse read counts
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

S7::method(autosomal_matrix, SeparatedStrandsSample) <-
  function(x) x@auto_fwd + x@auto_rev

S7::method(sex_matrix, SeparatedStrandsSample) <-
  function(x) x@sex_fwd + x@sex_rev

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

.sample_autosomal_reads <- function(x) {
  if (S7::S7_inherits(x, CombinedStrandsSample)) {
    return(list(x@auto_matrix))
  }
  if (S7::S7_inherits(x, SeparatedStrandsSample)) {
    return(list(x@auto_fwd, x@auto_rev))
  }
  .legacy_field(x, "autosomal_chromosome_reads")
}

.sample_sex_reads <- function(x) {
  if (S7::S7_inherits(x, CombinedStrandsSample)) {
    return(list(x@sex_matrix_))
  }
  if (S7::S7_inherits(x, SeparatedStrandsSample)) {
    return(list(x@sex_fwd, x@sex_rev))
  }
  .legacy_field(x, "sex_chromosome_reads")
}

.sample_name <- function(x) {
  if (S7::S7_inherits(x, NIPTSample)) x@sample_name
  else .legacy_field(x, "sample_name")
}

.sample_binsize <- function(x) {
  if (S7::S7_inherits(x, NIPTSample)) x@binsize else NA_integer_
}

.sample_correction_status <- function(x,
                                      scope = c("autosomal", "sex")) {
  scope <- match.arg(scope)
  if (S7::S7_inherits(x, NIPTSample)) {
    if (scope == "autosomal") x@correction@autosomal else x@correction@sex
  } else {
    if (scope == "autosomal") .legacy_field(x, "correction_status_autosomal")
    else .legacy_field(x, "correction_status_sex")
  }
}

.sample_with_reads <- function(x, autosomal = NULL, sex = NULL) {
  if (!S7::S7_inherits(x, NIPTSample)) {
    if (!is.null(autosomal)) {
      x <- .legacy_set_field(x, "autosomal_chromosome_reads", autosomal)
    }
    if (!is.null(sex)) {
      x <- .legacy_set_field(x, "sex_chromosome_reads", sex)
    }
    return(x)
  }

  if (!is.null(autosomal)) {
    if (S7::S7_inherits(x, CombinedStrandsSample)) {
      x@auto_matrix <- autosomal[[1L]]
    } else {
      x@auto_fwd <- autosomal[[1L]]
      x@auto_rev <- autosomal[[2L]]
    }
  }
  if (!is.null(sex)) {
    if (S7::S7_inherits(x, CombinedStrandsSample)) {
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
  if (!S7::S7_inherits(x, NIPTSample)) {
    if (!is.null(autosomal)) {
      x <- .legacy_set_field(x, "correction_status_autosomal", autosomal)
    }
    if (!is.null(sex)) {
      x <- .legacy_set_field(x, "correction_status_sex", sex)
    }
    return(x)
  }

  correction <- x@correction
  if (!is.null(autosomal)) correction@autosomal <- autosomal
  if (!is.null(sex)) correction@sex <- sex
  x@correction <- correction
  x
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
      x <- .sample_with_correction_status(x, autosomal = value)
      x
    },
    correction_status_sex = {
      x <- .sample_with_correction_status(x, sex = value)
      x
    },
    sample_name = {
      x@sample_name <- value
      x
    },
    binsize = {
      x@binsize <- value
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
#' @return A 22 x N numeric matrix of chromosome fractions.
#' @export
fractions_auto <- S7::new_generic("fractions_auto", "cg")

#' Extract fractions suitable for regression
#'
#' For \code{CombinedControlGroup}: 22 x N matrix.
#' For \code{SeparatedControlGroup}: 44 x N matrix (rows "1F".."22F","1R".."22R").
#'
#' @param cg A \code{NIPTControlGroup} object.
#' @return A numeric matrix.
#' @export
fractions_for_regression <- S7::new_generic("fractions_for_regression", "cg")

#' Number of samples in a control group
#' @param cg A \code{NIPTControlGroup} object.
#' @return An integer.
#' @export
n_controls <- S7::new_generic("n_controls", "cg")

#' Sample names of a control group
#' @param cg A \code{NIPTControlGroup} object.
#' @return A character vector.
#' @export
control_names <- S7::new_generic("control_names", "cg")


#' Abstract base class for NIPTeR control groups
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

    ok <- vapply(self@samples, function(s) {
      inherits(s, "NIPTeRSample") || S7::S7_inherits(s, NIPTSample)
    }, logical(1L))
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
    if (S7::S7_inherits(s, NIPTSample)) s@sample_name else .legacy_field(s, "sample_name")
  }, character(1L))


#' NIPTeR control group for CombinedStrands samples
#' @export
CombinedControlGroup <- S7::new_class(
  "CombinedControlGroup",
  parent = NIPTControlGroup
)

#' NIPTeR control group for SeparatedStrands samples
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
    if (S7::S7_inherits(s, SeparatedStrandsSample)) s@auto_fwd
    else s$autosomal_chromosome_reads[[1L]]
  })
  mats_rev <- lapply(cg@samples, function(s) {
    if (S7::S7_inherits(s, SeparatedStrandsSample)) s@auto_rev
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
  if (S7::S7_inherits(cg, NIPTControlGroup)) cg@sample_sex else NULL
}

.control_sex_source <- function(cg) {
  if (S7::S7_inherits(cg, NIPTControlGroup)) cg@sex_source else NULL
}

.control_with_sample_sex <- function(cg, sample_sex = NULL, sex_source = NULL) {
  if (!S7::S7_inherits(cg, NIPTControlGroup)) {
    return(cg)
  }
  cg@sample_sex <- .normalize_sample_sex(sample_sex, cg@samples)
  if (!is.null(sex_source)) cg@sex_source <- sex_source
  cg
}

.nipt_control_group_dollar <- function(x, name) {
  switch(
    name,
    samples = x@samples,
    Description = x@description,
    description = x@description,
    sample_sex = x@sample_sex,
    sex_source = x@sex_source,
    correction_status_autosomal = unique(vapply(
      x@samples, .sample_correction_status, character(1L), scope = "autosomal"
    )),
    correction_status_sex = unique(vapply(
      x@samples, .sample_correction_status, character(1L), scope = "sex"
    )),
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


# ---- S7 <-> S3 conversion helpers -------------------------------------------

# Convert S7 NIPTSample to legacy S3 NIPTeRSample list (for internal function compat)
.s7_to_s3_sample <- function(s) {
  if (S7::S7_inherits(s, CombinedStrandsSample)) {
    structure(
      list(
        autosomal_chromosome_reads  = list(s@auto_matrix),
        sex_chromosome_reads        = list(s@sex_matrix_),
        correction_status_autosomal = s@correction@autosomal,
        correction_status_sex       = s@correction@sex,
        sample_name                 = s@sample_name
      ),
      class = c("NIPTeRSample", "CombinedStrands")
    )
  } else {
    structure(
      list(
        autosomal_chromosome_reads  = list(s@auto_fwd, s@auto_rev),
        sex_chromosome_reads        = list(s@sex_fwd, s@sex_rev),
        correction_status_autosomal = s@correction@autosomal,
        correction_status_sex       = s@correction@sex,
        sample_name                 = s@sample_name
      ),
      class = c("NIPTeRSample", "SeparatedStrands")
    )
  }
}

# Convert legacy S3 NIPTeRSample list back to S7 NIPTSample
# `original_s7` is the original S7 object (used to know binsize)
.s3_to_s7_sample <- function(s3, original_s7) {
  .d <- function(m) m
  if (inherits(s3, "SeparatedStrands")) {
    SeparatedStrandsSample(
      sample_name = s3$sample_name,
      binsize     = original_s7@binsize,
      auto_fwd    = .d(s3$autosomal_chromosome_reads[[1L]]),
      auto_rev    = .d(s3$autosomal_chromosome_reads[[2L]]),
      sex_fwd     = .d(s3$sex_chromosome_reads[[1L]]),
      sex_rev     = .d(s3$sex_chromosome_reads[[2L]]),
      correction  = NIPTCorrectionRecord(
        autosomal = s3$correction_status_autosomal,
        sex       = s3$correction_status_sex
      )
    )
  } else {
    CombinedStrandsSample(
      sample_name = s3$sample_name,
      binsize     = original_s7@binsize,
      auto_matrix = .d(s3$autosomal_chromosome_reads[[1L]]),
      sex_matrix_ = .d(s3$sex_chromosome_reads[[1L]]),
      correction  = NIPTCorrectionRecord(
        autosomal = s3$correction_status_autosomal,
        sex       = s3$correction_status_sex
      )
    )
  }
}

# Convert S7 NIPTControlGroup to legacy S3 NIPTeRControlGroup list
.s7_cg_to_s3 <- function(cg) {
  s3_samples <- lapply(cg@samples, function(s) {
    if (S7::S7_inherits(s, NIPTSample)) .s7_to_s3_sample(s) else s
  })
  structure(
    list(
      samples = s3_samples,
      correction_status_autosomal_chromosomes = "Uncorrected",
      correction_status_sex_chromosomes = "Uncorrected",
      Description = cg@description
    ),
    class = c("NIPTeRControlGroup",
              if (S7::S7_inherits(cg, SeparatedControlGroup)) "SeparatedStrands" else "CombinedStrands")
  )
}

# Convert legacy S3 NIPTeRControlGroup back to S7 NIPTControlGroup
.s3_cg_to_s7 <- function(s3_cg, original_s7_cg) {
  s7_samples <- lapply(seq_along(s3_cg$samples), function(i) {
    orig <- if (i <= length(original_s7_cg@samples)) original_s7_cg@samples[[i]] else NULL
    s3s <- s3_cg$samples[[i]]
    if (!is.null(orig) && S7::S7_inherits(orig, NIPTSample)) {
      .s3_to_s7_sample(s3s, orig)
    } else {
      s3s
    }
  })
  if (S7::S7_inherits(original_s7_cg, SeparatedControlGroup)) {
    SeparatedControlGroup(samples = s7_samples, description = original_s7_cg@description)
  } else {
    CombinedControlGroup(samples = s7_samples, description = original_s7_cg@description)
  }
}


# ---- Strand-type helper used by guards in scoring functions -----------------

# Returns "separated" or "combined" for S7 NIPTSample/NIPTControlGroup or
# S3 NIPTeRSample/NIPTeRControlGroup objects.
.strand_type_of <- function(x) {
  if (S7::S7_inherits(x, NIPTSample)) return(strand_type(x))
  if (S7::S7_inherits(x, SeparatedControlGroup)) return("separated")
  if (S7::S7_inherits(x, NIPTControlGroup)) return("combined")
  if (inherits(x, "SeparatedStrands")) "separated" else "combined"
}
