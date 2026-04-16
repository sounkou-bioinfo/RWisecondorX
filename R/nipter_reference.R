#' Build a typed NIPT reference model
#'
#' Packages the control group, chromosome-level training frame, and optional
#' sex-prediction models into one validated reference object. The embedded
#' \code{reference_frame} is enriched with the derived sex-scoring columns used
#' by the production pipeline: \code{ConsensusGender}, \code{RR_X},
#' \code{RR_Y}, \code{RR_X_SexClassMAD}, \code{RR_Y_SexClassMAD}, and
#' \code{IsRefSexOutlier}. When \code{y_unique_ratios} are supplied, the frame
#' also carries \code{YUniqueRatio}.
#'
#' @param control_group A \code{NIPTControlGroup}.
#' @param sample_sex Optional character vector overriding the sex labels stored
#'   on \code{control_group}. Accepted values are \code{"female"},
#'   \code{"male"}, \code{"ambiguous"}, and \code{"unknown"}.
#' @param sex_source Optional scalar string describing where \code{sample_sex}
#'   came from when supplied here, for example \code{"explicit"} or
#'   \code{"laboratory_lims"}.
#' @param sex_methods Character vector of fraction-based sex model methods to
#'   build. Allowed values are \code{"y_fraction"} and \code{"xy_fraction"}.
#'   Default builds both.
#' @param y_unique_ratios Optional named numeric vector of Y-unique ratios for
#'   building an additional \code{"y_unique"} sex model and storing
#'   \code{YUniqueRatio} in the reference frame.
#' @param build_params Optional named list of provenance parameters to attach to
#'   the returned reference model.
#'
#' @return A list-like \code{NIPTReferenceModel} S7 object containing:
#' \describe{
#'   \item{control_group}{The input control group, optionally re-annotated with
#'     explicit sample sex labels.}
#'   \item{reference_frame}{A typed \code{NIPTReferenceFrame} with one row per
#'     control sample, including chromosome counts/fractions and the derived
#'     sex-aware reference columns used for X/Y scoring.}
#'   \item{sex_models}{Named list of \code{NIPTeRSexModel} objects keyed by
#'     method.}
#'   \item{sample_sex_source}{The provenance string for the control-group sex
#'     labels, if available.}
#'   \item{build_date}{UTC timestamp for the build.}
#'   \item{build_params}{Caller-supplied provenance metadata.}
#' }
#'
#' @export
nipter_build_reference <- function(control_group,
                                   sample_sex = NULL,
                                   sex_source = NULL,
                                   sex_methods = c("y_fraction", "xy_fraction"),
                                   y_unique_ratios = NULL,
                                   build_params = list()) {
  stopifnot(.is_nipt_control_group_object(control_group))
  stopifnot(is.list(build_params))

  sex_methods <- unique(as.character(sex_methods))
  allowed_methods <- c("y_fraction", "xy_fraction")
  if (!all(sex_methods %in% allowed_methods)) {
    stop(
      "'sex_methods' must contain only: ",
      paste(allowed_methods, collapse = ", "),
      call. = FALSE
    )
  }

  cg <- control_group
  if (!is.null(sample_sex) || !is.null(sex_source)) {
    cg <- .control_with_sample_sex(
      cg,
      sample_sex = if (is.null(sample_sex)) .control_sample_sex(cg) else sample_sex,
      sex_source = if (is.null(sex_source)) .control_sex_source(cg) else sex_source
    )
  }

  y_unique_ratios <- .normalize_named_numeric_by_sample(
    y_unique_ratios,
    cg$samples,
    arg = "y_unique_ratios"
  )

  reference_frame <- nipter_reference_frame(cg)

  sex_models <- list()
  for (method in sex_methods) {
    sex_models[[method]] <- nipter_sex_model(cg, method = method)
  }
  if (!is.null(y_unique_ratios)) {
    sex_models[["y_unique"]] <- nipter_sex_model_y_unique(y_unique_ratios)
  }

  reference_frame <- .augment_reference_frame_for_sex(
    reference_frame,
    cg,
    sex_models = sex_models,
    y_unique_ratios = y_unique_ratios
  )

  .as_nipt_reference_model(list(
    control_group = cg,
    reference_frame = reference_frame,
    sex_models = sex_models,
    sample_sex_source = .control_sex_source(cg),
    build_date = format(Sys.time(), tz = "UTC", usetz = TRUE),
    build_params = build_params
  ))
}

.normalize_named_numeric_by_sample <- function(values, samples, arg) {
  if (is.null(values)) {
    return(NULL)
  }
  if (!is.numeric(values) || !length(values)) {
    stop(sprintf("'%s' must be NULL or a non-empty numeric vector.", arg),
         call. = FALSE)
  }

  sample_names <- vapply(samples, .sample_name, character(1L))
  if (is.null(names(values))) {
    if (length(values) != length(sample_names)) {
      stop(
        sprintf(
          "Unnamed '%s' vectors must have length equal to the number of samples.",
          arg
        ),
        call. = FALSE
      )
    }
    names(values) <- sample_names
  } else {
    if (!all(sample_names %in% names(values))) {
      stop(
        sprintf("'%s' names must cover every control-group sample name.", arg),
        call. = FALSE
      )
    }
    values <- values[sample_names]
  }

  names(values) <- sample_names
  values
}

.zscore_minus_self <- function(x) {
  out <- rep(NA_real_, length(x))
  if (!length(x)) {
    return(out)
  }

  for (i in seq_along(x)) {
    ref <- x[-i]
    ref <- ref[is.finite(ref)]
    if (length(ref) < 2L) {
      next
    }
    ref_sd <- stats::sd(ref)
    if (!is.finite(ref_sd) || ref_sd <= .Machine$double.eps) {
      next
    }
    out[i] <- (x[[i]] - mean(ref)) / ref_sd
  }

  out
}

.predict_reference_sample_sex <- function(control_group,
                                          sex_models,
                                          y_unique_ratios = NULL) {
  if (!length(sex_models)) {
    return(NULL)
  }

  samples <- control_group$samples
  out <- stats::setNames(character(length(samples)),
                         vapply(samples, .sample_name, character(1L)))

  for (i in seq_along(samples)) {
    samp <- samples[[i]]
    samp_name <- .sample_name(samp)
    yu_ratio <- if (is.null(y_unique_ratios)) NULL else y_unique_ratios[[samp_name]]
    out[[samp_name]] <- nipter_predict_sex(
      samp,
      sex_models,
      y_unique_ratio = yu_ratio
    )$prediction
  }

  out
}

.resolve_reference_consensus_gender <- function(sample_sex, predicted_gender) {
  out <- predicted_gender

  if (!is.null(sample_sex)) {
    explicit <- ifelse(sample_sex %in% c("female", "male"),
                       sample_sex, NA_character_)
    if (is.null(out)) {
      out <- explicit
    } else {
      use_explicit <- !is.na(explicit)
      out[use_explicit] <- explicit[use_explicit]
    }
  }

  if (is.null(out) || anyNA(out) || !all(out %in% c("female", "male"))) {
    return(NULL)
  }

  out
}

.augment_reference_frame_for_sex <- function(reference_frame,
                                             control_group,
                                             sex_models,
                                             y_unique_ratios = NULL,
                                             outlier_threshold = 3) {
  frame <- as.data.frame(reference_frame, stringsAsFactors = FALSE)
  frame$Sample_name <- as.character(frame$Sample_name)
  auto_totals <- rowSums(frame[, .nipt_reference_count_cols[1:22], drop = FALSE])
  denom <- pmax(auto_totals, .Machine$double.eps)
  frame$RR_X <- frame$NChrReads_X / denom
  frame$RR_Y <- frame$NChrReads_Y / denom

  if (!is.null(y_unique_ratios)) {
    frame$YUniqueRatio <- unname(y_unique_ratios[frame$Sample_name])
  }

  predicted_gender <- .predict_reference_sample_sex(
    control_group,
    sex_models = sex_models,
    y_unique_ratios = y_unique_ratios
  )
  consensus_gender <- .resolve_reference_consensus_gender(
    sample_sex = .control_sample_sex(control_group),
    predicted_gender = predicted_gender
  )

  if (!is.null(consensus_gender)) {
    frame$ConsensusGender <- unname(consensus_gender[frame$Sample_name])
    frame$RR_X_SexClassMAD <- NA_real_
    frame$RR_Y_SexClassMAD <- NA_real_
    for (sex in c("female", "male")) {
      idx <- which(frame$ConsensusGender == sex)
      if (!length(idx)) {
        next
      }
      frame$RR_X_SexClassMAD[idx] <- .zscore_minus_self(frame$RR_X[idx])
      frame$RR_Y_SexClassMAD[idx] <- .zscore_minus_self(frame$RR_Y[idx])
    }
    frame$IsRefSexOutlier <- (
      is.finite(frame$RR_X_SexClassMAD) &
        abs(frame$RR_X_SexClassMAD) >= outlier_threshold
    ) | (
      is.finite(frame$RR_Y_SexClassMAD) &
        abs(frame$RR_Y_SexClassMAD) >= outlier_threshold
    )
  }

  .as_nipt_reference_frame(frame)
}
