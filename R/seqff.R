#' Predict Fetal Fraction with SeqFF
#'
#' @param input BAM path, SAM path, counts path, or a data frame with `binName`
#'   and `counts` columns.
#' @param input_type One of `bam`, `sam`, or `counts`.
#' @param samtools_bin Samtools executable to use for BAM input.
#' @param samtools_exclude_flags Optional integer flag passed to `samtools view -F`.
#' @param samtools_min_mapq Optional integer mapping-quality threshold passed to
#'   `samtools view -q`.
#'
#' @return Named numeric vector with `SeqFF`, `Enet`, and `WRSC`.
#' @export
seqff_predict <- function(input,
                          input_type = c("bam", "sam", "counts"),
                          samtools_bin = "samtools",
                          samtools_exclude_flags = NULL,
                          samtools_min_mapq = NULL) {
  input_type <- match.arg(input_type)
  support <- seqff_load_support_data()
  counts <- seqff_prepare_counts(
    input = input,
    input_type = input_type,
    samtools_bin = samtools_bin,
    samtools_exclude_flags = samtools_exclude_flags,
    samtools_min_mapq = samtools_min_mapq
  )

  seqff_predict_from_counts(
    counts = counts,
    bininfo = support$bininfo,
    model = support$model
  )
}

seqff_prepare_counts <- function(input,
                                 input_type,
                                 samtools_bin,
                                 samtools_exclude_flags,
                                 samtools_min_mapq) {
  if (input_type == "counts") {
    if (is.data.frame(input)) {
      counts <- input
    } else {
      counts <- utils::read.table(
        input,
        header = FALSE,
        colClasses = c("character", "integer")
      )
      colnames(counts) <- c("binName", "counts")
    }
    return(counts)
  }

  positions <- seqff_read_positions(
    input = input,
    input_type = input_type,
    samtools_bin = samtools_bin,
    samtools_exclude_flags = samtools_exclude_flags,
    samtools_min_mapq = samtools_min_mapq
  )

  seqff_bin_counts(positions)
}

seqff_read_positions <- function(input,
                                 input_type,
                                 samtools_bin,
                                 samtools_exclude_flags,
                                 samtools_min_mapq) {
  command <- switch(
    input_type,
    bam = seqff_samtools_view_command(
      bam = input,
      samtools_bin = samtools_bin,
      samtools_exclude_flags = samtools_exclude_flags,
      samtools_min_mapq = samtools_min_mapq
    ),
    sam = sprintf("cut -f3,4 %s", shQuote(normalizePath(input, winslash = "/", mustWork = TRUE)))
  )

  positions <- data.table::fread(
    cmd = command,
    sep = "\t",
    header = FALSE,
    select = 1:2,
    col.names = c("refChr", "begin"),
    colClasses = c("character", "integer"),
    data.table = FALSE
  )

  if (nrow(positions) == 0L) {
    stop("No alignments were read for SeqFF input.", call. = FALSE)
  }

  positions$refChr <- seqff_normalize_chr_names(positions$refChr)
  mitochondrial <- positions$refChr %in% c("*", "chrM", "chrMT", "M", "MT")
  positions[!mitochondrial, , drop = FALSE]
}

seqff_samtools_view_command <- function(bam,
                                        samtools_bin,
                                        samtools_exclude_flags,
                                        samtools_min_mapq) {
  args <- c("view")

  if (!is.null(samtools_exclude_flags)) {
    args <- c(args, "-F", as.character(samtools_exclude_flags))
  }
  if (!is.null(samtools_min_mapq)) {
    args <- c(args, "-q", as.character(samtools_min_mapq))
  }

  args <- c(args, shQuote(normalizePath(bam, winslash = "/", mustWork = TRUE)))
  paste(shQuote(samtools_bin), paste(args, collapse = " "), "| cut -f3,4")
}

seqff_normalize_chr_names <- function(ref_chr) {
  ref_chr <- as.character(ref_chr)
  has_chr_prefix <- grepl("^chr", ref_chr)
  ref_chr[!has_chr_prefix] <- paste0("chr", ref_chr[!has_chr_prefix])
  ref_chr
}

seqff_bin_counts <- function(positions) {
  binindex <- ifelse(
    (positions$begin %% 50000) == 0,
    floor(positions$begin / 50000) - 1L,
    floor(positions$begin / 50000)
  )
  bincounts <- table(paste(positions$refChr, binindex, sep = "_"))

  data.frame(
    binName = names(bincounts),
    counts = as.numeric(bincounts),
    stringsAsFactors = FALSE
  )
}

seqff_predict_from_counts <- function(counts, bininfo, model) {
  merged <- merge(bininfo, counts, by = "binName", all.x = TRUE)
  merged <- merged[order(merged$binorder), ]

  autosomes <- merged$BinFilterFlag == 1 &
    !merged$CHR %in% c("chrX", "chrY", "X", "Y")
  usable <- merged$BinFilterFlag == 1

  autoscaled <- merged$counts[autosomes] / sum(merged$counts[autosomes], na.rm = TRUE)
  allscaled <- merged$counts[usable] / sum(merged$counts[autosomes], na.rm = TRUE)
  mediancountpergc <- tapply(
    autoscaled,
    merged$GC[autosomes],
    median,
    na.rm = TRUE
  )

  loess_fitted <- stats::predict(
    stats::loess(mediancountpergc ~ as.numeric(names(mediancountpergc))),
    merged$GC[usable]
  )
  normalized <- allscaled + (stats::median(autoscaled, na.rm = TRUE) - loess_fitted)

  bincounts <- rep(0, nrow(bininfo))
  names(bincounts) <- merged$binName
  bincounts[usable] <- (normalized / sum(normalized, na.rm = TRUE)) * length(normalized)

  wrsc <- seqff_ff_pred(
    gc_norm_bc_61927 = bincounts,
    B = model$B,
    mu = model$mu,
    parameter_1 = model$parameter.1,
    parameter_2 = model$parameter.2
  )
  enet <- as.numeric(bincounts %*% model$elnetbeta + model$elnetintercept)
  ff <- c((wrsc + enet) / 2, enet, wrsc)
  names(ff) <- c("SeqFF", "Enet", "WRSC")
  ff
}

seqff_ff_pred <- function(gc_norm_bc_61927, B, mu, parameter_1, parameter_2) {
  gc_norm_bc_61927[is.na(gc_norm_bc_61927)] <- 0
  gc_norm_bc <- gc_norm_bc_61927[grepl("chr[0-9]", names(gc_norm_bc_61927))]
  gc_norm_bc_centered <- gc_norm_bc - mu
  y_hat <- matrix(c(1, gc_norm_bc_centered), nrow = 1) %*% B
  y_hat_rep <- sum(y_hat, na.rm = TRUE) / sum(gc_norm_bc)
  (y_hat_rep + parameter_1) * parameter_2
}

seqff_load_support_data <- function() {
  bininfo_path <- seqff_support_path("pd4615-sup-0010-table2.csv")
  model_path <- seqff_support_path("pd4615-sup-0008-file1.rdata")

  bininfo <- utils::read.csv(bininfo_path, stringsAsFactors = FALSE)
  colnames(bininfo)[1] <- "binName"
  bininfo$binorder <- seq_len(nrow(bininfo))

  load_env <- new.env(parent = emptyenv())
  load(model_path, envir = load_env)

  list(
    bininfo = bininfo,
    model = list(
      B = load_env$B,
      mu = load_env$mu,
      parameter.1 = load_env$parameter.1,
      parameter.2 = load_env$parameter.2,
      elnetbeta = load_env$elnetbeta,
      elnetintercept = load_env$elnetintercept
    )
  )
}

seqff_support_path <- function(filename) {
  installed_path <- system.file("extdata", "seqff", filename, package = "RWisecondorX")
  if (nzchar(installed_path)) {
    return(installed_path)
  }

  root <- .seqff_find_package_root()
  path <- file.path(root, "inst", "extdata", "seqff", filename)
  if (!file.exists(path)) {
    stop("SeqFF support file not found: ", filename, call. = FALSE)
  }
  path
}

.seqff_find_package_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate package root for SeqFF support files.", call. = FALSE)
    }
    current <- parent
  }
}
