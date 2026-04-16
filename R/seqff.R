#' Predict Fetal Fraction with SeqFF
#'
#' @param input BAM path, SAM path, counts path, or a data frame with `binName`
#'   and `counts` columns.
#' @param input_type One of `bam`, `sam`, or `counts`.
#' @param mapq Optional integer mapping-quality threshold for BAM input.
#' @param require_flags Optional integer bitmask; only reads with all bits set
#'   are kept for BAM input.
#' @param exclude_flags Optional integer bitmask; reads with any bit set are
#'   dropped for BAM input.
#' @param con Optional open DBI connection with duckhts already loaded.
#' @param reference Optional FASTA reference path for CRAM input.
#'
#' @return Named numeric vector with `SeqFF`, `Enet`, and `WRSC`.
#' @export
seqff_predict <- function(input,
                          input_type = c("bam", "sam", "counts"),
                          mapq = NULL,
                          require_flags = NULL,
                          exclude_flags = NULL,
                          con = NULL,
                          reference = NULL) {
  input_type <- match.arg(input_type)
  support <- seqff_load_support_data()
  counts <- seqff_prepare_counts(
    input = input,
    input_type = input_type,
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    con = con,
    reference = reference
  )

  seqff_predict_from_counts(
    counts = counts,
    bininfo = support$bininfo,
    model = support$model
  )
}

seqff_prepare_counts <- function(input,
                                 input_type,
                                 mapq,
                                 require_flags,
                                 exclude_flags,
                                 con,
                                 reference) {
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
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    con = con,
    reference = reference
  )

  seqff_bin_counts(positions)
}

seqff_read_positions <- function(input,
                                 input_type,
                                 mapq,
                                 require_flags,
                                 exclude_flags,
                                 con,
                                 reference) {
  positions <- switch(
    input_type,
    bam = seqff_read_positions_bam(
      bam = input,
      mapq = mapq,
      require_flags = require_flags,
      exclude_flags = exclude_flags,
      con = con,
      reference = reference
    ),
    sam = seqff_read_positions_sam(input)
  )

  if (nrow(positions) == 0L) {
    stop("No alignments were read for SeqFF input.", call. = FALSE)
  }

  positions$refChr <- seqff_normalize_chr_names(positions$refChr)
  mitochondrial <- positions$refChr %in% c("*", "chrM", "chrMT", "M", "MT")
  positions[!mitochondrial, , drop = FALSE]
}

seqff_read_positions_bam <- function(bam,
                                     mapq,
                                     require_flags,
                                     exclude_flags,
                                     con,
                                     reference) {
  stopifnot(is.character(bam), length(bam) == 1L, nzchar(bam))
  stopifnot(file.exists(bam))
  if (!is.null(reference)) {
    stopifnot(is.character(reference), length(reference) == 1L, nzchar(reference))
    stopifnot(file.exists(reference))
  }

  mapq <- if (is.null(mapq)) 0L else as.integer(mapq)
  require_flags <- if (is.null(require_flags)) 0L else as.integer(require_flags)
  exclude_flags <- if (is.null(exclude_flags)) 0L else as.integer(exclude_flags)

  own_con <- is.null(con)
  if (own_con) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
    Rduckhts::rduckhts_load(con)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  rows <- Rduckhts::rduckhts_bam_bin_counts(
    con,
    path = bam,
    bin_width = 50000L,
    reference = reference,
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    rmdup = "none"
  )

  if (!nrow(rows)) {
    return(data.frame(refChr = character(0), begin = integer(0), stringsAsFactors = FALSE))
  }

  chrom_col <- grep("^chrom$", names(rows), ignore.case = TRUE, value = TRUE)
  if (!length(chrom_col)) {
    stop("SeqFF BAM bin counts missing chrom column.", call. = FALSE)
  }

  data.frame(
    refChr = rows[[chrom_col[[1L]]]],
    begin = as.integer(rows$bin_id) * 50000L + 1L,
    stringsAsFactors = FALSE
  )
}

seqff_read_positions_sam <- function(input) {
  data.table::fread(
    cmd = sprintf("cut -f3,4 %s", shQuote(normalizePath(input, winslash = "/", mustWork = TRUE))),
    sep = "\t",
    header = FALSE,
    select = 1:2,
    col.names = c("refChr", "begin"),
    colClasses = c("character", "integer"),
    data.table = FALSE
  )
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
  stop("SeqFF support file not found in installed package: ", filename, call. = FALSE)
}
