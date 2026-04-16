#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package 'optparse' is required for this script.", call. = FALSE)
}
if (!requireNamespace("RWisecondorX", quietly = TRUE)) {
  stop("Package 'RWisecondorX' must be installed to use this script.",
       call. = FALSE)
}

.split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    character()
  } else {
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  }
}

.parse_logical_csv <- function(x) {
  vals <- tolower(.split_csv(x))
  if (length(vals) == 0L) {
    return(logical())
  }
  if (any(!vals %in% c("true", "false", "1", "0"))) {
    stop("`--mosaic-values` must contain only TRUE/FALSE or 1/0 values.",
         call. = FALSE)
  }
  vals %in% c("true", "1")
}

option_list <- list(
  optparse::make_option("--bam-dir", type = "character", dest = "bam_dir",
                        default = NULL,
                        help = "Directory of donor BAM/CRAM files"),
  optparse::make_option("--bams", type = "character", dest = "bams",
                        default = NULL,
                        help = "Comma-separated donor BAM/CRAM paths"),
  optparse::make_option("--out-dir", type = "character", dest = "out_dir",
                        help = "Output directory for simulated BAMs"),
  optparse::make_option("--donor-ids", type = "character",
                        dest = "donor_ids", default = NULL,
                        help = "Optional comma-separated donor IDs"),
  optparse::make_option("--chroms", type = "character", dest = "chroms",
                        default = "21,18,13",
                        help = "Comma-separated trisomy chromosomes; default %default"),
  optparse::make_option("--fetal-fractions", type = "character",
                        dest = "fetal_fractions", default = "0.10",
                        help = "Comma-separated fetal fractions; default %default"),
  optparse::make_option("--mosaic-values", type = "character",
                        dest = "mosaic_values", default = "FALSE",
                        help = paste("Comma-separated mosaic flags",
                                     "(TRUE/FALSE or 1/0); default %default")),
  optparse::make_option("--include-donors", action = "store_true",
                        dest = "include_donors", default = FALSE,
                        help = "Include donor negatives in the manifest"),
  optparse::make_option("--seed", type = "integer", dest = "seed",
                        default = 1L,
                        help = "Base RNG seed; default %default"),
  optparse::make_option("--threads", type = "integer", dest = "threads",
                        default = 1L,
                        help = "htslib I/O and indexing threads; default %default"),
  optparse::make_option("--reference", type = "character",
                        dest = "reference", default = NULL,
                        help = "Reference FASTA required for CRAM input"),
  optparse::make_option("--pattern", type = "character", dest = "pattern",
                        default = "\\.(bam|cram)$",
                        help = "Regex used with --bam-dir; default %default"),
  optparse::make_option("--no-index", action = "store_true",
                        dest = "no_index", default = FALSE,
                        help = "Do not create .bai indexes for simulated BAMs"),
  optparse::make_option("--overwrite", action = "store_true",
                        dest = "overwrite", default = FALSE,
                        help = "Overwrite existing simulated BAMs and manifest")
)

parser <- optparse::OptionParser(
  usage = paste(
    "%prog --bam-dir donors/ --out-dir simulated/ [options]",
    "or %prog --bams a.bam,b.bam --out-dir simulated/ [options]"
  ),
  option_list = option_list
)
opt <- optparse::parse_args(parser)

if (is.null(opt$bam_dir) == is.null(opt$bams)) {
  optparse::print_help(parser)
  stop("Supply exactly one of --bam-dir or --bams.", call. = FALSE)
}
if (is.null(opt$out_dir) || !nzchar(opt$out_dir)) {
  optparse::print_help(parser)
  stop("Missing required option: --out-dir", call. = FALSE)
}

bams <- if (is.null(opt$bams)) NULL else .split_csv(opt$bams)
donor_ids <- .split_csv(opt$donor_ids)
if (length(donor_ids) == 0L) {
  donor_ids <- NULL
}

chroms <- .split_csv(opt$chroms)
fetal_fractions <- as.numeric(.split_csv(opt$fetal_fractions))
mosaic <- .parse_logical_csv(opt$mosaic_values)

manifest <- RWisecondorX::simulate_trisomy_cohort(
  bams = bams,
  bam_dir = opt$bam_dir,
  out_dir = opt$out_dir,
  donor_ids = donor_ids,
  trisomy_chromosomes = chroms,
  fetal_fraction = fetal_fractions,
  mosaic = mosaic,
  seed = opt$seed,
  threads = opt$threads,
  reference = opt$reference,
  index = !isTRUE(opt$no_index),
  include_donors = isTRUE(opt$include_donors),
  pattern = opt$pattern,
  overwrite = isTRUE(opt$overwrite)
)

cat(sprintf("Wrote manifest: %s\n", attr(manifest, "manifest_path")))
cat(sprintf("Rows: %d\n", nrow(manifest)))
cat(sprintf("Simulated positives: %d\n",
            sum(manifest$cohort_role == "simulated_positive")))
