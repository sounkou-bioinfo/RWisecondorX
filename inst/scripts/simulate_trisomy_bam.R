#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package 'optparse' is required for this script.", call. = FALSE)
}
if (!requireNamespace("RWisecondorX", quietly = TRUE)) {
  stop("Package 'RWisecondorX' must be installed to use this script.",
       call. = FALSE)
}

option_list <- list(
  optparse::make_option("--bam", type = "character", dest = "bam",
                        help = "Input donor BAM/CRAM"),
  optparse::make_option("--out", type = "character", dest = "out",
                        help = "Output BAM path"),
  optparse::make_option("--chrom", type = "character", dest = "chrom",
                        help = "Target chromosome, e.g. 21 or chr21"),
  optparse::make_option("--fetal-fraction", type = "double",
                        dest = "fetal_fraction", default = 0.10,
                        help = "Fetal fraction in (0,1]; default %default"),
  optparse::make_option("--mosaic", action = "store_true",
                        dest = "mosaic", default = FALSE,
                        help = "Use the 50%% mosaic model"),
  optparse::make_option("--seed", type = "integer", dest = "seed",
                        default = 1L,
                        help = "RNG seed; default %default"),
  optparse::make_option("--threads", type = "integer", dest = "threads",
                        default = 1L,
                        help = "htslib I/O and indexing threads; default %default"),
  optparse::make_option("--reference", type = "character",
                        dest = "reference", default = NULL,
                        help = "Reference FASTA required for CRAM input"),
  optparse::make_option("--no-index", action = "store_true",
                        dest = "no_index", default = FALSE,
                        help = "Do not create a .bai index"),
  optparse::make_option("--overwrite", action = "store_true",
                        dest = "overwrite", default = FALSE,
                        help = "Overwrite existing output BAM and BAI")
)

parser <- optparse::OptionParser(
  usage = "%prog --bam donor.bam --out sim_t21.bam --chrom 21 [options]",
  option_list = option_list
)
opt <- optparse::parse_args(parser)

required <- c("bam", "out", "chrom")
missing_required <- required[!nzchar(vapply(required, function(x) {
  val <- opt[[x]]
  if (is.null(val)) "" else as.character(val)
}, character(1)))]
if (length(missing_required) > 0L) {
  optparse::print_help(parser)
  stop("Missing required options: ", paste(missing_required, collapse = ", "),
       call. = FALSE)
}

res <- RWisecondorX::simulate_trisomy_bam(
  bam = opt$bam,
  out_bam = opt$out,
  trisomy_chr = opt$chrom,
  fetal_fraction = opt$fetal_fraction,
  mosaic = isTRUE(opt$mosaic),
  seed = opt$seed,
  threads = opt$threads,
  reference = opt$reference,
  index = !isTRUE(opt$no_index),
  overwrite = isTRUE(opt$overwrite)
)

cat(sprintf("Wrote %s\n", res$output_bam))
cat(sprintf("Target chromosome: %s\n", res$target_chr))
cat(sprintf("Input records: %.0f\n", res$input_records))
cat(sprintf("Output records: %.0f\n", res$output_records))
cat(sprintf("Removal probability: %.6f\n", res$removal_prob))
