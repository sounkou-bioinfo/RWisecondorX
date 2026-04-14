#!/usr/bin/env Rscript
# inst/scripts/convert_sample.R
#
# CLI script for converting a single BAM/CRAM into a binned BED.gz or NPZ file.
#
# Supports both WisecondorX and NIPTeR binning modes.
#
# Usage:
#   Rscript convert_sample.R --mode wisecondorx --bam sample.bam --out sample.bed.gz
#   Rscript convert_sample.R --mode wisecondorx --bam sample.bam --out sample.npz --npz
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz --gc-table hg38_gc.tsv.bgz
#   Rscript convert_sample.R --mode nipter --bam sample.bam --out sample.bed.gz --fasta hg38.fa

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required. Install it with: install.packages('optparse')",
       call. = FALSE)
}

library(optparse)

option_list <- list(
  make_option("--mode", type = "character", default = "wisecondorx",
              help = "Binning mode: 'wisecondorx' or 'nipter' [default: %default]"),

  # --- Input / output ---
  make_option("--bam", type = "character", default = NULL,
              help = "Input BAM or CRAM file [required]"),
  make_option("--out", type = "character", default = NULL,
              help = "Output path (.bed.gz or .npz) [required]"),

  # --- Output format ---
  make_option("--npz", action = "store_true", default = FALSE,
              help = paste0("Write WisecondorX-compatible NPZ instead of BED.gz. ",
                            "Only valid in wisecondorx mode [default: %default]")),

  # --- Binning parameters ---
  make_option("--binsize", type = "integer", default = NULL,
              help = paste0("Bin size in bp. Defaults: 5000 (wisecondorx), ",
                            "50000 (nipter) [default: mode-dependent]")),
  make_option("--mapq", type = "integer", default = NULL,
              help = "Minimum MAPQ [default: 1 wisecondorx, 0 nipter]"),
  make_option("--rmdup", type = "character", default = NULL,
              help = paste0("Dedup mode: 'streaming', 'flag', 'none'. ",
                            "Defaults: 'streaming' (wisecondorx), 'none' (nipter)")),
  make_option("--exclude-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask to exclude (samtools -F) [default: %default]"),
  make_option("--require-flags", type = "integer", default = 0L,
              help = "SAM flag bitmask to require (samtools -f) [default: %default]"),
  make_option("--reference", type = "character", default = NULL,
              help = "FASTA reference for CRAM decoding"),

  # --- NIPTeR-specific ---
  make_option("--separate-strands", action = "store_true", default = FALSE,
              help = "NIPTeR mode: produce SeparatedStrands 9-column BED [default: %default]"),
  make_option("--gc-table", type = "character", default = NULL,
              help = paste0("Precomputed GC table (TSV.bgz from nipter_gc_precompute). ",
                            "NIPTeR mode: GC-correct and embed corrected counts in BED")),
  make_option("--fasta", type = "character", default = NULL,
              help = paste0("FASTA reference for on-the-fly GC correction (alternative to --gc-table). ",
                            "NIPTeR mode only. Slower than --gc-table for multiple samples"))
)

parser <- OptionParser(
  usage = "usage: %prog [options]",
  description = paste(
    "Convert a single BAM/CRAM file to a binned BED.gz or NPZ file.",
    "",
    "Supports WisecondorX and NIPTeR binning modes with all filtering options.",
    "",
    "Examples:",
    "  # WisecondorX BED output (default)",
    "  %prog --bam sample.bam --out sample.bed.gz",
    "",
    "  # WisecondorX NPZ output (for Python interop)",
    "  %prog --bam sample.bam --out sample.npz --npz",
    "",
    "  # NIPTeR BED output with GC correction (precomputed table)",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz --gc-table hg38_gc.tsv.bgz",
    "",
    "  # NIPTeR BED output with GC correction (on-the-fly from FASTA)",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz --fasta hg38.fa",
    "",
    "  # NIPTeR with pre-filtering (MAPQ 40, exclude duplicates)",
    "  %prog --mode nipter --bam sample.bam --out sample.bed.gz --mapq 40 --exclude-flags 1024",
    sep = "\n"
  ),
  option_list = option_list
)

opts <- parse_args(parser)

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

mode <- tolower(opts$mode)
if (!mode %in% c("wisecondorx", "nipter")) {
  stop("--mode must be 'wisecondorx' or 'nipter', got: ", opts$mode, call. = FALSE)
}

if (is.null(opts$bam)) {
  stop("--bam is required (path to input BAM/CRAM file)", call. = FALSE)
}

if (!file.exists(opts$bam)) {
  stop("BAM file not found: ", opts$bam, call. = FALSE)
}

if (is.null(opts$out)) {
  stop("--out is required (output path for BED.gz or NPZ file)", call. = FALSE)
}

if (opts$npz && mode != "wisecondorx") {
  stop("--npz is only valid in wisecondorx mode", call. = FALSE)
}

if (!is.null(opts$`gc-table`) && mode != "nipter") {
  stop("--gc-table is only valid in nipter mode", call. = FALSE)
}

if (!is.null(opts$fasta) && mode != "nipter") {
  stop("--fasta GC correction is only valid in nipter mode", call. = FALSE)
}

if (!is.null(opts$`gc-table`) && !is.null(opts$fasta)) {
  stop("--gc-table and --fasta are mutually exclusive; provide one or the other", call. = FALSE)
}

if (opts$`separate-strands` && mode != "nipter") {
  stop("--separate-strands is only valid in nipter mode", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Set mode-dependent defaults
# ---------------------------------------------------------------------------

if (mode == "wisecondorx") {
  binsize <- opts$binsize %||% 5000L
  mapq    <- opts$mapq    %||% 1L
  rmdup   <- opts$rmdup   %||% "streaming"
} else {
  binsize <- opts$binsize %||% 50000L
  mapq    <- opts$mapq    %||% 0L
  rmdup   <- opts$rmdup   %||% "none"
}

library(RWisecondorX)

# ---------------------------------------------------------------------------
# Convert
# ---------------------------------------------------------------------------

bam <- opts$bam
out <- opts$out

if (mode == "wisecondorx") {
  if (opts$npz) {
    cat(sprintf("Converting %s → %s (WisecondorX NPZ, binsize=%d, rmdup=%s)\n",
                basename(bam), basename(out), binsize, rmdup))
    bam_convert_npz(
      bam       = bam,
      npz       = out,
      binsize   = binsize,
      rmdup     = rmdup,
      reference = opts$reference
    )
  } else {
    cat(sprintf("Converting %s → %s (WisecondorX BED, binsize=%d, mapq=%d, rmdup=%s)\n",
                basename(bam), basename(out), binsize, mapq, rmdup))
    bam_convert_bed(
      bam           = bam,
      bed           = out,
      binsize       = binsize,
      mapq          = mapq,
      require_flags = opts$`require-flags`,
      exclude_flags = opts$`exclude-flags`,
      rmdup         = rmdup,
      reference     = opts$reference
    )
  }
} else {
  # NIPTeR mode
  gc_table <- opts$`gc-table`
  gc_fasta <- opts$fasta
  corrected_sample <- NULL

  if (!is.null(gc_table) || !is.null(gc_fasta)) {
    cat(sprintf("Binning %s (NIPTeR, binsize=%d, mapq=%d, rmdup=%s) ...\n",
                basename(bam), binsize, mapq, rmdup))
    raw_sample <- nipter_bin_bam(
      bam              = bam,
      binsize          = binsize,
      mapq             = mapq,
      require_flags    = opts$`require-flags`,
      exclude_flags    = opts$`exclude-flags`,
      rmdup            = rmdup,
      separate_strands = opts$`separate-strands`,
      reference        = opts$reference
    )
    cat("Applying GC correction ...\n")
    if (!is.null(gc_table)) {
      corrected_sample <- nipter_gc_correct(raw_sample, gc_table = gc_table)
    } else {
      corrected_sample <- nipter_gc_correct(raw_sample, fasta = gc_fasta,
                                            binsize = binsize)
    }
  }

  cat(sprintf("Converting %s → %s (NIPTeR BED%s, binsize=%d, mapq=%d, rmdup=%s)\n",
              basename(bam), basename(out),
              if (opts$`separate-strands`) " 9-col" else " 5-col",
              binsize, mapq, rmdup))

  nipter_bin_bam_bed(
    bam              = bam,
    bed              = out,
    binsize          = binsize,
    mapq             = mapq,
    require_flags    = opts$`require-flags`,
    exclude_flags    = opts$`exclude-flags`,
    rmdup            = rmdup,
    separate_strands = opts$`separate-strands`,
    corrected        = corrected_sample,
    reference        = opts$reference
  )
}

cat("Done.\n")
