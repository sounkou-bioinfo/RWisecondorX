library(tinytest)
library(RWisecondorX)

.helper_nipter <- system.file("tinytest", "helper_nipter.R", package = "RWisecondorX")
if (nzchar(.helper_nipter)) {
  sys.source(.helper_nipter, envir = environment())
} else {
  sys.source("inst/tinytest/helper_nipter.R", envir = environment())
}

.with_duckhts_con <- function(fun) {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  Rduckhts::rduckhts_load(con)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fun(con)
}

.check_generic_metadata_roundtrip <- function() {
  plain_bed <- tempfile(fileext = ".bed")
  bed_gz <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(c(plain_bed, bed_gz, paste0(bed_gz, ".tbi"))), add = TRUE)

  writeLines(
    c(
      "##RWX_format=rwisecondorx_bed",
      "##RWX_schema=count_v1",
      "1\t0\t50000\t3",
      "1\t50000\t100000\t4"
    ),
    con = plain_bed
  )

  .with_duckhts_con(function(con) {
    Rduckhts::rduckhts_bgzip(con, plain_bed, output_path = bed_gz,
                             overwrite = TRUE, threads = 1L)
    Rduckhts::rduckhts_tabix_index(con, bed_gz, preset = "bed",
                                   threads = 1L, comment_char = "#")
  })

  meta <- .with_duckhts_con(function(con) tabix_metadata(bed_gz, con = con))
  expect_identical(
    meta[c("format", "schema")],
    c(format = "rwisecondorx_bed", schema = "count_v1"),
    info = "tabix_metadata reads RWX key/value comment lines"
  )

  bins <- .with_duckhts_con(function(con) bed_to_sample(bed_gz, binsize = 50000L, con = con))
  expect_identical(
    as.integer(bins[["1"]][1:2]),
    c(3L, 4L),
    info = "bed_to_sample ignores optional metadata comment lines"
  )
}

.check_nipter_metadata_roundtrip <- function() {
  sample <- .sim_nipter_sample(
    c("1" = 40L, "2" = 20L),
    name = "meta_sample",
    n_bins = 4L
  )

  nipter_bed <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(c(nipter_bed, paste0(nipter_bed, ".tbi"))), add = TRUE)

  .with_duckhts_con(function(con) {
    nipter_sample_to_bed(
      sample = sample,
      bed = nipter_bed,
      binsize = 50000L,
      con = con,
      metadata = list(
        format = "nipter_bed",
        schema = "combined_v1",
        gc_include_sex = TRUE
      )
    )
  })

  nipter_meta <- .with_duckhts_con(function(con) tabix_metadata(nipter_bed, con = con))
  expect_identical(
    nipter_meta[c("format", "schema", "gc_include_sex")],
    c(format = "nipter_bed", schema = "combined_v1", gc_include_sex = "true"),
    info = "nipter_sample_to_bed writes optional RWX metadata lines"
  )

  nipter_rt <- .with_duckhts_con(function(con) bed_to_nipter_sample(nipter_bed, con = con))
  expect_true(
    S7::S7_inherits(nipter_rt, CombinedStrandsSample),
    info = "bed_to_nipter_sample ignores metadata comment lines"
  )
}

.check_bam_convert_native_stats_metadata <- function() {
  mixed_bam <- system.file("extdata", "fixture_mixed.bam", package = "Rduckhts")
  if (!nzchar(mixed_bam)) {
    return(invisible(NULL))
  }

  bed_gz <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(c(bed_gz, paste0(bed_gz, ".tbi"))), add = TRUE)

  .with_duckhts_con(function(con) {
    bam_convert_bed(
      bam = mixed_bam,
      bed = bed_gz,
      binsize = 5000L,
      rmdup = "streaming",
      con = con,
      metadata = list(format = "rwisecondorx_bed")
    )
  })

  meta <- .with_duckhts_con(function(con) tabix_metadata(bed_gz, con = con))
  expect_identical(
    meta[["format"]],
    "rwisecondorx_bed",
    info = "bam_convert_bed preserves explicit caller metadata"
  )
  expect_identical(
    meta[["native_bin_stats"]],
    "gc,mq",
    info = "bam_convert_bed appends native stats metadata when metadata output is enabled"
  )
  expect_true(
    "gc_perc_pre_weighted_mean" %in% names(meta),
    info = "bam_convert_bed metadata includes aggregated prefilter GC"
  )
  expect_true(
    "mean_mapq_post_weighted_mean" %in% names(meta),
    info = "bam_convert_bed metadata includes aggregated postfilter MAPQ"
  )

  bins <- .with_duckhts_con(function(con) bed_to_sample(bed_gz, binsize = 5000L, con = con))
  total_from_body <- sum(vapply(
    bins,
    function(x) if (is.null(x)) 0L else sum(as.integer(x)),
    integer(1L)
  ))
  expect_true(
    total_from_body > 0L,
    info = "bed_to_sample still reads a metadata-bearing bam_convert_bed body"
  )
  expect_true(
    as.numeric(meta[["count_total_sum"]]) >= total_from_body,
    info = "native metadata count summary is compatible with the BED body totals"
  )
}

.check_sample_metrics_metadata_roundtrip <- function() {
  sample <- .sim_nipter_sample(
    c("1" = 40L, "2" = 20L),
    name = "metric_sample",
    n_bins = 4L
  )

  regions_tsv <- tempfile(fileext = ".tsv")
  bed_gz <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(c(regions_tsv, bed_gz, paste0(bed_gz, ".tbi"))), add = TRUE)

  regions <- data.frame(
    Chromosome = c("Y", "Y"),
    Start = c(100L, 300L),
    End = c(200L, 400L),
    GeneName = c("GENE1", "GENE2"),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    regions,
    file = regions_tsv,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  seqff <- RWisecondorX:::.seqff_metric_record(
    pre = c(SeqFF = 0.07, Enet = 0.06, WRSC = 0.08),
    post = c(SeqFF = 0.06, Enet = 0.05, WRSC = 0.07),
    source = "computed"
  )
  y_unique <- RWisecondorX:::.y_unique_metric_record(
    pre = list(ratio = 0.001, y_unique_reads = 10L, total_nuclear_reads = 10000L),
    post = list(ratio = 0.0008, y_unique_reads = 8L, total_nuclear_reads = 10000L),
    regions_file = regions_tsv,
    source = "computed"
  )

  .with_duckhts_con(function(con) {
    nipter_sample_to_bed(
      sample = sample,
      bed = bed_gz,
      binsize = 50000L,
      con = con,
      metadata = RWisecondorX:::.sample_metrics_tabix_metadata(
        seqff = seqff,
        y_unique = y_unique,
        filters_pre = list(mapq = 0L, require_flags = 0L, exclude_flags = 0L),
        filters_post = list(mapq = 40L, require_flags = 0L, exclude_flags = 1024L)
      )
    )
  })

  meta <- .with_duckhts_con(function(con) tabix_metadata(bed_gz, con = con))
  expect_identical(meta[["ff_seqff_pre"]], "0.07",
                   info = "sample metrics metadata includes nonfiltered SeqFF")
  expect_identical(meta[["ff_seqff_post"]], "0.06",
                   info = "sample metrics metadata includes filtered SeqFF")
  expect_identical(meta[["y_unique_ratio_pre"]], "0.001",
                   info = "sample metrics metadata includes nonfiltered Y-unique ratio")
  expect_identical(meta[["y_unique_ratio_post"]], "8e-04",
                   info = "sample metrics metadata includes filtered Y-unique ratio")
  expect_identical(meta[["metrics_post_exclude_flags"]], "1024",
                   info = "sample metrics metadata includes filtered flag provenance")
  expect_true(
    identical(meta[["y_unique_regions_file"]], regions_tsv),
    info = "sample metrics metadata keeps the Y-region file provenance"
  )
  expect_true(
    !"y_unique_genes" %in% names(meta),
    info = "sample metrics metadata no longer embeds Y-unique gene labels"
  )
  expect_true(
    !"y_unique_regions" %in% names(meta),
    info = "sample metrics metadata no longer embeds explicit Y-unique interval payloads"
  )
}

.check_generic_metadata_roundtrip()
.check_nipter_metadata_roundtrip()
.check_bam_convert_native_stats_metadata()
.check_sample_metrics_metadata_roundtrip()
