library(tinytest)
library(RWisecondorX)
library(DBI)
library(duckdb)
library(Rduckhts)

local({
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  Rduckhts::rduckhts_load(con)

  fa <- tempfile(fileext = ".fa")
  writeLines(
    c(
      ">21",
      "AACCNNNNGG",
      ">22",
      "NNNNNNNNNN",
      ">X",
      "ACGTACGTAC",
      ">Y",
      "GGGGCCCCAA"
    ),
    con = fa
  )
  Rduckhts::rduckhts_fasta_index(con, fa)

  gc_tbl <- RWisecondorX:::.get_gc_table(fa, binsize = 10L, con = con)

  expect_equal(
    gc_tbl[["21"]][1L],
    4 / 6,
    tolerance = 1e-10,
    info = ".get_gc_table excludes N bases from the GC denominator"
  )
  expect_true(
    is.na(gc_tbl[["22"]][1L]),
    info = ".get_gc_table returns NA for all-N bins"
  )

  gc_out <- tempfile(fileext = ".tsv.bgz")
  nipter_gc_precompute(fa, binsize = 10L, out = gc_out, con = con)
  gc_bgz <- RWisecondorX:::.load_gc_table(gc_out, con = con)

  expect_equal(
    gc_bgz[["21"]][1L],
    4 / 6,
    tolerance = 1e-10,
    info = "nipter_gc_precompute writes non-N-normalized GC fractions"
  )
  expect_true(
    is.na(gc_bgz[["22"]][1L]),
    info = "nipter_gc_precompute preserves NA for all-N bins"
  )
})
