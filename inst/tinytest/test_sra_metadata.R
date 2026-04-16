library(tinytest)
library(RWisecondorX)

expect_identical(
  sra_runinfo_url("PRJNA400134"),
  "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA400134"
)

expect_identical(
  sra_runinfo_url("PRJNA400134", mode = "term"),
  "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?term=PRJNA400134"
)

tmp_csv <- tempfile(fileext = ".csv")

fake_download <- function(url, destfile, quiet = TRUE) {
  writeLines(
    c(
      "Run,BioProject,SampleName",
      "SRR000001,PRJNA400134,control_001"
    ),
    con = destfile
  )
  0L
}

returned_path <- download_sra_runinfo(
  accession = "PRJNA400134",
  dest = tmp_csv,
  downloader = fake_download
)

expect_identical(returned_path, tmp_csv)

runinfo <- read_sra_runinfo("PRJNA400134", path = tmp_csv)

expect_identical(nrow(runinfo), 1L)
expect_identical(runinfo$Run[[1]], "SRR000001")
expect_identical(runinfo$BioProject[[1]], "PRJNA400134")
