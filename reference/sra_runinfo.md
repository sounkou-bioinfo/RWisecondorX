# SRA Run Metadata Utilities

These helpers standardize how `RWisecondorX` stages SRA run metadata for
conformance fixtures. They use the NCBI SRA backend `runinfo` endpoint
and default to writing project metadata into `inst/extdata/` in the
source tree.

## Usage

``` r
sra_runinfo_url(accession, mode = c("acc", "term"))

download_sra_runinfo(
  accession,
  dest = file.path("inst", "extdata", paste0(accession, "_runinfo.csv")),
  overwrite = FALSE,
  quiet = TRUE,
  mode = c("acc", "term"),
  downloader = utils::download.file
)

read_sra_runinfo(
  accession,
  path = system.file("extdata", paste0(accession, "_runinfo.csv"), package =
    "RWisecondorX")
)
```

## Arguments

- accession:

  A study or project accession such as `PRJNA400134`.

- mode:

  Retrieval mode. `"acc"` queries the backend with `acc=<accession>`;
  `"term"` queries with `term=<accession>`.

- dest:

  Output CSV path. Defaults to `inst/extdata/<accession>_runinfo.csv` in
  the current source checkout.

- overwrite:

  Logical; overwrite an existing file.

- quiet:

  Logical; passed through to the downloader.

- downloader:

  Function used to retrieve the URL. Defaults to
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html).
  This is injectable so tests can avoid network calls.

- path:

  Path to the run metadata CSV. Defaults to the bundled file under
  `inst/extdata/<accession>_runinfo.csv`.

## Value

`sra_runinfo_url()` returns a length-1 character URL.
`download_sra_runinfo()` returns the output path, invisibly.
`read_sra_runinfo()` returns a `data.frame`.

## Details

The SRA Run Selector and backend `runinfo` endpoints can expose slightly
different columns. For reproducible package fixtures, `RWisecondorX`
standardizes on the backend `runinfo` CSV output.

## Examples

``` r
sra_runinfo_url("PRJNA400134")
#> [1] "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA400134"

if (FALSE) { # \dontrun{
download_sra_runinfo("PRJNA400134")
metadata <- read_sra_runinfo("PRJNA400134")
} # }
```
