# Read RWisecondorX tabix metadata headers

Reads optional provenance metadata stored in leading comment lines of a
generic tabix-indexed TSV/BED file. Metadata lines use the form
`##RWX_<key>=<value>`.

## Usage

``` r
tabix_metadata(path, prefix = .rwx_tabix_metadata_prefix, con = NULL)
```

## Arguments

- path:

  Path to a tabix-indexed TSV/BED file.

- prefix:

  Metadata line prefix to match. Defaults to `##RWX_`.

- con:

  Optional open DBI connection with duckhts already loaded.

## Value

A named character vector of parsed metadata values. Returns an empty
named character vector when no matching metadata lines are present.

## Details

This is intended for our bgzipped, tabix-indexed TSV artifacts such as
RWisecondorX and NIPTeR BED outputs. The data rows remain headerless;
this function only inspects the optional metadata comment block.
