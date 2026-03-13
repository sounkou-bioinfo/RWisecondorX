# RWisecondorX extdata

This directory holds small packaged fixtures and staged metadata snapshots used
for conformance and integration tests.

SRA run metadata can be staged here during development with:

```r
RWisecondorX::download_sra_runinfo("PRJNA400134")
```

Large raw sequencing inputs should not be committed here.
