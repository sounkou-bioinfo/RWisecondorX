# RWisecondorX extdata

This directory holds small packaged fixtures and staged metadata snapshots used
for conformance and integration tests.

Packaged alignment fixtures:

- `hg00106_chr11_fixture.bam` / `.bai`: real chromosome 11 subset used for
  convert-step conformance checks against upstream WisecondorX
- `fixture_paired.bam` / `.bai`: tiny proper-pair plus duplicate-pair fixture
- `fixture_single.bam` / `.bai`: tiny single-end plus duplicate single-end fixture
- `fixture_mixed.bam` / `.bai`: mixed paired, single-end, duplicate-flagged, and
  improper-pair fixture for unit tests
- `fixture_mixed.cram` / `.crai`: CRAM rendering of `fixture_mixed.bam`
- `fixture_ref.fa` / `.fai`: reference FASTA for the CRAM fixture

Regenerate the synthetic fixtures with:

```sh
make fixtures
```

The generator script is `scripts/make_fixtures.sh`. It uses hand-written,
coordinate-sorted SAM records plus `samtools` to produce deterministic BAM/CRAM
files with known duplicate-handling behaviour.

SRA run metadata can be staged here during development with:

```r
RWisecondorX::download_sra_runinfo("PRJNA400134")
```

Large raw sequencing inputs should not be committed here.
