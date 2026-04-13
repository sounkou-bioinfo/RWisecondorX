# RWisecondorX Agent Guidelines

This document provides guidance for AI agents working on `RWisecondorX`.

## Project Goal

Build an R toolkit for copy number analysis and trisomy prediction in
non-invasive prenatal testing (NIPT), comparable in scope to `NIPTUtils`
but built entirely on `Rduckhts`/`DuckDB` rather than
`Rsamtools`/Python. The package ports two upstream implementations:

- **WisecondorX** — copy number variation detection for cfDNA
  (Python/pysam upstream)
- **NIPTeR** — fast trisomy prediction via chromosomal fraction
  Z-scores, NCV scores, regression, and chi-squared correction
  (R/Rsamtools upstream)

Design priorities:

- **Exact conformance** with upstream implementations where documented
  and testable.
- **Performance** — DuckDB SQL replaces per-read Python loops and
  Rsamtools scans; large BAMs run in seconds, not minutes.
- **Interoperable file formats** — bgzipped, tabix-indexed BED files as
  the intermediate layer between binning and analysis, consumable by
  DuckDB, R, Python, or any tabix-aware tool.
- **No Python runtime dependency** — `Rduckhts` replaces `pysam` for all
  HTS operations; `reticulate` and `condathis` are optional and used
  only for conformance testing.

------------------------------------------------------------------------

## Implementation Status

### Completed

**Shared binning engine — `R/convert.R`** -
[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md):
DuckDB/SQL read-counting core. Supports `binsize`, `mapq`,
`require_flags` (samtools `-f`), `exclude_flags` (samtools `-F`),
`rmdup` (`"streaming"` / `"flag"` / `"none"`), CRAM via `reference`.
Achieves exact bin-for-bin conformance with WisecondorX on
HG00106.chrom11 (25,115 non-zero bins, 0 mismatches). -
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md):
bgzipped + tabix-indexed 4-column BED output (`chrom`, `start`, `end`,
`count`). Uses
[`Rduckhts::rduckhts_bgzip()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_bgzip.html)
and `rduckhts_tabix_index()` — no external tools.

**WisecondorX layer** - `R/wisecondorx_cli.R`:
[`wisecondorx_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_convert.md),
[`wisecondorx_newref()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_newref.md),
[`wisecondorx_predict()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/wisecondorx_predict.md)
— thin `condathis` wrappers delegating to the official bioconda
package. - `R/npz.R`:
[`bam_convert_npz()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_npz.md)
— WisecondorX-compatible NPZ output via `reticulate`.

**NIPTeR binning layer — `R/nipter_bin.R`** -
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md):
produces `NIPTeRSample` objects (class
`c("NIPTeRSample", "CombinedStrands")`). Exposes `mapq`,
`require_flags`, `exclude_flags`, `rmdup` for pre-filtering matching
real-world NIPT pipelines (e.g. `mapq=40L, exclude_flags=1024L`). -
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md):
5-column BED (`chrom`, `start`, `end`, `count`, `corrected_count`).
`corrected_count` is `NA` until GC correction is ported; column is
present so downstream tools can rely on the schema now.

**Tests** - `inst/tinytest/test_fixtures.R` (41 assertions): synthetic
BAM/CRAM fixtures, all three `rmdup` modes, CRAM reference round-trip. -
`inst/tinytest/test_nipter.R` (47 assertions): `NIPTeRSample` structure,
matrix dimensions and rownames, MAPQ filter, `exclude_flags=1024` dedup,
BED.gz output shape, optional bin-for-bin conformance against
[`NIPTeR::bin_bam_sample()`](https://rdrr.io/pkg/NIPTeR/man/bin_bam_sample.html)
via `NIPTER_CONFORMANCE_BAM` env var. - `inst/tinytest/test_cli_args.R`:
CLI argument builder unit tests. - `inst/tinytest/test_npz.R`: NPZ
round-trip tests. - `inst/tinytest/test_integration.R`: WisecondorX
bin-for-bin conformance via condathis. -
`inst/tinytest/test_sra_metadata.R`: SRA URL helpers.

**Infrastructure** - `R/aaa.R`: SRA metadata helpers. -
`scripts/make_fixtures.sh`: reproducible synthetic BAM/CRAM fixture
generation. - `DESCRIPTION`: NIPTeR authors (Dirk de Weerd, Lennart
Johansson) listed as `cph`. NIPTeR added to `Suggests`.

**NIPTeR statistical layer**

Each file is strictly separate; never mix NIPTeR and WisecondorX code.

- `R/nipter_control.R` —
  [`nipter_as_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_as_control_group.md):
  constructs a `NIPTeRControlGroup` from a list of `NIPTeRSample`
  objects with validation and dedup.
  [`nipter_diagnose_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_diagnose_control_group.md):
  per-chromosome Z-scores and Shapiro-Wilk normality tests across the
  control group.
  [`nipter_match_control_group()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_match_control_group.md):
  selects the best-fit control samples for a test sample by
  sum-of-squares chromosomal fraction distance.
- `R/nipter_gc.R` —
  [`nipter_gc_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_gc_correct.md):
  LOESS and bin-weight GC correction. GC content derived from the FASTA
  reference via `rduckhts_fasta_nuc()` — **no bundled correction
  tables** (unlike the original NIPTeR `sysdata.rda`). Supports both
  single samples and control groups. Sex chromosome correction via
  nearest-neighbour (LOESS) or same-bucket weights (bin-weight).
- `R/nipter_chi.R` —
  [`nipter_chi_correct()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_chi_correct.md):
  chi-squared overdispersion correction applied simultaneously to sample
  and control group. Uses normalised chi-squared threshold (default
  3.5). Sex chromosome correction uses mean autosomal correction factor
  per bin position.
- `R/nipter_score.R` —
  [`nipter_z_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_z_score.md):
  chromosomal fraction Z-score with Shapiro-Wilk normality test.
  [`nipter_ncv_score()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_ncv_score.md):
  normalised chromosome value with brute-force denominator search using
  [`utils::combn()`](https://rdrr.io/r/utils/combn.html) (replaces
  [`sets::set_combn()`](https://rdrr.io/pkg/sets/man/set.html)
  dependency).
- `R/nipter_regression.R` —
  [`nipter_regression()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md):
  forward stepwise regression Z-score with train/test split, practical
  vs theoretical CV selection. Currently CombinedStrands only;
  SeparatedStrands requires `separate_strands = TRUE` support.
- `inst/tinytest/test_nipter_stats.R` — 69 assertions covering all
  statistical functions.

### Open architectural questions

- **`separate_strands = TRUE`**: currently stubbed with
  [`stop()`](https://rdrr.io/r/base/stop.html) in
  [`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md).
  Requires `GROUP BY rname, bin, strand` in the SQL layer and returning
  two matrices (`forward_reads_matrix` and `reverse_reads_matrix`).
  Needed for the SeparatedStrands variant of
  [`nipter_regression()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_regression.md),
  which is already implemented structurally but gated on the input
  class.
- **Gender/sex prediction**: investigate whether to use NIPTeR’s
  Y-fraction method, WisecondorX’s `--yfrac` approach, or a unified
  strategy.
- **DNACopy replacement**: the segmentation step in
  `wisecondorx_predict` uses DNACopy internally; evaluate whether to
  expose or replace it.
- **Multi-chromosome NIPTeR conformance fixture**: `make fixtures`
  currently produces a chr11-only BAM. A multi-chromosome synthetic BAM
  (all 24 chroms, no unmapped reads, no same-position collisions) would
  allow `NIPTER_CONFORMANCE_BAM` to be populated automatically in CI
  without a real patient BAM.

------------------------------------------------------------------------

## Source Layout Expectations

- `R/convert.R` — shared BAM/CRAM binning engine; used by both
  WisecondorX and NIPTeR layers.
- `R/nipter_bin.R` — NIPTeR binning layer.
- `R/nipter_gc.R` — GC correction (LOESS and bin-weight).
- `R/nipter_score.R` — Z-score and NCV scoring.
- `R/nipter_chi.R` — chi-squared overdispersion correction.
- `R/nipter_control.R` — control group construction, diagnostics, and
  matching.
- `R/nipter_regression.R` — forward stepwise regression Z-score.
- `R/wisecondorx_cli.R` — CLI wrappers.
- `R/npz.R` — NPZ output.
- `R/aaa.R` — SRA metadata helpers.
- `inst/tinytest/` — unit tests (one file per feature family).
- `inst/extdata/` — synthetic BAM/CRAM fixtures.
- The WisecondorX upstream algorithm reference is
  `../../duckhts/.sync/WisecondorX/`.
- The NIPTeR upstream algorithm reference is
  `../../duckhts/.sync/NIPTeR/`.
- The WisecondorX conformance script is
  `../../duckhts/scripts/wisecondorx_convert_conformance.py`.
- There is NO `inst/python/` directory.

------------------------------------------------------------------------

## Interoperable File Formats

The BED.gz format is the language-agnostic handoff between binning and
downstream analysis:

- 4-column WisecondorX BED: `chrom`, `start`, `end`, `count` (written by
  [`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md)).
- 5-column NIPTeR BED: `chrom`, `start`, `end`, `count`,
  `corrected_count` (written by
  [`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md)).
- Coordinates are 0-based half-open intervals (BED convention).
  Chromosomes use no `chr` prefix.
- All files are bgzipped (BGZF) and tabix-indexed via
  [`Rduckhts::rduckhts_bgzip()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_bgzip.html)
  and
  [`Rduckhts::rduckhts_tabix_index()`](https://rgenomicsetl.r-universe.dev/Rduckhts/reference/rduckhts_tabix_index.html).
  Do not use [`gzfile()`](https://rdrr.io/r/base/connections.html) or
  external tools.

------------------------------------------------------------------------

## Convert Step: Exact SQL Replication of pysam larp/larp2

The WisecondorX convert step uses a streaming dedup state machine
(`larp`/`larp2`). The exact SQL equivalent uses DuckDB window functions
over `FILE_OFFSET`:

``` sql
LAG(pos) OVER (ORDER BY file_offset)                                          -- prev_pos (larp)
LAST_VALUE(CASE WHEN is_paired != 0 THEN pnext END IGNORE NULLS)              -- prev_pnext (larp2)
    OVER (ORDER BY file_offset ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING)
```

Key subtleties reproduced exactly:

- Improper pairs (`is_paired AND NOT is_proper`) are invisible to both
  larp and larp2. This filter is internal to `rmdup = "streaming"` — it
  is intrinsic to the WisecondorX algorithm, not a user flag option.
- Unpaired reads update larp but NOT larp2 — hence the `IGNORE NULLS`
  trick.
- `larp` is never reset between chromosomes — hence no `PARTITION BY`.
- Bin assignment uses integer division `pos // binsize` matching
  Python’s `int(pos / binsize)`.
- `pos` is 0-based: subtract 1 from duckhts 1-based POS before dividing.

------------------------------------------------------------------------

## Flag Filtering API

[`bam_convert()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert.md),
[`bam_convert_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/bam_convert_bed.md),
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md),
and
[`nipter_bin_bam_bed()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam_bed.md)
expose samtools-style `-f`/`-F` flag filtering:

- `require_flags` (integer bitmask, default `0L`): only reads where
  `(flag & require_flags) == require_flags` are kept. Equivalent to
  `samtools view -f`.
- `exclude_flags` (integer bitmask, default `0L`): reads where
  `(flag & exclude_flags) != 0` are dropped. Equivalent to
  `samtools view -F`.

Common real-world NIPT pre-filter pattern (matches
`samtools view --min-MQ 40 -F 1024`):

``` r
nipter_bin_bam("sample.dm.bam", mapq = 40L, exclude_flags = 1024L)
```

------------------------------------------------------------------------

## NIPTeR Layer: Known Issues and Design Decisions

### NIPTeR’s scanBam approach is “concat and pretend single-ended”

NIPTeR’s `bin_bam_sample()` calls `scanBam()` without filtering
parameters, then splits all read records by strand (`+` / `-`),
discarding only reads with strand `"*"` (unmapped). It then applies
[`unique()`](https://rdrr.io/r/base/unique.html) on positions per
chromosome per strand and sums forward and reverse counts into a single
matrix.

The consequence is that **paired-end data is treated as if you
concatenated both mates into a single-ended read pool**: R1 (forward)
and R2 (reverse) are each counted independently at their respective
start positions. A paired-end library produces roughly 2× the bin counts
of a single-end library for the same sequencing depth, because both
mates of every fragment are counted.

This creates a fundamental incompatibility: **a NIPTeR reference
(control group) built from single-end samples cannot be used to test a
paired-end sample, and vice versa.** In practice this goes unnoticed
when all samples in a cohort use the same library preparation, but
mixing single-end and paired-end data produces wildly incorrect Z-scores
and NCV values.

Our
[`nipter_bin_bam()`](https://sounkou-bioinfo.github.io/RWisecondorX/reference/nipter_bin_bam.md)
replicates this behaviour (both mates counted) as the default
(`rmdup = "none"`, no flag filters), keeping it compatible with NIPTeR
control groups built the same way. Users who need fragment-level counts
(one count per fragment, not per read) should pre-filter to R1 only
using `require_flags = 0x41L` (paired + read1) or count only unique
fragments — this is a future enhancement.

### NIPTeR’s implicit positional dedup is similar to WisecondorX streaming dedup, but different

NIPTeR’s `bin_reads()` applies
[`unique()`](https://rdrr.io/r/base/unique.html) to read positions per
chromosome per strand before binning. The effect is that two reads with
the same start position on the same chromosome and strand are counted
once, not twice. This is position-based deduplication.

WisecondorX’s `larp` also collapses consecutive reads at the same
position. The critical differences:

|                          | NIPTeR [`unique()`](https://rdrr.io/r/base/unique.html) | WisecondorX larp/larp2                     |
|--------------------------|---------------------------------------------------------|--------------------------------------------|
| Scope                    | Per-chromosome, per-strand                              | Global file order, no reset between chroms |
| Pairing awareness        | None — operates on positions only                       | larp2 uses pnext for paired reads          |
| Improper pair handling   | None                                                    | Improper pairs excluded entirely           |
| Equivalent samtools flag | Not equivalent to any flag                              | Not expressible as a single flag           |

Both collapse same-position reads, but NIPTeR’s
[`unique()`](https://rdrr.io/r/base/unique.html) also silently merges
any two reads that happen to share a start position, including distinct
reads that are not PCR duplicates. On the `hg00106_chr11_fixture.bam`
fixture (1,018 paired-end reads) NIPTeR’s approach would count 1,003 vs
1,018 from `rmdup = "none"` — a difference of 15 reads, all at positions
where the duplicate flag (`0x400`) is set on one mate.

### Two additional bugs in NIPTeR’s scanBam call

**Bug 1 — Split length mismatch on BAMs with unmapped reads.**

``` r
# NIPTeR source — bin_sample.R
splitted_reads <- split(x = reads, f = droplevels(strands[strands != "*"]))
```

`strands[strands != "*"]` has length `N − U` (U = unmapped reads) but
`reads` has length `N`. R’s
[`split()`](https://rdrr.io/r/base/split.html) requires equal lengths,
so this errors or silently misbehaves on any BAM containing unmapped
reads. NIPTeR implicitly requires `samtools view -F 4` pre-filtering.

**Bug 2 — Positional dedup via
[`unique()`](https://rdrr.io/r/base/unique.html) is flag-ignorant.**

``` r
# NIPTeR source — bin_sample.R / bin_reads()
reads <- sapply(X = unique(reads_data_frame[min_read:max_read]),
                FUN = getbins, bin_size = bin_size)
bins <- tabulate(reads, nbins = n_bins)
```

[`tabulate()`](https://rdrr.io/r/base/tabulate.html) here counts unique
*positions per bin*, not reads. A PCR duplicate at a unique position
survives; two distinct reads at the same position are silently merged.
This is not equivalent to `-F 1024`.

**Consequence for conformance testing.** Exact bin-for-bin match against
[`NIPTeR::bin_bam_sample()`](https://rdrr.io/pkg/NIPTeR/man/bin_bam_sample.html)
requires a BAM with no unmapped reads and no two reads sharing a start
position per strand. The `NIPTER_CONFORMANCE_BAM` env var must point to
such a BAM. Our `rmdup = "none"` is the closest mode but counts each
read independently rather than unique positions.

------------------------------------------------------------------------

## Agent Working Instructions

Always read the existing implementation before changing it:

1.  Read the relevant R files, package metadata, and tests before
    proposing changes.
2.  When upstream WisecondorX behaviour is unclear, consult
    `../../duckhts/.sync/WisecondorX/`. Do NOT look for `inst/python/` —
    it has been removed.
3.  When upstream NIPTeR behaviour is unclear, consult
    `../../duckhts/.sync/NIPTeR/`.
4.  For conformance testing, use the official `wisecondorx` bioconda
    package via `condathis`. Do NOT vendor Python source.
5.  Preserve CRAN-friendly package behaviour. Do not make package
    loading depend on a preconfigured external Python installation.
6.  Keep package runtime logic and conformance tooling separate.
7.  Keep WisecondorX and NIPTeR layers in **separate files**
    (`nipter_*.R` vs `wisecondorx_*.R`). The shared engine is
    `convert.R`.

------------------------------------------------------------------------

## Documentation Conventions

- `README.Rmd` is the primary editable README source. Run `make readme`
  to regenerate `README.md`.
- R documentation and `NAMESPACE` are generated from roxygen2 tags. Run
  `make rd` after changing `@export`, `@param`, or `@seealso` tags.
- Do not manually edit `NAMESPACE` or `README.md`.

## Build and Workflow Rules

- Use make targets: `make rd`, `make test`, `make readme`,
  `make fixtures`.
- Keep build steps deterministic and non-interactive.

## Conformance Testing Rules

- Conformance tests skip cleanly when `condathis`, `reticulate`, or
  required BAMs are unavailable.
- `WISECONDORX_TEST_BAM` — a human BAM for WisecondorX integration
  tests.
- `NIPTER_CONFORMANCE_BAM` — a coordinate-sorted, pre-filtered (no
  unmapped reads, no same-position strand collisions) whole-genome BAM
  for NIPTeR conformance.

## Testing Rules

- Use `tinytest`. One file per feature family under `inst/tinytest/`.
- Update tests when R-facing behaviour changes.
- Tests requiring Python or `reticulate` must skip cleanly when
  unavailable.

## Copyright and Attribution

- NIPTeR authors Dirk de Weerd and Lennart Johansson are listed as `cph`
  in `DESCRIPTION`.
- WisecondorX authors are listed as `ctb`.
- Any new upstream algorithm ported must have its original authors
  credited with an appropriate role.
