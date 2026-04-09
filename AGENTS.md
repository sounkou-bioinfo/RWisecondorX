# RWisecondorX Agent Guidelines

This document provides guidance for AI agents working on `RWisecondorX`.

## Project Goal

Build an R package that exposes WisecondorX-style copy number analysis
workflows on top of `Rduckhts`, while maintaining exact conformance with
the upstream WisecondorX Python implementation.

## Current Direction

- Runtime R package functionality avoids requiring Python by default.
- `condathis` is optional and used only for reproducible conformance
  test environments; it installs the official `wisecondorx` package from
  bioconda, not any vendored copy.
- The vendored Python implementation has been **removed** from
  `inst/python/wisecondorx/` to avoid licence drift and maintenance
  burden. The canonical upstream reference is
  `../../duckhts/.sync/WisecondorX/` (algorithm reference only) and the
  official bioconda package (execution reference).
- `reticulate` remains an optional dependency used only for conformance
  tests and optional helpers.
- Conformance work invokes the official `wisecondorx` CLI via
  `condathis` and compares bin count vectors against the DuckDB SQL
  implementation.

## Convert Step: Exact SQL Replication of pysam larp/larp2

The WisecondorX convert step uses a streaming dedup state machine
(`larp`/`larp2`). The exact SQL equivalent uses DuckDB window functions
over `FILE_OFFSET` (a column added to `read_bam()` by the duckhts
extension that captures the BGZF virtual offset after each read, giving
monotonic BAM file order):

``` sql
LAG(pos) OVER (ORDER BY file_offset)                                          -- prev_pos (larp)
LAST_VALUE(CASE WHEN is_paired != 0 THEN pnext END IGNORE NULLS)              -- prev_pnext (larp2)
    OVER (ORDER BY file_offset ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING)
```

Key subtleties reproduced exactly: - Improper pairs
(`is_paired AND NOT is_proper`) are invisible to both larp and larp2. -
Unpaired reads update larp but NOT larp2 — hence the `IGNORE NULLS`
trick. - `larp` is never reset between chromosomes — hence no
`PARTITION BY`. - Bin assignment uses integer division `pos // binsize`
(not `/`), matching Python’s `int(pos / binsize)`. - `pos` is 0-based:
subtract 1 from duckhts 1-based POS before dividing.

The conformance script is at
`../../duckhts/scripts/wisecondorx_convert_conformance.py` and achieves
exact bin-for-bin match on HG00106.chrom11 (25,115 non-zero bins, 0
mismatches).

## Agent Working Instructions

Always read the existing implementation before changing it: 1. Read the
relevant R files, package metadata, and tests. 2. When upstream
WisecondorX behavior is unclear, consult the local upstream mirror at
`../../duckhts/.sync/WisecondorX` (algorithm reference). Do NOT look for
a `inst/python/` directory — it has been removed. 3. For conformance,
use the official `wisecondorx` bioconda package via `condathis`, not any
vendored copy. 4. Preserve CRAN-friendly package behavior. Do not make
package loading depend on a preconfigured external Python installation.
5. Keep package runtime logic and conformance tooling separate: use
`reticulate` only for optional Python-backed checks or helpers, and
`condathis` only for optional reproducible test environments.

## Source Layout Expectations

- R package code lives under `R/`.
- Tinytest files live under `inst/tinytest/`.
- There is NO `inst/python/` directory — vendored Python has been
  removed.
- The canonical upstream algorithm reference is
  `../../duckhts/.sync/WisecondorX/`.
- The conformance script is
  `../../duckhts/scripts/wisecondorx_convert_conformance.py`.

## Documentation Conventions

- `README.Rmd` is the primary editable README source.
- Keep package metadata in `DESCRIPTION` aligned with the actual runtime
  strategy.
- If Python dependency requirements change, update test helpers and
  `DESCRIPTION` together.
- R documentation and `NAMESPACE` should be generated from roxygen2 tags
  in `R/` files.
- Do not manually edit generated artifacts when a source file exists
  upstream of them.

## Generated Files Rules

- Do not hand-edit `NAMESPACE`; regenerate it from roxygen2.
- Do not hand-edit `README.md`; edit `README.Rmd` and render it.
- Treat generated files as outputs of the package workflow, not primary
  editing targets.

## Build And Workflow Rules

- Prefer the repository `Makefile` over ad hoc command sequences when an
  appropriate target exists.
- Use make targets for common tasks such as building, roxygenizing,
  testing, and package checks.
- Keep build steps deterministic and non-interactive.
- Keep CRAN-facing package behavior intact when changing build or test
  workflows.

## Python Integration Rules

- Do not make package load depend on `reticulate` or Python
  availability.
- Use `reticulate` primarily in tests, conformance helpers, or
  explicitly optional workflows.
- If optional Python-backed helpers are added, guard them with
  [`requireNamespace("reticulate", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html).
- Avoid hard-coding machine-local Python paths in package code.
- Keep version constraints as loose as practical; pin only when
  conformance requires it.

## Conformance Testing Rules

- Conformance tests are optional and skip cleanly when `condathis` or
  `reticulate` are unavailable.
- Install the official `wisecondorx` bioconda package via `condathis`
  for reproducible reference runs. Do NOT vendor Python source.
- The `rmdup` parameter in the R convert function must have at least
  three modes:
  - `"streaming"` — exact larp/larp2 replication via FILE_OFFSET LAG
    (WisecondorX default)
  - `"none"` — no deduplication (recommended for NIPT per WisecondorX
    docs: `--normdup`)
  - `"flag"` — use SAM 0x400 flag (for Picard/sambamba pre-marked BAMs)
- Reuse ideas from `../../duckhts/scripts/vendor_conformance_data.sh`:
  deterministic setup, local staging, and checksum manifests where
  relevant.
- The conformance reference for the convert step is
  `../../duckhts/scripts/wisecondorx_convert_conformance.py`, which
  achieves exact bin-for-bin match.

## Testing Rules

- Use `tinytest` as the default R package test framework.
- Add package tests under `inst/tinytest/`, preferably one file per
  feature family or workflow.
- When R-facing behavior changes, update or add `tinytest` coverage
  rather than relying only on manual checks.
- Prefer running tests through the repository `Makefile` when a suitable
  target exists.
- Tests that depend on Python or `reticulate` must skip cleanly when
  those dependencies are unavailable.

## Build And Packaging Rules

- The R package tarball should exclude repo-maintenance files such as
  Python packaging metadata, container files, and agent-only
  documentation.
- Do not exclude `inst/python/wisecondorx`; it is part of the package
  payload.
- Keep `.Rbuildignore` aligned with the repository layout whenever new
  top-level maintenance files are added.
