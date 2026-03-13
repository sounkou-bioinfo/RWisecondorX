# RWisecondorX Agent Guidelines

This document provides guidance for AI agents working on `RWisecondorX`.

## Project Goal
Build an R package that exposes WisecondorX-style copy number analysis workflows on top of `Rduckhts`, while keeping conformance with the upstream Python implementation bundled under `inst/python/wisecondorx`.

## Current Direction
- Runtime Python integration should use `reticulate`.
- `condathis` is optional and should be used for reproducible conformance environments, not as the primary runtime path.
- The package vendors the upstream Python reference implementation inside `inst/python/wisecondorx`.
- Conformance work should compare package behavior against the bundled Python implementation or an explicitly managed reference environment.

## Agent Working Instructions
Always read the existing implementation before changing it:
1. Read the relevant R files, package metadata, and tests.
2. Inspect the bundled Python sources under `inst/python/wisecondorx` before changing behavior that mirrors upstream WisecondorX.
3. When upstream behavior is unclear, consult the local upstream mirror at `../../duckhts/.sync/WisecondorX` before relying on secondary summaries.
4. Preserve CRAN-friendly package behavior. Do not make package loading depend on a preconfigured external Python installation.
5. Keep runtime integration and conformance tooling separate: `reticulate` for package code, `condathis` only for optional reproducible tests.

## Source Layout Expectations
- R package code lives under `R/`.
- Tinytest files live under `inst/tinytest/`.
- Bundled Python sources live under `inst/python/wisecondorx/`.
- The canonical upstream reference mirror is `../../duckhts/.sync/WisecondorX/`.

## Documentation Conventions
- `README.Rmd` is the primary editable README source.
- Keep package metadata in `DESCRIPTION` aligned with the actual runtime strategy.
- If Python dependency requirements change, update both `DESCRIPTION` and `R/zzz.R`.
- R documentation and `NAMESPACE` should be generated from roxygen2 tags in `R/` files.
- Do not manually edit generated artifacts when a source file exists upstream of them.

## Generated Files Rules
- Do not hand-edit `NAMESPACE`; regenerate it from roxygen2.
- Do not hand-edit `README.md`; edit `README.Rmd` and render it.
- Treat generated files as outputs of the package workflow, not primary editing targets.

## Build And Workflow Rules
- Prefer the repository `Makefile` over ad hoc command sequences when an appropriate target exists.
- Use make targets for common tasks such as building, roxygenizing, testing, and package checks.
- Keep build steps deterministic and non-interactive.
- Keep CRAN-facing package behavior intact when changing build or test workflows.

## Python Integration Rules
- Declare Python requirements with `reticulate::py_require()` in `.onLoad()`.
- Prefer delay-loaded imports when adding Python module bindings from R.
- Avoid hard-coding machine-local Python paths in package code.
- Keep version constraints as loose as practical; pin only when conformance requires it.

## Conformance Testing Rules
- Conformance tests should be optional and skip cleanly when Python dependencies are unavailable.
- If you need reproducible external environments or CLI execution, follow the `duckhts` pattern of using scripts for reproducible setup rather than ad hoc commands.
- Reuse ideas from `../../duckhts/scripts/vendor_conformance_data.sh`: deterministic setup, local staging, and checksum manifests where relevant.

## Testing Rules
- Use `tinytest` as the default R package test framework.
- Add package tests under `inst/tinytest/`, preferably one file per feature family or workflow.
- When R-facing behavior changes, update or add `tinytest` coverage rather than relying only on manual checks.
- Prefer running tests through the repository `Makefile` when a suitable target exists.

## Build And Packaging Rules
- The R package tarball should exclude repo-maintenance files such as Python packaging metadata, container files, and agent-only documentation.
- Do not exclude `inst/python/wisecondorx`; it is part of the package payload.
- Keep `.Rbuildignore` aligned with the repository layout whenever new top-level maintenance files are added.
