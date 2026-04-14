# Makefile for RWisecondorX package

# Package name
PKG_NAME := RWisecondorX
PKG_VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

# Default target
.PHONY: all
all: build

# Install package dependencies (if needed)
.PHONY: deps
deps:
	R -e "for (pkg in c('tinytest', 'roxygen2', 'BiocManager')) if (!require(pkg, quietly = TRUE, character.only = TRUE)) install.packages(pkg, repos = 'https://cran.r-project.org')"
	R -e "if (!require('DNAcopy', quietly = TRUE)) BiocManager::install('DNAcopy', ask = FALSE, update = FALSE)"
	R -e "if (!require('ParDNAcopy', quietly = TRUE)) install.packages('ParDNAcopy', repos = c('https://gsun2018.r-universe.dev', 'https://cloud.r-project.org'))"

# Build the package
.PHONY: build
build: deps
	R CMD build .

.PHONY: install check build rd test test-source fixtures cohort conformance nipter-fixture clean help

install: build
	THREADS=4 R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz

.PHONY: check

check: build
	R CMD check --as-cran  $(PKG_NAME)_$(PKG_VERSION).tar.gz
readme:
	R -e "rmarkdown::render('README.Rmd', output_format = rmarkdown::github_document(html_preview = FALSE))"
fixtures:
	bash scripts/make_fixtures.sh
cohort:
	R -e "library(RWisecondorX); generate_cohort('cohort_out')"
rd:
	R -e 'roxygen2::roxygenize(".",load_code=NULL)'

test: install
	time R -e "tinytest::test_package('$(PKG_NAME)')"

test-source:
	R -e "tinytest::run_test_dir('inst/tinytest')"

conformance:
	R -e "tinytest::run_test_file('inst/tinytest/test_wisecondorx_e2e.R')"

nipter-fixture:
	Rscript inst/scripts/make_nipter_fixture.R

# Clean build artifacts
.PHONY: clean
clean:
	rm -f $(PKG_NAME)_*.tar.gz
	rm -rf $(PKG_NAME).Rcheck/
	rm -f src/*.o src/*.so

# Help target
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all        - Build the package (default)"
	@echo "  deps       - Install package dependencies"
	@echo "  build      - Build source package"
	@echo "  install    - Build and install package locally"
	@echo "  check      - Run R CMD check --as-cran"
	@echo "  rd         - Regenerate NAMESPACE and man/ from roxygen2"
	@echo "  test       - Run installed-package tinytest tests"
	@echo "  test-source - Run tinytest directly from the source tree"
	@echo "  fixtures   - Regenerate packaged BAM/CRAM fixtures under inst/extdata"
	@echo "  cohort     - Generate synthetic BAM cohort into cohort_out/"
	@echo "  conformance - Run E2E WisecondorX conformance test (Python arm conditional on condathis)"
	@echo "  nipter-fixture - Build multi-chromosome NIPTeR conformance BAM fixture"
	@echo "  clean      - Clean build artifacts"
	@echo "  help       - Show this help message"
