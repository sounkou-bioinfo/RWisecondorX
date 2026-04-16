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

.PHONY: install check build rd test test-source test-fast test-io test-nipter test-rwisecondorx test-real-nipter test-seqff clean help

TEST_RUN = R -e

install: build
	THREADS=12 R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz

.PHONY: check

check: build
	R CMD check --as-cran  $(PKG_NAME)_$(PKG_VERSION).tar.gz
readme:
	R -e "rmarkdown::render('README.Rmd', output_format = rmarkdown::github_document(html_preview = FALSE))"
rd:
	R -e 'roxygen2::roxygenize(".",load_code=NULL)'

test: install
	time R -e "tinytest::test_package('$(PKG_NAME)')"

test-source:
	R -e "tinytest::run_test_dir('inst/tinytest')"

test-fast:
	$(TEST_RUN) "tinytest::run_test_file('inst/tinytest/test_cli_args.R'); tinytest::run_test_file('inst/tinytest/test_seqff.R'); tinytest::run_test_file('inst/tinytest/test_sra_metadata.R')"

test-io:
	$(TEST_RUN) "tinytest::run_test_file('inst/tinytest/test_bed_reader.R'); tinytest::run_test_file('inst/tinytest/test_npz.R'); tinytest::run_test_file('inst/tinytest/test_integration.R'); tinytest::run_test_file('inst/tinytest/test_newref_bed_dir.R')"

test-nipter:
	$(TEST_RUN) "tinytest::run_test_file('inst/tinytest/test_nipter.R'); tinytest::run_test_file('inst/tinytest/test_nipter_stats.R'); tinytest::run_test_file('inst/tinytest/test_nipter_matching.R'); tinytest::run_test_file('inst/tinytest/test_nipter_sex.R'); tinytest::run_test_file('inst/tinytest/test_nipter_reference_frame.R'); tinytest::run_test_file('inst/tinytest/test_nipter_control_qc.R'); tinytest::run_test_file('inst/tinytest/test_nipter_conformance.R')"

test-rwisecondorx:
	$(TEST_RUN) "tinytest::run_test_file('inst/tinytest/test_rwisecondorx.R')"

test-seqff:
	$(TEST_RUN) "tinytest::run_test_file('inst/tinytest/test_seqff.R')"

test-real-nipter: install
	THREADS=$${THREADS:-20} NIPTER_REAL_BAM_ENABLE=1 NIPTER_REAL_BAM_LIST=$${NIPTER_REAL_BAM_LIST:-/mnt/data/BixCTF/NiptSeqNeo/all_bam_list_sample_500.txt} NIPTER_REAL_BAM_LIMIT=$${NIPTER_REAL_BAM_LIMIT:-50} R -e "tinytest::run_test_file('inst/tinytest/test_nipter_real_manifest.R')"

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
	@echo "  test-fast  - Run small fast source tests"
	@echo "  test-io    - Run convert/BED/NPZ/integration source tests"
	@echo "  test-nipter - Run NIPTeR source tests"
	@echo "  test-rwisecondorx - Run native WisecondorX source tests"
	@echo "  test-seqff - Run SeqFF source tests"
	@echo "  test-real-nipter - Run opt-in internal NIPTeR real-BAM structural validation"
	@echo "  clean      - Clean build artifacts"
	@echo "  help       - Show this help message"
