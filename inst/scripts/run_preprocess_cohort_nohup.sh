#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage:
  run_preprocess_cohort_nohup.sh [threads] [source_manifest] [jobs] [sample_count] [sample_seed]

Arguments:
  threads          Per-tool thread budget passed to preprocess_cohort.R (default: 20)
  source_manifest  Source BAM/CRAM manifest to sample from
                   (default: /mnt/data/BixCTF/NiptSeqNeo/all_bam_list_last_100_days.txt)
  jobs             Number of samples to process in parallel per stage (default: 4)
  sample_count     Number of unique BAM/CRAM paths to sample into the generated manifest (default: 500)
  sample_seed      RNG seed used for manifest sampling (default: 1995)

Environment:
  SOURCE_MANIFEST              Source BAM/CRAM manifest override
  GENERATED_MANIFEST           Path for the sampled manifest written by this wrapper
  FASTA                        Reference FASTA override
  SEQFF_TSV                    Optional precomputed SeqFF TSV/CSV passed through to preprocess_cohort.R
  Y_UNIQUE_TSV                 Optional precomputed Y-unique TSV/CSV passed through to preprocess_cohort.R
  Y_REGIONS_FILE               Explicit Y-unique regions BED/TSV path
  OUT_ROOT                     Output root override
  DRY_RUN                      1 to stop after manifest generation and command printing [default: 0]
  OVERWRITE                    1 to force regeneration of existing outputs [default: 0]
  NIPTER_GC_INCLUDE_SEX        1 to GC-correct X/Y before NIPTeR BED export [default: 1]
  NIPTER_GC_LOESS_SPAN         LOESS span for NIPTeR GC correction [default: 0.75]
  NIPTER_GC_CURVES             1 to write per-sample GC plots [default: 1]
  NIPTER_SEPARATE_STRANDS      1 to write 9-column NIPTeR BEDs [default: 1]
  SEQFF                        1 to compute SeqFF metrics [default: 1]
  WCX_WRITE_NPZ                1 to also write upstream WisecondorX NPZ files [default: 1]
  QC_PLOT_THEME                QC plot theme passed to preprocess_cohort.R [default: minimal]
  QC_PLOT_BASE_SIZE            QC plot base size [default: 11]
  WCX_BINSIZE                  WisecondorX binsize [default: 100000]
  WCX_MAPQ                     WisecondorX MAPQ filter [default: 1]
  WCX_REQUIRE_FLAGS            WisecondorX required SAM flags [default: 0]
  WCX_EXCLUDE_FLAGS            WisecondorX excluded SAM flags [default: 0]
  WCX_RMDUP                    WisecondorX duplicate mode [default: streaming]
  NIPTER_BINSIZE               NIPTeR binsize [default: 50000]
  NIPTER_MAPQ                  NIPTeR MAPQ filter [default: 40]
  NIPTER_EXCLUDE_FLAGS         NIPTeR excluded SAM flags [default: 1024]
  METRICS_PRE_MAPQ             MAPQ for nonfiltered SeqFF/Y-unique metrics [default: 0]
  METRICS_PRE_REQUIRE_FLAGS    Required flags for nonfiltered metrics [default: 0]
  METRICS_PRE_EXCLUDE_FLAGS    Excluded flags for nonfiltered metrics [default: 0]
  METRICS_POST_MAPQ            MAPQ for filtered SeqFF/Y-unique metrics [default: NIPTER_MAPQ]
  METRICS_POST_REQUIRE_FLAGS   Required flags for filtered metrics [default: 0]
  METRICS_POST_EXCLUDE_FLAGS   Excluded flags for filtered metrics [default: NIPTER_EXCLUDE_FLAGS]
  MOSDEPTH_EXCLUDE_FLAGS_BASE  Base excluded flags for mosdepth-style outputs [default: 1796]
  MOSDEPTH_FILTERED_MAPQ       MAPQ for filtered mosdepth outputs [default: NIPTER_MAPQ]
  MOSDEPTH_FILTERED_EXCLUDE_FLAGS
                               Excluded flags for filtered mosdepth outputs [default: MOSDEPTH_EXCLUDE_FLAGS_BASE]
EOF
  exit 0
fi

bool_to_r_logical() {
  case "${1}" in
    1|true|TRUE|yes|YES) echo "TRUE" ;;
    0|false|FALSE|no|NO) echo "FALSE" ;;
    *)
      echo "Expected boolean 0/1/TRUE/FALSE, got '${1}'" >&2
      exit 1
      ;;
  esac
}

THREADS="${1:-30}"
DEFAULT_SOURCE_MANIFEST="/mnt/data/BixCTF/NiptSeqNeo/all_bam_list_last_500.txt"
SOURCE_MANIFEST="${2:-${SOURCE_MANIFEST:-${DEFAULT_SOURCE_MANIFEST}}}"
JOBS="${3:-6}"
SAMPLE_COUNT="${4:-${SAMPLE_COUNT:-500}}"
SAMPLE_SEED="${5:-${SAMPLE_SEED:-1995}}"

DEFAULT_FASTA="/mnt/data/niptseq/resources/genomes/human_g1k_v37.fasta"
DEFAULT_OUT_ROOT_BASE="/mnt/data/niptseq/testdata"
DEFAULT_Y_REGIONS_FILE="${REPO_ROOT}/inst/extdata/grch37_Y_chrom_blacklist.bed"

SOURCE_STEM="$(basename "${SOURCE_MANIFEST}" .txt)"
GENERATED_MANIFEST="${GENERATED_MANIFEST:-${DEFAULT_OUT_ROOT_BASE}/manifests/${SOURCE_STEM}_sample_${SAMPLE_COUNT}_seed${SAMPLE_SEED}.txt}"
OUT_ROOT="${OUT_ROOT:-${DEFAULT_OUT_ROOT_BASE}/${SOURCE_STEM}_sample_${SAMPLE_COUNT}_seed${SAMPLE_SEED}}"
FASTA="${FASTA:-${DEFAULT_FASTA}}"
SEQFF_TSV="${SEQFF_TSV:-}"
Y_UNIQUE_TSV="${Y_UNIQUE_TSV:-}"
Y_REGIONS_FILE="${Y_REGIONS_FILE:-${DEFAULT_Y_REGIONS_FILE}}"
DRY_RUN="${DRY_RUN:-0}"

OVERWRITE="${OVERWRITE:-0}"
NIPTER_GC_INCLUDE_SEX="${NIPTER_GC_INCLUDE_SEX:-1}"
NIPTER_GC_LOESS_SPAN="${NIPTER_GC_LOESS_SPAN:-0.75}"
NIPTER_GC_CURVES="${NIPTER_GC_CURVES:-1}"
NIPTER_SEPARATE_STRANDS="${NIPTER_SEPARATE_STRANDS:-1}"
SEQFF="${SEQFF:-1}"
WCX_WRITE_NPZ="${WCX_WRITE_NPZ:-1}"
QC_PLOT_THEME="${QC_PLOT_THEME:-minimal}"
QC_PLOT_BASE_SIZE="${QC_PLOT_BASE_SIZE:-11}"

WCX_BINSIZE="${WCX_BINSIZE:-100000}"
WCX_MAPQ="${WCX_MAPQ:-1}"
WCX_REQUIRE_FLAGS="${WCX_REQUIRE_FLAGS:-0}"
WCX_EXCLUDE_FLAGS="${WCX_EXCLUDE_FLAGS:-0}"
WCX_RMDUP="${WCX_RMDUP:-streaming}"

NIPTER_BINSIZE="${NIPTER_BINSIZE:-50000}"
NIPTER_MAPQ="${NIPTER_MAPQ:-40}"
NIPTER_EXCLUDE_FLAGS="${NIPTER_EXCLUDE_FLAGS:-1024}"

METRICS_PRE_MAPQ="${METRICS_PRE_MAPQ:-0}"
METRICS_PRE_REQUIRE_FLAGS="${METRICS_PRE_REQUIRE_FLAGS:-0}"
METRICS_PRE_EXCLUDE_FLAGS="${METRICS_PRE_EXCLUDE_FLAGS:-0}"
METRICS_POST_MAPQ="${METRICS_POST_MAPQ:-${NIPTER_MAPQ}}"
METRICS_POST_REQUIRE_FLAGS="${METRICS_POST_REQUIRE_FLAGS:-0}"
METRICS_POST_EXCLUDE_FLAGS="${METRICS_POST_EXCLUDE_FLAGS:-${NIPTER_EXCLUDE_FLAGS}}"

MOSDEPTH_EXCLUDE_FLAGS_BASE="${MOSDEPTH_EXCLUDE_FLAGS_BASE:-1796}"
MOSDEPTH_FILTERED_MAPQ="${MOSDEPTH_FILTERED_MAPQ:-${NIPTER_MAPQ}}"
MOSDEPTH_FILTERED_EXCLUDE_FLAGS="${MOSDEPTH_FILTERED_EXCLUDE_FLAGS:-${MOSDEPTH_EXCLUDE_FLAGS_BASE}}"

SEQFF_LOGICAL="$(bool_to_r_logical "${SEQFF}")"
OVERWRITE_LOGICAL="$(bool_to_r_logical "${OVERWRITE}")"
NIPTER_GC_INCLUDE_SEX_LOGICAL="$(bool_to_r_logical "${NIPTER_GC_INCLUDE_SEX}")"
NIPTER_GC_CURVES_LOGICAL="$(bool_to_r_logical "${NIPTER_GC_CURVES}")"
NIPTER_SEPARATE_STRANDS_LOGICAL="$(bool_to_r_logical "${NIPTER_SEPARATE_STRANDS}")"
WCX_WRITE_NPZ_LOGICAL="$(bool_to_r_logical "${WCX_WRITE_NPZ}")"
if [[ ! -f "${SOURCE_MANIFEST}" ]]; then
  echo "Source manifest does not exist: ${SOURCE_MANIFEST}" >&2
  exit 1
fi
if [[ ! -f "${FASTA}" ]]; then
  echo "Reference FASTA does not exist: ${FASTA}" >&2
  exit 1
fi
if [[ -n "${SEQFF_TSV}" && ! -f "${SEQFF_TSV}" ]]; then
  echo "SEQFF_TSV does not exist: ${SEQFF_TSV}" >&2
  exit 1
fi
if [[ -n "${Y_UNIQUE_TSV}" && ! -f "${Y_UNIQUE_TSV}" ]]; then
  echo "Y_UNIQUE_TSV does not exist: ${Y_UNIQUE_TSV}" >&2
  exit 1
fi
if [[ ! -f "${Y_REGIONS_FILE}" ]]; then
  echo "Y-unique regions file does not exist: ${Y_REGIONS_FILE}" >&2
  exit 1
fi


mkdir -p "$(dirname "${GENERATED_MANIFEST}")"
mkdir -p "${OUT_ROOT}"
mkdir -p "${OUT_ROOT}/logs"

cd "${REPO_ROOT}"


Rscript - "${SOURCE_MANIFEST}" "${GENERATED_MANIFEST}" "${SAMPLE_COUNT}" "${SAMPLE_SEED}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
source_manifest <- args[[1L]]
generated_manifest <- args[[2L]]
sample_count <- as.integer(args[[3L]])
sample_seed <- as.integer(args[[4L]])

lines <- readLines(source_manifest, warn = FALSE)
lines <- trimws(lines)
lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
if (length(lines) >= 1L && identical(lines[[1L]], "bam")) {
  lines <- lines[-1L]
}
lines <- unique(lines)
if (length(lines) < sample_count) {
  stop(
    "Source manifest contains only ", length(lines),
    " unique BAM/CRAM paths, fewer than requested sample_count=",
    sample_count,
    call. = FALSE
  )
}
set.seed(sample_seed)
writeLines(sample(lines, sample_count), generated_manifest)
RS

echo "Running preprocess_cohort.R with PID $$"
echo "Configuration:"
echo "  REPO_ROOT=${REPO_ROOT}"
echo "  R_LIBS_USER=${R_LIBS_USER}"
echo "  SOURCE_MANIFEST=${SOURCE_MANIFEST}"
echo "  GENERATED_MANIFEST=${GENERATED_MANIFEST}"
echo "  OUT_ROOT=${OUT_ROOT}"
echo "  FASTA=${FASTA}"
echo "  SEQFF_TSV=${SEQFF_TSV:-<computed>}"
echo "  Y_UNIQUE_TSV=${Y_UNIQUE_TSV:-<computed>}"
echo "  Y_REGIONS_FILE=${Y_REGIONS_FILE}"
echo "  DRY_RUN=${DRY_RUN}"
echo "  THREADS=${THREADS}"
echo "  JOBS=${JOBS}"
echo "  SAMPLE_COUNT=${SAMPLE_COUNT}"
echo "  SAMPLE_SEED=${SAMPLE_SEED}"
echo "  WCX_BINSIZE=${WCX_BINSIZE}"
echo "  WCX_MAPQ=${WCX_MAPQ}"
echo "  WCX_REQUIRE_FLAGS=${WCX_REQUIRE_FLAGS}"
echo "  WCX_EXCLUDE_FLAGS=${WCX_EXCLUDE_FLAGS}"
echo "  WCX_RMDUP=${WCX_RMDUP}"
echo "  NIPTER_BINSIZE=${NIPTER_BINSIZE}"
echo "  NIPTER_MAPQ=${NIPTER_MAPQ}"
echo "  NIPTER_EXCLUDE_FLAGS=${NIPTER_EXCLUDE_FLAGS}"
echo "  NIPTER_GC_INCLUDE_SEX=${NIPTER_GC_INCLUDE_SEX_LOGICAL}"
echo "  NIPTER_GC_LOESS_SPAN=${NIPTER_GC_LOESS_SPAN}"
echo "  NIPTER_GC_CURVES=${NIPTER_GC_CURVES_LOGICAL}"
echo "  NIPTER_SEPARATE_STRANDS=${NIPTER_SEPARATE_STRANDS_LOGICAL}"
echo "  QC_PLOT_THEME=${QC_PLOT_THEME}"
echo "  QC_PLOT_BASE_SIZE=${QC_PLOT_BASE_SIZE}"
echo "  SEQFF=${SEQFF_LOGICAL}"
echo "  METRICS_PRE_MAPQ=${METRICS_PRE_MAPQ}"
echo "  METRICS_PRE_REQUIRE_FLAGS=${METRICS_PRE_REQUIRE_FLAGS}"
echo "  METRICS_PRE_EXCLUDE_FLAGS=${METRICS_PRE_EXCLUDE_FLAGS}"
echo "  METRICS_POST_MAPQ=${METRICS_POST_MAPQ}"
echo "  METRICS_POST_REQUIRE_FLAGS=${METRICS_POST_REQUIRE_FLAGS}"
echo "  METRICS_POST_EXCLUDE_FLAGS=${METRICS_POST_EXCLUDE_FLAGS}"
echo "  MOSDEPTH_EXCLUDE_FLAGS_BASE=${MOSDEPTH_EXCLUDE_FLAGS_BASE}"
echo "  MOSDEPTH_FILTERED_MAPQ=${MOSDEPTH_FILTERED_MAPQ}"
echo "  MOSDEPTH_FILTERED_EXCLUDE_FLAGS=${MOSDEPTH_FILTERED_EXCLUDE_FLAGS}"
echo "  WCX_WRITE_NPZ=${WCX_WRITE_NPZ_LOGICAL}"
echo "  OVERWRITE=${OVERWRITE_LOGICAL}"
echo "  generated_manifest_lines=$(wc -l < "${GENERATED_MANIFEST}")"

cmd=(
  Rscript "${REPO_ROOT}/inst/scripts/preprocess_cohort.R"
  --bam-list "${GENERATED_MANIFEST}"
  --out-root "${OUT_ROOT}"
  --fasta "${FASTA}"
  --threads "${THREADS}"
  --jobs "${JOBS}"
  --wcx-binsize "${WCX_BINSIZE}"
  --wcx-mapq "${WCX_MAPQ}"
  --wcx-require-flags "${WCX_REQUIRE_FLAGS}"
  --wcx-exclude-flags "${WCX_EXCLUDE_FLAGS}"
  --wcx-rmdup "${WCX_RMDUP}"
  --nipter-binsize "${NIPTER_BINSIZE}"
  --nipter-mapq "${NIPTER_MAPQ}"
  --nipter-exclude-flags "${NIPTER_EXCLUDE_FLAGS}"
  --nipter-gc-include-sex="${NIPTER_GC_INCLUDE_SEX_LOGICAL}"
  --nipter-gc-loess-span "${NIPTER_GC_LOESS_SPAN}"
  --nipter-gc-curves="${NIPTER_GC_CURVES_LOGICAL}"
  --nipter-separate-strands="${NIPTER_SEPARATE_STRANDS_LOGICAL}"
  --qc-plot-theme "${QC_PLOT_THEME}"
  --qc-plot-base-size "${QC_PLOT_BASE_SIZE}"
  --seqff="${SEQFF_LOGICAL}"
  --y-regions-file "${Y_REGIONS_FILE}"
  --metrics-pre-mapq "${METRICS_PRE_MAPQ}"
  --metrics-pre-require-flags "${METRICS_PRE_REQUIRE_FLAGS}"
  --metrics-pre-exclude-flags "${METRICS_PRE_EXCLUDE_FLAGS}"
  --metrics-post-mapq "${METRICS_POST_MAPQ}"
  --metrics-post-require-flags "${METRICS_POST_REQUIRE_FLAGS}"
  --metrics-post-exclude-flags "${METRICS_POST_EXCLUDE_FLAGS}"
  --mosdepth-exclude-flags-base "${MOSDEPTH_EXCLUDE_FLAGS_BASE}"
  --mosdepth-filtered-mapq "${MOSDEPTH_FILTERED_MAPQ}"
  --mosdepth-filtered-exclude-flags "${MOSDEPTH_FILTERED_EXCLUDE_FLAGS}"
  --wcx-write-npz="${WCX_WRITE_NPZ_LOGICAL}"
  --overwrite="${OVERWRITE_LOGICAL}"
)

if [[ -n "${SEQFF_TSV}" ]]; then
  cmd+=(--seqff-tsv "${SEQFF_TSV}")
fi

if [[ -n "${Y_UNIQUE_TSV}" ]]; then
  cmd+=(--y-unique-tsv "${Y_UNIQUE_TSV}")
fi

printf 'Command:\n  '
printf '%q ' "${cmd[@]}"
printf '\n'

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=1, stopping before preprocess_cohort.R execution."
  exit 0
fi

time "${cmd[@]}"
