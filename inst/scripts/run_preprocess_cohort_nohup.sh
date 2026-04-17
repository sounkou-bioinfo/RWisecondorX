#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage:
  run_preprocess_cohort_nohup.sh [threads] [manifest] [jobs]

Arguments:
  threads   Per-tool thread budget passed to preprocess_cohort.R (default: 20)
  manifest  BAM/CRAM manifest path (default: /mnt/data/BixCTF/NiptSeqNeo/all_bam_list_sample_500.txt)
  jobs      Number of samples to process in parallel per stage (default: 4)

Environment:
  FASTA     Reference FASTA (default: /mnt/data/niptseq/resources/genomes/human_g1k_v37.fasta)
  OUT_ROOT  Output root (default: /mnt/data/niptseq/testdata/<manifest-basename>)
  OVERWRITE Set to 1 to force regeneration of existing outputs
  NIPTER_GC_INCLUDE_SEX Set to 1 to GC-correct X/Y before NIPTeR BED export
  TABIX_METADATA Set to 1 to write optional ##RWX_<key>=<value> BED metadata headers
EOF
  exit 0
fi

THREADS="${1:-20}"
MANIFEST="${2:-/mnt/data/BixCTF/NiptSeqNeo/all_bam_list_sample_500.txt}"
JOBS="${3:-4}"
FASTA="${FASTA:-/mnt/data/niptseq/resources/genomes/human_g1k_v37.fasta}"
OUT_ROOT="${OUT_ROOT:-/mnt/data/niptseq/testdata/$(basename "${MANIFEST}" .txt)}"
OVERWRITE="${OVERWRITE:-0}"
NIPTER_GC_INCLUDE_SEX="${NIPTER_GC_INCLUDE_SEX:-0}"
TABIX_METADATA="${TABIX_METADATA:-0}"

mkdir -p "${OUT_ROOT}"
mkdir -p "${OUT_ROOT}/logs"

cd "${REPO_ROOT}"

# echo current process pid
echo "Running preprocess_cohort.R with PID $$"

cmd=(
  Rscript "${REPO_ROOT}/inst/scripts/preprocess_cohort.R"
  --bam-list "${MANIFEST}" \
  --out-root "${OUT_ROOT}" \
  --fasta "${FASTA}" \
  --threads "${THREADS}" \
  --jobs "${JOBS}" \
  --wcx-binsize 100000 \
  --nipter-binsize 50000 \
  --nipter-mapq 40 \
  --nipter-exclude-flags 1024 \
  --wcx-write-npz
)

if [[ "${OVERWRITE}" == "1" ]]; then
  cmd+=(--overwrite)
fi

if [[ "${NIPTER_GC_INCLUDE_SEX}" == "1" ]]; then
  cmd+=(--nipter-gc-include-sex)
fi

if [[ "${TABIX_METADATA}" == "1" ]]; then
  cmd+=(--tabix-metadata)
fi

time "${cmd[@]}"
