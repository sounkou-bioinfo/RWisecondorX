/*  samtools_stats_summary.cpp -- minimal samtools-stats summary port

    Portions of the record-processing logic below are copied and adapted
    from samtools stats.c.

    Copyright (C) 2012-2025 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>
    Author: Sam Nicholls <sam@samnicholls.net>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdint>
#include <string>

using namespace Rcpp;

// These macros are copied from samtools stats.c so the summary semantics stay
// aligned with samtools' own flag interpretation.
#define IS_PAIRED(bam) ((bam)->core.flag&BAM_FPAIRED)
#define IS_PAIRED_AND_MAPPED(bam) (((bam)->core.flag&BAM_FPAIRED) && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))
#define IS_PROPERLYPAIRED(bam) (((bam)->core.flag&(BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR) && !((bam)->core.flag&BAM_FUNMAP))
#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#define IS_DUP(bam) ((bam)->core.flag&BAM_FDUP)

#define READ_ORDER_FIRST 1
#define READ_ORDER_LAST 2

namespace {

struct SamtoolsSummary {
  std::uint64_t nreads_filtered = 0;
  std::uint64_t nreads_1st = 0;
  std::uint64_t nreads_2nd = 0;
  std::uint64_t nreads_other = 0;
  std::uint64_t nreads_unmapped = 0;
  std::uint64_t nreads_single_mapped = 0;
  std::uint64_t nreads_paired_and_mapped = 0;
  std::uint64_t nreads_properly_paired = 0;
  std::uint64_t nreads_paired_tech = 0;
  std::uint64_t nreads_anomalous = 0;
  std::uint64_t nreads_dup = 0;
  std::uint64_t nreads_mq0 = 0;
  std::uint64_t nreads_qcfail = 0;
  std::uint64_t nreads_secondary = 0;
  std::uint64_t nreads_supplementary = 0;
  std::uint64_t total_len = 0;
  std::uint64_t total_len_1st = 0;
  std::uint64_t total_len_2nd = 0;
  std::uint64_t total_len_dup = 0;
  std::uint64_t nbases_mapped = 0;
  std::uint64_t nbases_mapped_cigar = 0;
  std::uint64_t nmismatches = 0;
};

inline double to_double(std::uint64_t x) {
  return static_cast<double>(x);
}

void collect_orig_read_stats_minimal(const bam1_t* bam_line, SamtoolsSummary& stats) {
  const int seq_len = bam_line->core.l_qseq;
  stats.total_len += static_cast<std::uint64_t>(seq_len);

  if (bam_line->core.flag & BAM_FQCFAIL) {
    stats.nreads_qcfail++;
  }
  if (bam_line->core.flag & BAM_FPAIRED) {
    stats.nreads_paired_tech++;
  }

  const std::uint32_t order = IS_PAIRED(bam_line)
    ? (IS_READ1(bam_line) ? READ_ORDER_FIRST : 0) +
        (IS_READ2(bam_line) ? READ_ORDER_LAST : 0)
    : READ_ORDER_FIRST;

  switch (order) {
    case READ_ORDER_FIRST:
      stats.nreads_1st++;
      stats.total_len_1st += static_cast<std::uint64_t>(seq_len);
      break;
    case READ_ORDER_LAST:
      stats.nreads_2nd++;
      stats.total_len_2nd += static_cast<std::uint64_t>(seq_len);
      break;
    default:
      stats.nreads_other++;
      break;
  }

  if (IS_UNMAPPED(bam_line)) {
    stats.nreads_unmapped++;
    return;
  }

  stats.nbases_mapped += static_cast<std::uint64_t>(seq_len);

  if (!bam_line->core.qual) {
    stats.nreads_mq0++;
  }

  if (!IS_PAIRED_AND_MAPPED(bam_line)) {
    stats.nreads_single_mapped++;
  } else {
    stats.nreads_paired_and_mapped++;
    if (IS_PROPERLYPAIRED(bam_line)) {
      stats.nreads_properly_paired++;
    }
    if (bam_line->core.tid != bam_line->core.mtid) {
      stats.nreads_anomalous++;
    }
  }
}

void collect_nm_and_cigar_metrics(const bam1_t* bam_line, SamtoolsSummary& stats) {
  // Number of mismatches
  uint8_t* nm = bam_aux_get(const_cast<bam1_t*>(bam_line), "NM");
  if (nm) {
    stats.nmismatches += static_cast<std::uint64_t>(bam_aux2i(nm));
  }

  // Number of mapped bases from cigar
  if (bam_line->core.n_cigar == 0) {
    stop("Encountered mapped read with no CIGAR while computing samtools-style summary.");
  }

  const std::uint32_t* cigar = bam_get_cigar(bam_line);
  for (std::uint32_t i = 0; i < bam_line->core.n_cigar; ++i) {
    const int cig = bam_cigar_op(cigar[i]);
    if (cig == BAM_CMATCH || cig == BAM_CINS || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
      stats.nbases_mapped_cigar += static_cast<std::uint64_t>(bam_cigar_oplen(cigar[i]));
    }
  }
}

}  // namespace

// [[Rcpp::export]]
List samtools_stats_summary_cpp(std::string path,
                                std::string reference = "",
                                int min_mapq = 0,
                                int require_flags = 0,
                                int exclude_flags = 0,
                                int decompression_threads = 0,
                                bool report_filtered_stream_bookkeeping = false) {
  if (min_mapq < 0) {
    stop("min_mapq must be >= 0.");
  }
  if (decompression_threads < 0) {
    stop("decompression_threads must be >= 0.");
  }

  samFile* sam = sam_open(path.c_str(), "r");
  if (sam == NULL) {
    stop(std::string("Failed to open alignment file: ") + path);
  }

  if (!reference.empty()) {
    if (hts_set_opt(sam, CRAM_OPT_REFERENCE, reference.c_str()) != 0) {
      sam_close(sam);
      stop(std::string("Failed to set CRAM reference for: ") + path);
    }
  }

  if (decompression_threads > 0) {
    if (hts_set_threads(sam, decompression_threads) != 0) {
      sam_close(sam);
      stop(std::string("Failed to configure htslib decompression threads for: ") + path);
    }
  }

  bam_hdr_t* header = sam_hdr_read(sam);
  if (header == NULL) {
    sam_close(sam);
    stop(std::string("Failed to read BAM/CRAM header for: ") + path);
  }

  bam1_t* bam_line = bam_init1();
  if (bam_line == NULL) {
    bam_hdr_destroy(header);
    sam_close(sam);
    stop("Failed to allocate bam1_t buffer.");
  }

  SamtoolsSummary stats;
  int ret = 0;
  while ((ret = sam_read1(sam, header, bam_line)) >= 0) {
    if (require_flags && (bam_line->core.flag & require_flags) != require_flags) {
      stats.nreads_filtered++;
      continue;
    }
    if (exclude_flags && (bam_line->core.flag & exclude_flags)) {
      stats.nreads_filtered++;
      continue;
    }
    if (min_mapq > 0 && bam_line->core.qual < min_mapq) {
      stats.nreads_filtered++;
      continue;
    }

    // Secondary reads don't count for most stats purposes
    if (bam_line->core.flag & BAM_FSECONDARY) {
      stats.nreads_secondary++;
      continue;
    }

    if (bam_line->core.flag & BAM_FSUPPLEMENTARY) {
      stats.nreads_supplementary++;
    }

    // If line has no sequence cannot continue
    const int seq_len = bam_line->core.l_qseq;
    if (!seq_len) {
      continue;
    }

    if (IS_DUP(bam_line)) {
      stats.total_len_dup += static_cast<std::uint64_t>(seq_len);
      stats.nreads_dup++;
    }

    if ((bam_line->core.flag & BAM_FSUPPLEMENTARY) == 0) {
      collect_orig_read_stats_minimal(bam_line, stats);
    }

    if (IS_UNMAPPED(bam_line)) {
      continue;
    }

    collect_nm_and_cigar_metrics(bam_line, stats);
  }

  bam_destroy1(bam_line);
  bam_hdr_destroy(header);
  if (sam_close(sam) != 0) {
    stop(std::string("Failed to close alignment file cleanly: ") + path);
  }

  if (ret < -1) {
    stop(std::string("htslib failed while reading alignment records from: ") + path);
  }

  const std::uint64_t sequences = stats.nreads_1st + stats.nreads_2nd + stats.nreads_other;
  const std::uint64_t raw_total_sequences = report_filtered_stream_bookkeeping
    ? sequences
    : stats.nreads_filtered + sequences;
  const std::uint64_t filtered_sequences = report_filtered_stream_bookkeeping
    ? 0
    : stats.nreads_filtered;
  const std::uint64_t reads_mapped = stats.nreads_paired_and_mapped + stats.nreads_single_mapped;
  const double error_rate = stats.nbases_mapped_cigar
    ? static_cast<double>(stats.nmismatches) / static_cast<double>(stats.nbases_mapped_cigar)
    : 0.0;
  const double duplicated_read_fraction = sequences
    ? static_cast<double>(stats.nreads_dup) / static_cast<double>(sequences)
    : 0.0;
  const double duplicated_base_fraction = stats.total_len
    ? static_cast<double>(stats.total_len_dup) / static_cast<double>(stats.total_len)
    : 0.0;

  List out(34);
  CharacterVector names(34);
  int i = 0;

  out[i] = path;
  names[i++] = "path";
  out[i] = reference.empty() ? CharacterVector::create(NA_STRING) : CharacterVector::create(reference);
  names[i++] = "reference";
  out[i] = min_mapq;
  names[i++] = "min_mapq";
  out[i] = require_flags;
  names[i++] = "require_flags";
  out[i] = exclude_flags;
  names[i++] = "exclude_flags";
  out[i] = decompression_threads;
  names[i++] = "decompression_threads";
  out[i] = report_filtered_stream_bookkeeping;
  names[i++] = "report_filtered_stream_bookkeeping";
  out[i] = to_double(raw_total_sequences);
  names[i++] = "raw_total_sequences";
  out[i] = to_double(filtered_sequences);
  names[i++] = "filtered_sequences";
  out[i] = to_double(sequences);
  names[i++] = "sequences";
  out[i] = to_double(stats.nreads_1st);
  names[i++] = "first_fragments";
  out[i] = to_double(stats.nreads_2nd);
  names[i++] = "last_fragments";
  out[i] = to_double(stats.nreads_other);
  names[i++] = "other_fragments";
  out[i] = to_double(reads_mapped);
  names[i++] = "reads_mapped";
  out[i] = to_double(stats.nreads_paired_and_mapped);
  names[i++] = "reads_mapped_and_paired";
  out[i] = to_double(stats.nreads_unmapped);
  names[i++] = "reads_unmapped";
  out[i] = to_double(stats.nreads_properly_paired);
  names[i++] = "reads_properly_paired";
  out[i] = to_double(stats.nreads_paired_tech);
  names[i++] = "reads_paired";
  out[i] = to_double(stats.nreads_dup);
  names[i++] = "reads_duplicated";
  out[i] = to_double(stats.nreads_mq0);
  names[i++] = "reads_mq0";
  out[i] = to_double(stats.nreads_qcfail);
  names[i++] = "reads_qc_failed";
  out[i] = to_double(stats.nreads_secondary);
  names[i++] = "non_primary_alignments";
  out[i] = to_double(stats.nreads_supplementary);
  names[i++] = "supplementary_alignments";
  out[i] = to_double(stats.nreads_anomalous / 2);
  names[i++] = "pairs_on_different_chromosomes";
  out[i] = to_double(stats.total_len);
  names[i++] = "total_length";
  out[i] = to_double(stats.total_len_1st);
  names[i++] = "total_first_fragment_length";
  out[i] = to_double(stats.total_len_2nd);
  names[i++] = "total_last_fragment_length";
  out[i] = to_double(stats.nbases_mapped);
  names[i++] = "bases_mapped";
  out[i] = to_double(stats.nbases_mapped_cigar);
  names[i++] = "bases_mapped_cigar";
  out[i] = to_double(stats.total_len_dup);
  names[i++] = "bases_duplicated";
  out[i] = to_double(stats.nmismatches);
  names[i++] = "mismatches_from_nm";
  out[i] = error_rate;
  names[i++] = "error_rate";
  out[i] = duplicated_read_fraction;
  names[i++] = "duplicated_read_fraction";
  out[i] = duplicated_base_fraction;
  names[i++] = "duplicated_base_fraction";

  out.attr("names") = names;
  return out;
}
