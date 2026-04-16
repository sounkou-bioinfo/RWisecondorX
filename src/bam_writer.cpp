#include <Rcpp.h>

#include <htslib/sam.h>

#include <sstream>
#include <string>

using namespace Rcpp;

namespace {

std::string build_sam_header(const CharacterVector& chr_names,
                             const IntegerVector& chr_lengths) {
  std::ostringstream out;
  out << "@HD\tVN:1.6\tSO:coordinate\n";
  for (R_xlen_t i = 0; i < chr_names.size(); ++i) {
    out << "@SQ\tSN:" << Rcpp::as<std::string>(chr_names[i])
        << "\tLN:" << chr_lengths[i] << "\n";
  }
  return out.str();
}

void stop_errno(const std::string& msg) {
  Rcpp::stop("%s (errno=%d)", msg, errno);
}

} // namespace

// [[Rcpp::export]]
double write_synthetic_bam_cpp(const std::string& path,
                               const CharacterVector& chr_names,
                               const IntegerVector& chr_lengths,
                               const List& positions_by_chr,
                               const int read_len = 36,
                               const int mapq = 60,
                               const int threads = 1) {
  if (chr_names.size() != chr_lengths.size() ||
      chr_names.size() != positions_by_chr.size()) {
    Rcpp::stop("chr_names, chr_lengths, and positions_by_chr must have equal length");
  }
  if (read_len < 1) {
    Rcpp::stop("read_len must be >= 1");
  }
  if (mapq < 0 || mapq > 255) {
    Rcpp::stop("mapq must be in [0, 255]");
  }

  const std::string header_text = build_sam_header(chr_names, chr_lengths);
  sam_hdr_t* header = sam_hdr_parse(header_text.size(), header_text.c_str());
  if (header == NULL) {
    Rcpp::stop("Failed to build SAM header");
  }

  samFile* fp = sam_open(path.c_str(), "wb");
  if (fp == NULL) {
    sam_hdr_destroy(header);
    stop_errno("Failed to open BAM for writing: " + path);
  }

  if (threads > 1 && hts_set_threads(fp, threads) != 0) {
    sam_close(fp);
    sam_hdr_destroy(header);
    Rcpp::stop("Failed to enable BAM writer threads");
  }

  if (sam_hdr_write(fp, header) < 0) {
    sam_close(fp);
    sam_hdr_destroy(header);
    stop_errno("Failed to write BAM header: " + path);
  }

  bam1_t* record = bam_init1();
  if (record == NULL) {
    sam_close(fp);
    sam_hdr_destroy(header);
    Rcpp::stop("Failed to allocate BAM record");
  }

  const std::string seq(read_len, 'A');
  const uint32_t cigar = (static_cast<uint32_t>(read_len) << BAM_CIGAR_SHIFT) |
    BAM_CMATCH;

  double n_written = 0.0;

  for (R_xlen_t chr_i = 0; chr_i < positions_by_chr.size(); ++chr_i) {
    IntegerVector positions = positions_by_chr[chr_i];
    const int chr_len = chr_lengths[chr_i];
    int prev_pos = 0;

    for (R_xlen_t j = 0; j < positions.size(); ++j) {
      const int pos1 = positions[j];
      if (IntegerVector::is_na(pos1)) {
        bam_destroy1(record);
        sam_close(fp);
        sam_hdr_destroy(header);
        Rcpp::stop("Encountered NA position in synthetic read positions");
      }
      if (pos1 < 1 || pos1 + read_len - 1 > chr_len) {
        bam_destroy1(record);
        sam_close(fp);
        sam_hdr_destroy(header);
        Rcpp::stop("Synthetic read position out of range for chromosome");
      }
      if (pos1 < prev_pos) {
        bam_destroy1(record);
        sam_close(fp);
        sam_hdr_destroy(header);
        Rcpp::stop("Synthetic read positions must be sorted within chromosome");
      }
      prev_pos = pos1;

      const std::string qname = "r" + std::to_string(static_cast<long long>(n_written) + 1LL);
      if (bam_set1(record,
                   qname.size(), qname.c_str(),
                   0, static_cast<int32_t>(chr_i), pos1 - 1, static_cast<uint8_t>(mapq),
                   1, &cigar,
                   -1, -1, 0,
                   static_cast<size_t>(read_len), seq.c_str(), NULL,
                   0) < 0) {
        bam_destroy1(record);
        sam_close(fp);
        sam_hdr_destroy(header);
        stop_errno("Failed to populate BAM record");
      }

      if (sam_write1(fp, header, record) < 0) {
        bam_destroy1(record);
        sam_close(fp);
        sam_hdr_destroy(header);
        stop_errno("Failed to write BAM record");
      }

      n_written += 1.0;
    }
  }

  bam_destroy1(record);
  if (sam_close(fp) < 0) {
    sam_hdr_destroy(header);
    stop_errno("Failed to close BAM writer");
  }
  sam_hdr_destroy(header);

  return n_written;
}
