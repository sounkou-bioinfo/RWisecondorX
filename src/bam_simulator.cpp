#include <Rcpp.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include <algorithm>
#include <cctype>
#include <random>
#include <string>

using namespace Rcpp;

namespace {

std::string normalize_chr_name_cpp(const std::string& x) {
  std::string out = x;
  if (out.size() >= 3 &&
      (out[0] == 'c' || out[0] == 'C') &&
      (out[1] == 'h' || out[1] == 'H') &&
      (out[2] == 'r' || out[2] == 'R')) {
    out = out.substr(3);
  }
  return out;
}

int find_tid_by_normalized_name(const sam_hdr_t* header,
                                const std::string& target_chr) {
  const std::string want = normalize_chr_name_cpp(target_chr);
  for (int tid = 0; tid < header->n_targets; ++tid) {
    const char* name = sam_hdr_tid2name(header, tid);
    if (name == NULL) {
      continue;
    }
    if (normalize_chr_name_cpp(name) == want) {
      return tid;
    }
  }
  return -1;
}

void stop_errno(const std::string& msg) {
  Rcpp::stop("%s (errno=%d)", msg, errno);
}

} // namespace

// [[Rcpp::export]]
List simulate_trisomy_bam_cpp(const std::string& input_path,
                              const std::string& output_path,
                              const std::string& trisomy_chr,
                              const double other_keep_prob,
                              const int seed = 1,
                              const int threads = 1,
                              const std::string& reference = "") {
  if (other_keep_prob < 0.0 || other_keep_prob > 1.0) {
    Rcpp::stop("other_keep_prob must be in [0, 1]");
  }
  if (threads < 1) {
    Rcpp::stop("threads must be >= 1");
  }

  samFile* in = sam_open(input_path.c_str(), "r");
  if (in == NULL) {
    stop_errno("Failed to open input alignment: " + input_path);
  }

  if (!reference.empty() &&
      hts_set_opt(in, CRAM_OPT_REFERENCE, reference.c_str()) < 0) {
    sam_close(in);
    Rcpp::stop("Failed to set CRAM reference on input alignment");
  }

  if (threads > 1 && hts_set_threads(in, threads) != 0) {
    sam_close(in);
    Rcpp::stop("Failed to enable input alignment threads");
  }

  sam_hdr_t* header = sam_hdr_read(in);
  if (header == NULL) {
    sam_close(in);
    Rcpp::stop("Failed to read input alignment header");
  }

  const int target_tid = find_tid_by_normalized_name(header, trisomy_chr);
  if (target_tid < 0) {
    sam_hdr_destroy(header);
    sam_close(in);
    Rcpp::stop("Chromosome '%s' not found in input alignment header", trisomy_chr);
  }

  samFile* out = sam_open(output_path.c_str(), "wb");
  if (out == NULL) {
    sam_hdr_destroy(header);
    sam_close(in);
    stop_errno("Failed to open output BAM: " + output_path);
  }

  if (threads > 1 && hts_set_threads(out, threads) != 0) {
    sam_close(out);
    sam_hdr_destroy(header);
    sam_close(in);
    Rcpp::stop("Failed to enable output BAM threads");
  }

  if (sam_hdr_write(out, header) < 0) {
    sam_close(out);
    sam_hdr_destroy(header);
    sam_close(in);
    stop_errno("Failed to write output BAM header");
  }

  bam1_t* record = bam_init1();
  if (record == NULL) {
    sam_close(out);
    sam_hdr_destroy(header);
    sam_close(in);
    Rcpp::stop("Failed to allocate BAM record");
  }

  std::mt19937_64 rng(static_cast<uint64_t>(seed));
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  double input_records = 0.0;
  double output_records = 0.0;
  double target_input = 0.0;
  double target_output = 0.0;
  double other_input = 0.0;
  double other_output = 0.0;
  double unmapped_input = 0.0;
  double unmapped_output = 0.0;

  int status = 0;
  while ((status = sam_read1(in, header, record)) >= 0) {
    input_records += 1.0;

    const int tid = record->core.tid;
    double keep_prob = 1.0;
    if (tid < 0) {
      unmapped_input += 1.0;
      keep_prob = 1.0;
    } else if (tid == target_tid) {
      target_input += 1.0;
      keep_prob = 1.0;
    } else {
      other_input += 1.0;
      keep_prob = other_keep_prob;
    }

    if (keep_prob >= 1.0 || unif(rng) < keep_prob) {
      if (sam_write1(out, header, record) < 0) {
        bam_destroy1(record);
        sam_close(out);
        sam_hdr_destroy(header);
        sam_close(in);
        stop_errno("Failed to write output BAM record");
      }

      output_records += 1.0;
      if (tid < 0) {
        unmapped_output += 1.0;
      } else if (tid == target_tid) {
        target_output += 1.0;
      } else {
        other_output += 1.0;
      }
    }
  }

  if (status < -1) {
    bam_destroy1(record);
    sam_close(out);
    sam_hdr_destroy(header);
    sam_close(in);
    Rcpp::stop("Error while reading input alignment");
  }

  bam_destroy1(record);
  if (sam_close(out) < 0) {
    sam_hdr_destroy(header);
    sam_close(in);
    stop_errno("Failed to close output BAM");
  }
  sam_hdr_destroy(header);
  if (sam_close(in) < 0) {
    stop_errno("Failed to close input alignment");
  }

  return List::create(
    Named("input_records") = input_records,
    Named("output_records") = output_records,
    Named("target_input") = target_input,
    Named("target_output") = target_output,
    Named("other_input") = other_input,
    Named("other_output") = other_output,
    Named("unmapped_input") = unmapped_input,
    Named("unmapped_output") = unmapped_output,
    Named("target_chr") = normalize_chr_name_cpp(trisomy_chr),
    Named("other_keep_prob") = other_keep_prob
  );
}
