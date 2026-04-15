// KNN reference bin finding for WisecondorX newref
//
// Computes the K-nearest-neighbour reference bins for each target bin,
// measured by squared Euclidean distance in PCA-corrected space. Uses
// leave-one-chromosome-out exclusion (a bin's references must come from
// other chromosomes).
//
// Parallelized with OpenMP when cpus > 1.
//
// Performance note: R stores matrices column-major (n_bins × n_samples),
// so row i's features are at offsets i, i + n_bins, i + 2*n_bins, ...
// — a stride of n_bins doubles (~225KB for 28K bins).  We transpose to
// row-major at the start so each bin's feature vector is contiguous
// (n_samples doubles ≈ 400 bytes for 50 samples), giving ~10-20× speedup
// from cache locality.

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
List knn_reference_cpp(NumericMatrix pca_corrected,
                       IntegerVector masked_bins_per_chr,
                       IntegerVector masked_bins_per_chr_cum,
                       int refsize,
                       std::string gender,
                       int cpus) {
    int n_bins    = pca_corrected.nrow();
    int n_samples = pca_corrected.ncol();
    int n_chrs    = masked_bins_per_chr.size();

    // Output matrices
    IntegerMatrix indexes(n_bins, refsize);
    NumericMatrix distances(n_bins, refsize);
    // Initialize: indexes = 0, distances = 1.0
    std::fill(indexes.begin(), indexes.end(), 0);
    std::fill(distances.begin(), distances.end(), 1.0);

    bool is_gonosomal = (gender == "F" || gender == "M");

    // --- Transpose pca_corrected from column-major (n_bins × n_samples)
    //     to row-major layout: trans[bin * n_samples + s]
    //     This makes each bin's feature vector contiguous in memory.
    const double* pca_ptr = REAL(pca_corrected);
    std::vector<double> trans(static_cast<size_t>(n_bins) * n_samples);
    for (int s = 0; s < n_samples; ++s) {
        const double* col = pca_ptr + static_cast<size_t>(s) * n_bins;
        for (int b = 0; b < n_bins; ++b) {
            trans[static_cast<size_t>(b) * n_samples + s] = col[b];
        }
    }

    // Pre-compute chromosome start/end for each chromosome (1-based R indexing
    // converted to 0-based C indexing)
    std::vector<int> chr_start_0(n_chrs), chr_end_0(n_chrs);
    for (int c = 0; c < n_chrs; ++c) {
        chr_start_0[c] = (c == 0) ? 0 : masked_bins_per_chr_cum[c - 1];
        chr_end_0[c]   = masked_bins_per_chr_cum[c] - 1;  // inclusive
    }

    // Process each chromosome
    for (int chr_idx = 0; chr_idx < n_chrs; ++chr_idx) {
        int cs = chr_start_0[chr_idx];
        int ce = chr_end_0[chr_idx];
        int n_chr_bins = ce - cs + 1;
        if (n_chr_bins <= 0) continue;

        // For gonosomal references, skip autosomal chromosomes (0-20, i.e. chr 1-21)
        if (is_gonosomal && chr_idx < 21) continue;

        // Build the "other rows" index — all bins NOT on this chromosome
        int n_other = n_bins - n_chr_bins;
        if (n_other <= 0) continue;

        std::vector<int> other_rows(n_other);
        {
            int k = 0;
            for (int i = 0; i < cs; ++i)        other_rows[k++] = i;
            for (int i = ce + 1; i < n_bins; ++i) other_rows[k++] = i;
        }

        int k_actual = std::min(refsize, n_other);

        // Parallel over bins within this chromosome
        #pragma omp parallel for schedule(dynamic, 64) num_threads(cpus) if(cpus > 1)
        for (int bin_i = cs; bin_i <= ce; ++bin_i) {
            // Pointer to this bin's contiguous feature vector in transposed layout
            const double* this_row = trans.data() +
                                     static_cast<size_t>(bin_i) * n_samples;

            // Compute squared distances to all candidate bins
            std::vector<std::pair<double, int>> dist_idx(n_other);
            for (int j = 0; j < n_other; ++j) {
                const double* cand_row = trans.data() +
                    static_cast<size_t>(other_rows[j]) * n_samples;
                double d = 0.0;
                for (int s = 0; s < n_samples; ++s) {
                    double diff = this_row[s] - cand_row[s];
                    d += diff * diff;
                }
                dist_idx[j] = std::make_pair(d, other_rows[j]);
            }

            // Partial sort to get K nearest
            std::partial_sort(dist_idx.begin(),
                              dist_idx.begin() + k_actual,
                              dist_idx.end());

            // Store results (convert to 1-based R indexing)
            for (int j = 0; j < k_actual; ++j) {
                indexes(bin_i, j)   = dist_idx[j].second + 1;  // 1-based
                distances(bin_i, j) = dist_idx[j].first;
            }
        }
    }

    return List::create(Named("indexes") = indexes,
                        Named("distances") = distances);
}


// Compute null ratios — vectorized version
// For each null sample and each bin, compute log2(sample[bin] / median(sample[ref_bins]))
// [[Rcpp::export]]
NumericMatrix null_ratios_cpp(NumericMatrix pca_corrected,
                              IntegerMatrix indexes,
                              IntegerVector null_sample_idx,
                              int cpus) {
    int n_bins    = pca_corrected.nrow();
    int n_null    = null_sample_idx.size();
    int refsize   = indexes.ncol();
    (void)pca_corrected.ncol();  // unused but kept for clarity

    NumericMatrix null_ratios(n_bins, n_null);
    // Initialize to 0
    std::fill(null_ratios.begin(), null_ratios.end(), 0.0);

    const double* pca_ptr = REAL(pca_corrected);

    #pragma omp parallel for schedule(dynamic, 1) num_threads(cpus) if(cpus > 1)
    for (int null_i = 0; null_i < n_null; ++null_i) {
        int case_i = null_sample_idx[null_i] - 1;  // Convert to 0-based

        // Extract this sample's values (column case_i of pca_corrected)
        std::vector<double> sample_vec(n_bins);
        for (int b = 0; b < n_bins; ++b) {
            sample_vec[b] = pca_ptr[case_i * n_bins + b];
        }

        // Temporary buffer for reference values
        std::vector<double> ref_vals(refsize);

        for (int bin_i = 0; bin_i < n_bins; ++bin_i) {
            // Gather reference values
            int n_pos = 0;
            for (int j = 0; j < refsize; ++j) {
                int ref_idx = indexes(bin_i, j) - 1;  // Convert to 0-based
                if (ref_idx >= 0 && ref_idx < n_bins) {
                    double v = sample_vec[ref_idx];
                    if (v > 0.0) {
                        ref_vals[n_pos++] = v;
                    }
                }
            }

            if (n_pos == 0) continue;

            // Compute median of positive reference values
            std::sort(ref_vals.begin(), ref_vals.begin() + n_pos);
            double med;
            if (n_pos % 2 == 1) {
                med = ref_vals[n_pos / 2];
            } else {
                med = (ref_vals[n_pos / 2 - 1] + ref_vals[n_pos / 2]) / 2.0;
            }

            if (med > 0.0 && sample_vec[bin_i] > 0.0) {
                null_ratios(bin_i, null_i) = std::log2(sample_vec[bin_i] / med);
            }
        }
    }

    return null_ratios;
}
