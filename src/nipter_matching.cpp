// NIPTeR control matching — Rcpp+OpenMP SSD computation
//
// Computes sum-of-squared-differences (SSD) distances between a query
// chromosomal fraction vector and every column of a fractions matrix.
// This replaces the R-level for-loop in nipter_match_control_group() that
// becomes the bottleneck when iterating over every sample in a large control
// group (the production pipeline computes N×N pairwise SSDs).
//
// nipter_ssd_scores_cpp: one query vs N controls, returns N SSDs.
// nipter_ssd_matrix_cpp: N controls vs N controls, returns N×N SSD matrix
//   (symmetric, diagonal zero). Used to compute the mean SSD of each control
//   against all others (the "matching round" QC statistic in production code).

#include <Rcpp.h>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector nipter_ssd_scores_cpp(NumericMatrix fracs,
                                    NumericVector query,
                                    IntegerVector compare_idx,
                                    int cpus = 1) {
    // fracs:       n_chr × n_controls (column-major, one column per sample)
    // query:       n_chr-element vector for the test sample
    // compare_idx: 0-based row indices (chromosome subset to use)
    // cpus:        OpenMP thread count

    int n_controls = fracs.ncol();
    int n_cmp      = compare_idx.size();

    NumericVector scores(n_controls, 0.0);

    // Pre-extract query values for the comparison chromosomes
    std::vector<double> q(n_cmp);
    for (int i = 0; i < n_cmp; ++i) {
        q[i] = query[compare_idx[i]];
    }

    // Raw pointer for column-major access
    const double* fp = REAL(fracs);
    int n_rows = fracs.nrow();

#ifdef _OPENMP
    if (cpus > 1) omp_set_num_threads(cpus);
#pragma omp parallel for schedule(static)
#endif
    for (int s = 0; s < n_controls; ++s) {
        double ssd = 0.0;
        const double* col = fp + s * n_rows;
        for (int i = 0; i < n_cmp; ++i) {
            double diff = col[compare_idx[i]] - q[i];
            ssd += diff * diff;
        }
        scores[s] = ssd;
    }

    return scores;
}


// [[Rcpp::export]]
NumericMatrix nipter_ssd_matrix_cpp(NumericMatrix fracs,
                                    IntegerVector compare_idx,
                                    int cpus = 1) {
    // Returns the symmetric N×N SSD matrix.
    // fracs: n_chr × n_controls
    // compare_idx: 0-based row indices to compare

    int n_controls = fracs.ncol();
    int n_cmp      = compare_idx.size();
    int n_rows     = fracs.nrow();
    const double* fp = REAL(fracs);

    NumericMatrix out(n_controls, n_controls);
    std::fill(out.begin(), out.end(), 0.0);

#ifdef _OPENMP
    if (cpus > 1) omp_set_num_threads(cpus);
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n_controls; ++i) {
        const double* ci = fp + i * n_rows;
        for (int j = i + 1; j < n_controls; ++j) {
            const double* cj = fp + j * n_rows;
            double ssd = 0.0;
            for (int k = 0; k < n_cmp; ++k) {
                int r = compare_idx[k];
                double diff = ci[r] - cj[r];
                ssd += diff * diff;
            }
            out(i, j) = ssd;
            out(j, i) = ssd;
        }
    }

    return out;
}
