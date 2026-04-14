// NIPTeR NCV denominator search — Rcpp implementation
//
// nipter_ncv_search_cpp: brute-force search over all combinations of up to
// max_elements denominator chromosomes that minimise the coefficient of
// variation (CV = sd/mean) of the NCV ratio across control samples.
//
// The NCV ratio for a given denominator set D is:
//   ratio_i = reads_focus_i / sum_{d in D}(reads_d_i)  for each control i
// CV = sd(ratios) / |mean(ratios)|
//
// Returns: integer vector of 0-based indices into 'candidates' for the
// best denominator set, plus the achieved minimum CV.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Knuth's iterative combination generator:
// Given n candidates and combination size k, iterate over all C(n,k) combos.
// Each combo is a vector of k indices (0-based) into 'candidates'.
// ---------------------------------------------------------------------------

static double compute_cv(const std::vector<double>& focus_reads,
                          const NumericMatrix&        ctrl_reads,  // n_chr × n_ctrl
                          const std::vector<int>&     denom_rows,  // row indices
                          int                         n_ctrl) {
    // Compute denominator sums and ratios across controls
    std::vector<double> ratios(n_ctrl);
    double sum_r = 0.0, sum_r2 = 0.0;

    for (int c = 0; c < n_ctrl; ++c) {
        double denom = 0.0;
        for (int r : denom_rows) {
            denom += ctrl_reads(r, c);
        }
        if (denom <= 0.0) return std::numeric_limits<double>::infinity();
        double ratio = focus_reads[c] / denom;
        ratios[c] = ratio;
        sum_r  += ratio;
        sum_r2 += ratio * ratio;
    }

    double mean_r = sum_r / n_ctrl;
    if (mean_r == 0.0) return std::numeric_limits<double>::infinity();

    // NIPTeR uses stats::sd (sample SD, n-1 divisor) — replicate that:
    double sd_r = (n_ctrl > 1)
        ? std::sqrt(sum_r2 / (n_ctrl - 1) -
                    (sum_r * sum_r) / ((double)n_ctrl * (n_ctrl - 1)))
        : 0.0;

    return sd_r / std::abs(mean_r);
}


// [[Rcpp::export]]
List nipter_ncv_search_cpp(NumericMatrix ctrl_reads,
                           IntegerVector candidates,
                           int           focus_row,
                           int           max_elements) {
    // ctrl_reads:   n_chr × n_ctrl (0-based row indices)
    // candidates:   0-based row indices of denominator candidates
    // focus_row:    0-based row index of focus chromosome
    // max_elements: maximum denominator set size

    int n_ctrl       = ctrl_reads.ncol();
    int n_candidates = candidates.size();
    max_elements     = std::min(max_elements, n_candidates);

    // Extract focus chromosome reads across all controls
    std::vector<double> focus_reads(n_ctrl);
    for (int c = 0; c < n_ctrl; ++c) {
        focus_reads[c] = ctrl_reads(focus_row, c);
    }

    // Convert candidates to std::vector for convenience
    std::vector<int> cands(candidates.begin(), candidates.end());

    double best_cv = std::numeric_limits<double>::infinity();
    std::vector<int> best_denom;

    // Iterate over combination sizes 1..max_elements
    for (int k = 1; k <= max_elements; ++k) {
        // Initialise index array [0, 1, ..., k-1]
        std::vector<int> idx(k);
        std::iota(idx.begin(), idx.end(), 0);

        while (true) {
            // Build current denominator row indices
            std::vector<int> denom_rows(k);
            for (int i = 0; i < k; ++i) denom_rows[i] = cands[idx[i]];

            double cv = compute_cv(focus_reads, ctrl_reads, denom_rows, n_ctrl);
            if (cv < best_cv) {
                best_cv    = cv;
                best_denom = denom_rows;
            }

            // Advance to next combination (Knuth Algorithm L)
            int i = k - 1;
            while (i >= 0 && idx[i] == n_candidates - k + i) --i;
            if (i < 0) break;  // all combinations exhausted for this k
            ++idx[i];
            for (int j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
        }
    }

    // Return: best_denom as 0-based row indices, best_cv
    return List::create(
        Named("denom_rows") = IntegerVector(best_denom.begin(), best_denom.end()),
        Named("best_cv")    = best_cv
    );
}
