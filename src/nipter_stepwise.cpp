// NIPTeR forward stepwise OLS predictor selection — Rcpp implementation
//
// Replaces the R-level loop that calls lm() for every candidate at every
// step of the greedy forward stepwise regression in nipter_regression().
//
// Algorithm: incremental Gram-Schmidt orthogonalisation.
//   - Centre y and all candidate columns (removes the implicit intercept).
//   - Maintain an orthonormal basis for the column space of selected predictors.
//   - At each step, for every remaining candidate:
//       * Compute its component orthogonal to the current basis (v_perp).
//       * Compute RSS_new = RSS_current - dot(y_resid, v_perp)^2 / ||v_perp||^2.
//       * Compute adj.R^2 = 1 - (RSS_new / (n-k-1)) / (TSS / (n-1)).
//       * Select the candidate with the highest adj.R^2.
//   - Update y_resid and the perp components of remaining candidates.
//
// This is O(n_train * n_candidates * n_step) — same asymptotic as the R
// version, but ~50-100x faster because it avoids formula parsing, S3
// dispatch, data.frame allocation, and full QR re-decomposition at each step.
//
// Input:
//   fracs      n_chr × n_train   (training subset, column-major)
//   focus_row  0-based row index of focus chromosome (response)
//   cand_rows  0-based row indices of candidate predictors
//   n_step     maximum number of predictors to select
//
// Returns: IntegerVector of 0-based indices INTO cand_rows, in selection order.
//   Length ≤ n_step (shorter if candidates become collinear).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector nipter_stepwise_cpp(NumericMatrix fracs,
                                  int           focus_row,
                                  IntegerVector cand_rows,
                                  int           n_step) {
    int n_train  = fracs.ncol();
    int n_cands  = cand_rows.size();
    n_step = std::min(n_step, n_cands);

    if (n_train < 2 || n_cands == 0 || n_step <= 0)
        return IntegerVector(0);

    // -----------------------------------------------------------------------
    // Extract and centre the response vector y
    // -----------------------------------------------------------------------
    std::vector<double> y(n_train);
    {
        double mean_y = 0.0;
        for (int i = 0; i < n_train; ++i) mean_y += fracs(focus_row, i);
        mean_y /= n_train;
        for (int i = 0; i < n_train; ++i) y[i] = fracs(focus_row, i) - mean_y;
    }

    // TSS (total sum of squares): fixed throughout
    double TSS = 0.0;
    for (double v : y) TSS += v * v;

    if (TSS < 1e-15) return IntegerVector(0);   // degenerate: y is constant

    // -----------------------------------------------------------------------
    // Extract and centre candidate columns.
    // perp[c] accumulates the component of candidate c that is orthogonal
    // to all already-selected basis vectors.  Initialised to centred column.
    // -----------------------------------------------------------------------
    std::vector<std::vector<double>> perp(n_cands, std::vector<double>(n_train));
    for (int c = 0; c < n_cands; ++c) {
        int row = cand_rows[c];
        double col_mean = 0.0;
        for (int i = 0; i < n_train; ++i) col_mean += fracs(row, i);
        col_mean /= n_train;
        for (int i = 0; i < n_train; ++i)
            perp[c][i] = fracs(row, i) - col_mean;
    }

    // -----------------------------------------------------------------------
    // Running residual of y (starts as y_centred)
    // -----------------------------------------------------------------------
    std::vector<double> y_resid = y;
    double RSS = TSS;

    std::vector<bool> available(n_cands, true);
    std::vector<int>  sel;
    sel.reserve(n_step);

    for (int step = 0; step < n_step; ++step) {
        // After adding this predictor the model has (step+1) predictors + intercept
        // df_residual = n - (step+1) - 1 = n - step - 2
        int df_res = n_train - step - 2;
        if (df_res <= 0) break;

        double best_adj_r2 = -std::numeric_limits<double>::infinity();
        int    best_c      = -1;

        for (int c = 0; c < n_cands; ++c) {
            if (!available[c]) continue;

            double norm_sq = 0.0;
            for (double v : perp[c]) norm_sq += v * v;

            if (norm_sq < 1e-12) continue;   // numerically collinear — skip

            double proj = 0.0;
            for (int i = 0; i < n_train; ++i) proj += y_resid[i] * perp[c][i];

            double rss_new = RSS - proj * proj / norm_sq;
            // adj.R^2 = 1 - (RSS_new / df_res) / (TSS / (n-1))
            double adj_r2 = 1.0 - (rss_new / df_res) / (TSS / (n_train - 1));

            if (adj_r2 > best_adj_r2) {
                best_adj_r2 = adj_r2;
                best_c      = c;
            }
        }

        if (best_c < 0) break;   // no usable candidate found

        sel.push_back(best_c);
        available[best_c] = false;

        // -------------------------------------------------------------------
        // Build the new orthonormal basis vector q from perp[best_c]
        // -------------------------------------------------------------------
        double norm_sq = 0.0;
        for (double v : perp[best_c]) norm_sq += v * v;
        double norm = std::sqrt(norm_sq);

        std::vector<double> q(n_train);
        for (int i = 0; i < n_train; ++i) q[i] = perp[best_c][i] / norm;

        // Update y_resid: remove component along q
        double proj_y = 0.0;
        for (int i = 0; i < n_train; ++i) proj_y += y_resid[i] * q[i];
        for (int i = 0; i < n_train; ++i) y_resid[i] -= proj_y * q[i];
        RSS -= proj_y * proj_y;
        if (RSS < 0.0) RSS = 0.0;   // numerical guard

        // Update perp for all remaining candidates: remove component along q
        for (int c = 0; c < n_cands; ++c) {
            if (!available[c]) continue;
            double proj_c = 0.0;
            for (int i = 0; i < n_train; ++i) proj_c += perp[c][i] * q[i];
            for (int i = 0; i < n_train; ++i) perp[c][i] -= proj_c * q[i];
        }
    }

    return IntegerVector(sel.begin(), sel.end());
}
