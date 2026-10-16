#include "inference.h"
#include <cmath>
#include "sampling_safety.h"

using namespace Rcpp;

namespace sceptre {
namespace inference {

// Shared score shape guard: a/w length n and D columns n. The Rcpp entry
// separately checks y/mu/observed cardinality. D rows may be zero; this does not
// revalidate a fitted model's rank or scan producer-owned bank coordinates.
void check_score_dimensions(const NumericVector& a, const NumericVector& w,
                            const NumericMatrix& D) {
  if (a.size() != w.size() || a.size() != D.ncol())
    stop("Score arrays a/w and D columns must have aligned cell dimensions.");
}

// Allocate an owned vector preserving observed order/repetitions and subtract
// one once. Positivity/NA is checked here; the Rcpp entry checks upper bounds.
// The caller's R vector is never decremented in place.
IntegerVector observed_positions_zero_based(const IntegerVector& positions) {
  IntegerVector zero_based(positions.size());
  for (R_xlen_t i = 0; i < positions.size(); ++i) {
    if (positions[i] < 1) stop("Observed treatment positions must be positive and not NA.");
    zero_based[i] = positions[i] - 1;
  }
  return zero_based;
}

// Ratio of treated sums y/mu and its existing residual-based standard error,
// both on the fold-change (not log) scale. y/mu share the fit-cell view;
// trt_idxs are zero-based, n_trt entries, in caller order (repetitions count).
// Degenerate sums/zero n_trt keep ordinary floating NaN/Inf. No mutation or RNG.
EffectEstimate estimate_effect(const NumericVector& y, const NumericVector& mu, const IntegerVector& trt_idxs, int n_trt) {
  double sum_y = 0, sum_mu = 0, sum_top = 0, curr_top, n_trt_double = (double) n_trt;
  for (int i = 0; i < n_trt; i ++) {
    sum_y += y[trt_idxs[i]];
    sum_mu += mu[trt_idxs[i]];
  }
  double fc = (sum_y/n_trt_double)/(sum_mu/n_trt_double);
  for (int i = 0; i < n_trt; i ++) {
    curr_top = y[trt_idxs[i]] - fc * mu[trt_idxs[i]];
    sum_top += curr_top * curr_top;
  }
  double se = sqrt(sum_top/(sum_mu * sum_mu));
  return EffectEstimate {fc, se};
}


// Studentized efficient score: sum(a[T]) / sqrt(sum(w[T])-||sum(D[,T])||^2).
// D is eigenbasis-rows (eigen() order) x fit-cells; its row count equals the
// design coefficient count. T is zero-based and already validated by
// the owning boundary. Preserve summation order, including repeated positions;
// no clamping/new fallback for a zero or negative studentizer. Inputs read only.
double observed_score(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, const IntegerVector& trt_idxs) {
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow();

  // iterate over the rows of D
  for (int i = 0; i < D_nrow; i ++) {
    inner_sum = 0;
    for (int j = 0; j < trt_idxs.size(); j ++) {
      inner_sum += D(i, trt_idxs[j]);
    }
    lower_right += inner_sum * inner_sum;
  }

  // second, compute the lower-left hand of the denominator; also, compute the top
  for (int j = 0; j < trt_idxs.size(); j ++) {
    top += a[trt_idxs[j]];
    lower_left += w[trt_idxs[j]];
  }

  // finally, compute the z-score
  return(top/sqrt(lower_left - lower_right));
}


// Evaluate only bank rows [start_pos,start_pos+B), returning B ordered scores.
// Bank entries are already zero-based fit positions. Fixed-cardinality views
// read exactly n_trt entries per row; use_all_cells (CRT) instead uses the whole
// current row, whose length can differ from observed n_trt or be zero. Producers
// establish coordinate validity; check only selected row lengths/segment bounds,
// never a hybrid suffix that belongs to a larger fit. The complete bank remains
// const and reusable across responses/stages. The owning entry checks aligned
// shapes before native ownership; do not repeat D.ncol()/its Rcpp protection
// inside this compute-only scope. No draws, input mutation or RNG.
std::vector<double> null_score_segment(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, int start_pos, int B, int n_trt, bool use_all_cells, const DrawBank& synthetic_idxs) {
  sceptre::sampling::nonnegative(n_trt, "n_trt");
  const int end_pos = sceptre::sampling::sum(start_pos, B, "Draw segment end");
  if (static_cast<size_t>(end_pos) > synthetic_idxs.size())
    stop("Selected draw segment exceeds the bank length.");
  // initialize variables
  const int* curr_vect;
  // begin() returns Rcpp's cached pointer without allocation or materialization.
  const double* D_data = D.begin();
  std::vector<double> out(B);
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), idx;

  // iterate over the B tests
  for (int k = start_pos; k < end_pos; k ++) {
    if (synthetic_idxs[k].size() > static_cast<size_t>(std::numeric_limits<int>::max()))
      stop("Draw row exceeds native integer capacity.");
    if (use_all_cells) n_trt = synthetic_idxs[k].size();
    else if (synthetic_idxs[k].size() < static_cast<size_t>(n_trt))
      stop("Selected draw row is shorter than the requested treatment prefix.");
    curr_vect = synthetic_idxs[k].data();
    lower_right = 0;
    inner_sum = 0;

    // iterate over the rows of D
    for (int i = 0; i < D_nrow; i ++) {
      inner_sum = 0;
      for (int j = 0; j < n_trt; j ++) {
        inner_sum += D_data[i + static_cast<R_xlen_t>(D_nrow) * curr_vect[j]];
      }
      lower_right += inner_sum * inner_sum;
    }

    // second, compute the lower-left hand of the denominator; also, compute the top
    lower_left = 0;
    top = 0;
    for (int j = 0; j < n_trt; j ++) {
      idx = curr_vect[j];
      top += a[idx];
      lower_left += w[idx];
    }

    // compute the z-score
    out[k - start_pos] = top/sqrt(lower_left - lower_right);
  }

  return(out);
}

} // namespace inference
} // namespace sceptre
