#include "inference.h"
#include <cmath>

using namespace Rcpp;

namespace sceptre {
namespace inference {

EffectEstimate estimate_effect(const NumericVector& y, const NumericVector& mu, const IntegerVector& trt_idxs, int n_trt) {
  double sum_y = 0, sum_mu = 0, sum_top = 0, curr_top, n_trt_double = (double) n_trt;
  for (int i = 0; i < n_trt; i ++) {
    sum_y += y[trt_idxs[i] - 1];
    sum_mu += mu[trt_idxs[i] - 1];
  }
  double fc = (sum_y/n_trt_double)/(sum_mu/n_trt_double);
  for (int i = 0; i < n_trt; i ++) {
    curr_top = y[trt_idxs[i] - 1] - fc * mu[trt_idxs[i] - 1];
    sum_top += curr_top * curr_top;
  }
  double se = sqrt(sum_top/(sum_mu * sum_mu));
  return EffectEstimate {fc, se};
}


double observed_score(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, const IntegerVector& trt_idxs) {
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow();

  // iterate over the rows of D
  for (int i = 0; i < D_nrow; i ++) {
    inner_sum = 0;
    for (int j = 0; j < trt_idxs.size(); j ++) {
      inner_sum += D(i, trt_idxs[j] - 1);
    }
    lower_right += inner_sum * inner_sum;
  }

  // second, compute the lower-left hand of the denominator; also, compute the top
  for (int j = 0; j < trt_idxs.size(); j ++) {
    top += a[trt_idxs[j] - 1];
    lower_left += w[trt_idxs[j] - 1];
  }

  // finally, compute the z-score
  return(top/sqrt(lower_left - lower_right));
}


std::vector<double> null_score_segment(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, int start_pos, int B, int n_trt, bool use_all_cells, const DrawBank& synthetic_idxs) {
  // initialize variables
  const int* curr_vect;
  std::vector<double> out(B);
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), idx;

  // iterate over the B tests
  for (int k = start_pos; k < B + start_pos; k ++) {
    if (use_all_cells) n_trt = synthetic_idxs[k].size();
    curr_vect = synthetic_idxs[k].data();
    lower_right = 0;
    inner_sum = 0;

    // iterate over the rows of D
    for (int i = 0; i < D_nrow; i ++) {
      inner_sum = 0;
      for (int j = 0; j < n_trt; j ++) {
        inner_sum += D(i, curr_vect[j]);
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

// Retain the existing native helper signature for callers in the package.
std::vector<double> estimate_log_fold_change_v2(NumericVector y, NumericVector mu, IntegerVector trt_idxs, int n_trt) {
  sceptre::inference::EffectEstimate effect = sceptre::inference::estimate_effect(y, mu, trt_idxs, n_trt);
  return std::vector<double> {effect.fold_change, effect.standard_error};
}
