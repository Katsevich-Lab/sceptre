#include <Rcpp.h>
#include "inference.h"

using namespace Rcpp;

// Keep the existing native helper signatures while ordinary numeric components
// own the arithmetic. Observed R indices are converted once at these adapters;
// synthetic draw indices are already zero-based and are never converted again.
double compute_observed_full_statistic_v2(NumericVector a, NumericVector w,
                                          NumericMatrix D, IntegerVector trt_idxs) {
  return sceptre::inference::observed_score(a, w, D,
    sceptre::inference::observed_positions_zero_based(trt_idxs));
}

std::vector<double> compute_null_full_statistics(const NumericVector& a,
                                                const NumericVector& w,
                                                const NumericMatrix& D,
                                                int start_pos, int B, int n_trt,
                                                bool use_all_cells,
                                                SEXP synthetic_idxs) {
  XPtr<sceptre::inference::DrawBank> draws(synthetic_idxs);
  return sceptre::inference::null_score_segment(a, w, D, start_pos, B, n_trt,
                                               use_all_cells, *draws);
}

// [[Rcpp::export]]
SEXP run_low_level_test_full_v4(NumericVector y,
                                NumericVector mu,
                                NumericVector a,
                                NumericVector w,
                                NumericMatrix D,
                                IntegerVector trt_idxs,
                                int n_trt,
                                bool use_all_cells,
                                SEXP synthetic_idxs,
                                int B1,
                                int B2,
                                int B3,
                                bool fit_parametric_curve,
                                bool return_resampling_dist,
                                int side_code) {
  // The compatibility boundary alone unwraps R's external pointer. Numeric
  // components receive a const bank reference with worker-local lifetime.
  const IntegerVector observed_positions =
    sceptre::inference::observed_positions_zero_based(trt_idxs);
  sceptre::inference::EffectEstimate effect =
    sceptre::inference::estimate_effect(y, mu, observed_positions, n_trt);
  double z_orig = sceptre::inference::observed_score(a, w, D, observed_positions);
  XPtr<sceptre::inference::DrawBank> draws(synthetic_idxs);
  sceptre::inference::TailResult tail = sceptre::inference::evaluate_staged_tail(
    z_orig, a, w, D, n_trt, use_all_cells, *draws, B1, B2, B3,
    fit_parametric_curve, side_code);
  return sceptre::inference::legacy_result_record(effect, z_orig, tail,
                                                 return_resampling_dist);
}
