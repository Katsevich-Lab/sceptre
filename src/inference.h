#ifndef SCEPTRE_INFERENCE_H
#define SCEPTRE_INFERENCE_H

#include <Rcpp.h>
#include <array>
#include <vector>

// Private numeric components of the existing full-score test. Observed
// treatment indices are R one-based positions; draw-bank entries are already
// native zero-based positions. Only the compatibility wrapper accepts SEXP.
namespace sceptre {
namespace inference {

typedef std::vector<std::vector<int>> DrawBank;

struct EffectEstimate {
  double fold_change;
  double standard_error;
};

struct TailResult {
  double p;
  int stage;
  std::array<double, 3> sn_params;
  // The last selected segment, including the legacy sorting performed by a
  // skew-normal fit. This owns its storage and never aliases the draw bank.
  std::vector<double> null_statistics;
};

EffectEstimate estimate_effect(const Rcpp::NumericVector& y,
                               const Rcpp::NumericVector& mu,
                               const Rcpp::IntegerVector& trt_idxs, int n_trt);

double observed_score(const Rcpp::NumericVector& a,
                      const Rcpp::NumericVector& w,
                      const Rcpp::NumericMatrix& D,
                      const Rcpp::IntegerVector& trt_idxs);

std::vector<double> null_score_segment(const Rcpp::NumericVector& a,
                                       const Rcpp::NumericVector& w,
                                       const Rcpp::NumericMatrix& D,
                                       int start_pos, int B, int n_trt,
                                       bool use_all_cells,
                                       const DrawBank& synthetic_idxs);

TailResult evaluate_staged_tail(double z_orig,
                                const Rcpp::NumericVector& a,
                                const Rcpp::NumericVector& w,
                                const Rcpp::NumericMatrix& D,
                                int n_trt, bool use_all_cells,
                                const DrawBank& synthetic_idxs,
                                int B1, int B2, int B3,
                                bool fit_parametric_curve, int side_code);

Rcpp::List legacy_result_record(const EffectEstimate& effect, double z_orig,
                               const TailResult& tail,
                               bool return_resampling_dist);

} // namespace inference
} // namespace sceptre

#endif
