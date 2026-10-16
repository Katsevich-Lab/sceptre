#ifndef SCEPTRE_INFERENCE_H
#define SCEPTRE_INFERENCE_H

#include <Rcpp.h>
#include <array>
#include <vector>

// Private numeric components of the existing full-score test. Both observed
// treatment indices and draw-bank entries are native zero-based positions.
// The sole Rcpp association entry converts R positions and borrows its XPtr.
namespace sceptre {
namespace inference {

// B ordered rows; permutation consumers select fixed-cardinality prefixes,
// whereas CRT consumers select entire variable-length rows. A hybrid suffix
// can belong to a larger fit view than the current prefix. Reuse this bank
// const across response tests and disjoint stage segments; no source-cell IDs.
typedef std::vector<std::vector<int>> DrawBank;

// Adapters check aligned cell arrays and observed positions once. Draw-bank
// producers establish coordinate validity: fixed-size consumers read only the
// requested prefix, while CRT consumers read the complete current row.
void check_score_dimensions(const Rcpp::NumericVector& a,
                            const Rcpp::NumericVector& w,
                            const Rcpp::NumericMatrix& D);

// Own a converted copy of positive R one-based observed positions, preserving
// order/repeats. The owning Rcpp entry checks the current fit upper bound.
Rcpp::IntegerVector observed_positions_zero_based(const Rcpp::IntegerVector& positions);

// Fold-change ratio and its residual-based SE, not log-transformed.
struct EffectEstimate {
  double fold_change;
  double standard_error;
};

struct TailResult {
  double p;
  int stage;
  std::array<double, 3> sn_params;
  // The last selected segment, including observable sorting performed by a
  // skew-normal fit. This owns its storage and never aliases the draw bank.
  std::vector<double> null_statistics;
};

// y/mu are n-cell arrays; trt_idxs are n_trt zero-based ordered positions.
// Preserve degenerate floating results and summation order; inputs read only.
EffectEstimate estimate_effect(const Rcpp::NumericVector& y,
                               const Rcpp::NumericVector& mu,
                               const Rcpp::IntegerVector& trt_idxs, int n_trt);

// Efficient studentized score from aligned a/w (length n), D (q x n) and
// zero-based treatment positions. No hidden fit, draw or numeric fallback.
double observed_score(const Rcpp::NumericVector& a,
                      const Rcpp::NumericVector& w,
                      const Rcpp::NumericMatrix& D,
                      const Rcpp::IntegerVector& trt_idxs);

// B ordered null scores from the bank segment [start_pos,start_pos+B).
// Fixed n_trt prefix unless use_all_cells selects the entire CRT row. Caller
// validates aligned score shapes before native result ownership begins.
std::vector<double> null_score_segment(const Rcpp::NumericVector& a,
                                       const Rcpp::NumericVector& w,
                                       const Rcpp::NumericMatrix& D,
                                       int start_pos, int B, int n_trt,
                                       bool use_all_cells,
                                       const DrawBank& synthetic_idxs);

// Empirical inclusive tails with add-one correction; two-sided values are
// twice the smaller tail, deliberately not clamped to one.
double compute_empirical_p_value(const std::vector<double>& null_statistics,
                                 double z_orig, int side);

// Fits/checks an owned selected segment (>=100 scores), sorting it and possibly
// leaving descending order after a left-tail check. Returns xi/omega/alpha/p;
// p=-1 signals the existing empirical fallback. Never mutates the draw bank.
std::vector<double> fit_and_evaluate_skew_normal(double z_orig,
                                               std::vector<double>& null_statistics,
                                               int side_code);

// Stage threshold .02; optional accepted SN fit, otherwise empirical fallback.
// B3=0 retains the previously selected segment rather than concatenating draws.
TailResult evaluate_staged_tail(double z_orig,
                                const Rcpp::NumericVector& a,
                                const Rcpp::NumericVector& w,
                                const Rcpp::NumericMatrix& D,
                                int n_trt, bool use_all_cells,
                                const DrawBank& synthetic_idxs,
                                int B1, int B2, int B3,
                                bool fit_parametric_curve, int side_code);

// Consume a raw owned tail only after its compute/unique_ptr scope has ended.
// continuation is already PROTECTed before native ownership; the raw-C-API leaf
// copies the exact selected segment/schema, then deletes tail on success or R
// jump. No R operation may intervene between ownership release and this call.
// Caller UNPROTECTs its continuation once after the returned SEXP is available.
SEXP association_result_record(const EffectEstimate& effect, double z_orig,
                               TailResult* tail, bool return_resampling_dist,
                               SEXP continuation);

} // namespace inference
} // namespace sceptre

#endif
