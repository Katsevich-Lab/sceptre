#include <Rcpp.h>
#include "inference.h"
#include "sampling_safety.h"

using namespace Rcpp;


// Sole Rcpp association boundary. y/mu/a/w have n ordered fit cells and D is
// nuisance-columns x n; observed trt_idxs are R one-based positions and n_trt
// is their observed cardinality. Validate dimensions/positions/budget arithmetic
// before converting observed positions exactly once. Empty/degenerate numeric
// fits retain their existing NaN/Inf arithmetic rather than a new fallback.
// synthetic_idxs is a live worker-local XPtr owning a zero-based full DrawBank;
// validate its address, then borrow it const. Fixed-size tests use n_trt prefixes;
// CRT tests use each selected row's entire cardinality, including empty rows.
// Return the unchanged effect/score/tail schema, optionally the last selected
// null segment. Components allocate their own result and never mutate inputs
// or the bank; no RNG or file effects. This is not a persisted pointer format.
// [[Rcpp::export]]
SEXP run_association_test(NumericVector y,
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
  sceptre::inference::check_score_dimensions(a, w, D);
  if (y.size() != a.size() || mu.size() != a.size())
    stop("Response, fitted means and score arrays must have aligned cell dimensions.");
  if (n_trt < 0 || trt_idxs.size() != n_trt)
    stop("Observed n_trt must equal the observed treatment position count.");
  for (int position : trt_idxs) {
    if (position < 1 || position > a.size())
      stop("Observed treatment position is outside the aligned cell range.");
  }
  sceptre::sampling::sum(sceptre::sampling::sum(B1, B2, "B1 + B2"), B3,
                         "B1 + B2 + B3");
  const IntegerVector observed_positions =
    sceptre::inference::observed_positions_zero_based(trt_idxs);
  sceptre::inference::EffectEstimate effect =
    sceptre::inference::estimate_effect(y, mu, observed_positions, n_trt);
  double z_orig = sceptre::inference::observed_score(a, w, D, observed_positions);
  XPtr<sceptre::inference::DrawBank> draws(synthetic_idxs);
  if (draws.get() == nullptr) stop("Synthetic draw bank pointer is null.");
  sceptre::inference::TailResult tail = sceptre::inference::evaluate_staged_tail(
    z_orig, a, w, D, n_trt, use_all_cells, *draws, B1, B2, B3,
    fit_parametric_curve, side_code);
  return sceptre::inference::association_result_record(effect, z_orig, tail,
                                                       return_resampling_dist);
}
