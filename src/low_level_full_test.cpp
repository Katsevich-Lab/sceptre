#include <Rcpp.h>
#include <cmath>
#include "inference.h"
#include "sampling_safety.h"
#include "result_lifetime.h"

using namespace Rcpp;


// Entry called from R for one gene-guide association test.
// y, mu, a and w contain the selected cells in the same order; D has one column
// per cell and eigenvectors in rows. Observed trt_idxs start at 1, and n_trt is
// their count. Check lengths, positions and resampling sizes before converting
// observed positions once to C++ numbering. Nonfinite observed scores return
// p = 1 with stage = 0, without ranking against null scores.
// synthetic_idxs points to this worker's full set of draws, with positions
// already starting at 0. Check the pointer, then read without changing it.
// Permutations use the first n_trt entries of each draw; CRT uses the full draw,
// including empty ones. Return effects, score, tail diagnostics and, if requested,
// the final stage's null statistics. This pointer is temporary, not a saved format.
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
  SEXP continuation;
  sceptre::inference::EffectEstimate effect;
  double z_orig;
  sceptre::inference::TailResult* owned_tail;
  {
    sceptre::ResultContinuation token;
    const sceptre::inference::DrawBank* bank;
    {
      // Finish temporary Rcpp objects before allocating the C++ tail result.
      // An R error can skip C++ destructors, so keep those lifetimes separate.
      const IntegerVector observed_positions =
        sceptre::inference::observed_positions_zero_based(trt_idxs);
      effect = sceptre::inference::estimate_effect(y, mu, observed_positions, n_trt);
      z_orig = sceptre::inference::observed_score(a, w, D, observed_positions);
      XPtr<sceptre::inference::DrawBank> draws(synthetic_idxs);
      if (draws.get() == nullptr) stop("Synthetic draw bank pointer is null.");
      bank = draws.get();
    }
    std::unique_ptr<sceptre::inference::TailResult> tail(new sceptre::inference::TailResult);
    if (!std::isfinite(z_orig)) {
      // Keep the requested hypothesis in the family, conservatively. Stage 0
      // identifies this policy fallback; no null statistics were evaluated.
      tail->p = 1.0;
      tail->stage = 0;
      tail->sn_params = {{NA_REAL, NA_REAL, NA_REAL}};
    } else {
      *tail = sceptre::inference::evaluate_staged_tail(
        z_orig, a, w, D, n_trt, use_all_cells, *bank, B1, B2, B3,
        fit_parametric_curve, side_code);
    }
    owned_tail = tail.release();
    continuation = token.release();
  }
  // Transfer directly to protected conversion: no intervening R call may
  // raise an error before the conversion function takes responsibility for cleanup.
  SEXP out = sceptre::inference::association_result_record(
    effect, z_orig, owned_tail, return_resampling_dist, continuation);
  UNPROTECT(1);
  return out;
}
