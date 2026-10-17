// Expose internal C++ functions for direct mathematical tests. The R runner adds
// a generated prefix identifying the source directory. Include current package
// code, rather than a frozen reference implementation. These test functions are
// compiled separately and are not installed with the package.
// [[Rcpp::depends(BH)]]
#include SCEPTRE_COMPONENT_SCORES
#include SCEPTRE_COMPONENT_TAIL
#include SCEPTRE_COMPONENT_STAGES
#include SCEPTRE_COMPONENT_LIFETIME

using namespace sceptre::inference;

// [[Rcpp::export]]
Rcpp::NumericVector component_effect(Rcpp::NumericVector y,
                                    Rcpp::NumericVector mu,
                                    Rcpp::IntegerVector treatment, int n_trt) {
  EffectEstimate effect = estimate_effect(y, mu, treatment, n_trt);
  return Rcpp::NumericVector::create(effect.fold_change, effect.standard_error);
}

// [[Rcpp::export]]
double component_observed(Rcpp::NumericVector a, Rcpp::NumericVector w,
                          Rcpp::NumericMatrix D, Rcpp::IntegerVector treatment) {
  return observed_score(a, w, D, treatment);
}

// [[Rcpp::export]]
std::vector<double> component_segment(Rcpp::NumericVector a, Rcpp::NumericVector w,
                                      Rcpp::NumericMatrix D, Rcpp::List draws,
                                      int start, int count, int n_trt,
                                      bool use_all_cells) {
  check_score_dimensions(a, w, D);
  const DrawBank bank = Rcpp::as<DrawBank>(draws);
  return null_score_segment(a, w, D, start, count, n_trt, use_all_cells, bank);
}

// [[Rcpp::export]]
SEXP component_tail_record(double z_orig, Rcpp::NumericVector a,
                                 Rcpp::NumericVector w, Rcpp::NumericMatrix D,
                                 Rcpp::List draws, int n_trt, bool use_all_cells,
                                 int B1, int B2, int B3,
                                 bool fit_parametric_curve, int side_code,
                                 double fold_change, double standard_error,
                                 bool return_resampling_dist) {
  check_score_dimensions(a, w, D);
  SEXP continuation;
  TailResult* owned_result;
  {
    sceptre::ResultContinuation token;
    std::unique_ptr<TailResult> result;
    {
      const DrawBank bank = Rcpp::as<DrawBank>(draws);
      const DrawBank before = bank;
      result.reset(new TailResult);
      *result = evaluate_staged_tail(z_orig, a, w, D, n_trt, use_all_cells,
                                     bank, B1, B2, B3,
                                     fit_parametric_curve, side_code);
      if (bank != before) Rcpp::stop("Staged evaluation mutated the draw bank.");
    }
    owned_result = result.release();
    continuation = token.release();
  }
  EffectEstimate effect = {fold_change, standard_error};
  SEXP out = association_result_record(effect, z_orig, owned_result,
                                        return_resampling_dist, continuation);
  UNPROTECT(1);
  return out;
}

// [[Rcpp::export]]
double component_empirical(std::vector<double> scores, double observed, int side) {
  return compute_empirical_p_value(scores, observed, side);
}
