#include "inference.h"
#include "shared_low_level_functions.h"
#include "sampling_safety.h"

namespace sceptre {
namespace inference {

TailResult evaluate_staged_tail(double z_orig,
                                const Rcpp::NumericVector& a,
                                const Rcpp::NumericVector& w,
                                const Rcpp::NumericMatrix& D,
                                int n_trt, bool use_all_cells,
                                const DrawBank& synthetic_idxs,
                                int B1, int B2, int B3,
                                bool fit_parametric_curve, int side_code) {
  sceptre::sampling::sum(sceptre::sampling::sum(B1, B2, "B1 + B2"), B3,
                         "B1 + B2 + B3");
  const double P_THRESH = 0.02;
  bool sn_fit_used = false;
  TailResult out;
  out.sn_params = {{NA_REAL, NA_REAL, NA_REAL}};
  out.stage = 1;

  out.null_statistics = null_score_segment(a, w, D, 0, B1, n_trt,
                                          use_all_cells, synthetic_idxs);
  out.p = compute_empirical_p_value(out.null_statistics, z_orig, side_code);

  if (out.p <= P_THRESH) {
    if (fit_parametric_curve) {
      out.null_statistics = null_score_segment(a, w, D, B1, B2, n_trt,
                                              use_all_cells, synthetic_idxs);
      // The existing fitter may sort/reverse this owned vector. Preserve that
      // selected distribution, including when B3 is zero and the fit fails.
      std::vector<double> fit_sn_out =
        fit_and_evaluate_skew_normal(z_orig, out.null_statistics, side_code);
      out.p = fit_sn_out[3];
      sn_fit_used = out.p > -0.5;
      if (sn_fit_used) for (int i = 0; i < 3; i ++) out.sn_params[i] = fit_sn_out[i];
      out.stage = 2;
    }

    if (!fit_parametric_curve || !sn_fit_used) {
      if (B3 > 0) out.null_statistics = null_score_segment(a, w, D, B1 + B2,
                                                          B3, n_trt,
                                                          use_all_cells,
                                                          synthetic_idxs);
      out.p = compute_empirical_p_value(out.null_statistics, z_orig, side_code);
      out.stage = 3;
    }
  }
  return out;
}

Rcpp::List legacy_result_record(const EffectEstimate& effect, double z_orig,
                               const TailResult& tail,
                               bool return_resampling_dist) {
  Rcpp::NumericVector sn_params = Rcpp::NumericVector::create(
    tail.sn_params[0], tail.sn_params[1], tail.sn_params[2]);
  if (return_resampling_dist) {
    return Rcpp::List::create(Rcpp::Named("p") = tail.p,
                              Rcpp::Named("z_orig") = z_orig,
                              Rcpp::Named("fc") = effect.fold_change,
                              Rcpp::Named("se") = effect.standard_error,
                              Rcpp::Named("stage") = tail.stage,
                              Rcpp::Named("sn_params") = sn_params,
                              Rcpp::Named("resampling_dist") = tail.null_statistics);
  }
  return Rcpp::List::create(Rcpp::Named("p") = tail.p,
                            Rcpp::Named("z_orig") = z_orig,
                            Rcpp::Named("fc") = effect.fold_change,
                            Rcpp::Named("se") = effect.standard_error,
                            Rcpp::Named("stage") = tail.stage,
                            Rcpp::Named("sn_params") = sn_params);
}

} // namespace inference
} // namespace sceptre
