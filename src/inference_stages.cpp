#include "inference.h"
#include "sampling_safety.h"

namespace sceptre {
namespace inference {

// Consume disjoint segments of the already-generated full bank: B1 at offset
// 0, B2 at B1, B3 at B1+B2. Stage 1 uses add-one empirical tails; p<=.02 reaches
// optional skew-normal stage 2. A declined/failed fit uses empirical stage 3.
// B3=0 deliberately retains the previous selected segment (B2 if fitting was
// attempted, otherwise B1), including its fitter-produced sorting/reversal.
// Only reached B2 needs the fitter's >=100 observations. No stage redraws or
// concatenates earlier samples. Output owns its last segment; bank stays const,
// missing SN parameters stay NA and existing sidedness/floor conventions stay.
// No RNG or file effects; selected-stage validation errors propagate unchanged.
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

// Wrap numeric components into the real R result schema p/z_orig/fc/se/stage/
// sn_params, optionally resampling_dist. Copy only the selected owned segment;
// never serialize/copy a bank pointer. No numeric transformation or mutation.
Rcpp::List association_result_record(const EffectEstimate& effect, double z_orig,
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
