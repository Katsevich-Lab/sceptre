#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
#include "shared_low_level_functions.h"


double compute_observed_resid_statistic(const NumericVector& resids, int s, IntegerVector trt_idxs) {
  double tot_sum = 0, trt_sum = 0, control_sum = 0, stat = 0;
  int n_control = resids.size() - s;
  for (int i = 0; i < resids.size(); i ++) tot_sum += resids[i];
  for (int i = 0; i < s; i ++) trt_sum += resids[trt_idxs[i] - 1];
  control_sum = tot_sum - trt_sum;
  stat = trt_sum/s - control_sum/n_control;
  return(stat);
}


std::vector<double> compute_null_resid_statistics(const NumericVector& resids, int start_pos, int n_stats_to_compute, int s, SEXP synthetic_idxs) {
  // dereference the synthetic treatment idxs
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idxs);

  // initialize variables
  int* curr_vect;
  std::vector<double> out(n_stats_to_compute);
  double tot_sum = 0, trt_sum = 0, control_sum = 0, stat = 0;

  // compute control_sum and n_control
  int n_control = resids.size() - s;
  for (int i = 0; i < resids.size(); i ++) tot_sum += resids[i];

  for (int k = start_pos; k < n_stats_to_compute + start_pos; k ++) {
    trt_sum = 0;
    curr_vect = &(*synth_idx_list)[k][0];
    for (int j = 0; j < s; j ++) {
      trt_sum += resids[curr_vect[j]];
    }
    control_sum = tot_sum - trt_sum;
    stat = trt_sum/s - control_sum/n_control;
    out[k - start_pos] = stat;
  }
  return(out);
}


// [[Rcpp::export]]
SEXP run_low_level_test_resid(NumericVector y,
                              NumericVector mu,
                              NumericVector resids,
                              IntegerVector trt_idxs,
                              int n_trt,
                              SEXP synthetic_idxs,
                              int B1,
                              int B2,
                              int B3,
                              bool fit_parametric_curve,
                              bool return_resampling_dist,
                              int side_code) {
  double P_THRESH = 0.02, p;
  bool sn_fit_used = false;
  NumericVector sn_params = NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  std::vector<double> fit_sn_out;
  int stage = 1;
  List out;

  // estimate the log fold change
  double lfc = estimate_log_fold_change_v2(y, mu, trt_idxs, n_trt);

  // compute the original statistic
  double z_orig = compute_observed_resid_statistic(resids, n_trt, trt_idxs);

  // compute the stage 1 vector of null statistics and p-value
  std::vector<double> null_statistics = compute_null_resid_statistics(resids, 0, B1, n_trt, synthetic_idxs);
  p = compute_empirical_p_value(null_statistics, z_orig, side_code);

  if (p <= P_THRESH) {
    // stage 2: if fit_parametric_curve true, draw stage 2 null statistics and get the SN p-value
    if (fit_parametric_curve) {
      // compute the stage 2 vector of null statistics
      null_statistics = compute_null_resid_statistics(resids, B1, B2, n_trt, synthetic_idxs);

      // compute the skew-normal p-value (set to -1 if fit bad)
      fit_sn_out = fit_and_evaluate_skew_normal(z_orig, null_statistics, side_code);
      p = fit_sn_out[3];
      sn_fit_used = p > -0.5;
      if (sn_fit_used) for (int i = 0; i < 3; i ++) sn_params[i] = fit_sn_out[i];
      stage = 2;
    }

    // stage 3: if skew normal fit failed, or if fit_parametric_curve false, draw stage 3 statistics and compute empirical p-value
    if (!fit_parametric_curve || !sn_fit_used) {
      if (B3 > 0) null_statistics = compute_null_resid_statistics(resids, B1 + B2, B3, n_trt, synthetic_idxs);
      p = compute_empirical_p_value(null_statistics, z_orig, side_code);
      stage = 3;
    }
  }

  // construct output
  if (return_resampling_dist) {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("lfc") = lfc, Named("stage") = stage, Named("sn_params") = sn_params, Named("resampling_dist") = null_statistics);
  } else {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("lfc") = lfc, Named("stage") = stage, Named("sn_params") = sn_params);
  }
  return out;
}
