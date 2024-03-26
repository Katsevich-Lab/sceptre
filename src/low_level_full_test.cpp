#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>
#include "shared_low_level_functions.h"

// [[Rcpp::export]]
double compute_observed_full_statistic_v2(NumericVector a, NumericVector w, NumericMatrix D, IntegerVector trt_idxs) {
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


std::vector<double> compute_null_full_statistics(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, int start_pos, int B, int n_trt, bool use_all_cells, SEXP synthetic_idxs) {
  // dereference the synthetic treatment idxs
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idxs);

  // initialize variables
  int* curr_vect;
  std::vector<double> out(B);
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), idx;

  // iterate over the B tests
  for (int k = start_pos; k < B + start_pos; k ++) {
    if (use_all_cells) n_trt = (*synth_idx_list)[k].size();
    curr_vect = &(*synth_idx_list)[k][0];
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
                                bool asymptotic_normality,
                                int side_code) {
  // initialize the output list and p-value
  List out;
  double p, P_THRESH = 0.02;
  bool sn_fit_used = false;
  NumericVector sn_params = NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  std::vector<double> fit_sn_out;
  int stage = 1;

  // estimate the log fold change
  double lfc = estimate_log_fold_change_v2(y, mu, trt_idxs, n_trt);

  // compute the original statistic
  double z_orig = compute_observed_full_statistic_v2(a, w, D, trt_idxs);

  // if p-value based on asymptotic normality is requested, compute and return it
  if (asymptotic_normality) {
    p = compute_normal_p_value(z_orig, side_code);
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("lfc") = lfc, Named("stage") = 0, Named("sn_params") = sn_params);
    return(out);
  }

  // compute the stage 1 vector of null statistics and p-value
  std::vector<double> null_statistics = compute_null_full_statistics(a, w, D, 0, B1, n_trt, use_all_cells, synthetic_idxs);
  p = compute_empirical_p_value(null_statistics, z_orig, side_code);

  if (p <= P_THRESH) {
    // stage 2: if fit_parametric_curve true, draw stage 2 null statistics and get the SN p-value
    if (fit_parametric_curve) {
      // compute the stage 2 vector of null statistics
      null_statistics = compute_null_full_statistics(a, w, D, B1, B2, n_trt, use_all_cells, synthetic_idxs);

      // compute the skew-normal p-value (set to -1 if fit bad)
      fit_sn_out = fit_and_evaluate_skew_normal(z_orig, null_statistics, side_code);
      p = fit_sn_out[3];
      sn_fit_used = p > -0.5;
      if (sn_fit_used) for (int i = 0; i < 3; i ++) sn_params[i] = fit_sn_out[i];
      stage = 2;
    }

    // stage 3: if skew normal fit failed, or if fit_parametric_curve false, draw stage 3 statistics and compute empirical p-value
    if (!fit_parametric_curve || !sn_fit_used) {
      if (B3 > 0) null_statistics = compute_null_full_statistics(a, w, D, B1 + B2, B3, n_trt, use_all_cells, synthetic_idxs);
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

  return(out);
}
