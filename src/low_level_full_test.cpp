#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>
#include "shared_low_level_functions.h"


double compute_observed_full_statistic(NumericVector a, NumericVector w, NumericMatrix D, int n_cntrl) {
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), D_ncol = D.ncol();
  
  // iterate over the rows of D
  for (int i = 0; i < D_nrow; i ++) {
    inner_sum = 0;
    for (int j = n_cntrl; j < D_ncol; j ++) {
      inner_sum += D(i, j);
    }
    lower_right += inner_sum * inner_sum;
  }
  
  // second, compute the lower-left hand of the denominator; also, compute the top
  for (int j = n_cntrl; j < D_ncol; j ++) {
    top += a[j];
    lower_left += w[j];
  }
  
  // finally, compute the z-score
  return(top/sqrt(lower_left - lower_right));
}


std::vector<double> compute_null_full_statistics(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, int start_pos, int B, int n_trt, SEXP synthetic_idxs) {
  // dereference the synthetic treatment idxs
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idxs);
  
  // initialize variables
  int* curr_vect;
  std::vector<double> out(B);
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), D_ncol = D.ncol(), idx;
  
  // iterate over the B tests
  for (int k = start_pos; k < B + start_pos; k ++) {
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
List run_low_level_test_full_v2(NumericVector y,
                                NumericVector mu,
                                NumericVector a,
                                NumericVector w,
                                NumericMatrix D,
                                int n_cntrl,
                                int n_trt,
                                SEXP synthetic_idxs,
                                int B1,
                                int B2,
                                int B3,
                                bool fit_skew_normal,
                                bool return_resampling_dist) {
  double P_THRESH = 0.02, p;
  bool sn_fit_used = false;
  int round = 1;
  List out;
  
  // estimate the log fold change
  double lfc = estimate_log_fold_change(y, mu, n_cntrl, n_trt);
  
  // compute the original statistic
  double z_orig = compute_observed_full_statistic(a, w, D, n_cntrl);
  
  // compute the round 1 vector of null statistics and p-value
  std::vector<double> null_statistics = compute_null_full_statistics(a, w, D, 0, B1, n_trt, synthetic_idxs);
  p = compute_empirical_p_value(null_statistics, z_orig, 0);
  
  if ((p <= P_THRESH) & !return_resampling_dist) {
    // round 2: if fit_skew_normal true, draw round 2 null statistics and get the SN p-value
    if (fit_skew_normal) {
      // compute the round 2 vector of null statistics
      null_statistics = compute_null_full_statistics(a, w, D, B1, B2, n_trt, synthetic_idxs);
      
      // compute the skew-normal p-value (set to -1 if fit bad)
      p = fit_and_evaluate_skew_normal(z_orig, null_statistics);
      sn_fit_used = p > -0.5;
      round = 2;
    }
    
    // round 3: if skew normal fit failed, or if fit_skew_normal false, draw round 3 statistics and compute empirical p-value
    if (!fit_skew_normal || !sn_fit_used) {
      null_statistics = compute_null_full_statistics(a, w, D, B1 + B2, B3, n_trt, synthetic_idxs);
      p = compute_empirical_p_value(null_statistics, z_orig, 0);
      round = 3;
    }
  }
  
  // construct output
  if (return_resampling_dist) {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("lfc") = lfc, Named("sn_fit_used") = sn_fit_used, Named("round") = round, Named("resampling_dist") = null_statistics);
  } else {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("lfc") = lfc, Named("sn_fit_used") = sn_fit_used, Named("round") = round);
  }
  
  return(out);
}