#include <Rcpp.h>
#include <math.h>
#include "shared_low_level_functions.h"
using namespace Rcpp;


// [[Rcpp::export]]
double compute_observed_distilled_statistic(NumericVector a, NumericVector b, int n_cntrl) {
  double top = 0, bottom = 0, out;
  for (int i = n_cntrl; i < a.size(); i++) {
      top += a[i];
      bottom += b[i];
    }
  return top/sqrt(bottom);
}


// [[Rcpp::export]]
std::vector<double> compute_null_distilled_statistics(const NumericVector& a, const NumericVector& b, int start_pos, int B, int n_trt, SEXP synthetic_idxs) {
  // dereference the synthetic treatment idxs
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idxs);
  
  // initialize variables
  std::vector<double> out(B);
  double top, bottom;
  int idx;
  int* curr_vect;

  // compute the null statistics
  for (int i = start_pos; i < B + start_pos; i++) {
    top = 0;
    bottom = 0;
    curr_vect = &(*synth_idx_list)[i][0];
    
    for (int j = 0; j < n_trt; j++) {
      idx = curr_vect[j];
      top += a[idx];
      bottom += b[idx];
    }
    out[i - start_pos] = top/sqrt(bottom);
  }
  
  return out;
}


// [[Rcpp::export]]
List run_low_level_test_distilled(NumericVector y,
                                  NumericVector mu,
                                  NumericVector a,
                                  NumericVector b,
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
  double z_orig = compute_observed_distilled_statistic(a, b, n_cntrl);
  
  // compute the round 1 vector of null statistics and p-value
  std::vector<double> null_statistics = compute_null_distilled_statistics(a, b, 0, B1, n_trt, synthetic_idxs);
  p = compute_empirical_p_value(null_statistics, z_orig, 0);
  
  if ((p <= P_THRESH) & !return_resampling_dist) {
    // round 2: if fit_skew_normal true, draw round 2 null statistics and get the SN p-value
    if (fit_skew_normal) {
      // compute the round 2 vector of null statistics
      null_statistics = compute_null_distilled_statistics(a, b, B1, B2, n_trt, synthetic_idxs);
      
      // compute the skew-normal p-value (set to -1 if fit bad)
      p = fit_and_evaluate_skew_normal(z_orig, null_statistics);
      sn_fit_used = p > -0.5;
      round = 2;
    }
    
    // round 3: if skew normal fit failed, or if fit_skew_normal false, draw round 3 statistics and compute empirical p-value
    if (!fit_skew_normal || !sn_fit_used) {
      null_statistics = compute_null_distilled_statistics(a, b, B1 + B2, B3, n_trt, synthetic_idxs);
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
