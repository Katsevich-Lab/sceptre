#include <Rcpp.h>
using namespace Rcpp;

#ifndef ESTIMATE_LOG_FOLD_CHANGE_V2
#define ESTIMATE_LOG_FOLD_CHANGE_V2
double estimate_log_fold_change_v2(NumericVector y, NumericVector mu, IntegerVector trt_idxs, int n_trt);
#endif

#ifndef COMPUTE_EMPIRICAL_P_VALUE
#define COMPUTE_EMPIRICAL_P_VALUE
double compute_empirical_p_value(const std::vector<double>& null_statistics, double z_orig, int side);
#endif


#ifndef FIT_SKEW_NORMAL_FUNCT
#define FIT_SKEW_NORMAL_FUNCT
std::vector<double> fit_skew_normal_funct(const std::vector<double>& y);
#endif


#ifndef CHECK_SN_TAIL
#define CHECK_SN_TAIL
bool check_sn_tail (const std::vector<double>& y, double xi_hat, double omega_hat, double alpha_hat);
#endif


#ifndef FIT_AND_EVALUATE_SKEW_NORMAL
#define FIT_AND_EVALUATE_SKEW_NORMAL
double fit_and_evaluate_skew_normal(double z_orig, std::vector<double>& null_statistics);
#endif
