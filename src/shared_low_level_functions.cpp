// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
using namespace Rcpp;
#include <boost/math/distributions/skew_normal.hpp>
using boost::math::skew_normal;
#include <cmath>
#include <math.h>
#include <algorithm>

double estimate_log_fold_change(NumericVector y, NumericVector mu, int n_cntrl, int n_trt) {
  double sum_y = 0, sum_mu = 0, n_trt_double = (double) n_trt;
  for (int i = n_cntrl; i < y.length(); i ++) {
    sum_y += y[i];
    sum_mu += mu[i];
  }
  return(log2(sum_y/n_trt_double) - log2(sum_mu/n_trt_double));
}


// side: -1 is left, 0 is both, 1 is right
// [[Rcpp::export]]
double compute_empirical_p_value(const std::vector<double>& null_statistics, double z_orig, int side) {
  double p;
  if (side == -1 || side == 1) { // left or right
    double counter = 0;
    double B = (double) null_statistics.size();
    if (side == -1) { // left
      for (int i = 0; i < null_statistics.size(); i ++) if (z_orig >= null_statistics[i]) counter ++;
    } else { // right
      for (int i = 0; i < null_statistics.size(); i ++) if (z_orig <= null_statistics[i]) counter ++;
    }
    p = (1.0 + counter)/(1.0 + B);
  } else { // two-sided
    p = 2 * std::min(compute_empirical_p_value(null_statistics, z_orig, -1),
                     compute_empirical_p_value(null_statistics, z_orig, 1));
  }
  return p;
}


// [[Rcpp::export]]
std::vector<double> fit_skew_normal_funct(const std::vector<double>& y) {
  // initialize variables
  int n = y.size();
  double MAX_GAMMA_1 = 0.995;
  double n_doub = (double) n, s_1 = 0, s_2 = 0, s_3 = 0;

  // compute mean and standard deviation
  for (int i = 0; i < n; i ++) {
    s_1 += y[i];
    s_2 += y[i] * y[i];
  }
  double m_y = s_1/n_doub;
  double sd_y = sqrt(s_2/n_doub - m_y * m_y);

  // compute gamma1
  for (int i = 0; i < n; i ++) {
    s_3 += pow(y[i] - m_y, 3);
  }
  double gamma1 = s_3/(n_doub * pow(sd_y, 3));
  if (gamma1 > MAX_GAMMA_1) gamma1 = 0.9 * MAX_GAMMA_1;

  // reparameterize solution
  double b = sqrt(2.0/M_PI);
  double r = std::copysign(1.0, gamma1) * pow(2 * std::abs(gamma1)/(4 - M_PI), 1.0/3.0);
  double delta = r/(b * sqrt(1 + r * r));
  double alpha = delta/sqrt(1 - delta * delta);
  double mu_z = b * delta;
  double sd_z = sqrt(1 - mu_z * mu_z);
  double omega = sd_y/sd_z;
  double xi = m_y - omega * mu_z;

  return std::vector<double> {xi, omega, alpha};
}


// [[Rcpp::export]]
bool check_sn_tail (const std::vector<double>& y, double xi_hat, double omega_hat, double alpha_hat) {
  // define variables
  double n = y.size(), ratio, quantile, p, sn_tail_prob;
  double RATIO_THRESH = 2.0;
  int idx;
  bool good_fit = true;

  // initialize the fitted skew normal distribution
  skew_normal dist(xi_hat, omega_hat, alpha_hat);

  // loop over the probabilities
  for (int i = 180; i < 199; i ++) {
    p = ((double) i)/200.0;
    idx = ceil(n * p);
    quantile = y[idx];
    sn_tail_prob = cdf(complement(dist, quantile));
    ratio = (1.0 - p)/sn_tail_prob;
    if (ratio > RATIO_THRESH) {
      good_fit = false;
      break;
    }
  }
  return good_fit;
}


// [[Rcpp::export]]
double fit_and_evaluate_skew_normal(double z_orig, std::vector<double>& null_statistics) {
  // 1. fit the skew normal
  std::vector<double> fitted_params = fit_skew_normal_funct(null_statistics);
  double p = -1.0;

  // 2. sort the vector of null statistics
  sort(null_statistics.begin(), null_statistics.end(), std::less<double>());

  // 3. compute the median of the null statistics
  int median_idx = (null_statistics.size() - 1)/2;
  double median = null_statistics[median_idx];

  // 4. determine the tail to check; if z_orig >= median, right; else, left
  bool check_right_tail = (z_orig >= median);
  bool tail_ok;

  // 5. check the appropriate tail and, if applicable, compute p-value
  if (check_right_tail) { // right tail check
    tail_ok = check_sn_tail(null_statistics, fitted_params[0], fitted_params[1], fitted_params[2]);
  } else { // left tail check
    std::reverse(null_statistics.begin(), null_statistics.end());
    for (int i = 0; i < null_statistics.size(); i ++) null_statistics[i] *= -1.0;
    tail_ok = check_sn_tail(null_statistics, -fitted_params[0], fitted_params[1], -fitted_params[2]);
  }

  // 6. if the tail is OK, compute the SN p-value (using the appropriate tail)
  if (tail_ok) {
    skew_normal dist(fitted_params[0], fitted_params[1], fitted_params[2]);
    p = 2.0 * (check_right_tail ? cdf(complement(dist, z_orig)) : cdf(dist, z_orig));
    if (p <= 1.0e-250) p = 1.0e-250;
  }

  // 7. return p
  return p;
}
