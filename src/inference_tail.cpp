#include "inference.h"
#include <boost/math/distributions/skew_normal.hpp>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace Rcpp;
using boost::math::skew_normal;

namespace sceptre {
namespace inference {

// Compute empirical p-values with one added to both tail counts and draw count.
// side -1 selects the left tail, +1 the right, and other values twice the smaller
// tail, capped at 1. Ties count in both tails; empty draws use the same
// formula. Leave scores unchanged.
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
    p = std::min(1.0, 2 * std::min(compute_empirical_p_value(null_statistics, z_orig, -1),
                                  compute_empirical_p_value(null_statistics, z_orig, 1)));
  }
  return p;
}


namespace {
// Fit skew-normal location, scale and shape from the selected null scores'
// moments; also return the sample mean and population standard deviation.
// Keep the calculation order and current skewness rule: values above .995 become
// .9*.995. Do not add a matching negative cap or repair degenerate moments.
// The calling fitter handles non-finite parameters by using empirical p-values.
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

  return std::vector<double> {xi, omega, alpha, m_y, sd_y};
}


// Check how well the skew-normal fit describes the upper tail of sorted scores.
// At probabilities .90, .905, ..., .99, compare the observed order statistic's
// fitted tail probability with its expected probability. Decline the fit if the
// ratio exceeds 2. The caller sorts (and may negate) the scores and checks that
// there are at least 100, so these order-statistic lookups are valid.
bool check_sn_tail(const std::vector<double>& y, double xi_hat, double omega_hat, double alpha_hat) {
  // The owning fitter checks n >= 100 before this internal .99 quantile lookup.
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


// Check for extreme scores in either tail using the current normal-envelope
// ratio cutoff of 1.5. Input must be nonempty and sorted in increasing order.
// Do not change the scores. Keep the existing divisions, including non-finite
// results in degenerate cases.
bool check_for_outliers (const std::vector<double>& null_statistics, double mu, double sd) {
  double min_z = null_statistics[0], max_z = null_statistics[null_statistics.size() - 1];
  double RATIO_THRESH = 1.5;
  double B = (double) null_statistics.size();
  double R_max = max_z/(mu + sd * sqrt(2 * log(B)));
  double R_min = min_z/(mu - sd * sqrt(2 * log(B)));
  bool ok = (R_max <= RATIO_THRESH && R_min <= RATIO_THRESH);
  return ok;
}

} // namespace

// Try a skew-normal p-value on a separate vector of at least 100 null scores.
// Return xi, omega, alpha and p. Use p = -1 to request empirical p-values when
// moments are non-finite, outliers are present or the tail fit is poor. Accepted
// p-values use the requested tail and have a minimum of 1e-250.
// This function sorts its score vector in place. A left-tail check reverses and
// negates it, then restores signs but not order; the returned scores may remain
// descending. Keep this observable ordering and leave the full random draws
// untouched. Error clearly if a reached fit has fewer than 100 scores.
std::vector<double> fit_and_evaluate_skew_normal(double z_orig, std::vector<double>& null_statistics, int side_code) {
  if (null_statistics.size() < 100)
    stop("Reached skew-normal fitting requires at least 100 null statistics.");
  // 0. define variables
  double p = -1.0;
  bool finite_params = true;

  // 1. fit the skew normal
  std::vector<double> fitted_params = fit_skew_normal_funct(null_statistics);

  // 2. verify that the fitted parameters are finite
  for (int i = 0; i < fitted_params.size(); i ++) {
    if (!std::isfinite(fitted_params[i])) finite_params = false;
  }

  if (finite_params) {
    // 3. sort the vector of null statistics
    sort(null_statistics.begin(), null_statistics.end(), std::less<double>());

    // 4. compute the median of the null statistics
    int median_idx = (null_statistics.size() - 1)/2;
    double median = null_statistics[median_idx];

    // 5. determine the tail to check; if z_orig >= median, right; else, left
    bool check_right_tail = (z_orig >= median);
    bool outlier_ok, fit_ok = false, use_sn;

    // 6. check for outliers in both tails
    outlier_ok = check_for_outliers(null_statistics, fitted_params[3], fitted_params[4]);

    // 7. check for goodness of fit in the appropriate tail
    if (outlier_ok) {
      if (check_right_tail) { // right tail check
        fit_ok = check_sn_tail(null_statistics, fitted_params[0], fitted_params[1], fitted_params[2]);
      } else { // left tail check
        std::reverse(null_statistics.begin(), null_statistics.end());
        for (int i = 0; i < null_statistics.size(); i ++) null_statistics[i] *= -1.0; // flip sign
        fit_ok = check_sn_tail(null_statistics, -fitted_params[0], fitted_params[1], -fitted_params[2]);
        for (int i = 0; i < null_statistics.size(); i ++) null_statistics[i] *= -1.0; // flip sign back
      }
    }

    // 8. compute the SN p-value (using the appropriate tail) if there are no outliers and the fit is OK
    use_sn = outlier_ok && fit_ok;
    if (use_sn) {
      skew_normal dist(fitted_params[0], fitted_params[1], fitted_params[2]);
      if (side_code == 0) { // two-tailed
        p = 2.0 * (check_right_tail ? cdf(complement(dist, z_orig)) : cdf(dist, z_orig));
      } else if (side_code == 1) { // right-tailed
        p = cdf(complement(dist, z_orig));
      } else { // left-tailed
        p = cdf(dist, z_orig);
      }
      if (p <= 1.0e-250) p = 1.0e-250;
    }
  }
  // 9. return p, xi, omega, alpha
  return std::vector<double> {fitted_params[0], fitted_params[1], fitted_params[2], p};
}

} // namespace inference
} // namespace sceptre
