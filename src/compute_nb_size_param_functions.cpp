// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <math.h>
using namespace Rcpp;


double nb_score(double th, NumericVector mu, NumericVector y) {
  double sum = 0;
  for (int i = 0; i < mu.size(); i ++) {
    sum += boost::math::digamma(th + y[i]) - boost::math::digamma(th) + std::log(th) + 1 - std::log(th + mu[i]) - (y[i] + th)/(mu[i] + th);
  }
  return(sum);
}


double nb_info(double th, NumericVector mu, NumericVector y) {
 double sum = 0;
  for (int i = 0; i < mu.size(); i ++) {
    sum += boost::math::trigamma(th) - boost::math::trigamma(th + y[i]) - 1/th + 2/(mu[i] + th) - (y[i] + th)/((mu[i] + th) * (mu[i] + th));
  }
  return(sum);
}


double nb_theta_pilot_est(NumericVector y, NumericVector mu) {
  double n = (double) y.size();
  double denom = Rcpp::sum((y/mu - 1) * (y/mu - 1));
  return(n/denom);
}


std::vector<double> nb_theta_mm(double t0, NumericVector y, NumericVector mu, double dfr, int limit, double eps) {
  int it = 0;
  double del = 1;
  while (++it < limit && fabs(del) > eps) {
    t0 = fabs(t0);
    del = (Rcpp::sum(Rcpp::pow(y - mu, 2)/(mu + Rcpp::pow(mu, 2)/t0)) - dfr)/Rcpp::sum(Rcpp::pow(y - mu, 2)/Rcpp::pow(mu + t0, 2));
    t0 -= del;
  }

  double warning = 0.0;
  if (t0 < 0 || it == limit || !std::isfinite(t0)) warning = 1.0;
  std::vector<double> out { t0, warning };
  return(out);
}


std::vector<double> nb_theta_mle(double t0, NumericVector y, NumericVector mu, int limit, double eps) {
  int it = 0;
  double del = 1;

  while (++it < limit && fabs(del) > eps) {
    t0 = fabs(t0);
    del = nb_score(t0, mu, y)/nb_info(t0, mu, y);
    t0 += del;
  }

  double warning = 0.0;
  if (t0 < 0 || it == limit || !std::isfinite(t0)) warning = 1.0;
  std::vector<double> out { t0, warning };
  return(out);
}


// [[Rcpp::export]]
List estimate_theta(NumericVector y, NumericVector mu, double dfr, int limit, double eps) {
  // first, attempt to estimate theta via MLE
  double t0 = nb_theta_pilot_est(y, mu);
  double estimate = t0;
  int method = 3;
  try {
    std::vector<double> out = nb_theta_mle(t0, y, mu, limit, eps);
    estimate = out[0];
    method = 1;
    // if there is a warning, swtich to MM
    if (out[1] > 0.5) {
      out = nb_theta_mm(t0, y, mu, dfr, limit, eps);
      estimate = out[0];
      method = 2;
      if (out[1] > 0.5) {
        estimate = t0;
        method = 3;
      }
    }
    // return both the estimate and method indicator
  } catch (...) {}
  return(List::create(estimate, method));
}
