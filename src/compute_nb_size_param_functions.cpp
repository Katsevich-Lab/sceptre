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


double nb_theta_mm(double t0, NumericVector y, NumericVector mu, double dfr, int limit, double eps) {
  int it = 0;
  double del = 1;
  while (++it < limit && fabs(del) > eps) {
    t0 = fabs(t0);
    del = (Rcpp::sum(Rcpp::pow(y - mu, 2)/(mu + Rcpp::pow(mu, 2)/t0)) - dfr)/Rcpp::sum(Rcpp::pow(y - mu, 2)/Rcpp::pow(mu + t0, 2));
    t0 -= del;
  }

  return(t0);
}


double nb_theta_mle(double t0, NumericVector y, NumericVector mu, int limit, double eps, int* warning) {
  int it = 0;
  double del = 1;

  while (++it < limit && fabs(del) > eps) {
    t0 = fabs(t0);
    del = nb_score(t0, mu, y)/nb_info(t0, mu, y);
    t0 += del;
  }

  if (t0 < 0 || it == limit) *warning = 1;
  return(t0);
}


//' Estimate theta
//'
//' This function estimates the negative binomial size parameter theta using the fitted means of a Poisson GLM.
//'
//' @param y a vector of expressions
//' @param mu a vector of fitted means from the Poisson regression
//' @param dfr the residual degrees of freedom of the Poisson regression
//' @param limit iteration limit
//' @param eps convergence threshold
//' @return the estimated theta
//' @noRd
// [[Rcpp::export]]
List estimate_theta(NumericVector y, NumericVector mu, double dfr, int limit, double eps) {
  // first, attempt to estimate theta via MLE
  double t0 = nb_theta_pilot_est(y, mu);
  double estimate = t0;
  int warning = 0;
  int method = 3;
  try {
    estimate = nb_theta_mle(t0, y, mu, limit, eps, &warning);
    method = 1;
    // if there is a warning, swtich to MM
    if (warning) {
      method = 2;
      estimate = nb_theta_mm(t0, y, mu, dfr, limit, eps);
    }
    // return both the estimate and method indicator
  } catch (...) {}
  return(List::create(estimate, method));
}
