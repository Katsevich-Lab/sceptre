// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <math.h>
#include <stdexcept>
using namespace Rcpp;

// Estimate negative-binomial dispersion theta, where Var(Y | mu) = mu + mu^2/theta.
// R supplies counts y and fixed Poisson fitted means mu in the same cell order.
// We estimate dispersion here, rather than refitting regression coefficients.
// R later restricts the returned estimate to its allowed range.

// Evaluate the negative-binomial dispersion score at trial theta th.
// y and mu must describe the same cells in the same order. Use the current Boost
// special-function calculations; invalid domains or overflow can throw numeric
// exceptions for the estimator below to handle.
double nb_score(double th, NumericVector mu, NumericVector y) {
  double sum = 0;
  for (int i = 0; i < mu.size(); i ++) {
    sum += boost::math::digamma(th + y[i]) - boost::math::digamma(th) + std::log(th) + 1 - std::log(th + mu[i]) - (y[i] + th)/(mu[i] + th);
  }
  return(sum);
}


// Evaluate observed information (minus the score derivative) at th, using the
// same cells as nb_score(). Preserve sum and special-function evaluation order
// to keep the fitted numerical results repeatable.
double nb_info(double th, NumericVector mu, NumericVector y) {
 double sum = 0;
  for (int i = 0; i < mu.size(); i ++) {
    sum += boost::math::trigamma(th) - boost::math::trigamma(th + y[i]) - 1/th + 2/(mu[i] + th) - (y[i] + th)/((mu[i] + th) * (mu[i] + th));
  }
  return(sum);
}


// Compute the starting dispersion n / sum((y/mu - 1)^2). Keep zero/non-finite
// arithmetic as is. This initial estimate is calculated before the try block
// that handles numerical failures in the later iterative estimators.
double nb_theta_pilot_est(NumericVector y, NumericVector mu) {
  double n = (double) y.size();
  double denom = Rcpp::sum((y/mu - 1) * (y/mu - 1));
  return(n/denom);
}


// Estimate dispersion by matching moments, starting from the supplied estimate
// and residual degrees of freedom. Return the last theta and a warning flag.
// Flag negative/non-finite theta or reaching the iteration limit. The caller
// decides whether to fall back to the starting estimate.
std::vector<double> nb_theta_mm(double t0, NumericVector y, NumericVector mu, double dfr, int limit, double eps) {
  int it = 0;
  double del = 1;
  while (++it < limit && fabs(del) > eps) {
    t0 = fabs(t0);
    // Match the Pearson residual sum to the supplied residual degrees of
    // freedom. Keep the existing step denominator and absolute-value restart;
    // neither a warning nor an exceptional step selects a new estimator here.
    del = (Rcpp::sum(Rcpp::pow(y - mu, 2)/(mu + Rcpp::pow(mu, 2)/t0)) - dfr)/Rcpp::sum(Rcpp::pow(y - mu, 2)/Rcpp::pow(mu + t0, 2));
    t0 -= del;
  }

  double warning = 0.0;
  if (t0 < 0 || it == limit || !std::isfinite(t0)) warning = 1.0;
  std::vector<double> out { t0, warning };
  return(out);
}


// Estimate dispersion by Newton updates of the score/information ratio.
// Return theta and a warning flag, with the same warning rules as moment fitting.
// Numeric exceptions pass to the main estimator below; inputs stay unchanged.
std::vector<double> nb_theta_mle(double t0, NumericVector y, NumericVector mu, int limit, double eps) {
  int it = 0;
  double del = 1;

  while (++it < limit && fabs(del) > eps) {
    t0 = fabs(t0);
    // Observed information is minus the score derivative, hence the plus
    // Newton step. Do not reorder its cell sums or change the stopping test.
    del = nb_score(t0, mu, y)/nb_info(t0, mu, y);
    t0 += del;
  }

  double warning = 0.0;
  if (t0 < 0 || it == limit || !std::isfinite(t0)) warning = 1.0;
  std::vector<double> out { t0, warning };
  return(out);
}


// Estimate dispersion from counts and fixed means for the same cells, with the
// supplied degrees of freedom, iteration limit and tolerance. Return the value
// and method code: 1 = maximum likelihood, 2 = moments, 3 = starting estimate.
// The code records the last assigned method, not a promise of convergence.
// Keep the assignment order within the try block: if moment fitting throws,
// the preceding likelihood result and method 1 remain, even if it warned.
// Catch only the expected Boost domain/overflow errors. Memory and system errors
// must propagate. Starting-estimate calculation and R-result allocation happen
// outside this numerical error handler.
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
    // A numerical warning selects MM from the original pilot (not the last
    // MLE iterate). A second warning restores that pilot and method 3.
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
  } catch (const std::domain_error&) {
  } catch (const std::overflow_error&) {}
  return(List::create(estimate, method));
}
