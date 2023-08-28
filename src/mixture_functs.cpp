#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>
#include <algorithm>


// [[Rcpp::export]]
double compute_tolerance_cpp(double curr_log_lik, double prev_log_lik) {
  double tol;
  if (curr_log_lik == -std::numeric_limits<double>::infinity() || prev_log_lik == -std::numeric_limits<double>::infinity()) {
    tol = 1.0;
  } else {
    tol = std::abs(curr_log_lik - prev_log_lik)/std::min(std::abs(curr_log_lik), std::abs(prev_log_lik));
  }
  return tol;
}


// [[Rcpp::export]]
List run_reduced_em_algo_cpp(NumericVector pi_guesses, NumericVector g_pert_guesses, NumericVector g,
                             NumericVector g_mus_pert0, NumericVector log_g_factorial) {
  // define global options
  int B = pi_guesses.size(), n = g.size(), outer_i = 0;
  double ep_tol = 0.5 * 1e-4, n_doub = (double) g.size(), outer_log_lik = -std::numeric_limits<double>::infinity();
  int max_it = 50, min_it = 3;

  // define "outer" Ti1s
  NumericVector outer_Ti1s = NumericVector(n);
  bool outer_converged = false;

  // define "inner" variables
  bool converged, all_zero, any_na;
  double prev_log_lik, curr_log_lik, curr_g_pert, curr_pi, tol, e1, e2;
  int iteration;
  NumericVector g_mus_pert1 = NumericVector(n), p0 = NumericVector(n), p1 = NumericVector(n),
    s = NumericVector(n), quotient = NumericVector(n), Ti1s = NumericVector(n);

  for (int i = 0; i < B; i++) {
    converged = false;
    prev_log_lik = -std::numeric_limits<double>::infinity();
    curr_g_pert = g_pert_guesses[i];
    curr_pi = pi_guesses[i];
    iteration = 1;

    while(true) {
      g_mus_pert1 = g_mus_pert0 + curr_g_pert;

      // compute log likelihood
      p0 = exp(log(1 - curr_pi) + g * log(g_mus_pert0) - g_mus_pert0 - log_g_factorial);
      p1 = exp(log(curr_pi) + g * log(g_mus_pert1) - g_mus_pert1 - log_g_factorial);
      s = p0 + p1;
      for (int j = 0; j < s.size(); j++) if (s[j] < 1e-100) s[j] = 1e-100;
      curr_log_lik = Rcpp::sum(log(s));

      // compute membership probabilities
      quotient = log(1 - curr_pi) - log(curr_pi) + g * (log(g_mus_pert0) - log(g_mus_pert1)) + g_mus_pert1 - g_mus_pert0;
      Ti1s = 1/(exp(quotient) + 1);

      // verify that Ti1s are not NA or uniformly 0
      all_zero = true;
      any_na = false;
      for (int j = 0; j < Ti1s.size(); j ++) {
        if (Ti1s[j] > 1e-100) all_zero = false;
        if (!std::isfinite(Ti1s[j])) any_na = true;
      }
      if (all_zero || any_na) {
        curr_log_lik = -std::numeric_limits<double>::infinity();
        break;
      }

      // estimate new_pi and ensure new_pi is less than 0.5
      curr_pi = Rcpp::sum(Ti1s)/n_doub;
      if (curr_pi > 0.5) {
        Ti1s = 1 - Ti1s;
        curr_pi = 1 - curr_pi;
      }

      // check for convergence
      tol = compute_tolerance_cpp(curr_log_lik, prev_log_lik);
      if (tol < ep_tol && iteration >= min_it) {
        converged = true;
        break;
      } else {
        prev_log_lik = curr_log_lik;
        iteration ++;
        if (iteration >= max_it) break;
      }

      // M step
      e1 = Rcpp::sum(Ti1s * g);
      e2 = Rcpp::sum(Ti1s * g_mus_pert0);
      curr_g_pert = log(e1) - log(e2);
    }
    if (curr_log_lik > outer_log_lik && converged) {
      for (int j = 0; j < n; j ++) outer_Ti1s[j] = Ti1s[j];
      outer_i = i;
      outer_log_lik = curr_log_lik;
      outer_converged = true;
    }
  }

  return List::create(Named("outer_Ti1s") = outer_Ti1s, Named("outer_i") = outer_i,
                      Named("outer_converged") = outer_converged, Named("outer_log_lik") = outer_log_lik);
}
