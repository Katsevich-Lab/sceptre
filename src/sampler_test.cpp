#include <Rcpp.h>
#include <random>
#include <boost/random.hpp>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector test_boost_sampler(int n) {
  NumericVector x(n);
  boost::random::mt19937 rng(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);
  for (int i = 0; i < n; i ++) {
    x[i] = distribution(rng);
  }
  return(x);
}


// [[Rcpp::export]]
NumericVector test_std_sampler(int n) {
  NumericVector x(n);
  std::mt19937 rng(4);
  std::uniform_real_distribution<double> distribution(0, 1);
  for (int i = 0; i < n; i ++) {
    x[i] = distribution(rng);
  }
  return(x);
}
