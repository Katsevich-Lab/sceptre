// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/policies/error_handling.hpp>
#include <stdexcept>
#include <typeinfo>
std::string recorded_type, recorded_message;
std::string exception_class(const std::exception& e) {
 if (dynamic_cast<const std::domain_error*>(&e)) return "std::domain_error";
 if (dynamic_cast<const std::overflow_error*>(&e)) return "std::overflow_error";
 if (dynamic_cast<const boost::math::evaluation_error*>(&e)) return "boost::math::evaluation_error";
 if (dynamic_cast<const boost::math::rounding_error*>(&e)) return "boost::math::rounding_error";
 return typeid(e).name();
}
#define estimate_theta instrumented_estimate_theta
#include "instrumented.inc"
#undef estimate_theta
// [[Rcpp::export]]
Rcpp::List probe_estimate(Rcpp::NumericVector y, Rcpp::NumericVector mu, double dfr, int limit=50, double eps=0.0001220703125) {
 recorded_type="none";recorded_message="";
 auto result=instrumented_estimate_theta(y,mu,dfr,limit,eps);
 return Rcpp::List::create(Rcpp::_["result"]=result,Rcpp::_["type"]=recorded_type,Rcpp::_["message"]=recorded_message);
}
// [[Rcpp::export]]
Rcpp::List probe_special(double x, bool trigamma) {
 try {
  double result=trigamma ? boost::math::trigamma(x) : boost::math::digamma(x);
  return Rcpp::List::create(Rcpp::_["result"]=result,Rcpp::_["type"]="none");
 } catch(const std::exception& e) {
  return Rcpp::List::create(Rcpp::_["type"]=exception_class(e),Rcpp::_["message"]=e.what());
 }
}
