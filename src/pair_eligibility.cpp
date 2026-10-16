#include <Rcpp.h>
#include "pair_eligibility.h"

// [[Rcpp::export]]
Rcpp::LogicalVector evaluate_pair_eligibility_cpp(
    Rcpp::IntegerVector treatment, Rcpp::IntegerVector control,
    int treatment_threshold, int control_threshold) {
  if (treatment.size() != control.size()) Rcpp::stop("QC counts must align.");
  Rcpp::LogicalVector result(treatment.size());
  for (R_xlen_t i = 0; i < treatment.size(); ++i) {
    result[i] = legacy_pair_eligibility_value(treatment[i], control[i],
        treatment_threshold, control_threshold, NA_INTEGER);
  }
  return result;
}
