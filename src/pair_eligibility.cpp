#include <Rcpp.h>
#include "pair_eligibility.h"

// Apply the shared QC rule to treatment/control count vectors in pair order.
// The caller supplies aligned integer counts and validated scalar thresholds.
// Keep repeats and R's missing-value rule: a known failure gives FALSE even if
// the other count is NA. Leave the supplied vectors unchanged.
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
