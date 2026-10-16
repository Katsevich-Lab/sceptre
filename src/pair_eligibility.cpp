#include <Rcpp.h>
#include "pair_eligibility.h"

// [[Rcpp::export]]
Rcpp::LogicalVector evaluate_pair_eligibility_cpp(
    Rcpp::IntegerVector treatment, Rcpp::IntegerVector control,
    int treatment_threshold, int control_threshold) {
  if (treatment.size() != control.size()) Rcpp::stop("QC counts must align.");
  Rcpp::LogicalVector result(treatment.size());
  for (R_xlen_t i = 0; i < treatment.size(); ++i) {
    // Preserve R's three-valued vector predicate for missing measurements.
    const bool treatment_missing = treatment[i] == NA_INTEGER;
    const bool control_missing = control[i] == NA_INTEGER;
    if ((!treatment_missing && treatment[i] < treatment_threshold) ||
        (!control_missing && control[i] < control_threshold)) {
      result[i] = false;
    } else if (treatment_missing || control_missing) {
      result[i] = NA_LOGICAL;
    } else {
      result[i] = legacy_pair_eligible(treatment[i], control[i],
                                      treatment_threshold, control_threshold);
    }
  }
  return result;
}
