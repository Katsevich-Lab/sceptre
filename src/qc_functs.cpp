#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector compute_n_nonzero_trt_vector(NumericVector expression_vector, List grna_group_idxs, IntegerVector grna_group_posits) {
  IntegerVector out(grna_group_posits.size());
  int counter = 0;
  IntegerVector curr_idxs;
  for (int i = 0; i < grna_group_posits.size(); i ++) {
    counter = 0;
    curr_idxs = grna_group_idxs[grna_group_posits[i] - 1];
    for (int j = 0; j < curr_idxs.size(); j ++) {
      if (expression_vector[curr_idxs[j] - 1] > 0.5) counter ++;
    }
    out[i] = counter;
  }
  return(out);
}
