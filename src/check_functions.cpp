#include <Rcpp.h>
using namespace Rcpp;

/*************************************************************************************
 * The functions below are helpful for checking the output of the other C++ functions.
 *************************************************************************************/

// [[Rcpp::export]]
IntegerVector synth_idx_list_to_matrix(SEXP synthetic_idx_ptr) {
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idx_ptr);
  // create long vector containing all entries
  int B = (*synth_idx_list).size();
  int M = (*synth_idx_list)[0].size();
  int counter = 0;
  int n_entries = B * M;
  IntegerVector v(n_entries);
  for (int i = 0; i < B; i ++) {
    for (int j = 0; j < M; j++) {
      v[counter] = (*synth_idx_list)[i][j];
      counter ++;
    }
  }
  v.attr("dim") = Dimension(M, B);
  return v;
}


// [[Rcpp::export]]
// List synth_idx_list_to_r_list(SEXP synthetic_idx_ptr) {

// }
