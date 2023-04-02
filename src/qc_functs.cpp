#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void load_nonzero_posits(IntegerVector j, IntegerVector p, int column_idx, std::vector<bool>& y) {
  int start = p[column_idx];
  int end = p[column_idx + 1];
  for (int k = 0; k < y.size(); k ++) y[k] = false;
  for (int k = start; k < end; k ++) y[j[k]] = true;
  return;
}


// [[Rcpp::export]]
List compute_nt_nonzero_matrix_and_n_ok_pairs(IntegerVector j, IntegerVector p, int n_cells, List grna_group_idxs, List indiv_nt_grna_idxs, IntegerVector all_nt_idxs, IntegerVector to_analyze_response_idxs, IntegerVector to_analyze_grna_idxs, int n_nonzero_trt_thresh, int n_nonzero_cntrl_thresh, bool compute_n_ok_pairs) {
  // 0. initialize variables and objects
  std::vector<bool> y(n_cells);
  int n_nt_grnas = indiv_nt_grna_idxs.size(), n_genes = p.size() - 1, n_ok_pairs = 0, pair_pointer = 0, n_nonzero = 0;
  IntegerVector curr_idxs;
  bool n_cntrl_cells_ok;
  IntegerMatrix M(n_nt_grnas, n_genes);
  
  // 1. iterate over genes
  for (int column_idx = 0; column_idx < n_genes; column_idx++) {
    load_nonzero_posits(j, p, column_idx, y);
    // 1.1 iterate over nt grnas, adding n nonzero trt to M matrix
    for (int grna_idx = 0; grna_idx < n_nt_grnas; grna_idx ++) {
      n_nonzero = 0;
      curr_idxs = indiv_nt_grna_idxs[grna_idx];
      for (int k = 0; k < curr_idxs.size(); k ++) {
        if (y[all_nt_idxs[curr_idxs[k] - 1] - 1]) n_nonzero ++;
      }
      M(grna_idx, column_idx) = n_nonzero;
    }
    
    // 1.2 if compute n_ok_pairs
    if (compute_n_ok_pairs) {
      // determine if n control cells exceeds threshold
      n_nonzero = 0;
      for (int k = 0; k < n_nt_grnas; k ++) n_nonzero += M(k, column_idx);
      n_cntrl_cells_ok = n_nonzero >= n_nonzero_cntrl_thresh;
      while (to_analyze_response_idxs[pair_pointer] - 1 == column_idx) {
        if (n_cntrl_cells_ok) {
          curr_idxs = grna_group_idxs[to_analyze_grna_idxs[pair_pointer] - 1];
          n_nonzero = 0;
          for (int k = 0; k < curr_idxs.size(); k ++) if (y[curr_idxs[k] - 1]) n_nonzero ++;
          if (n_nonzero >= n_nonzero_trt_thresh) n_ok_pairs ++;  
        }
        pair_pointer ++;
      }
    }
  }
  return List::create(Named("n_nonzero_mat") = M, Named("n_ok_pairs") = n_ok_pairs);
}


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
