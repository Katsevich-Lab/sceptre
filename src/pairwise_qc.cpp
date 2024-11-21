#include <Rcpp.h>
using namespace Rcpp;
#include "shared_low_level_functions.h"

// this function outputs four pieces:
// 1. N nonzero mat; this is the number of nonzero cells for each gene-NT gRNA pair (same regardless of control group)
// 2. n_nonzero_tot; this is the vector giving the number of nonzero cells per gene; when using the NT cells as the control group, we restrict our attention to the NT cells; when using the complement set as the control group, by contrast, we sum over all cells.
// 3,4. n nonzero trt, n nonzero cntrl vectors
// [[Rcpp::export]]
List compute_nt_nonzero_matrix_and_n_ok_pairs_v3(IntegerVector j, IntegerVector p, int n_cells_orig, int n_cells_sub,
                                                 List grna_group_idxs, List indiv_nt_grna_idxs, IntegerVector all_nt_idxs,
                                                 IntegerVector to_analyze_response_idxs, IntegerVector to_analyze_grna_idxs,
                                                 bool control_group_complement, IntegerVector cells_in_use) {
  // 0. initialize variables and objects
  int n_nt_grnas = indiv_nt_grna_idxs.size(), n_genes = p.size() - 1, n_pairs = to_analyze_response_idxs.size(), pair_pointer = 0, n_nonzero = 0, curr_n_nonzero_cntrl = 0, curr_n_nonzero_trt = 0;
  IntegerVector curr_idxs, n_nonzero_tot(n_genes), n_nonzero_trt(n_pairs), n_nonzero_cntrl(n_pairs);
  IntegerMatrix M(n_nt_grnas, n_genes);

  std::vector<bool> y_sub(n_cells_sub), y_orig(n_cells_orig);
  std::vector<int> cells_in_use_zero_idx(n_cells_sub);
  for (int i = 0; i < cells_in_use_zero_idx.size(); i ++) cells_in_use_zero_idx[i] = cells_in_use[i] - 1;

  // 1. iterate over genes
  for (int row_idx = 0; row_idx < n_genes; row_idx++) {
    // load nonzero positions into the boolean vector y_sub
    load_nonzero_posits(j, p, row_idx, y_orig, y_sub, cells_in_use_zero_idx);
    // 1.2 iterate over nt grnas, adding n nonzero trt to M matrix
    for (int grna_idx = 0; grna_idx < n_nt_grnas; grna_idx ++) {
      n_nonzero = 0;
      curr_idxs = indiv_nt_grna_idxs[grna_idx];
      for (int k = 0; k < curr_idxs.size(); k ++) {
        if (!control_group_complement) { // NT cells control group
          if (y_sub[all_nt_idxs[curr_idxs[k] - 1] - 1]) n_nonzero ++; // redirection
        } else { // complement set cells control group
          if (y_sub[curr_idxs[k] - 1]) n_nonzero ++; // no redirection
        }
      }
      M(grna_idx, row_idx) = n_nonzero;
    }

    // 1.2 update n_nonzero_tot for this gene
    n_nonzero_tot[row_idx] = 0;
    if (!control_group_complement) {
      for (int k = 0; k < n_nt_grnas; k ++) n_nonzero_tot[row_idx] += M(k, row_idx);
    } else {
      for (int k = 0; k < y_sub.size(); k++) if (y_sub[k]) n_nonzero_tot[row_idx] ++;
    }

    // iterate through the discovery pairs containing this gene
    while (pair_pointer < n_pairs && to_analyze_response_idxs[pair_pointer] - 1 == row_idx) {
      curr_idxs = grna_group_idxs[to_analyze_grna_idxs[pair_pointer] - 1];
      curr_n_nonzero_trt = 0;
      // first, get n nonzero trt
      for (int k = 0; k < curr_idxs.size(); k ++) if (y_sub[curr_idxs[k] - 1]) curr_n_nonzero_trt ++;
      // second, get n nonzero cntrl
      if (!control_group_complement) { // NT cells
        curr_n_nonzero_cntrl = n_nonzero_tot[row_idx];
      } else { // complement set
        curr_n_nonzero_cntrl = n_nonzero_tot[row_idx] - curr_n_nonzero_trt;
      }
      n_nonzero_trt[pair_pointer] = curr_n_nonzero_trt;
      n_nonzero_cntrl[pair_pointer] = curr_n_nonzero_cntrl;
      pair_pointer ++;
    }
  }

  return List::create(Named("n_nonzero_mat") = M, Named("n_nonzero_tot") = n_nonzero_tot,
                      Named("n_nonzero_trt") = n_nonzero_trt, Named("n_nonzero_cntrl") = n_nonzero_cntrl);
}
