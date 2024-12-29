#include <Rcpp.h>
using namespace Rcpp;
#include "shared_low_level_functions.h"

// this function outputs four pieces:
// 1. N nonzero mat; this is the number of nonzero cells for each gene-NT gRNA pair (same regardless of control group)
// 2. n_nonzero_tot; this is the vector giving the number of nonzero cells per gene; when using the NT cells as the control group, we restrict our attention to the NT cells; when using the complement set as the control group, by contrast, we sum over all cells.
// 3,4. n nonzero trt, n nonzero cntrl vectors
// [[Rcpp::export]]
List compute_nt_nonzero_matrix_and_n_ok_pairs_v3(
    IntegerVector j,
    IntegerVector p,
    int n_cells_orig,
    int n_cells_sub,
    List grna_group_idxs,
    List indiv_nt_grna_idxs,
    IntegerVector all_nt_idxs,
    IntegerVector to_analyze_response_idxs,
    IntegerVector to_analyze_grna_idxs,
    bool control_group_complement,
    bool treatment_group_inclusive,
    IntegerVector cells_in_use
) {
  // 0. initialize variables and objects
  int n_nt_grnas = indiv_nt_grna_idxs.size();
  int n_genes    = p.size() - 1;
  int n_pairs    = to_analyze_response_idxs.size();
  int pair_pointer = 0, n_nonzero = 0, curr_n_nonzero_cntrl = 0, curr_n_nonzero_trt = 0;

  IntegerVector curr_idxs, n_nonzero_tot(n_genes), n_nonzero_trt(n_pairs), n_nonzero_cntrl(n_pairs);

  // M will remain all zero if treatment_group_inclusive == true && !control_group_complement;
  // the calibration check is currently not supported for this combination of parameters
  IntegerMatrix M(n_nt_grnas, n_genes);

  std::vector<bool> y_sub(n_cells_sub), y_orig(n_cells_orig);
  std::vector<int> cells_in_use_zero_idx(n_cells_sub);
  for (int i = 0; i < cells_in_use_zero_idx.size(); i++) {
    cells_in_use_zero_idx[i] = cells_in_use[i] - 1;
  }

  // 1. iterate over genes
  for (int row_idx = 0; row_idx < n_genes; row_idx++) {
    // load nonzero positions into y_sub
    load_nonzero_posits(j, p, row_idx, y_orig, y_sub, cells_in_use_zero_idx);

    // -------------------------------------------------------------
    // 1.2 Optionally fill M, unless (treatment_group_inclusive && !control_group_complement)
    // -------------------------------------------------------------
    if (!(treatment_group_inclusive && !control_group_complement)) {
      // Normal path: compute # of nonzero cells in each NT gRNA
      for (int grna_idx = 0; grna_idx < n_nt_grnas; grna_idx++) {
        n_nonzero = 0;
        curr_idxs = indiv_nt_grna_idxs[grna_idx];

        if (!control_group_complement) {
          // NT cells control: redirect through all_nt_idxs
          for (int k = 0; k < curr_idxs.size(); k++) {
            int cell0 = all_nt_idxs[curr_idxs[k] - 1] - 1;
            if (y_sub[cell0]) {
              n_nonzero++;
            }
          }
        } else {
          // complement set: just use curr_idxs directly
          for (int k = 0; k < curr_idxs.size(); k++) {
            int cell0 = curr_idxs[k] - 1;
            if (y_sub[cell0]) {
              n_nonzero++;
            }
          }
        }
        M(grna_idx, row_idx) = n_nonzero;
      }
    }
    // else: do nothing, M remains zero-initialized for this row

    // -------------------------------------------------------------
    // 1.3 update n_nonzero_tot for this gene
    // -------------------------------------------------------------
    n_nonzero_tot[row_idx] = 0;
    if (!control_group_complement) {
      // directly count how many of the all_nt_idxs cells are nonzero
      for (int idx = 0; idx < all_nt_idxs.size(); idx++) {
        int cell0 = all_nt_idxs[idx] - 1;
        if (y_sub[cell0]) {
          n_nonzero_tot[row_idx]++;
        }
      }
    } else {
      // complement scenario: sum over all y_sub
      for (int k = 0; k < (int)y_sub.size(); k++) {
        if (y_sub[k]) {
          n_nonzero_tot[row_idx]++;
        }
      }
    }

    // 1.4 iterate through the discovery pairs containing this gene
    while (pair_pointer < n_pairs &&
           (to_analyze_response_idxs[pair_pointer] - 1 == row_idx))
    {
      curr_idxs = grna_group_idxs[to_analyze_grna_idxs[pair_pointer] - 1];
      curr_n_nonzero_trt = 0;

      // compute treatment
      for (int k = 0; k < curr_idxs.size(); k++) {
        int cell0 = curr_idxs[k] - 1;
        if (y_sub[cell0]) {
          curr_n_nonzero_trt++;
        }
      }

      // compute control
      if (!control_group_complement) {
        curr_n_nonzero_cntrl = n_nonzero_tot[row_idx];
      } else {
        // complement set
        curr_n_nonzero_cntrl = n_nonzero_tot[row_idx] - curr_n_nonzero_trt;
      }

      n_nonzero_trt[pair_pointer]   = curr_n_nonzero_trt;
      n_nonzero_cntrl[pair_pointer] = curr_n_nonzero_cntrl;
      pair_pointer++;
    }
  }

  return List::create(
    Named("n_nonzero_mat")  = M,
    Named("n_nonzero_tot")  = n_nonzero_tot,
    Named("n_nonzero_trt")  = n_nonzero_trt,
    Named("n_nonzero_cntrl")= n_nonzero_cntrl
  );
}
