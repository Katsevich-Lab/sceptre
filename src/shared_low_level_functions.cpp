#include <Rcpp.h>
using namespace Rcpp;

// output: a vector y whose length is equal to the number of cells, containing a true/false indicating presence/absence of a gene expression
void load_nonzero_posits(IntegerVector j, IntegerVector p, int row_idx, std::vector<bool>& y_orig,
                         std::vector<bool>& y_sub, std::vector<int>& cells_in_use_zero_idx) {
  int start = p[row_idx];
  int end = p[row_idx + 1];
  for (int k = 0; k < y_orig.size(); k ++) y_orig[k] = false;
  for (int k = start; k < end; k ++) y_orig[j[k]] = true;
  for (int k = 0; k < y_sub.size(); k++) y_sub[k] = y_orig[cells_in_use_zero_idx[k]];
  return;
}
