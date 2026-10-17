#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// Load one sparse row or column into the supplied dense output vector, clearing
// it first. The caller supplies valid Matrix column/row numbers and compressed
// pointers, all starting at 0, and enough output space for the full dimension.
// Load actual values, including stored zeros and fractions. Valid sparse rows/
// columns have unique positions; otherwise the last value wins. Only out changes.
void load_sparse_vector(IntegerVector l, IntegerVector p, NumericVector x, int idx, int dim, NumericVector out) {
  int start = p[idx], end = p[idx + 1];
  for (int k = 0; k < dim; k ++) out[k] = 0.0;
  for (int k = start; k < end; k ++) out[l[k]] = x[k];
}


// Read one complete gene row into a new dense vector in original cell order.
// row_idx starts at 1 in R and is converted once for load_sparse_vector().
// n_cells is the full matrix width; unstored entries become zero. R storage
// helpers select and reorder the requested cells afterward. Source slots stay
// unchanged.
// [[Rcpp::export]]
NumericVector load_csr_row(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, int n_cells) {
  NumericVector out(n_cells);
  load_sparse_vector(j, p, x, row_idx - 1, n_cells, out);
  return out;
}


// Build cumulative row/column pointers from sorted row/column numbers starting
// at 0. Numbers must be in [0, dim); the caller sorts and checks them first.
// Return dim+1 offsets, repeating offsets for empty rows/columns. With dim = 1,
// all entries belong to that single row/column. Repeated positions remain
// separate stored entries here; this helper does not combine them.
// [[Rcpp::export]]
IntegerVector obtain_pointer_vector(IntegerVector i, int dim) {
  int i_size = i.size();
  IntegerVector p(dim + 1);
  p[0] = 0;
  // special case of dim = 1
  if (dim == 1) {
    p[1] = i_size;
  } else {
    int curr_idx = 0, counter = 0;
    for (curr_idx = 0; curr_idx < dim; curr_idx ++) {
      while (counter < i_size && i[counter] == curr_idx) counter ++;
      p[curr_idx + 1] = counter;
    }
  }
  return p;
}


// Calculate per-cell summaries from a valid column-compressed feature-by-cell
// matrix. Matrix's i/p/x slots use positions starting at 0; n_genes and n_cells
// give its dimensions. Outputs follow original cell order.
// Here nonzero means a loaded value > 0.5; total UMIs sum actual values, including
// fractions. Mitochondrial rows arrive from R starting at 1 and are copied before
// conversion. The largest feature uses positions starting at 0 and takes the
// first maximum. All-zero columns keep feature 0 and a NaN fraction. Disabled
// optional summaries remain zero.
// Prepare all R data and named outputs first. Then use cached data pointers in
// a C++ calculation scope that ends before returning. Original slots and vectors
// must stay unchanged, including the mitochondrial-row vector.
// [[Rcpp::export]]
List compute_cell_covariates_cpp(IntegerVector i, IntegerVector p, NumericVector x, int n_genes,
                                 int n_cells, IntegerVector mt_gene_idxs, bool compute_p_mito, bool compute_max_feature) {
  IntegerVector n_nonzero(n_cells), max_feature(n_cells);
  NumericVector n_umi(n_cells), p_mito(n_cells), y(n_genes), frac_umis_max_feature(n_cells);
  List out = List::create(Named("n_nonzero") = n_nonzero, Named("n_umi") = n_umi, Named("p_mito") = p_mito,
                          Named("max_feature") = max_feature, Named("frac_umis_max_feature") = frac_umis_max_feature);
  const int* mt_rows = nullptr;
  R_xlen_t n_mt = 0;
  if (compute_p_mito) {
    mt_rows = INTEGER(mt_gene_idxs);
    n_mt = mt_gene_idxs.size();
    for (R_xlen_t row = 0; row < n_mt; ++row)
      if (mt_rows[row] < 1 || mt_rows[row] > n_genes)
        stop("Mitochondrial row is outside the one-based feature range.");
  }
  const int* indices = INTEGER(i);
  const int* pointers = INTEGER(p);
  const double* values = REAL(x);
  int* nonzero_data = INTEGER(n_nonzero);
  int* max_feature_data = INTEGER(max_feature);
  double* umi_data = REAL(n_umi);
  double* mito_data = REAL(p_mito);
  double* loaded = REAL(y);
  double* fraction_data = REAL(frac_umis_max_feature);
  {
    std::vector<int> owned_mt(n_mt);
    for (R_xlen_t row = 0; row < n_mt; ++row) owned_mt[row] = mt_rows[row] - 1;
    for (int k = 0; k < n_cells; k ++) {
      int curr_max_feature = 0, n_nonzero_count = 0;
      double umi_count_curr_max_feature = 0, n_umi_count = 0;
      for (int row = 0; row < n_genes; ++row) loaded[row] = 0;
      for (int entry = pointers[k]; entry < pointers[k + 1]; ++entry)
        loaded[indices[entry]] = values[entry];
      for (int r = 0; r < n_genes; r ++) {
        n_umi_count += loaded[r];
        if (loaded[r] > 0.5) n_nonzero_count ++;
      }

      if (compute_p_mito) {
        double n_umi_count_mito = 0;
        for (int row : owned_mt) n_umi_count_mito += loaded[row];
        mito_data[k] = n_umi_count_mito/n_umi_count;
      }

      if (compute_max_feature) {
        for (int r = 0; r < n_genes; r ++) {
          if (loaded[r] > umi_count_curr_max_feature) {
            curr_max_feature = r;
            umi_count_curr_max_feature = loaded[r];
          }
        }
        max_feature_data[k] = curr_max_feature;
        fraction_data[k] = umi_count_curr_max_feature/n_umi_count;
      }

      nonzero_data[k] = n_nonzero_count;
      umi_data[k] = n_umi_count;
    }
  }
  return out;
}


// Count how many guide memberships each cell has. Each group lists cell numbers
// starting at 1 among n_cells; the assignment caller has already checked them.
// Repeated memberships contribute repeatedly rather than being removed.
// Return a new integer vector, leaving input group lists unchanged.
// [[Rcpp::export]]
IntegerVector compute_n_grnas_per_cell_vector(List grna_assignments, int n_cells) {
  IntegerVector out(n_cells), current_assignments;
  for (int i = 0; i < grna_assignments.size(); i ++) {
    current_assignments = grna_assignments[i];
    for (int j = 0; j < current_assignments.size(); j ++) {
      out[current_assignments[j] - 1] ++;
    }
  }
  return out;
}


// Find cells whose stored value in one sparse feature row is >= threshold.
// row_idx and returned original cell-column numbers start at 1. Keep storage
// order and compare fractional values directly. Explicitly stored zeros can
// qualify when threshold <= 0; unstored zeros are not visited. Count matches
// first, then allocate the exact output size. Leave all source slots unchanged.
// [[Rcpp::export]]
IntegerVector threshold_count_matrix(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, double threshold) {
  IntegerVector out;
  int n_gte_threshold = 0, counter = 0, start = p[row_idx - 1], end = p[row_idx];
  for (int k = start; k < end; k ++) if (x[k] >= threshold) n_gte_threshold ++;
  out = IntegerVector(n_gte_threshold);
  for (int k = start; k < end; k ++) if (x[k] >= threshold) out[counter ++] = j[k] + 1;
  return out;
}
