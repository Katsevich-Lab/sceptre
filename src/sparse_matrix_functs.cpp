#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

void load_sparse_vector(IntegerVector l, IntegerVector p, NumericVector x, int idx, int dim, NumericVector out) {
  int start = p[idx], end = p[idx + 1];
  for (int k = 0; k < dim; k ++) out[k] = 0.0;
  for (int k = start; k < end; k ++) out[l[k]] = x[k];
}


// [[Rcpp::export]]
NumericVector load_csr_row(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, int n_cells) {
  NumericVector out(n_cells);
  load_sparse_vector(j, p, x, row_idx - 1, n_cells, out);
  return out;
}


// [[Rcpp::export]]
List compute_cell_covariates_cpp(IntegerVector i, IntegerVector p, NumericVector x, int n_genes,
                                 int n_cells, IntegerVector mt_gene_idxs, bool compute_p_mito, bool compute_max_feature) {
  IntegerVector n_nonzero(n_cells), max_feature(n_cells);
  NumericVector n_umi(n_cells), p_mito(n_cells), frac_umis_max_feature(n_cells);
  List out;

  // visit only the stored entries of each column, in increasing row order; as
  // the counts are nonnegative, the sums, the first maximum and so every
  // covariate equal those over the full (dense) column
  vector<bool> is_mt(n_genes, false);
  for (int k = 0; k < mt_gene_idxs.size(); k ++) is_mt[mt_gene_idxs[k] - 1] = true;

  int n_nonzero_count, curr_max_feature;
  double n_umi_count, n_umi_count_mito, umi_count_curr_max_feature, y;
  for (int k = 0; k < n_cells; k ++) {
    curr_max_feature = 0;
    umi_count_curr_max_feature = 0;
    n_umi_count = 0;
    n_umi_count_mito = 0;
    n_nonzero_count = 0;
    for (int l = p[k]; l < p[k + 1]; l ++) {
      y = x[l];
      n_umi_count += y;
      if (y > 0.5) n_nonzero_count ++;
      if (is_mt[i[l]]) n_umi_count_mito += y;
      if (y > umi_count_curr_max_feature) {
        curr_max_feature = i[l];
        umi_count_curr_max_feature = y;
      }
    }

    // compute p mito
    if (compute_p_mito) p_mito[k] = n_umi_count_mito/n_umi_count;

    // compute max feature, fraction expressed
    if (compute_max_feature) {
      max_feature[k] = curr_max_feature;
      frac_umis_max_feature[k] = umi_count_curr_max_feature/n_umi_count;
    }

    n_nonzero[k] = n_nonzero_count;
    n_umi[k] = n_umi_count;
  }

  out = List::create(Named("n_nonzero") = n_nonzero, Named("n_umi") = n_umi, Named("p_mito") = p_mito,
                     Named("max_feature") = max_feature, Named("frac_umis_max_feature") = frac_umis_max_feature);
  return out;
}


// compute_cell_covariates_cpp() (without the maximum feature) for counts held
// as cells x features, column-compressed: the features are the columns
// `feature_cols` (1-based, increasing) and `mt_feature_idxs` (1-based) indexes
// into them. Each cell's entries are visited in increasing feature order, as
// compute_cell_covariates_cpp() visits them, so the covariates are the same.
// [[Rcpp::export]]
List compute_cell_covariates_by_feature_cpp(IntegerVector i, IntegerVector p, NumericVector x, IntegerVector feature_cols,
                                            int n_cells, IntegerVector mt_feature_idxs, bool compute_p_mito) {
  IntegerVector n_nonzero(n_cells);
  NumericVector n_umi(n_cells), n_umi_mito(n_cells), p_mito(n_cells);
  vector<bool> is_mt(feature_cols.size(), false);
  for (int k = 0; k < mt_feature_idxs.size(); k ++) is_mt[mt_feature_idxs[k] - 1] = true;

  double y;
  for (int k = 0; k < feature_cols.size(); k ++) {
    int col = feature_cols[k] - 1;
    for (int l = p[col]; l < p[col + 1]; l ++) {
      y = x[l];
      n_umi[i[l]] += y;
      if (y > 0.5) n_nonzero[i[l]] ++;
      if (is_mt[k]) n_umi_mito[i[l]] += y;
    }
  }
  if (compute_p_mito) for (int c = 0; c < n_cells; c ++) p_mito[c] = n_umi_mito[c]/n_umi[c];
  return List::create(Named("n_nonzero") = n_nonzero, Named("n_umi") = n_umi, Named("p_mito") = p_mito);
}


List compute_colwise_max(IntegerVector i, IntegerVector p, NumericVector x, int n_cells, NumericVector grna_lib_size) {
  IntegerVector assignment_vect(n_cells);
  NumericVector frac_umis(n_cells);
  int p_start, p_end, curr_maximizer;
  double curr_max;

  for (int k = 0; k < n_cells; k ++) {
    p_start = p[k];
    p_end = p[k+1];
    curr_max = -1.0;
    curr_maximizer = -1;
    for (int l = p_start; l < p_end; l ++) {
      if (x[l] > curr_max) {
        curr_max = x[l];
        curr_maximizer = i[l];
      }
    }
    assignment_vect[k] = curr_maximizer + 1;
    frac_umis[k] = curr_max/grna_lib_size[k];
  }
  return(List::create(Named("assignment_vect") = assignment_vect, Named("frac_umis") = frac_umis));
}


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


// [[Rcpp::export]]
IntegerVector threshold_count_matrix(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, double threshold) {
  IntegerVector out;
  int n_gte_threshold = 0, counter = 0, start = p[row_idx - 1], end = p[row_idx];
  for (int k = start; k < end; k ++) if (x[k] >= threshold) n_gte_threshold ++;
  out = IntegerVector(n_gte_threshold);
  for (int k = start; k < end; k ++) if (x[k] >= threshold) out[counter ++] = j[k] + 1;
  return out;
}


// The columns of several column-compressed matrices with the same columns,
// stacked by row: matrix k (row indices i_list[[k]], values x_list[[k]],
// column pointers p_list[[k]]) supplies rows row_offsets[k] + 1, ...; within a
// column, the entries of each matrix follow those of the matrices before it.
// The values are integer counts.
// [[Rcpp::export]]
List stack_csc_rows(List i_list, List x_list, List p_list, IntegerVector row_offsets, int n_cols) {
  int n_matrices = i_list.size();
  std::vector<IntegerVector> is(n_matrices), xs(n_matrices), ps(n_matrices);
  long long n_entries = 0;
  for (int k = 0; k < n_matrices; k ++) {
    is[k] = i_list[k];
    xs[k] = x_list[k];
    ps[k] = p_list[k];
    n_entries += is[k].size();
  }
  IntegerVector i(n_entries), x(n_entries), p(n_cols + 1);
  int l = 0;
  for (int c = 0; c < n_cols; c ++) {
    for (int k = 0; k < n_matrices; k ++) {
      for (int m = ps[k][c]; m < ps[k][c + 1]; m ++) {
        i[l] = is[k][m] + row_offsets[k];
        x[l] = xs[k][m];
        l ++;
      }
    }
    p[c + 1] = l;
  }
  return List::create(Named("i") = i, Named("p") = p, Named("x") = x);
}


// Whether the values x[from], ..., x[from + count - 1] are counts, in one pass
// and without copying them: the first problem found, in the order that
// check_integer_counts() reports them (1: missing, 2: negative, 3: infinite,
// 4: not an integer, 0: none), the number of non-integer values and the first
// of them.
// [[Rcpp::export]]
List scan_count_values(NumericVector x, double from, double count) {
  bool missing = false, negative = false, infinite = false;
  double n_non_integer = 0, example = NA_REAL;
  R_xlen_t start = static_cast<R_xlen_t>(from), end = start + static_cast<R_xlen_t>(count);
  for (R_xlen_t k = start; k < end; k ++) {
    double v = x[k];
    if (std::isnan(v)) {
      missing = true;
    } else if (v < 0) {
      negative = true;
    } else if (std::isinf(v)) {
      infinite = true;
    } else if (v != std::round(v)) {
      if (n_non_integer == 0) example = v;
      n_non_integer ++;
    }
  }
  int problem = missing ? 1 : negative ? 2 : infinite ? 3 : n_non_integer > 0 ? 4 : 0;
  return List::create(Named("problem") = problem, Named("n_non_integer") = n_non_integer,
                      Named("example") = example);
}
