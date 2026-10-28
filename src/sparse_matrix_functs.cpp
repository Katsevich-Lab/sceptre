#include <Rcpp.h>
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


// [[Rcpp::export]]
List compute_n_nonzero_and_max_feature(IntegerVector i, IntegerVector p, NumericVector x) {
  // for each column of a CSC matrix: the number of entries above 0.5, and the
  // largest positive entry with its row (the first such row; row 0 if no entry
  // is positive)
  int n_cells = p.size() - 1;
  IntegerVector n_nonzero(n_cells), max_feature(n_cells);
  NumericVector max_value(n_cells);
  for (int k = 0; k < n_cells; k ++) {
    for (int l = p[k]; l < p[k + 1]; l ++) {
      if (x[l] > 0.5) n_nonzero[k] ++;
      if (x[l] > max_value[k]) {
        max_value[k] = x[l];
        max_feature[k] = i[l];
      }
    }
  }
  return List::create(Named("n_nonzero") = n_nonzero, Named("max_feature") = max_feature,
                      Named("max_value") = max_value);
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
void increment_vector(IntegerVector x, int value) {
  for (int i = 0; i < x.size(); i++) x[i] += value;
  return;
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
