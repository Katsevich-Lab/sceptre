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
  IntegerVector p(dim + 1);
  p[0] = 0;
  int curr_idx = 0, counter = 0, j = 0;
  for (curr_idx = 0; curr_idx < dim; curr_idx ++) {
    while (i[j] == curr_idx) {
      counter ++;
      j ++;
    }
    p[curr_idx + 1] = counter;
  }
  return p;
}


// [[Rcpp::export]]
List compute_cell_covariates_cpp(IntegerVector i, IntegerVector p, NumericVector x, int n_genes, int n_cells, IntegerVector mt_gene_idxs, bool compute_p_mito) {
  IntegerVector n_nonzero(n_cells);
  NumericVector n_umi(n_cells), p_mito(n_cells), y(n_genes);
  List out;

  int n_nonzero_count;
  double n_umi_count, n_umi_count_mito;
  for (int k = 0; k < mt_gene_idxs.size(); k ++) mt_gene_idxs[k] --;

  for (int k = 0; k < n_cells; k ++) {
    n_umi_count = 0;
    n_nonzero_count = 0;
    load_sparse_vector(i, p, x, k, n_genes, y);
    for (int r = 0; r < n_genes; r ++) {
      n_umi_count += y[r];
      if (y[r] > 0.5) n_nonzero_count ++;
    }

    if (compute_p_mito) {
      n_umi_count_mito = 0;
      for (int r = 0; r < mt_gene_idxs.size(); r ++) n_umi_count_mito += y[mt_gene_idxs[r]];
      p_mito[k] = n_umi_count_mito/n_umi_count;
    }

    n_nonzero[k] = n_nonzero_count;
    n_umi[k] = n_umi_count;
  }

  if (compute_p_mito) {
    out = List::create(Named("n_nonzero") = n_nonzero, Named("n_umi") = n_umi, Named("p_mito") = p_mito);
  } else {
    out = List::create(Named("n_nonzero") = n_nonzero, Named("n_umi") = n_umi);
  }

  return out;
}


// [[Rcpp::export]]
IntegerVector compute_colwise_max(IntegerVector i, IntegerVector p, NumericVector x, int n_cells) {
  IntegerVector out(n_cells);
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
    out[k] = curr_maximizer + 1;
  }
  return(out);
}


// [[Rcpp::export]]
IntegerVector group_and_threshold(IntegerVector j, IntegerVector p, NumericVector x, IntegerVector row_idxs, int threshold) {
  int row_idx, col_idx, start, end, counter = 0;
  double count;
  std::unordered_set<int> set;
  for (int i = 0; i < row_idxs.size(); i++) {
    row_idx = row_idxs[i] - 1;
    int start = p[row_idx], end = p[row_idx + 1];
    for (int k = start; k < end; k ++) {
      count = x[k];
      col_idx = j[k];
      if (count >= threshold) {
        set.insert(col_idx + 1);
      }
    }
  }
  IntegerVector out(set.size());
  std::unordered_set<int>::iterator itr;
  for (itr = set.begin(); itr != set.end(); itr++) {
    out[counter] = *itr;
    counter ++;
  }
  return(out);
}
